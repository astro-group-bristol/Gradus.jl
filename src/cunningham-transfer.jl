@with_kw struct CunninghamTransferFunction{T}
    rₑ::T
    gs::Vector{T}
    gstar::Vector{T}
    f::Vector{T}
end

"""
    integrate_single_geodesic(m::AbstractMetricParams, u, d::AbstractAccretionDisc, rₒ, θₒ; kwargs...)

Integrate a single geodesic with impact parameters calculated via

```math
\\begin{align}
    \\alpha &= r_\\text{o} \\cos (\\theta_\\text{o}), \\
    \\beta &= r_\\text{o} \\sin (\\theta_\\text{o}).
\\end{align}
```

Returns an [AbstractGeodesicPoint](@ref), depending on `m`.

"""
function integrate_single_geodesic(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionDisc,
    rₒ,
    θₒ;
    max_time = 2e3,
    kwargs...,
)
    α = rₒ * cos(θₒ)
    β = rₒ * sin(θₒ)
    v = map_impact_parameters(m, u, α, β)
    sol = tracegeodesics(m, u, v, d, (0.0, max_time); save_on = false, kwargs...)
    getgeodesicpoint(m, sol)
end

function find_offset_for_radius(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionDisc,
    radius,
    θₒ;
    zero_atol = 1e-7,
    offset_max = 20.0,
    kwargs...,
)
    f = r -> begin
        gp = integrate_single_geodesic(m, u, d, r, θₒ; kwargs...)
        r = if gp.retcode == :Intersected
            gp.u2[2] * sin(gp.u2[3])
        else
            0.0
        end
        radius - r
    end

    r = Roots.find_zero(f, (0.0, offset_max); atol = zero_atol)
    if !isapprox(f(r), 0.0, atol = 10 * zero_atol)
        return NaN
    end
    r
end

function impact_parameters_for_radius(
    m::AbstractMetricParams,
    u::AbstractVector{T},
    d::AbstractAccretionDisc,
    radius;
    N = 500,
    kwargs...,
) where {T}
    α = zeros(T, N)
    β = zeros(T, N)

    θs = range(0, 2π, N)
    ThreadsX.foreach(enumerate(θs)) do (i, θ)
        r = find_offset_for_radius(m, u, d, radius, θ; kwargs...)
        α[i] = r * cos(θ)
        β[i] = r * sin(θ)
    end
    (filter(!isnan, α), filter(!isnan, β))
end

function redshift_ratio(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    max_time,
    αs,
    βs;
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    kwargs...,
)
    velfunc(i) = map_impact_parameters(m, u, αs[i], βs[i])
    simsols = tracegeodesics(
        m,
        u,
        velfunc,
        d,
        (0.0, max_time);
        trajectories = length(αs),
        save_on = false,
        kwargs...,
    )
    map(simsols) do sol
        gp = getgeodesicpoint(m, sol)
        redshift_pf(m, gp, max_time)
    end
end

function jacobian_∂αβ_∂gr(
    m,
    u,
    d,
    max_time,
    gs,
    αs,
    βs;
    order = 5,
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    kwargs...,
)
    gmin, gmax = extrema(gs)
    gstar(g) = (g - gmin) / (gmax - gmin)

    f =
        ((α, β),) -> begin
            v = map_impact_parameters(m, u, α, β)
            sol =
                tracegeodesics(m, u, v, d, (0.0, max_time); save_on = false, kwargs...)
            gp = getgeodesicpoint(m, sol)
            g = redshift_pf(m, gp, 2000.0)
            # return r and g*
            @SVector [gp.u2[2], gstar(g)]
        end

    # choice between FiniteDifferences and FiniteDiff is tricky
    # since FiniteDiff is so much faster, but seems to yield really bad jacobians
    # for this specific problem, so instead stenciling with a given order
    cfdm = FiniteDifferences.central_fdm(order, 1)
    ThreadsX.map(zip(αs, βs)) do (α, β)
        x = @SVector [α, β]
        j = FiniteDifferences.jacobian(cfdm, f, x) |> first
        abs(inv(det(j)))
    end
end

gstar(g::AbstractArray) =
    gstar(g, extrema(g)...)
function gstar(g, gmin, gmax)
    Δg = gmax - gmin
    @. (g - gmin) / Δg
end

cunningham_transfer_function(
    rₑ::Number,
    gs::AbstractArray,
    gstars::AbstractArray,
    js::AbstractArray,
) = @. (1 / (π * rₑ)) * gs * √(gstars * (1 - gstars)) * js

function cunningham_transfer_function(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    rₑ,
    max_time;
    num_points = 1000,
    finite_diff_order = 6,
    redshift_pf::PointFunction = Gradus.ConstPointFunctions.redshift,
    offset_max = 20.0,
    zero_atol = 1e-7,
    tracer_kwargs...,
)
    αs, βs = Gradus.impact_parameters_for_radius(
        m,
        u,
        d,
        rₑ;
        N = num_points,
        offset_max = offset_max,
        zero_atol = zero_atol,
        tracer_kwargs...,
    )

    gs = redshift_ratio(
        m,
        u,
        d,
        max_time,
        αs,
        βs;
        redshift_pf = redshift_pf,
        tracer_kwargs...,
    )
    gstars = gstar(gs)

    js = jacobian_∂αβ_∂gr(
        m,
        u,
        d,
        max_time,
        gs,
        αs,
        βs;
        order = finite_diff_order,
        redshift_pf = redshift_pf,
        tracer_kwargs...,
    )

    f = cunningham_transfer_function(rₑ, gs, gstars, js)

    # package and return
    # todo: maybe already interpolate here?
    CunninghamTransferFunction(rₑ, gs, gstars, f)
end

function _split_branches(gstars::AbstractArray{T}, f::AbstractArray{T}) where {T}
    decreasing = true
    gprev = 1.0
    upper = Tuple{T,T}[]
    lower = Tuple{T,T}[]
    for (i, g) in enumerate(gstars)
        if gprev > g
            decreasing = true
        else
            decreasing = false
        end
        gprev = g
        if decreasing
            push!(lower, (g, f[i]))
        else
            push!(upper, (g, f[i]))
        end
    end
    upper, lower
end

struct InterpolatedCunninghamTransferBranches{T,I}
    lower::I
    upper::I
    g_extrema::Tuple{T,T}
    radius::T
end

function _make_sorted_interpolation(g, f)
    I = sortperm(g)
    g = @inbounds g[I]
    f = @inbounds f[I]
    Interpolations.deduplicate_knots!(g)
    linear_interpolation(g, f, extrapolation_bc = Line())
end

function _interpolate_branches(ctf::CunninghamTransferFunction; offset=1e-4)
    # avoid extrema
    mask = @. (ctf.gstar > offset) & (ctf.gstar < 1 - offset)

    gstars = @inbounds @views ctf.gstar[mask]
    f = @inbounds @views ctf.f[mask]

    lower, upper = _split_branches(gstars, f)

    # interpolate
    f1 = _make_sorted_interpolation(first.(lower), last.(lower))
    f2 = _make_sorted_interpolation(first.(upper), last.(upper))

    InterpolatedCunninghamTransferBranches(f1, f2, extrema(ctf.gs), ctf.rₑ)
end

 _cunning_integrand(f, g, rₑ, gs, gmin, gmax) =
    f * g^3 * π * rₑ / (√(gs * (1 - gs)) * (gmax - gmin))

function _integrate_tranfer_function_branches(
    ictb::InterpolatedCunninghamTransferBranches,
    g, a, b; offset = 1e-4
)
    gmin, gmax = ictb.g_extrema

    # limiting approximations
    if (g > gmax) || (g < gmin)
        return 0.0
    elseif (g > gmax - offset) || (g < gmin + offset)
        g_edge = g < gmax ? gmin + offset : gmax - offset

        gs_edge = gstar(g_edge, gmin, gmax)
        f1_edge = ictb.lower(gs_edge)
        f2_edge = ictb.upper(gs_edge)
        norm = _cunning_integrand(f1_edge + f2_edge, g_edge, ictb.radius, gs_edge, gmin, gmax)

        return 2 * norm * √offset * (√b - √a) * (gmax - gmin)
    end

    gs = gstar(g, gmin, gmax)
    f_lower = ictb.lower(gs)
    f_upper = ictb.upper(gs)

    _cunning_integrand(f_lower + f_upper, g, ictb.radius, gs, gmin, gmax)
end

function _integrate_tranfer_function_branches(
    ictbs::AbstractVector{<:InterpolatedCunninghamTransferBranches},
    g, a, b
)
    map(i -> _integrate_tranfer_function_branches(i, g, a, b), ictbs)
end

export impact_parameters_for_radius,
    redshift_ratio,
    jacobian_∂αβ_∂gr,
    gstar,
    cunningham_transfer_function,
    CunninghamTransferFunction,
    InterpolatedCunninghamTransferBranches
