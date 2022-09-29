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
    kwargs...,
)
    α = rₒ * cos(θₒ)
    β = rₒ * sin(θₒ)
    v = map_impact_parameters(m, u, α, β)
    sol = tracegeodesics(m, u, v, d, (0.0, 2000.0); save_on = false, kwargs...)
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
    f(r) = begin
        gp = integrate_single_geodesic(m, u, d, r, θₒ; kwargs...)
        gp.u2[2] - radius
    end

    Roots.find_zero(f, (0.0, offset_max); atol = zero_atol)
end

function impact_parameters_for_radius(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionDisc,
    radius;
    N = 500,
    kwargs...,
)
    θs = range(0, 2π, N)
    rs = ThreadsX.map(θs) do θ
        find_offset_for_radius(m, u, d, radius, θ; kwargs...)
    end
    α = @. rs * cos(θs)
    β = @. rs * sin(θs)
    (α, β)
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
    order = 3,
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    kwargs...,
)
    gmin, gmax = extrema(gs)
    gstar(g) = (g - gmin) / (gmax - gmin)

    f((α, β)) = begin
        v = map_impact_parameters(m, u, α, β)
        sol = tracegeodesics(m, u, v, d, (0.0, max_time); save_on = false, kwargs...)
        gp = getgeodesicpoint(m, sol)
        g = redshift_pf(m, gp, 2000.0)
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

function gstar(g::AbstractArray)
    gmin, gmax = extrema(g)
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
    finite_diff_order = 3,
    redshift_pf::PointFunction = Gradus.ConstPointFunctions.redshift,
    tracer_kwargs...,
)
    αs, βs =
        Gradus.impact_parameters_for_radius(m, u, d, rₑ; N = num_points, tracer_kwargs...)

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

export impact_parameters_for_radius,
    redshift_ratio,
    jacobian_∂αβ_∂gr,
    gstar,
    cunningham_transfer_function,
    CunninghamTransferFunction
