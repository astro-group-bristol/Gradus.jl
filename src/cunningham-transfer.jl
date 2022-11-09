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
        r = if gp.status == StatusCodes.IntersectedWithGeometry
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

function _find_extremal_redshift_with_guess(
    m,
    u,
    d,
    rₑ,
    g_guess,
    θ_guess;
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    δθ = π,
    minimal = true,
    zero_atol = 1e-7,
    offset_max = 20.0,
    kwargs...,
)
    f =
        θ -> begin
            r_offset = find_offset_for_radius(
                m,
                u,
                d,
                rₑ,
                θ;
                zero_atol = zero_atol,
                offset_max = offset_max,
                kwargs...,
            )
            if !isnan(r_offset)
                gp = integrate_single_geodesic(m, u, d, r_offset, θ; kwargs...)
                g = redshift_pf(m, gp, gp.t2)
                minimal ? g : -g
            else
                minimal ? 100.0 : -100.0 
                # minimal ? g_guess : -g_guess
            end
        end

    res = optimize(f, θ_guess - δθ, θ_guess + δθ, Optim.GoldenSection())
    θopt = Optim.minimizer(res)
    (
        res,
        find_offset_for_radius(
            m,
            u,
            d,
            rₑ,
            θopt;
            zero_atol = zero_atol,
            offset_max = offset_max,
            kwargs...,
        ),
    )
end

function impact_parameters_for_radius!(
    α,
    β,
    m::AbstractMetricParams,
    u::AbstractVector{T},
    d::AbstractAccretionDisc,
    radius;
    kwargs...,
) where {T}
    if size(α) != size(β)
        throw(DimensionMismatch("α, β must have the same dimensions and size."))
    end
    θs = range(0, 2π, length(α))
    @inbounds @threads for i in eachindex(θs)
        θ = θs[i]
        r = find_offset_for_radius(m, u, d, radius, θ; kwargs...)
        α[i] = r * cos(θ)
        β[i] = r * sin(θ)
    end
    (α, β)
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
    impact_parameters_for_radius!(α, β, m, u, d, radius; kwargs...)
    (filter(!isnan, α), filter(!isnan, β))
end

function redshift_ratio!(
    gs,
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    max_time,
    αs,
    βs;
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    kwargs...,
)
    if size(gs) != size(αs)
        throw(DimensionMismatch("αs, gs must have the same dimensions and size."))
    end

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

    @inbounds for (i, sol) in enumerate(simsols.u)
        gp = getgeodesicpoint(m, sol)
        g = redshift_pf(m, gp, max_time)
        gs[i] = g
    end
    gs
end

function redshift_ratio(
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    max_time,
    αs::AbstractVector{T},
    βs;
    kwargs...,
) where {T}
    gs = zeros(T, length(αs))
    redshift_ratio!(gs, m, u, d, max_time, αs, βs; kwargs...)
    gs
end

function jacobian_∂αβ_∂gr!(
    Js,
    m,
    u,
    d,
    max_time,
    gs,
    gmin,
    gmax,
    αs,
    βs;
    order = 5,
    redshift_pf = Gradus.ConstPointFunctions.redshift,
    kwargs...,
)
    if size(gs) != size(Js) || size(Js) != size(αs) || size(αs) != size(βs)
        throw(
            DimensionMismatch("αs, βs, gs, and Js must have the same dimensions and size."),
        )
    end

    f =
        ((α, β),) -> begin
            v = map_impact_parameters(m, u, α, β)
            sol =
                tracegeodesics(m, u, v, d, (0.0, max_time); save_on = false, kwargs...)
            gp = getgeodesicpoint(m, sol)
            g = redshift_pf(m, gp, 2000.0)
            # return r and g*
            @SVector [gp.u2[2], gstar(g, gmin, gmax)]
        end

    # choice between FiniteDifferences and FiniteDiff is tricky
    # since FiniteDiff is so much faster, but seems to yield really bad jacobians
    # for this specific problem, so instead stenciling with a given order
    cfdm = FiniteDifferences.central_fdm(order, 1)
    @inbounds @threads for i in eachindex(αs)
        x = @SVector [αs[i], βs[i]]
        j = FiniteDifferences.jacobian(cfdm, f, x) |> first
        Js[i] = abs(inv(det(j)))
    end
    Js
end

function jacobian_∂αβ_∂gr(
    m,
    u,
    d,
    max_time,
    gs::AbstractArray{T},
    gmin,
    gmax,
    αs::AbstractArray{T},
    βs;
    kwargs...,
) where {T}
    Js = zeros(T, size(αs))
    jacobian_∂αβ_∂gr!(Js, m, u, d, max_time, gs, gmin, gmax, αs, βs; kwargs...)
    Js
end

gstar(g::AbstractArray) = gstar(g, extrema(g)...)
function gstar(g, gmin, gmax)
    Δg = gmax - gmin
    @. (g - gmin) / Δg
end

g_to_gstar(args...) = gstar(args...)
function gstar_to_g(gs, gmin, gmax)
    Δg = gmax - gmin
    @. gs * Δg + gmin
end

@muladd cunningham_transfer_function(
    rₑ::Number,
    gs::AbstractArray,
    gstars::AbstractArray,
    js::AbstractArray,
) = @. (1 / (π * rₑ)) * gs * √(gstars * (1 - gstars)) * js

function cunningham_transfer_function!(
    αs,
    βs,
    Js::AbstractVector{T},
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    rₑ,
    max_time;
    finite_diff_order = 5,
    redshift_pf::PointFunction = ConstPointFunctions.redshift,
    offset_max = 20.0,
    zero_atol = 1e-7,
    tracer_kwargs...,
) where {T}
    impact_parameters_for_radius!(
        αs,
        βs,
        m,
        u,
        d,
        rₑ;
        offset_max = offset_max,
        zero_atol = zero_atol,
        max_time = max_time,
        tracer_kwargs...,
    )
    # filter and trim arrays
    mask = @. !(isnan(αs) | isnan(βs))
    _αs = @views αs[mask]
    _βs = @views βs[mask]

    N = length(_αs)
    @assert N > 0 "No valid impact parameters found (all are NaN)."
    _Js = @views Js[1:N]

    gs = zeros(T, N)
    redshift_ratio!(
        gs,
        m,
        u,
        d,
        max_time,
        _αs,
        _βs;
        redshift_pf = redshift_pf,
        tracer_kwargs...,
    )

    gmin, gmax = extrema(gs)
    jacobian_∂αβ_∂gr!(
        _Js,
        m,
        u,
        d,
        max_time,
        gs,
        gmin,
        gmax,
        _αs,
        _βs;
        order = finite_diff_order,
        redshift_pf = redshift_pf,
        tracer_kwargs...,
    )

    gstars = gstar(gs, gmin, gmax)
    f = cunningham_transfer_function(rₑ, gs, gstars, _Js)

    # package and return
    CunninghamTransferFunction(rₑ, gs, gstars, f)
end

function cunningham_transfer_function(
    m::AbstractMetricParams{T},
    u,
    d::AbstractAccretionGeometry,
    rₑ,
    max_time;
    num_points = 1000,
    kwargs...,
) where {T}
    α = zeros(T, num_points)
    β = zeros(T, num_points)
    Js = zeros(T, num_points)
    cunningham_transfer_function!(α, β, Js, m, u, d, rₑ, max_time; kwargs...)
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
        if decreasing
            push!(upper, (g, f[i]))
        else
            push!(lower, (g, f[i]))
        end
        gprev = g
    end
    lower, upper
end

struct InterpolatedCunninghamTransferBranches{T,I}
    lower::I
    upper::I
    radius::T
    g_extrema::Tuple{T,T}
    g_limits::Tuple{T,T}
end

_gs_in_limits(ictb, gs) = ictb.g_limits[2] > gs > ictb.g_limits[1]

function _make_sorted_interpolation(g, f)
    I = sortperm(g)
    _g = @inbounds g[I]
    _f = @inbounds f[I]
    # Interpolations.deduplicate_knots!(_g)
    # linear_interpolation(_g, _f, extrapolation_bc = NaN)
    DataInterpolations.LinearInterpolation(_f, _g)
end

function _interpolate_branches(ctf::CunninghamTransferFunction; offset = 1e-7)
    # avoid extrema
    mask = @. (ctf.gstar > offset) & (ctf.gstar < 1 - offset)

    gstars = @inbounds @views ctf.gstar[mask]
    f = @inbounds @views ctf.f[mask]

    lower, upper = _split_branches(gstars, f)

    gs_lower = first.(lower)
    gs_upper = first.(upper)
    # valid interval
    gs_min_limit = max(minimum(gs_lower), minimum(gs_upper))
    gs_max_limit = min(maximum(gs_lower), maximum(gs_upper))

    # interpolate
    f1 = _make_sorted_interpolation(gs_lower, last.(lower))
    f2 = _make_sorted_interpolation(gs_upper, last.(upper))

    InterpolatedCunninghamTransferBranches(
        f1,
        f2,
        ctf.rₑ,
        extrema(ctf.gs),
        (gs_min_limit, gs_max_limit),
    )
end

@muladd _cunning_integrand(f, g, rₑ, gs, gmin, gmax) =
    f * g^3 * π * rₑ / (√(gs * (1 - gs)) * (gmax - gmin))

function integrate_transfer_functions(
    ε,
    ictbs::Vector{<:InterpolatedCunninghamTransferBranches{T1}},
    bins,
) where {T1}
    # global min and max
    ggmin = maximum(i -> first(i.g_limits), ictbs)
    ggmax = minimum(i -> last(i.g_limits), ictbs)

    radii = map(i -> i.radius, ictbs)

    r_low, r_high = extrema(radii)

    segbuf = QuadGK.alloc_segbuf(T1, eltype(bins))
    y = map(eachindex(@view(bins[1:end-1]))) do index
        # bin edges
        bin_low = bins[index]
        bin_high = bins[index+1]
        # calculate area under bin
        integrated_bin_for_radii =
            _areas_under_transfer_functions(ε, ictbs, bin_low, bin_high, ggmin, ggmax)

        # interpolate across different radii
        intp = DataInterpolations.LinearInterpolation(integrated_bin_for_radii, radii)
        # use (smooth) interpolation to integrate fully
        res::T1, _ = quadgk(intp, r_low, r_high; segbuf = segbuf, order = 4)
        res
    end

    (bins, y)
end

@muladd function _cunningham_line_profile_integrand(ictb)
    gs -> begin
        if !Gradus._gs_in_limits(ictb, gs)
            0.0
        else
            gmin, gmax = ictb.g_extrema
            g = gstar_to_g(gs, gmin, gmax)
            w = g^3 / √(gs * (1 - gs))
            f = ictb.lower(gs) + ictb.upper(gs)
            w * f / (gmax - gmin)
        end
    end
end

function _calculate_interpolated_transfer_branches(
    m::AbstractMetricParams{T},
    u,
    d,
    radii;
    num_points = 100,
    verbose = false,
    offset = 1e-7,
    max_time = 2000.0,
    kwargs...,
) where {T}
    # pre-allocate arrays
    αs = zeros(T, num_points)
    βs = zeros(T, num_points)
    Js = zeros(T, num_points)

    progress_bar = init_progress_bar("Transfer functions:", length(radii), verbose)

    # IILF for calculating the interpolated branches
    𝔉 =
        rₑ -> begin
            # this feels like such a bad practice
            # calling GC on young objects to clear the temporarily allocated memory
            # but we get a significant speedup
            GC.gc(false)
            ctf = cunningham_transfer_function!(
                αs,
                βs,
                Js,
                m,
                u,
                d,
                rₑ,
                max_time;
                offset_max = 2rₑ + 20.0,
                kwargs...,
            )
            ProgressMeter.next!(progress_bar)
            _interpolate_branches(ctf; offset = offset)
        end

    # calculate interpolated transfer functions for each emission radius
    map(𝔉, radii)
end

function _areas_under_transfer_functions(
    ε,
    ictbs::Vector{<:InterpolatedCunninghamTransferBranches{T1}},
    bin_low::T2,
    bin_high::T2,
    ggmin,
    ggmax,
) where {T1,T2}
    segbuf = QuadGK.alloc_segbuf(T1, T2)
    map(ictbs) do ictb
        𝒮 = _cunningham_line_profile_integrand(ictb)
        Sg = g -> begin
            gs = gstar(g, ictb.g_extrema...)
            if ggmin < gs < ggmax
                𝒮(gs)
            else
                0.0
            end
        end
        res, _ = quadgk(Sg, bin_low, bin_high; order = 2, segbuf = segbuf)
        r = ictb.radius
        res * r * ε(r)
    end
end

export impact_parameters_for_radius,
    redshift_ratio,
    jacobian_∂αβ_∂gr,
    gstar,
    g_to_gstar,
    gstar_to_g,
    cunningham_transfer_function,
    CunninghamTransferFunction,
    InterpolatedCunninghamTransferBranches
