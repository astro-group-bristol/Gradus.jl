@with_kw struct CunninghamTransferFunction{T}
    "``g^\\ast`` values."
    g✶::Vector{T}
    "Transfer function data."
    f::Vector{T}
    gmin::T
    gmax::T
    "Emission radius."
    rₑ::T
end

@with_kw struct InterpolatedCunninghamTransferFunction{T,U,L}
    upper_f::U
    lower_f::L
    gmin::T
    gmax::T
    rₑ::T
end

function _make_sorted_interpolation(g, f)
    I = sortperm(g)
    _g = @inbounds g[I]
    _f = @inbounds f[I]
    DataInterpolations.LinearInterpolation(_f, _g)
end

function splitbranches(ctf::CunninghamTransferFunction{T}) where {T}
    upper_f = T[]
    lower_f = T[]
    upper_g✶ = T[]
    lower_g✶ = T[]
    N = (length(ctf.f) ÷ 2) + 1
    sizehint!(upper_f, N)
    sizehint!(lower_f, N)
    sizehint!(upper_g✶, N)
    sizehint!(lower_g✶, N)

    decreasing = true
    g✶previous = 0.0
    for (i, g✶) in enumerate(ctf.g✶)
        if g✶previous > g✶
            decreasing = true
        else
            decreasing = false
        end
        if decreasing
            push!(lower_f, ctf.f[i])
            push!(lower_g✶, g✶)
        else
            push!(upper_f, ctf.f[i])
            push!(upper_g✶, g✶)
        end
        g✶previous = g✶
    end

    (lower_g✶, lower_f, upper_g✶, upper_f)
end

function interpolate_transfer_function(ctf::CunninghamTransferFunction{T}) where {T}
    (lower_g✶, lower_f, upper_g✶, upper_f) = splitbranches(ctf)
    lower_branch = _make_sorted_interpolation(lower_g✶, lower_f)
    upper_branch = _make_sorted_interpolation(upper_g✶, upper_f)
    InterpolatedCunninghamTransferFunction(
        upper_branch,
        lower_branch,
        ctf.gmin,
        ctf.gmax,
        ctf.rₑ,
    )
end

@muladd function _calculate_transfer_function(rₑ, g, g✶, J)
    @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J
end

function cunningham_transfer_function(
    m::AbstractMetricParams{T},
    u,
    d,
    rₑ;
    max_time = 2e3,
    diff_order = 4,
    redshift_pf = ConstPointFunctions.redshift,
    offset_max = 20.0,
    zero_atol = 1e-7,
    N = 100,
    tracer_kwargs...,
) where {T}
    Js = zeros(T, N)
    gs = zeros(T, N)

    θs = range(π / 2, 2π + π / 2, N)
    @inbounds for i in eachindex(θs)
        θ = θs[i]
        r, gp = find_offset_for_radius(
            m,
            u,
            d,
            rₑ,
            θ;
            zero_atol = zero_atol,
            offset_max = offset_max,
            max_time = max_time,
            tracer_kwargs...,
        )
        α = r * cos(θ)
        β = r * sin(θ)
        g = redshift_pf(m, gp, max_time)
        gs[i] = g
        Js[i] = jacobian_∂αβ_∂gr(
            m,
            u,
            d,
            α,
            β,
            max_time;
            diff_order = diff_order,
            redshift_pf = redshift_pf,
            tracer_kwargs...,
        )
    end

    gmin, gmax = infer_extremal(gs, θs, π, 2π)
    # convert from ∂g to ∂g✶
    @. Js = (gmax - gmin) * Js
    @inbounds for i in eachindex(Js)
        g✶ = g_to_g✶(gs[i], gmin, gmax)
        # Js is now storing f
        Js[i] = _calculate_transfer_function(rₑ, gs[i], g✶, Js[i])
        # gs is now storing g✶
        gs[i] = g✶
    end

    CunninghamTransferFunction(gs, Js, gmin, gmax, rₑ)
end

function infer_extremal(y, x, x0, x1)
    _, y0 = interpolate_extremal(y, x, x0)
    _, y1 = interpolate_extremal(y, x, x1)
    if y0 > y1
        return y1, y0
    else
        return y0, y1
    end
end

function interpolate_extremal(y, x, x0)
    interp = DataInterpolations.CubicSpline(y, x)
    ∂(f) = x -> ForwardDiff.derivative(f, x)
    x̄ = find_zero(∂(interp), x0)
    x̄, interp(x̄)
end

function interpolated_transfer_branches(
    m::AbstractMetricParams{T},
    u,
    d,
    radii;
    verbose = false,
    kwargs...,
) where {T}
    # IILF for calculating the interpolated branches
    𝔉 =
        rₑ -> begin
            ctf = cunningham_transfer_function(
                m,
                u,
                d,
                rₑ,
                ;
                offset_max = 2rₑ + 20.0,
                kwargs...,
            )
            interpolate_transfer_function(ctf)
        end

    # calculate interpolated transfer functions for each emission radius
    ThreadsX.map(𝔉, radii)
end
