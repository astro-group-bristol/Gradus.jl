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

function splitbranches(ctf::CunninghamTransferFunction)
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
    g✶previous = 1.0
    for (i, g✶) in enumerate(ctf.g✶)
        if g✶previous > g✶
            decreasing = true
        else
            decreasing = false
        end
        if decreasing
            push!(upper_f, ctf.f[i])
            push!(upper_g✶, g✶)
        else
            push!(lower_f, ctf.f[i])
            push!(lower_g✶, g✶)
        end
        g✶previous = g✶
    end
    lower_branch = CunninghamTransferFunction(
        lower_g✶, lower_f, ctf.gmin, ctf.gmax, ctf.rₑ
    )
    upper_branch = CunninghamTransferFunction(
        lower_g✶, lower_f, ctf.gmin, ctf.gmax, ctf.rₑ
    )

    (lower_branch, upper_branch)
end

@muladd function _calculate_transfer_function(rₑ, g, g✶, J)
    @. (1 / (π * rₑ)) * g * √(g✶ * (1 - g✶)) * J
end

function cunningham_transfer_function(
    m::AbstractMetricParams{T}, u, d, rₑ; max_time = 2e3, diff_order = 4,
    redshift_pf = ConstPointFunctions.redshift, offset_max = 20.0,
    zero_atol=1e-7, N=100, tracer_kwargs...
) where {T}
    Js = zeros(T, N)
    gs = zeros(T, N)

    θs = range(0, 2π, N)
    @inbounds @threads for i in eachindex(θs)
        θ = θs[i]
        r, gp = find_offset_for_radius(m, u, d, rₑ, θ; zero_atol = zero_atol, offset_max = offset_max, max_time = max_time, tracer_kwargs...)
        α = r * cos(θ)
        β = r * sin(θ)
        g = redshift_pf(m, gp, max_time)
        gs[i] = g
        Js[i] = jacobian_∂αβ_∂gr(m, u, d, α, β, max_time; diff_order = diff_order, redshift_pf = redshift_pf, tracer_kwargs...)
    end

    gmin, gmax = extrema(gs)
    @. Js = (gmax - gmin) * Js
    @inbounds @threads for i in eachindex(Js)
        g✶ = (gs[i] - gmin) / (gmax - gmin)
        # Js is now storing f
        Js[i] = _calculate_transfer_function(rₑ, gs[i], g✶, Js[i])
        # gs is now storing g✶
        gs[i] = g✶
    end

    CunninghamTransferFunction(
       gs, Js, gmin, gmax, rₑ 
    )
end