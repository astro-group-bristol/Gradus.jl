module __MorrisThorneAD
using ..StaticArrays
@fastmath function metric_components(b, lθ)
    (l, θ) = lθ

    tt = -1
    rr = 1
    θθ = b^2 + l^2
    ϕϕ = (b^2 + l^2) * sin(θ)

    # no time-azimuth couple
    @SVector [tt, rr, θθ, ϕϕ, 0.0]
end

end # module

# new structure for our spacetime
"""
    struct MorrisThorneAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}

Morris-Thorne wormhole metric.

$(FIELDS)
"""
@with_kw struct MorrisThorneAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}
    @deftype T
    "Throat size."
    b = 1.0
end

# implementation
metric_components(m::MorrisThorneAD{T}, rθ) where {T} =
    __MorrisThorneAD.metric_components(m.b, rθ)
GradusBase.inner_radius(m::MorrisThorneAD{T}) where {T} = 0.0
