module __MorrisThorneAD
using ..StaticArrays

@fastmath function metric_components(b, lθ)
    (l, θ) = lθ

    tt = -1
    rr = 1
    θθ = b^2 + l^2
    ϕϕ = (b^2 + l^2) * sin(θ)

    # no time-azimuth couple
    @SVector [tt, rr, θθ, ϕϕ, 0.0]
end

end # module

# new structure for our spacetime
"""
    struct MorrisThorneWormhole{T} <: AbstractStaticAxisSymmetric{T}

Morris-Thorne wormhole metric.

- `b = 1.0`: Throat size.
"""
@with_kw struct MorrisThorneWormhole{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Throat size."
    b = 1.0
end

# implementation
metric_components(m::MorrisThorneWormhole{T}, rθ) where {T} =
    __MorrisThorneAD.metric_components(m.b, rθ)
inner_radius(m::MorrisThorneWormhole{T}) where {T} = 0.0

export MorrisThorneWormhole
