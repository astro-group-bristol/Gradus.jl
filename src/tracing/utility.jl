"""
    map_impact_parameters(m::AbstractMetricParams{T}, u, α, β)

Map impact parameters `α` and `β` to a velocity vector at some position `u` in the given metric `m`.
"""
function map_impact_parameters(
    m::AbstractMetricParams{T},
    u::SVector{S,T},
    α::N,
    β::N,
) where {S,T,N<:Number}
    SVector{S,T}(T(0.0), impact_parameters_to_vel(m, u, α, β)...)
end

function map_impact_parameters(
    m::AbstractMetricParams{T},
    u,
    α::AbstractVector{P},
    β::AbstractVector{P},
) where {T,P}
    [map_impact_parameters(m, u, a, b) for (a, b) in zip(α, β)]
end

@inline function impact_parameters_to_vel(m::AbstractMetricParams{T}, u, α, β) where {T}
    mcomp = metric_components(m, @view(u[2:3]))
    T(-1.0), -β / mcomp[3], -α / √(mcomp[3] * mcomp[4])
end

export map_impact_parameters
