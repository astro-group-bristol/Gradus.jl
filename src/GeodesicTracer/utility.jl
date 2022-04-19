"""
    map_impact_parameters(m::AbstractMetricParams{T}, u, α, β)

Map impact parameters `α` and `β` to a velocity vector at some position `u` in the given metric `m`.
"""
function map_impact_parameters(
    m::AbstractMetricParams{T},
    u::AbstractVector{T},
    α::N,
    β::N
) where {T,N<:Number}
    [0.0, alpha_beta_to_vel(m, u, α, β)...]
end

function map_impact_parameters(
    m::AbstractMetricParams{T},
    u::SVector{S,T},
    α::N,
    β::N
) where {S,T,N<:Number}
    SVector{S,T}(0.0, alpha_beta_to_vel(m, u, α, β)...)
end

function map_impact_parameters(
    m::AbstractMetricParams{T},
    u,
    α::AbstractVector{P},
    β::AbstractVector{P}
) where {T,P}
    curried(_α, _β) = map_impact_parameters(m, u, _α, _β)
    curried.(α, β)
end

function alpha_beta_to_vel(m::AbstractMetricParams{T}, u, α, β) where {T}
    reg = u[2]^2
    -1.0, -β / reg, -α / reg
end

export map_impact_parameters
