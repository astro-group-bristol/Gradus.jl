function split_options(::AbstractMetricParameters, opts)
    opts, (;)
end

"""
    map_impact_parameters(m::AbstractMetricParameters{T}, u, α, β)

Map impact parameters `α` and `β` to a velocity vector at some position `u` in the given metric `m`.
"""
function map_impact_parameters(
    m::AbstractMetricParameters{T},
    x::SVector{S,T},
    α::N,
    β::N,
) where {S,T,N<:Number}
    SVector{S,T}(T(0.0), impact_parameters_to_vel(m, x, α, β)...)
end

function map_impact_parameters(
    m::AbstractMetricParameters{T},
    x,
    α::AbstractVector{P},
    β::AbstractVector{P},
) where {T,P}
    [map_impact_parameters(m, x, a, b) for (a, b) in zip(α, β)]
end

function map_impact_parameters(m::AbstractMetricParameters, x, α::AbstractVector, β::Number)
    [map_impact_parameters(m, x, a, β) for a in α]
end
function map_impact_parameters(m::AbstractMetricParameters, x, α::Number, β::AbstractVector)
    [map_impact_parameters(m, x, α, b) for b in β]
end

@inline function impact_parameters_to_vel(m::AbstractMetricParameters{T}, u, α, β) where {T}
    mcomp = metric_components(m, @view(u[2:3]))
    T(-1.0), -β / mcomp[3], -α / √(mcomp[3] * mcomp[4])
end

function faraday_tensor(m::AbstractMetricParameters, x)
    ST = SVector{4,eltype(x)}
    dA = ForwardDiff.jacobian(t -> electromagnetic_potential(m, t), SVector(x[2], x[3]))
    ∂A = hcat(zeros(ST), dA, zeros(ST))
    g = inv(metric(m, x))
    # raise first index: F^μ_κ
    # @tullio F[μ, κ] := g[μ, σ] * (∂A[σ, κ] - ∂A[κ, σ])

    # faster
    g * (∂A - ∂A')
end

export map_impact_parameters
