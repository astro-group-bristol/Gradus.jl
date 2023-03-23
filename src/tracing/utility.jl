function split_options(::AbstractMetric, opts)
    opts, (;)
end

"""
    map_impact_parameters(m::AbstractMetric{T}, u, α, β)

Map impact parameters `α` and `β` to a velocity vector at some position `u` in the given metric `m`.
"""
function map_impact_parameters(
    m::AbstractMetric{T},
    x::SVector{S,T},
    α::N,
    β::N,
) where {S,T,N<:Number}
    SVector{S,T}(T(0.0), impact_parameters_to_three_velocity(m, x, α, β)...)
end

function map_impact_parameters(
    m::AbstractMetric{T},
    x,
    α::AbstractVector{P},
    β::AbstractVector{P},
) where {T,P}
    [map_impact_parameters(m, x, a, b) for (a, b) in zip(α, β)]
end

function map_impact_parameters(m::AbstractMetric, x, α::AbstractVector, β::Number)
    [map_impact_parameters(m, x, a, β) for a in α]
end
function map_impact_parameters(m::AbstractMetric, x, α::Number, β::AbstractVector)
    [map_impact_parameters(m, x, α, b) for b in β]
end


"""
    impact_parameters_to_three_velocity(m::AbstractMetric, x, α, β)

Return the three velocity components ``(v^r, v^\\theta, v^\\phi)`` corresponding to impact parameters ``\\alpha`` and
``\\beta``, at position `x` for metric `m`.

For static, axis-symmetric metrics, the map is defined
```math
(\\alpha, \\beta) \\mapsto \\vec{v} =
\\left(
\\begin{matrix}
    -1 \\
    \\frac{-\\beta}{g_{\\theta\\theta}} \\
    \\frac{-\\alpha}{\\sqrt{g_{\\theta\\theta} g_{\\phi\\phi}}} \\
\\end{matrix}
\\right).
```

The impact parameters are interpreted as follows:

- if the geodesic were a straight line path, the impact parameter in a given dimension is the distance to the origin from 
the closest point along the geodesic in ``r_\\text{g}``.

For example, for the Schwarzschild metric with ``M = 1``, the impact parameters ``(\\alpha = 2, \\beta = 0)`` would travel 
tangential to the event horizon (``r_\\text{s} = 2M``) if space were flat. 
"""
@inline function impact_parameters_to_three_velocity(
    m::AbstractStaticAxisSymmetric{T},
    x,
    α,
    β,
) where {T}
    mcomp = metric_components(m, @view(x[2:3]))
    T(-1.0), -β / mcomp[3], -α / √(mcomp[3] * mcomp[4])
end

function faraday_tensor(m::AbstractMetric, x)
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
