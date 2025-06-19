# vectors of vectors 
function constrain_all(
    m::AbstractMetric,
    xs::AbstractArray{<:SVector},
    vs::AbstractArray{<:SVector},
    μ,
)
    @inbounds for i in eachindex(vs)
        vs[i] = constrain_all(m, xs[i], vs[i], μ)
    end
    vs
end

constrain_all(m::AbstractMetric, x::SVector{4}, v::SVector{4,T}, μ) where {T<:Number} =
    SVector{4,T}(constrain(m, x, v, μ = μ), v[2], v[3], v[4])

function wrap_constraint(m::AbstractMetric, position, velfunc::Function, μ)
    (i) -> constrain_all(m, position, velfunc(i), μ)
end

"""
    constrain_normalize(m::AbstractMetric, x, v::SVector{4,T}; μ = zero(T))

Normalize and then constrain the vector. This is useful for when a certain proportionality
between e.g. `v[1]` and `v[2]` is desired.
"""
function constrain_normalize(m::AbstractMetric, x, v::SVector{4,T}; μ = zero(T)) where {T}
    v_norm = v ./ √(abs(propernorm(m, x, v)))
    constrain_all(m, x, v_norm, μ)
end
