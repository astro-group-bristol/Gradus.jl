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

@inline function constrain_all(
    m::AbstractMetric,
    x::SVector{4},
    v::SVector{4,T},
    μ,
) where {T<:Number}
    @inbounds SVector{4,T}(_constrain(m, x, v, μ = μ), v[2], v[3], v[4])
end

function wrap_constraint(m::AbstractMetric, position, velfunc::Function, μ)
    (i) -> constrain_all(m, position, velfunc(i), μ)
end
