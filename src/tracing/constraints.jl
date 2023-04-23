# vectors of vectors 
function constrain_all(
    m::AbstractMetric,
    us::AbstractArray{<:StaticVector},
    vs::AbstractArray{<:StaticVector},
    μ,
)
    @inbounds for i in eachindex(vs)
        vs[i] = constrain_all(m, us[i], vs[i], μ)
    end
    vs
end

@inline function constrain_all(
    m::AbstractMetric,
    u::StaticVector{S},
    v::StaticVector{S,T},
    μ,
) where {S,T<:Number}
    # mut = MVector{S,T}(v)
    # mut[1] = constrain(m, u, v, μ = μ)
    # SVector{S,T}(mut)
    # this results in no performance gain over the above, but i think it
    # reads nicer
    @inbounds SVector{S,T}(constrain(m, u, v, μ = μ), v[2], v[3], v[4])
end

function wrap_constraint(m::AbstractMetric, position, velfunc::Function, μ)
    (i) -> constrain_all(m, position, velfunc(i), μ)
end
