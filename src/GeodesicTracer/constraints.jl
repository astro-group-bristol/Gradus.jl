# vectors of vectors 
function constrain_all(m::AbstractMetricParams{T}, us, vs, μ) where {T}
    @inbounds for i in eachindex(vs)
        vs[i] = constrain_all(m, us[i], vs[i], μ)
    end
    vs
end

function constrain_all(
    m::AbstractMetricParams{T},
    u::AbstractVector{T},
    v::AbstractVector{T},
    μ,
) where {T<:Number}
    v[1] = constrain(m, u, v, μ = μ)
    v
end

function constrain_all(
    m::AbstractMetricParams{T},
    u::StaticVector{S,T},
    v::StaticVector{S,T},
    μ,
) where {S,T<:Number}
    # mut = MVector{S,T}(v)
    # mut[1] = constrain(m, u, v, μ = μ)
    # SVector{S,T}(mut)
    #
    # this results in no performance gain over the above, but i think it
    # reads nicer
    @inbounds SVector{S,T}(constrain(m, u, v, μ = μ), v[2], v[3], v[4])
end

function wrap_constraint(m::AbstractMetricParams{T}, position, velfunc::Function, μ) where {T}
    (i) -> constrain_all(m, position, velfunc(i), μ)
end