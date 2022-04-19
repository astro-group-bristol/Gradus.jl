# inplace
function constrain_all(
    m::AbstractMetricParams{T},
    us::AbstractArray{T,2},
    vs::AbstractArray{T,2},
    μ
) where {T<:Number}
    us_cols = eachcol(us)
    vs_cols = eachcol(vs)
    curried(u, v) = constrain(m, u, v, μ = μ)
    @. vs[1, :] = curried(us_cols, vs_cols)
    vs
end

# SMatrix -- do we actually want to support this?
function constrain_all(
    m::AbstractMetricParams{T},
    us::SMatrix{4,M,T,L},
    vs::SMatrix{4,M,T,L},
    μ
) where {M,T<:Number,L}
    alloc_vs = constrain_all(m, us, similar(vs), μ)
    SMatrix{4,M}(alloc_vs)
end

# vectors of vectors 
function constrain_all(m::AbstractMetricParams{T}, us, vs, μ) where {T<:Number}
    curried(u, v) = constrain_all(m, u, v, μ)
    curried.(us, vs)
end

function constrain_all(
    m::AbstractMetricParams{T},
    u::AbstractVector{T},
    v::AbstractVector{T},
    μ
) where {T<:Number}
    v[1] = constrain(m, u, v, μ = μ)
    v
end

function constrain_all(
    m::AbstractMetricParams{T},
    u::StaticVector{S,T},
    v::StaticVector{S,T},
    μ
) where {S,T<:Number}
    mut = MVector{S,T}(v)
    mut[1] = constrain(m, u, v, μ = μ)
    SVector{S,T}(mut)
end
