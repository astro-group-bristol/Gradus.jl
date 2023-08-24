
struct SphericalMetric{T} <: AbstractStaticSphericallySymmetric{T} end
SphericalMetric() = SphericalMetric{Float64}()

function Gradus.metric_components(::SphericalMetric, rθ)
    T = eltype(rθ)
    r, θ = rθ
    tt = -one(T)
    rr = one(T)
    θθ = r^2
    ϕϕ = r^2 * sin(θ)^2

    SVector(tt, rr, θθ, ϕϕ, zero(T))
end
isco(::SphericalMetric{T}) where {T} = zero(T)
inner_radius(::SphericalMetric{T}) where {T} = 4eps(T)

struct CartesianMetric{T} <: AbstractMetric{T,FlatCartesian{(:x, :y, :z)}} end
CartesianMetric() = CartesianMetric{Float64}()

function Gradus.metric_components(::CartesianMetric, x)
    T = eltype(x)
    SVector{4,T}(-1, 1, 1, 1)
end
isco(::CartesianMetric{T}) where {T} = zero(T)
inner_radius(::CartesianMetric{T}) where {T} = 4eps(T)

# trait based dispatch to get the minkowski metric depending on the coordinates of the spacetime
minkowski_with_coords(::M) where {M<:AbstractMetric} = minkowski_metric(M)
minkowski_with_coords(::Type{<:AbstractMetric{T,<:BoyerLindquist}}) where {T} =
    SphericalMetric{T}()
minkowski_with_coords(::Type{<:AbstractMetric{T,<:FlatCartesian}}) where {T} =
    CartesianMetric{T}()

minkowski_matrix(m::AbstractMetric, x) = metric(minkowski_metric(m), x)
minkowski_matrix() = @SMatrix [
    -1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 0 1
]

export minkowski_matrix
