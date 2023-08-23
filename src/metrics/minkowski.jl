
struct SphericalMetric{T} <: AbstractStaticAxisSymmetric{T} end

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

inner_radius(::SphericalMetric) = 4eps()
