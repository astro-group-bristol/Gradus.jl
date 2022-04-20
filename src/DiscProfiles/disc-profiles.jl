# work in progress code
# not actually used yet

abstract type AbstractAccretionProfile end
abstract type AbstractInterpolatedProfile <: AbstractAccretionProfile end
abstract type AbstractDelaunayProfile <: AbstractAccretionProfile end

evaluate(aap::AbstractAccretionProfile, p) =
    error("Not implemented for `$(typeof(aap))` yet.")

struct WeightedPoint{T} <: Point2D
    _x::Float64
    _y::Float64
    _weight::T
end
GeometricalPredicates.getx(p::WeightedPoint{T}) where {T} = p._x
GeometricalPredicates.gety(p::WeightedPoint{T}) where {T} = p._y
getw(p::WeightedPoint{T}) where {T} = p._weight

struct DelaunayDiscProfile{T} <: AbstractDelaunayProfile
    tesselation::DelaunayTessellation2D{WeightedPoint{T}}
    radius::Float64
end
xfm_forward(ddp::DelaunayDiscProfile{T}, p) where {T} =
    ((p .+ ddp.radius) ./ (2 .* ddp.radius)) .* 0.99 .+ 1.005
xfm_backward(ddp::DelaunayDiscProfile{T}, p) where {T} =
    (((points .- 1.005) ./ 0.99) .* 2 .* ddp.radius) .- ddp.radius


function __renderprofile(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T}N,
    time_domain;
    sampler,
    kwargs...,
) where {T}
    # TODO
    error("Not implemented for $(typeof(d)).")
end


function __renderprofile(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::GeometricThinDisc{T}N,
    time_domain;
    sampler,
    kwargs...,
) where {T}
    # TODO
    error("Not implemented for $(typeof(d)).")
end
