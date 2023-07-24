@with_kw struct LampPostModel{T} <: AbstractCoronaModel{T}
    @deftype T
    h = 5.0
    θ = 0.01
    ϕ = 0.0
end

function sample_position_velocity(m::AbstractMetric, model::LampPostModel{T}) where {T}
    x = SVector{4,T}(0, model.h, model.θ, model.ϕ)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-gcomp[1])) * SVector{4,T}(1, 0, 0, 0)
    x, v
end

# can exploit point source symmetry for certain disc models
emissivity_profile(
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    model::LampPostModel;
    kwargs...,
) = _point_source_symmetric_emissivity_profile(m, d, model; kwargs...)
