sample_position(m::AbstractMetricParams{T}, model::AbstractCoronaModel{T}, N) where {T} =
    error("Not implemented for $(typeof(model)).")

@with_kw struct LampPostModel{T} <: AbstractCoronaModel{T}
    @deftype T
    h = 5.0
    θ = 0.01
    ϕ = 0.0
end

function sample_position(::AbstractMetricParams{T}, model::LampPostModel{T}, N) where {T}
    u = @SVector [T(0.0), model.h, model.θ, model.ϕ]
    fill(u, N)
end

function sample_velocity(
    m::AbstractMetricParams{T},
    ::AbstractCoronaModel{T},
    sampler::AbstractDirectionSampler,
    us,
    N,
) where {T}
    map(1:N) do index
        # todo: sampler should have proper iterator interface
        i = geti(sampler, index, N)
        θ, ϕ = sample_angles(sampler, i, N)
        r, t, p = vector_to_local_sky(m, us[index], θ, ϕ)
        @SVector [T(0.0), r, t, p]
    end
end

function source_velocity(m::AbstractMetricParams, model::AbstractCoronaModel)
    error("Not implemented for $(typeof(model)).")
end

function source_velocity(m::AbstractMetricParams, model::LampPostModel)
    # stationary source
    rθ = @SVector[model.h, model.θ]
    gcomp = metric_components(m, rθ)
    inv(√(-gcomp[1])) * @SVector[1.0, 0.0, 0.0, 0.0]
end

export LampPostModel
