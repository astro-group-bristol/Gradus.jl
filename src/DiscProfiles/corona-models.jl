abstract type AbstractCoronaModel{T} end

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
    map(1:N) do i
        θ, ϕ = sample_angles(sampler, i, N)
        r, t, p = vector_to_local_sky(m, us[i], θ, ϕ)
        @SVector [T(0.0), r, t, p]
    end
end
