# todo: this should take into account position
function source_velocity(::AbstractMetricParams, model::AbstractCoronaModel)
    error("Not implemented for $(typeof(model)).")
end

function sample_position(::AbstractMetricParams, model::AbstractCoronaModel, N)
    error("Not implemented for $(typeof(model)).")
end

"""
    sample_velocity(
        m::AbstractMetricParams, 
        model::AbstractCoronaModel, 
        sampler::AbstractDirectionSampler, 
        us, 
        N
    )

Sample `N` initial (un-normalised) velocities in the local coordinates of the `model`, according to
the [`AbstractDirectionSampler`](@ref) algorithm and domain, given initial 4-vector positions `us`.

This function is metric generic and is implemented in the following way:
- Sample angles ``\\alpha`` and ``\\beta`` in the local sky.
- Map these angles to a spherical polar vector (interally first to Cartesian and then using a Jacobian
transformation to spherical polar).
- Calculate the local tetrad given the source velocity at `u = us[i]`.
- Use this tetrad to map the local vector to the global coordinates.
"""
function sample_velocity(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    sampler::AbstractDirectionSampler,
    us,
    N,
) where {T}
    v_source = source_velocity(m, model)
    @inbounds map(1:N) do index
        # todo: sampler should have proper iterator interface
        i = geti(sampler, index, N)
        θ, ϕ = sample_angles(sampler, i, N)

        u = us[index]
        sky_angles_to_velocity(m, u, v_source, θ, ϕ)
    end
end


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

function source_velocity(m::AbstractMetricParams, model::LampPostModel)
    # stationary source
    rθ = @SVector[model.h, model.θ]
    gcomp = metric_components(m, rθ)
    inv(√(-gcomp[1])) * @SVector[1.0, 0.0, 0.0, 0.0]
end

# this function is a placeholder: needs to be able to work
# for any disc geometry and model type
function energy_ratio(m, gps, model::LampPostModel)
    energy_ratio(m, gps, SVector(0.0, model.h, model.θ, 0.0), source_velocity(m, model))
end
function energy_ratio(m, gps, u_src, v_src)
    g_src = metric(m, u_src)
    map(gps) do gp
        @tullio e_src := g_src[i, j] * gp.v1[i] * v_src[j]
        # at the disc
        g_disc = metric(m, gp.u2)
        v_disc = CircularOrbits.fourvelocity(m, SVector(gp.u2[2], gp.u2[3]))
        @tullio e_disc := g_disc[i, j] * gp.v2[i] * v_disc[j]
        # ratio g = E_source / E_disc
        e_src / e_disc
    end
end

export LampPostModel
