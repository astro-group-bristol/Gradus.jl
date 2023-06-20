function sample_position_direction_velocity(m::AbstractMetric, model::AbstractCoronaModel{T}, sampler::AbstractDirectionSampler, N::Int) where {T}
    xs = Vector{SVector{4,T}}(undef, N)
    vs = Vector{SVector{4,T}}(undef, N)
    vs_source = Vector{SVector{4,T}}(undef, N)

    rmin = inner_radius(m)
    
    i = 1
    while (i <= N)
        x, v = sample_position_velocity(m, model, sampler, i, N)
        if x[2] < rmin * 1.1
            # skip those within the event horizon
            continue
        end

        xs[i] = x
        vs_source[i] = v
        vs[i] = sample_local_velocity(m, sampler, x, v, i, N)
        i += 1
    end

    xs, vs, vs_source
end

function sample_position_velocity(::AbstractMetric, model::AbstractCoronaModel, ::AbstractDirectionSampler, i, N,)
    error("This functions needs to be implemented for $(typeof(model)). See the documentation for this function for instructions.")
end

"""
function sample_local_velocity(
    m::AbstractMetric,
    sampler::AbstractDirectionSampler,
    x,
    v,
    index,
    N,
)

Sample a single initial (un-normalised) velocity vector in the local coordinates of at `x`, according to
the [`AbstractDirectionSampler`](@ref) algorithm and domain, given an initial 4-vector velocity `v`.

This function is metric generic and is implemented in the following way:
- Sample angles ``\\alpha`` and ``\\beta`` in the local sky.
- Map these angles to a spherical polar vector (interally first to Cartesian and then using a Jacobian
transformation to spherical polar).
- Calculate the tetrad given position and velocity `x` and `v`.
- Use this tetrad to map the local vector to the global coordinates.
"""
function sample_local_velocity(
    m::AbstractMetric{T},
    sampler::AbstractDirectionSampler,
    x,
    v,
    index,
    N,
) where {T}
    i = geti(sampler, index, N)
    θ, ϕ = sample_angles(sampler, i, N)
    sky_angles_to_velocity(m, x, v, θ, ϕ)
end

@with_kw struct LampPostModel{T} <: AbstractCoronaModel{T}
    @deftype T
    h = 5.0
    θ = 0.01
    ϕ = 0.0
end

function sample_position_velocity(m::AbstractMetric, model::LampPostModel{T}, ::AbstractDirectionSampler, i, N) where {T}
    x = SVector{4, T}(0, model.h, model.θ, model.ϕ)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-gcomp[1])) * SVector{4, T}(1, 0, 0, 0)
    x, v
end

export LampPostModel
