function sample_position_direction_velocity(
    m::AbstractMetric,
    model::AbstractCoronaModel{T},
    sampler::AbstractDirectionSampler,
    N::Int,
) where {T}
    xs = Vector{SVector{4,T}}(undef, N)
    vs = Vector{SVector{4,T}}(undef, N)
    vs_source = Vector{SVector{4,T}}(undef, N)

    rmin = inner_radius(m)
    Threads.@threads for i = 1:N
        x, v = sample_position_velocity(m, model, sampler, i, N)
        while x[2] < rmin * 1.9
            x, v = sample_position_velocity(m, model, sampler, i, N)
        end

        # avoid coordinate singularities
        if x[3] < 1e-3
            x = SVector(x[1], x[2], 1e-3, x[4])
        end
        if x[3] > π - 1e-3
            x = SVector(x[1], x[2], π - 1e-3, x[4])
        end

        xs[i] = x
        vs_source[i] = v
        vs[i] = sample_local_velocity(m, sampler, x, v, i, N)
    end

    xs, vs, vs_source
end

function sample_position_velocity(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    ::AbstractDirectionSampler,
    i,
    N,
)
    sample_position_velocity(m, model)
end

function sample_position_velocity(::AbstractMetric, model::AbstractCoronaModel)
    error(
        "This functions needs to be implemented for $(typeof(model)). See the documentation for this function for instructions.",
    )
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

# bootstrap tracing function for convenience
function tracegeodesics(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    args...;
    n_samples = 1024,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
    kwargs...,
)
    xs, vs, _ = sample_position_direction_velocity(m, model, sampler, n_samples)
    tracegeodesics(m, xs, vs, args...; kwargs...)
end

struct CoronaGeodesics{T,M,G,C,P,V}
    trace::T
    metric::M
    geometry::G
    model::C
    geodesic_points::P
    source_velocity::V
end

function tracecorona(
    m::AbstractMetric,
    g::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    λ_max = 10_000,
    n_samples = 1024,
    sampler = EvenSampler(domain = BothHemispheres(), generator = RandomGenerator()),
    trace = TraceGeodesic(),
    callback = domain_upper_hemisphere(),
    kwargs...,
)
    xs, vs, source_vels = sample_position_direction_velocity(m, model, sampler, n_samples)
    gps = tracegeodesics(
        m,
        xs,
        vs,
        g,
        λ_max;
        trace = trace,
        save_on = false,
        ensemble = EnsembleEndpointThreads(),
        callback = callback,
        kwargs...,
    )
    mask = [i.status == StatusCodes.IntersectedWithGeometry for i in gps]
    CoronaGeodesics(trace, m, g, model, gps[mask], source_vels[mask])
end

include("models/lamp-post.jl")

export LampPostModel
