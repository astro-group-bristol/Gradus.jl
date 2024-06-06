struct RingCorona{T} <: Gradus.AbstractCoronaModel{T}
    "Radius of the ring"
    r::T
    "Height of the base of the cylinder above the disc"
    h::T
end

function sample_position_velocity(m::AbstractMetric, model::RingCorona{T}) where {T}
    r = sqrt(model.r^2 + model.h^2)
    # (x, y) flipped because coordinates are off the Z axis
    θ = atan(model.r, model.h)
    x = SVector{4,T}(0, r, θ, 0)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-gcomp[1])) * SVector{4,T}(1, 0, 0, 0)
    x, v
end

struct LongitudalArms{T,A}
    β::T
    left::A
    right::A
end

function emissivity_profile(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric{T},
    d::AbstractAccretionGeometry,
    model::RingCorona;
    num_β_slices::Int = 20,
    β_min = 0,
    β_max = π,
    callback = domain_upper_hemisphere(),
    kwargs...,
) where {T}
    # TODO: assume co-rotating with the disc portion below it
    # for now, assume stationary corona
    x, v = Gradus.sample_position_velocity(m, model)

    # trace the left half of the local sky
    left_δs = deg2rad.(range(setup.δmin, setup.δmax, setup.n_samples))
    right_δs = left_δs .+ π
    δs = vcat(right_δs, left_δs)

    map(range(β_min, β_max, num_β_slices)) do β
        all_gps = _ring_corona_emissivity_slice(
            setup,
            m,
            d,
            model,
            δs,
            x,
            v,
            β;
            callback = callback,
            kwargs...,
        )

        # regular sorting
        mask = [i.status == StatusCodes.IntersectedWithGeometry for i in all_gps]
        gps = all_gps[mask]
        δs_filtered = @views δs[mask]

        # split the sky into a left and right side, such that each side has strictly
        # sortable and monotonic r as a function of θ on the local sky
        ρs = map(i -> _equatorial_project(i.x), gps)
        _, min_ρ_index = findmin(ρs)

        left = @views _process_ring_traces(
            setup,
            m,
            d,
            v,
            gps[1:min_ρ_index-1],
            ρs[1:min_ρ_index-1],
            δs_filtered[1:min_ρ_index-1],
        )
        right = @views _process_ring_traces(
            setup,
            m,
            d,
            v,
            gps[min_ρ_index:end],
            ρs[min_ρ_index:end],
            δs_filtered[min_ρ_index:end],
        )
        LongitudalArms(β, left, right)
    end
end

# one slice of the orange
function _ring_corona_emissivity_slice(
    setup::EmissivityProfileSetup,
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::RingCorona,
    δs,
    x,
    v,
    β;
    kwargs...,
)
    θ₀ = atan(model.r, model.h)
    velfunc = rotated_polar_angle_to_velfunc(m, x, v, δs, β; θ₀ = θ₀)
    gps = tracegeodesics(
        m,
        x,
        velfunc,
        d,
        setup.λmax;
        save_on = false,
        ensemble = Gradus.EnsembleEndpointThreads(),
        trajectories = length(δs),
        kwargs...,
    )

    gps
end

function _process_ring_traces(setup::EmissivityProfileSetup, m, d, v, gps, rs, δs)
    # no need to filter intersected, as we have already done that before calling this process function    J = sortperm(rs)
    J = sortperm(rs)
    points = gps[J]
    δs_sorted = δs[J]
    r, ε = _point_source_emissivity(m, d, setup.spectrum, v, rs[J], δs_sorted, points)
    t = [i.x[1] for i in @views(points[1:end-1])]
    (; t, r, ε)
end


export RingCorona
