struct DiscCorona{T} <: AbstractCoronaModel{T}
    "Radius of the disc"
    r::T
    "Height of the base of the cylinder above the disc"
    h::T
end

function sample_position_velocity(m::AbstractMetric, model::DiscCorona{T}) where {T}
    x = rand() * model.r
    r = sqrt(x^2 + model.h^2)
    # (x, y) flipped because coordinates are off the Z axis
    θ = atan(x, model.h)
    x = SVector{4,T}(0, r, θ, 0)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    v = inv(√(-gcomp[1])) * SVector{4,T}(1, 0, 0, 0)
    x, v
end

struct RingCorona{T} <: AbstractCoronaModel{T}
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

"""
    LongitudalArms{T}

For internal use. Used to pass around information about a particular slice of
geodesics through an offset point-like corona. The `β` field stores the rotation
angle in the local sky.
"""
struct LongitudalArms{T}
    β::T
    left_r::Vector{T}
    left_t::Vector{T}
    left_ε::Vector{T}
    right_r::Vector{T}
    right_t::Vector{T}
    right_ε::Vector{T}
end

"""
    struct RingCoronaProfile{T} <: AbstractDiscProfile

A specialised disc profile that can stores the various
[`TimeDependentRadialDiscProfile`](@ref) for the ring-like extended corona.
"""
struct RingCoronaProfile{T} <: AbstractDiscProfile
    left_arm::TimeDependentRadialDiscProfile{T}
    right_arm::TimeDependentRadialDiscProfile{T}
end

function emissivity_at(prof::RingCoronaProfile{T}, ρ) where {T}
    emissivity_at(prof.left_arm, ρ) + emissivity_at(prof.right_arm, ρ)
end

function emissivity_interp(prof::RingCoronaProfile{T}, ρ) where {T}
    left_arm = emissivity_interp(prof.left_arm, ρ)
    right_arm = emissivity_interp(prof.right_arm, ρ)
    function _add_arms(t)
        left = if t >= left_arm.t[1] && t <= left_arm.t[end]
            left_arm(t)
        else
            zero(T)
        end
        right = if t >= right_arm.t[1] && t <= right_arm.t[end]
            right_arm(t)
        else
            zero(T)
        end

        left + right
    end
end

function emissivity_interp_limits(prof::RingCoronaProfile, ρ)
    l_min, l_max = emissivity_interp_limits(prof.left_arm, ρ)
    r_min, r_max = emissivity_interp_limits(prof.right_arm, ρ)
    min(l_min, r_min), max(l_max, r_max)
end

function _make_time_dependent_radial_emissivity(arms::Vector{LongitudalArms{T}}) where {T}
    N = length(arms)
    left_radii = Vector{Vector{T}}(undef, N)
    left_times = Vector{Vector{T}}(undef, N)
    left_emissivities = Vector{Vector{T}}(undef, N)

    right_radii = Vector{Vector{T}}(undef, N)
    right_times = Vector{Vector{T}}(undef, N)
    right_emissivities = Vector{Vector{T}}(undef, N)

    for (i, arm) in enumerate(arms)
        left_radii[i] = arm.left_r
        left_times[i] = arm.left_t
        left_emissivities[i] = arm.left_ε

        right_radii[i] = arm.right_r
        right_times[i] = arm.right_t
        right_emissivities[i] = arm.right_ε
    end

    weights = zeros(T, N)
    left =
        TimeDependentRadialDiscProfile(weights, left_radii, left_times, left_emissivities)

    right = TimeDependentRadialDiscProfile(
        weights,
        right_radii,
        right_times,
        right_emissivities,
    )
    RingCoronaProfile(left, right)
end

function emissivity_profile(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::RingCorona;
    kwargs...,
)
    arms = _calculate_ring_arms(setup, m, d, model; kwargs...)
    _make_time_dependent_radial_emissivity(arms)
end

function _calculate_ring_arms(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric{T},
    d::AbstractAccretionGeometry,
    model::RingCorona;
    num_β_slices::Int = 13,
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
        LongitudalArms(β, left.r, left.t, left.ε, right.r, right.t, right.ε)
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
    tracegeodesics(
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
end

function _process_ring_traces(setup::EmissivityProfileSetup, m, d, v, gps, rs, δs)
    # no need to filter intersected, as we have already done that before calling this process function    J = sortperm(rs)
    J = sortperm(rs)
    points = gps[J]
    δs_sorted = δs[J]
    r, ε = _point_source_emissivity(m, d, setup.spectrum, v, rs[J], δs_sorted, points)
    t = [i.x[1] for i in @views(points[1:end-1])]
    (; t, r, ε = abs.(ε))
end

function integrate_lagtransfer(
    profile::RingCoronaProfile,
    itb::InterpolatingTransferBranches{T},
    radii,
    g_grid,
    t_grid;
    Nr = 1000,
    t0 = 0,
    kwargs...,
) where {T}
    # pre-allocate output
    flux = zeros(T, (length(g_grid), length(t_grid)))
    fine_rₑ_grid = Grids._geometric_grid(extrema(radii)..., Nr) |> collect
    setup =
        _integration_setup(T, _lineprofile_integrand, nothing; time = nothing, kwargs...)
    _integrate_transfer_problem!(
        flux,
        setup,
        itb,
        profile,
        fine_rₑ_grid,
        g_grid,
        t_grid;
        t0 = t0,
    )
    _normalize(flux, g_grid)
end

# TODO: refactor the time-integration to make things like below possible without copy pasting the function
function _integrate_transfer_problem!(
    output::AbstractMatrix,
    setup::IntegrationSetup,
    itb::InterpolatingTransferBranches{T},
    profile::RingCoronaProfile,
    radii_grid,
    g_grid,
    t_grid;
    t0 = 0,
) where {T}
    g_grid_view = @views g_grid[1:end-1]
    # build fine radial grid for trapezoidal integration
    @inbounds for (i, rₑ) in enumerate(radii_grid)
        closures = _IntegrationClosures(setup, itb(rₑ))
        Δrₑ = _trapezoidal_weight(radii_grid, rₑ, i)
        θ = Δrₑ * rₑ * π / (closures.branch.gmax - closures.branch.gmin)

        # interpolate the emissivity as function of time
        t_ε_interp = emissivity_interp(profile, rₑ)
        t_ε_limts = emissivity_interp_limits(profile, rₑ)
        δt = (t_ε_limts[2] - t_ε_limts[1]) / 100

        @inbounds for j in eachindex(g_grid_view)
            g_grid_low = clamp(g_grid[j], closures.branch.gmin, closures.branch.gmax)
            g_grid_hi = clamp(g_grid[j+1], closures.branch.gmin, closures.branch.gmax)
            # skip if bin not relevant
            if g_grid_low == g_grid_hi
                continue
            end
            for (glo, ghi) in
                _g_fine_grid_iterate(setup.g_grid_upscale, g_grid_low, g_grid_hi)
                k1 = integrate_bin(closures, _lower_branch(closures), glo, ghi)
                k2 = integrate_bin(closures, _upper_branch(closures), glo, ghi)
                # find which bin to dump in
                t_lower_branch, t_upper_branch = _time_bins(closures, glo, ghi)

                imax = lastindex(t_grid)
                # loop over all times and find the offsets to dump flux into
                for time in range(t_ε_limts..., 100)
                    tlower = t_lower_branch + time - t0
                    i1 = @views searchsortedfirst(t_grid, tlower)

                    tupper = t_upper_branch + time - t0
                    i2 = @views searchsortedfirst(t_grid, tupper)

                    em = t_ε_interp(time)
                    if i1 <= imax
                        output[j, i1] += k1 * θ * em * δt
                    end
                    if i2 <= imax
                        output[j, i2] += k2 * θ * em * δt
                    end
                end
            end
        end
    end
end

export RingCorona, DiscCorona
