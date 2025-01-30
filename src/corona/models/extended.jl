module SourceVelocities

using ..Gradus: AbstractMetric, CircularOrbits, SVector

"""
    co_rotating(m::AbstractMetric, x::SVector{4})

Calculates the a source velocity assuming the cylinder described by the point
`x` co-rotates with the accretion disc below. This assumes Keplerian orbits, and
uses [`CircularOrbits`](@ref) to perform the calculation.
"""
function co_rotating(m::AbstractMetric, x::SVector{4})
    CircularOrbits.fourvelocity(m, SVector(x[2], x[3]))
end

"""
    stationary(m::AbstractMetric, x::SVector{4})

Assumes the source is stationay, and calculates a velocity vector with
components
```math
v^\\mu = (v^t, 0, 0, 0),
```
where the time-component is determined by the metric to satisfy
```math
g_{\\mu\\nu} v^\\mu v^\\nu = -1.
```
"""
function stationary(m::AbstractMetric, x)
    gcomp = metric_components(m, SVector(x[2], x[3]))
    inv(√(-gcomp[1])) * SVector{4,T}(1, 0, 0, 0)
end

end # module

"""
    RingCorona

A ring-like corona, representing an infinitely thin thing at some radius and
height above the accretion disc.
"""
struct RingCorona{T,VelFunc} <: AbstractCoronaModel{T}
    "Source velocity function. May be any one of [`SourceVelocities`](@ref) or a
    custom implementation."
    vf::VelFunc
    "Radius of the ring"
    r::T
    "Height of the base of the cylinder above the disc"
    h::T
end

RingCorona(r, h) = RingCorona(SourceVelocities.co_rotating, r, h)
RingCorona(; r = 5.0, h = 5.0, vf = SourceVelocities.co_rotating) = RingCorona(vf, r, h)

function sample_position_velocity(m::AbstractMetric, model::RingCorona{T}) where {T}
    r = sqrt(model.r^2 + model.h^2)
    # (x, y) flipped because coordinates are off the Z axis
    θ = atan(model.r, model.h)
    x = SVector{4,T}(0, r, θ, 0)
    x, model.vf(m, x)
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
    x, v = Gradus.sample_position_velocity(m, model)

    # trace the left half of the local sky
    left_δs = deg2rad.(range(setup.δmin, setup.δmax, setup.n_samples ÷ (2 * num_β_slices)))
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
            gps[1:(min_ρ_index-1)],
            ρs[1:(min_ρ_index-1)],
            δs_filtered[1:(min_ρ_index-1)],
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
    # no need to filter intersected, as we have already done that before calling this process function
    J = sortperm(rs)
    points = gps[J]
    δs_sorted = δs[J]
    r, ε = _point_source_emissivity(m, d, setup.spectrum, v, rs[J], δs_sorted, points)
    t = [i.x[1] for i in @views(points[1:(end-1)])]
    (; t, r, ε = abs.(ε))
end

# TODO: refactor the time-integration to make things like below possible without copy pasting the function
function _integrate_transfer_problem!(
    output::AbstractMatrix,
    setup::IntegrationSetup{T,Profile},
    transfer_function_radial_interpolation,
    r_limits,
    g_grid,
    t_grid;
    g_scale = 1,
) where {T,Profile<:Union{<:RingCoronaProfile,DiscCoronaProfile}}
    g_grid_view = @views g_grid[1:(end-1)]

    r_itterator = collect(Grids._geometric_grid(r_limits..., setup.n_radii))
    r2 = first(iterate(r_itterator, 2))
    # prime the first r_prev so that the bin width is r2 - r1
    r_prev = r_limits[1] - (r2 - r_limits[1])

    N_t_steps = 100

    for rₑ in r_itterator
        branch = transfer_function_radial_interpolation(rₑ)
        S_lower = _lower_branch(setup, branch)
        S_upper = _upper_branch(setup, branch)

        Δrₑ = rₑ - r_prev
        # integration weight for this annulus
        θ = Δrₑ * rₑ * π / (branch.gmax - branch.gmin)

        # interpolate the emissivity as function of time
        t_ε_interp = emissivity_interp(setup.profile, rₑ)
        t_ε_limts = emissivity_interp_limits(setup.profile, rₑ)
        δt = (t_ε_limts[2] - t_ε_limts[1]) / N_t_steps

        @inbounds for j in eachindex(g_grid_view)
            glo = clamp(g_grid[j] / g_scale, branch.gmin, branch.gmax)
            ghi = clamp(g_grid[j+1] / g_scale, branch.gmin, branch.gmax)
            # skip if bin not relevant
            if glo == ghi
                continue
            end

            Δg = (ghi - glo) / setup.g_grid_upscale

            for i = 1:setup.g_grid_upscale
                g_fine_lo = glo + (i - 1) * Δg
                g_fine_hi = g_fine_lo + Δg

                k1 = integrate_bin(
                    setup,
                    S_lower,
                    g_fine_lo,
                    g_fine_hi,
                    branch.gmin,
                    branch.gmax,
                )
                k2 = integrate_bin(
                    setup,
                    S_upper,
                    g_fine_lo,
                    g_fine_hi,
                    branch.gmin,
                    branch.gmax,
                )

                # find which bin to dump in
                (tl1, tl2), (tu1, tu2) = _time_bins(setup, branch, g_fine_lo, g_fine_hi)
                t_lower_branch = (tl1 + tl2) / 2
                t_upper_branch = (tu1 + tu2) / 2

                imax = lastindex(t_grid)
                # loop over all times and find the offsets to dump flux into
                for time in range(t_ε_limts..., N_t_steps)
                    tlower = t_lower_branch + time - setup.t0
                    i1 = @views searchsortedfirst(t_grid, tlower)

                    tupper = t_upper_branch + time - setup.t0
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

        r_prev = rₑ
    end
end

"""
    DiscCorona

A disk-like corona with no height but some radial extent.

Depending on the algorithm chosen, emissivity profiles for this corona are
either calculated via Monte-Carlo sampling, or by treating the extended source
as many concentric [`RingCorona`](@ref).
"""
struct DiscCorona{T,VelFunc} <: AbstractCoronaModel{T}
    "Source velocity function. May be any one of [`SourceVelocities`](@ref) or a
    custom implementation."
    vf::VelFunc
    "Radius of the disc"
    r::T
    "Height of the base of the cylinder above the disc"
    h::T
end

DiscCorona(; vf = SourceVelocities.co_rotating, r = 5.0, h = 5.0) =
    DiscCorona{typeof(vf)}(vf, r, h)

function sample_position_velocity(m::AbstractMetric, model::DiscCorona{T}) where {T}
    x = rand() * model.r
    r = sqrt(x^2 + model.h^2)
    # (x, y) flipped because coordinates are off the Z axis
    θ = atan(x, model.h)
    x = SVector{4,T}(0, r, θ, 0)
    x, model.vf(m, x)
end

function emissivity_profile(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::DiscCorona;
    n_rings = 10,
    kwargs...,
)
    radii = range(1e-2, model.r, n_rings)
    ring_profiles = map(radii) do r
        emissivity_profile(setup, m, d, RingCorona(model.vf, r, model.h); kwargs...)
    end
    DiscCoronaProfile(collect(radii), ring_profiles, _ -> 0)
end


export RingCorona, DiscCorona
