module SourceVelocities

using ..Gradus:
    AbstractMetric,
    CircularOrbits,
    SVector,
    metric_components,
    propernorm,
    metric,
    constrain_all,
    isco

"""
    co_rotating(m::AbstractMetric, x::SVector{4})

Calculates the a source velocity assuming the cylinder described by the point
`x` co-rotates with the accretion disc below. This assumes Keplerian orbits, and
uses [`CircularOrbits`](@ref) to perform the calculation.
"""
function co_rotating(m::AbstractMetric, x::SVector{4})
    sinθ = sin(x[3])
    v = CircularOrbits.fourvelocity(m, max(isco(m), x[2] * sinθ)) .* sinθ
    v = v ./ sqrt(abs(propernorm(metric(m, x), v)))
    v = constrain_all(m, x, v, 1.0)
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
function stationary(m::AbstractMetric{T}, x) where {T}
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

function DEFAULT_β_ANGLES(; n_regular = 100, n_refined = 100)
    angles_coarse = range(0, π, n_regular)
    theta_angles = vcat(angles_coarse, range(0.62, 2.2, n_refined)) |> sort
    theta_angles
end

function emissivity_profile(
    setup::EmissivityProfileSetup{true},
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::RingCorona;
    βs = DEFAULT_β_ANGLES(),
    kwargs...,
)
    arms = corona_arms(setup, m, d, model, βs; kwargs...)
    _make_time_dependent_radial_emissivity(arms)
end

function _process_ring_traces(setup::EmissivityProfileSetup, m, d, v, gps, rs, δs)
    # no need to filter intersected, as we have already done that before calling this process function
    J = sortperm(δs)

    J = unique(i -> rs[i], J)
    # now need to find the cycle permutation to make the minimal radius the first element
    _, index = findmax(i -> rs[i], J)
    J = vcat(J[index+1:end], J[1:index-1])

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
