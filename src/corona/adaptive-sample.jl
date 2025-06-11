struct CoronaGridValues{T}
    "Time"
    t::T
    "Radius on the disc"
    r::T
    "Radial angle on the disc"
    ϕ::T
    "Redshift"
    g::T
    "|∂(r, ϕ) / ∂(theta, phi)| / sin(theta)"
    J::T
end

_to_vector(v::CoronaGridValues) = SVector{5}(v.t, v.r, v.ϕ, v.g, v.J)
_from_vector(v::AbstractVector) = CoronaGridValues(v[1], v[2], v[3], v[4], v[5])

make_null(::Type{T}) where {T<:CoronaGridValues} = T(NaN, NaN, NaN, NaN, NaN)

vector_average(
    weights::AbstractVector{<:Number},
    values::AbstractVector{<:CoronaGridValues},
) = _from_vector(sum(i -> i[1] * _to_vector(i[2]), zip(weights, values)))

const Tag_Make_Emissivity = Val{:_make_emissivity_tracer}

_make_emiss_Dual(x; d1 = zero(typeof(x)), d2 = zero(typeof(x))) =
    ForwardDiff.Dual{Tag_Make_Emissivity}(x, d1, d2)

function _emissivity_jacobian!(integ, parameters, th::T, ph::T) where {T}
    th_dual = _make_emiss_Dual(th; d1 = one(T))
    ph_dual = _make_emiss_Dual(ph; d2 = one(T))

    v = sky_angles_to_velocity(
        parameters.m,
        parameters.x_src,
        parameters.v_src,
        th_dual,
        ph_dual,
    )
    gp = _solve_reinit!(integ, vcat(parameters.xinit, v))

    if gp.status != StatusCodes.IntersectedWithGeometry
        return CoronaGridValues{T}(NaN, NaN, NaN, NaN, NaN)
    end

    v_disc = parameters.disc_velocity(gp.x)
    _redshift = energy_ratio(parameters.m, gp, parameters.v_src, v_disc)
    r = _equatorial_project(gp.x)

    y_dual = SVector{4,typeof(th_dual)}(r, gp.x[4], gp.x[1], _redshift)

    res = ForwardDiff.value.(Tag_Make_Emissivity, y_dual)
    jac = _extract_jacobian(
        Tag_Make_Emissivity,
        SVector{2}(y_dual[1], y_dual[2]),
        SVector{2}(th_dual, ph_dual),
    )

    CoronaGridValues{T}(res[3], res[1], res[2], res[4], abs(det(jac)) / sin(th))
end

function _make_emissivity_tracer(
    m::AbstractMetric{NumT},
    corona::AbstractCoronaModel,
    d::AbstractAccretionDisc;
    t_max = 1_000_000.0,
    solver_opts...,
) where {NumT}
    dual_zero = _make_emiss_Dual(zero(NumT))

    x_src, v_src = sample_position_velocity(m, corona)

    xinit = SVector{4,typeof(dual_zero)}(
        _make_emiss_Dual(x_src[1]),
        _make_emiss_Dual(x_src[2]),
        _make_emiss_Dual(x_src[3]),
        _make_emiss_Dual(x_src[4]),
    )

    # function for obtaining keplerian velocities
    disc_velocity = _keplerian_velocity_projector(m, d)

    _v_temp = sky_angles_to_velocity(m, x_src, v_src, _make_emiss_Dual(0.1), dual_zero)

    # init a reusable integrator
    integ = _init_integrator(
        m,
        xinit,
        _v_temp,
        d,
        _make_emiss_Dual(t_max);
        save_on = false,
        callback = domain_upper_hemisphere(),
        chart = chart_for_metric(m, 2 * t_max),
        integrator_verbose = false,
        solver_opts...,
    )

    integ, (; m, x_src, v_src, xinit, disc_velocity)
end

function check_refine(
    sky::AdaptiveSky{T,<:CoronaGridValues},
    i1::Int,
    i2::Int;
    percentage = 2,
) where {T}
    v1 = sky.values[i1]
    v2 = sky.values[i2]

    if isnan(v1.r) && isnan(v2.r)
        return false
    end

    g_too_coarse = !isapprox(v1.g, v2.g, atol = percentage / 100)
    J_too_coarse = !isapprox(v1.J, v2.J, atol = percentage / 100)

    g_too_coarse || J_too_coarse
end

function AdaptiveSky(
    m::AbstractMetric{T},
    corona::AbstractCoronaModel,
    d::AbstractAccretionDisc;
    kwargs...,
) where {T}
    integ, ps = _make_emissivity_tracer(m, corona, d; kwargs...)

    # one integrator copy for each thread
    integrators = typeof(integ)[deepcopy(integ) for _ = 1:(Threads.nthreads()-1)]
    push!(integrators, integ)

    # tracing function
    function _trace(th, ph)
        t_id = Threads.threadid()
        _emissivity_jacobian!(integrators[t_id], ps, th, ph)
    end

    AdaptiveSky(CoronaGridValues{T}, _trace, check_refine)
end

"""
    fill_sky_values(sky::AdaptiveSky{T,<:CoronaGridValues}, N)

Returns three dense arrays, two vectors representing the θ and ϕ, and a grid
with the averaged [`CoronaGridValues`](@ref) at that point, computed using
[`vector_average`](@ref). This effectively acts to fill in the sky with `N x N`
points using the information calculated by the adaptive sky.

    fill_sky_values(
        sky::AdaptiveSky{T,<:CoronaGridValues},
        phi_grid::AbstractVector,
        theta_grid::AbstractVector
    )

Returns only the grid.
"""
function fill_sky_values(sky::AdaptiveSky{T,V}, N) where {T,V<:CoronaGridValues}
    phi_grid = collect(range(-π, π, N))
    theta_grid = collect(range(0, π, N))
    out = fill_sky_values(sky, phi_grid, theta_grid)
    phi_grid, theta_grid, out
end

function fill_sky_values(
    sky::AdaptiveSky{T,V},
    phi_grid::AbstractVector,
    theta_grid::AbstractVector,
) where {T,V<:CoronaGridValues}
    output = V[]
    sizehint!(output, length(theta_grid) * length(phi_grid))

    index_set = Set{Int}()

    for ph in phi_grid
        empty!(index_set)
        for th in theta_grid
            push!(index_set, Grids.get_parent_index(sky.grid, cos(th), ph))
        end
        # remove all those without parents
        setdiff!(index_set, 0)

        if isempty(index_set)
            # fill with dummy values, then on to the next column
            append!(output, fill(make_null(V), length(theta_grid)))
            continue
        end

        # get the thetas of all rows in the column
        # TODO: use buffers
        column = T[acos(sky.grid.cells[i].pos[1]) for i in index_set]
        values = V[sky.values[i] for i in index_set]

        I = sortperm(column)
        permute!(column, I)
        permute!(values, I)

        for th in theta_grid
            # linear interpolate the values
            idx = clamp(searchsortedlast(column, th), 1, lastindex(column) - 1)
            x1, y1 = column[idx], values[idx]
            x2, y2 = column[idx+1], values[idx+1]

            w = (th - x1) / (x2 - x1)
            y = vector_average(SVector{2,T}(1 - w, w), SVector{2,V}(y1, y2))

            push!(output, y)
        end
    end

    reshape(output, (length(phi_grid), length(theta_grid)))
end

"""
    function bin_emissivity_grid(
        m::AbstractMetric,
        d::AbstractAccretionDisc,
        r_bins,
        ϕ_bins,
        sky::AdaptiveSky{T, <:CoronaGridValues}
    )

Bin the emissivity into a grid described by the axes `r_bins` and `phi_bins`,
returning the emissivity calcuated for each bin in a grid.

See also [`bin_emissivity_grid!`](@ref).
"""
function bin_emissivity_grid(
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    r_bins,
    ϕ_bins,
    sky::AdaptiveSky{T,<:CoronaGridValues},
) where {T}
    _, out = bin_emissivity_grid!(
        zeros(eltype(r_bins), (length(r_bins), length(ϕ_bins))),
        zeros(eltype(r_bins), (length(r_bins), length(ϕ_bins))),
        m,
        d,
        r_bins,
        ϕ_bins,
        sky,
    )
    out
end

"""
    bin_emissivity_grid!(
        output,
        solid_angle,
        m::AbstractMetric,
        d::AbstractAccretionDisc,
        r_bins,
        ϕ_bins,
        sky::AdaptiveSky{T, <:CoronaGridValues}
    )

Like [`bin_emissivity_grid`](@ref), but the output and temporary `solid_angle`
grid can be pre-allocated. Asserts the dimensions of both are `(lenght(r_bins),
length(ϕ_bins))`.
"""
function bin_emissivity_grid!(
    output,
    solid_angle,
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    r_bins,
    ϕ_bins,
    sky::AdaptiveSky{T,V},
) where {T,V<:CoronaGridValues}
    @assert size(output) == size(solid_angle)
    @assert size(output) == (length(r_bins), length(ϕ_bins))

    _disc_velocity = Gradus._keplerian_velocity_projector(m, d)
    # calculates area with all relativistic corrections
    function area(r)
        x = SVector(0.0, r, π/2, 0.0)
        v_disc = _disc_velocity(x)
        _A = Gradus._proper_area(m, r, π/2)
        _lf = Gradus.lorentz_factor(m, x, v_disc)
        _lf * _A
    end

    r_max = maximum(r_bins)

    for (index, v) in enumerate(sky.values)
        if Grids.has_children(sky.grid, index) || isnan(v.r) || v.r > r_max
            continue
        end

        th = acos(sky.grid.cells[index].pos[1])
        ΔΩ = prod(Grids.get_cell_width(sky.grid, index)) / sin(th)

        r_i = searchsortedlast(r_bins, v.r)
        ϕ_i = searchsortedlast(ϕ_bins, mod2pi(v.ϕ))
        if (r_i != 0) && (ϕ_i != 0)
            # TODO: allow generic spectrum
            output[r_i, ϕ_i] += ΔΩ * v.J * (v.g^2 * area(v.r))
            solid_angle[r_i, ϕ_i] += ΔΩ
        end
    end

    for i in eachindex(solid_angle)
        if solid_angle[i] > 0
            output[i] /= solid_angle[i]
        end
    end

    solid_angle, output
end

"""
    interpolate_emissivity_grid!(output::AbstractMatrix{T}, r_bins, ϕ_bins)

Fills in the missing values from the `output` grid of
[`bin_emissivity_grid`](@ref) by interpolating over `ϕ_bins`.
"""
function interpolate_emissivity_grid!(output::AbstractMatrix{T}, r_bins, ϕ_bins) where {T}
    @assert size(output) == (length(r_bins), length(ϕ_bins))
    _zero = zero(eltype(output))
    for (i, r) in enumerate(r_bins)
        col = @views output[i, :]
        I = @. col != _zero

        if count(I) < 2
            continue
        end

        x = @views ϕ_bins[I]
        y = @views col[I]

        itp = NaNLinearInterpolator(x, y, NaN)

        for (j, ϕ) in enumerate(ϕ_bins)
            output[i, j] = abs(itp(ϕ))
        end
    end
    output
end
