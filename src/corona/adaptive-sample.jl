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

    reshape(output, (length(theta_grid), length(phi_grid)))
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
    sky::AdaptiveSky{T,<:CoronaGridValues};
    kwargs...
) where {T}
    solid_angle, output = bin_emissivity_grid!(
        zeros(eltype(r_bins), (length(r_bins), length(ϕ_bins))),
        zeros(eltype(r_bins), (length(r_bins), length(ϕ_bins))),
        m,
        d,
        r_bins,
        ϕ_bins,
        sky;
        kwargs...
    )

    for i in eachindex(solid_angle)
        if solid_angle[i] > 0
            output[i] /= solid_angle[i]
        end
    end

    output
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
    sky::AdaptiveSky{T,V};
    Γ = 2,
) where {T,V<:CoronaGridValues}
    @assert size(output) == size(solid_angle)
    @assert size(output) == (length(r_bins), length(ϕ_bins))

    _disc_velocity = _keplerian_velocity_projector(m, d)
    # calculates area with all relativistic corrections
    function area(r)
        x = SVector{4,typeof(r)}(0, r, π/2, 0)
        v_disc = _disc_velocity(x)
        _A = _proper_area(m, r, π/2)
        _lf = lorentz_factor(m, x, v_disc)
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
            output[r_i, ϕ_i] += ΔΩ * v.J * (v.g^Γ * area(v.r))
            solid_angle[r_i, ϕ_i] += ΔΩ
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

"""
    empty_fraction(r, ϕ, block)

Count the number of empty entries in `block` (i.e. `==(0)`) and return the
fraction of empty entries to total entries.
"""
empty_fraction(r, ϕ, block) = count(!=(0), block) / length(block)

"""
    function evaluate_refinement_metric(
        m::AbstractMetric,
        d::AbstractAccretionGeometry,
        sky::AdaptiveSky,
        r_bins,
        phi_bins;
        split = 6,
        metric = empty_fraction,
    )

Bins the [`AdaptiveSky`](@ref) into the ``(r, \\phi)`` plane, whereupon a
stencil of dimensions `size(plane) / split` is used to split the plane into a
number of blocks. For each block, `metric` is evaluated.

A 2D matrix is returned with element type `@NamedTuple(;r, ϕ, score)`, where the
`r` and `ϕ` are the ranges of `r` and `ϕ` that the block applies to, and `score`
is the evaluated metric for that block.
"""
function evaluate_refinement_metric(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    sky::AdaptiveSky,
    r_bins,
    phi_bins;
    split = 6,
    metric = empty_fraction,
)
    N = length(r_bins)
    output = bin_emissivity_grid(m, d, r_bins, phi_bins, sky)
    stencil_size = div(N, split)

    function _eval_metric(i, j)
        r = i:(i+stencil_size)
        ϕ = j:(j+stencil_size)
        block = @views output[i:(i+stencil_size), j:(j+stencil_size)]
        (; r, ϕ, score = metric(r, ϕ, block))
    end

    axis_itt = 1:(stencil_size÷4):(N-stencil_size)

    [_eval_metric(i, j) for i in axis_itt, j in axis_itt]
end

"""
    function step_block!(
        m::AbstractMetric,
        d::AbstractAccretionGeometry,
        sky::AdaptiveSky;
        check_threshold = true,
        threshold = 0.2,
        N = 500,
        top = 3,
        split = 6,
        percentage = 10,
        metric = empty_fraction,
        kwargs...,
    )

Refine the sky using a block-based strategy, where the block that scores the
lowest `metric` is refined. The intended use is to fill in gaps on the ``(r,
\\phi)`` plane on the disc.


Grids the sky into `N` ``r`` and ``\\phi`` bins, then uses
[`evaluate_refinement_metric`](@ref) to determine a given `metric` on the sky.
Selects the `top` number of lowest scoring cells for refinement.

If `check_threshold = true`, filters to keep only those blocks who's score is
less than `threshold`. With the [`empty_fraction`](@ref), this will act to
insure a given fraction of the cell is covered.

Returns `true` if some cells did not meet the threshold, else `false`.

The cells are then refined to an accuracy of `percentage`. All other kwargs are
forwarded to [`trace_step!`](@ref).
"""
function step_block!(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    sky::AdaptiveSky;
    check_threshold = true,
    threshold = 0.2,
    N = 500,
    top = 3,
    split = 6,
    percentage = 10,
    metric = empty_fraction,
    kwargs...,
)
    r_bins = collect(Grids._geometric_grid(isco(m), 1000.0, N))
    phi_bins = range(0, 2π, N)

    _counts = evaluate_refinement_metric(
        m,
        d,
        sky,
        r_bins,
        phi_bins;
        split = split,
        metric = metric,
    )

    _counts = filter(i -> i.score != 0, vec(_counts))
    sort!(_counts; by = i -> i.score)

    if check_threshold
        _counts = filter(i -> i.score < threshold, _counts)
        if length(_counts) == 0
            return false
        end
    end

    lims = map(1:min(length(_counts), top)) do index
        r_min, r_max = @views extrema(r_bins[_counts[index].r])
        phi_min, phi_max = @views extrema(phi_bins[_counts[index].ϕ])
        (; r = (r_min, r_max), ϕ = (phi_min, phi_max))
    end

    refiner = fine_refine_function(; percentage = percentage) do v
        for lim in lims
            in_r = (v.r >= lim.r[1]) && (v.r <= lim.r[2])
            in_ϕ = (mod2pi(v.ϕ) >= lim.ϕ[1]) && (mod2pi(v.ϕ) <= lim.ϕ[2])
            if in_r && in_ϕ
                return true
            end
        end
        false
    end

    trace_step!(sky; check_refine = refiner, kwargs...)
    true
end

"""
    refine_function(f)

Used to define a new `check_refine` function for an [`AdaptiveSky`](@ref). This
can be passed e.g. to [`refine_all`](@ref).

The function `f` is given each cell's value and should return `true` if the
cell should be refined, else `false`.

# Warning

This function currently only works with `AdaptiveSky` when the value type is
[`CoronaGridValues`](@ref), as it first checks if the radii are `NaN`.
"""
function refine_function(f)
    function _refine(sky::AdaptiveSky{T,<:CoronaGridValues}, i1, i2) where {T}
        v1 = sky.values[i1]
        v2 = sky.values[i2]
        if isnan(v1.r) && isnan(v2.r)
            return false
        end
        f(v1) || f(v2)
    end
end

"""
    fine_refine_function(f)

Like [`refine_function`](@ref) but calls the original refine function in
[`AdaptiveSky`](@ref) first. This is "fine" in the sense that it applies
fine-grained refinement criteria.
"""
function fine_refine_function(f; kwargs...)
    function _refine(sky, i1, i2)
        if check_refine(sky, i1, i2; kwargs...)
            v1 = sky.values[i1]
            v2 = sky.values[i2]
            f(v1) || f(v2)
        else
            false
        end
    end
end

"""
    function adaptive_solve!(
        m::AbstractMetric,
        d::AbstractAccretionDisc,
        sky::AdaptiveSky;
        verbose = true,
        limit = 50,
        kwargs...,
    )

Adaptively solve an illumination profile on the accretion disc for the corona's
[`AdaptiveSky`](@ref) using repeated calls to [`step_block!`](@ref) until the
threshold is satisfied for all blocks.

This method traces the initial setup using [`trace_initial!`](@ref), followed
by `trace_calls + 1` calls to [`trace_step!`](@ref).

All `kwargs` are forwarded to [`step_block!`](@ref). The `limit` keyword can be
used to limit the number of iterations of [`step_block!`](@ref).

If `verbose=true` a progress indicator is printed to the terminal.

Returns the total number of refinement iterations taken.
"""
function adaptive_solve!(
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    sky::AdaptiveSky;
    verbose = true,
    limit = 50,
    trace_calls = 2,
    kwargs...,
)
    i::Int = 0
    progress = ProgressMeter.ProgressUnknown(
        desc = "Geodesics",
        color = :none,
        showspeed = true,
        enabled = verbose,
    )

    function showvals(N)
        () -> [("Refining", N), ("Iteration", i)]
    end

    # trace the initial sky
    Gradus.trace_initial!(sky)
    # this happens outside the loop as the first is so fast it tends to ruin
    # the output of the progress bar
    Gradus.trace_step!(sky)
    i += 1

    for _ = 1:trace_calls
        Gradus.trace_step!(sky; progress_bar = progress, showvalues = showvals)
        i += 1
    end

    # ensure the distant radii are well refined
    Gradus.trace_step!(
        sky;
        check_refine = fine_refine_function(v -> v.r > 800; percentage = 10),
        progress_bar = progress,
        showvalues = showvals,
    )
    i += 1

    while step_block!(
        m,
        d,
        sky;
        split = 5,
        top = 4,
        progress_bar = progress,
        showvalues = showvals,
        kwargs...,
    )
        i += 1
        if i >= limit
            break
        end
    end

    Gradus.ProgressMeter.finish!(progress)

    i
end
