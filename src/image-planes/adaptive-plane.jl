function _init_adaptive_tracer(
    m::AbstractMetric,
    x::SVector,
    d::AbstractAccretionGeometry;
    t_max = 2x[2],
    solver_opts...,
)
    v = map_impact_parameters(m, x, 1.0, 1.0)
    # init a reusable integrator
    integ = _init_integrator(
        m,
        x,
        v,
        d,
        t_max;
        save_on = false,
        # callback = domain_upper_hemisphere(),
        chart = chart_for_metric(m, 2 * 2x[2]),
        integrator_verbose = false,
        solver_opts...,
    )
end

function check_refine(
    sky::AdaptiveSky{T,<:GeodesicPoint},
    i1::Int,
    i2::Int;
    percentage = 20,
) where {T}
    v1 = sky.values[i1]
    v2 = sky.values[i2]

    hit1 = (v1.status == StatusCodes.IntersectedWithGeometry)
    hit2 = (v2.status == StatusCodes.IntersectedWithGeometry)

    if (!hit1) && (!hit2)
        return false
    end

    # if one hit but not the other
    if !hit1 || !hit2
        return true
    end

    r_too_coarse = !isapprox(log10(v1.x[2]), log10(v2.x[2]), rtol = percentage / 100)
    ph_too_coarse = !isapprox(mod2pi(v1.x[4]), mod2pi(v2.x[4]), rtol = percentage / 100)
    th_too_coarse = !isapprox(mod2pi(v1.x[3]), mod2pi(v2.x[3]), rtol = percentage / 100)

    r_too_coarse || th_too_coarse
end

"""
    function AdaptivePlane(
        m::AbstractMetric{T},
        x::SVector{4,T},
        d::AbstractAccretionGeometry;
        kwargs...,
    )
"""
function AdaptivePlane(
    m::AbstractMetric{T},
    x::SVector{4,T},
    d::AbstractAccretionGeometry;
    check_refine = check_refine,
    α_lims = (-10, 10),
    β_lims = (-10, 10),
    kwargs...,
) where {T}
    integ = _init_adaptive_tracer(m, x, d; kwargs...)

    # one integrator copy for each thread
    integrators = typeof(integ)[deepcopy(integ) for _ = 1:(Threads.nthreads()-1)]
    push!(integrators, integ)

    # tracing function
    function _trace(α, β)
        t_id = Threads.threadid()
        _integ = integrators[t_id]
        v = map_impact_parameters(m, x, -α, β)
        _solve_reinit!(_integ, vcat(x, v))
    end

    AdaptiveSky(
        Grids.AdaptiveGrid(α_lims..., β_lims...),
        GeodesicPoint{T,Nothing}[],
        _trace,
        check_refine,
    )
end

function fill_sky_values(sky::AdaptiveSky, N; kwargs...)
    x_grid = collect(range(sky.grid.limits[1]..., N))
    y_grid = collect(range(sky.grid.limits[2]..., N))
    out = fill_sky_values(sky, x_grid, y_grid; kwargs...)
    x_grid, y_grid, out
end

function fill_sky_values(
    sky::AdaptiveSky{T,V},
    x_grid::AbstractVector,
    y_grid::AbstractVector;
    metric,
    pf = PointFunction((m, gp, t) -> gp.x[1]),
) where {T,V<:GeodesicPoint}
    output = Union{Nothing,V}[]
    sizehint!(output, length(x_grid) * length(y_grid))

    column = T[]
    values = T[]
    sizehint!(column, length(y_grid))
    sizehint!(values, length(y_grid))

    for x in x_grid
        empty!(column)
        empty!(values)

        curr_index::Int = 0
        for y in y_grid
            index = Grids.get_parent_index(sky.grid, x, y)
            if index > 0 && (index != curr_index)
                curr_index = index

                cell_x = sky.grid.cells[index].pos[1]

                direction = (x < cell_x) ? 4 : 2
                # interpolate the value between closes neighbour
                neigh_index = sky.grid.neighbours[index][direction]
                neigh_x = sky.grid.cells[neigh_index].pos[2]

                v1 = pf(metric, sky.values[index], 0.0)
                v2 = pf(metric, sky.values[neigh_index], 0.0)

                if sky.values[neigh_index].status != StatusCodes.IntersectedWithGeometry
                    push!(values, v1)
                    push!(column, sky.grid.cells[index].pos[2])
                else
                    y = if direction == 4
                        w = (x - neigh_x) / (cell_x - neigh_x)
                        (1 - w) * v2 + w * v1
                    else
                        w = (x - cell_x) / (neigh_x - cell_x)
                        (1 - w) * v1 + w * v2
                    end
                    push!(values, y)
                    push!(column, sky.grid.cells[index].pos[2])
                end
            end
        end

        if isempty(column)
            # fill with dummy values, then on to the next column
            append!(output, fill(NaN, length(y_grid)))
            continue
        end

        # TODO: use buffers
        # column = T[sky.grid.cells[i].pos[2] for i in index_set]
        # values = V[sky.values[i] for i in index_set]

        # I = sortperm(column)
        # permute!(column, I)
        # permute!(values, I)

        for y in y_grid
            # linear interpolate the values
            idy = clamp(searchsortedlast(column, y), 1, lastindex(column) - 1)
            x1, y1 = column[idy], values[idy]
            x2, y2 = column[idy+1], values[idy+1]

            w = (y - x1) / (x2 - x1)

            q = (1 - w) * y1 + w * y2

            push!(output, q)
        end
    end

    reshape(output, (length(x_grid), length(y_grid)))
end
