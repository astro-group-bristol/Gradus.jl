"""
    AdaptiveSky(V::Type, calc_v::Function, check_refine::Function)

Create an adaptive grid for the local sky, where each point is assigned a value
`V`. Must also provide two functions with the following prototypes:

```julia
# Calculate the new grid value at `(theta, phi)`
calculate_value(theta, phi)::V

# Return true if the cells being compared need to be refined.
check_refine(sky::AdaptiveSky, i1::Int, i2::Int)::Bool
```

The values for a given cell can be access with `sky.values[i]`.

For the interpolation schemes to work, the following function must also be
defined for the value type:
```julia
vector_average(weight::AbstractVector{<:Number}, values::AbstractVector{V})::V where{V}
# Return something that represents a null / NaN of the given type
make_null(V::Type)::V
```

"""
struct AdaptiveSky{T,V,CalcV,RefineV}
    "Backing adaptive grid"
    grid::Grids.AdaptiveGrid{SVector{2,T}}
    "User defined values for each cell."
    values::Vector{V}
    "Constructs new values"
    calculate_value::CalcV
    "Should return `true` to refine cells"
    check_refine::RefineV
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(sky::AdaptiveSky))
    println(io, "AdaptiveSky[N=$(length(sky.values))]")
end

AdaptiveSky(V::Type, calc_v, check_refine; pole_offset = 1e-6) = AdaptiveSky(
    Grids.AdaptiveGrid(cos(pole_offset), cos(π - pole_offset), -π, π),
    V[],
    calc_v,
    check_refine,
)

"""
trace_all!(sky::AdaptiveSky, to_refine::Vector{Int})

Refine all cells in `to_refine` using a multi-threaded loop.
"""
function trace_all!(
    sky::AdaptiveSky,
    to_refine::Vector{Int};
    verbose = false,
    progress_bar = init_progress_bar(
        "Refining $(length(to_refine)) cell(s): ",
        length(to_refine),
        verbose,
    ),
    showvalues = nothing,
)
    prog = if isnothing(progress_bar)
        progress_bar = init_progress_bar(
            "Refining $(length(to_refine)) cell(s): ",
            length(to_refine),
            verbose,
        )
    else
        progress_bar
    end

    Threads.@threads for index in to_refine
        cell = sky.grid.cells[index]
        th = acos(cell.pos[1])
        value = sky.calculate_value(th, cell.pos[2])
        sky.values[index] = value
        # update progress
        if isnothing(showvalues)
            ProgressMeter.next!(prog)
        else
            ProgressMeter.next!(prog; showvalues = showvalues(length(to_refine)))
        end
    end
end

function Grids.refine!(sky::AdaptiveSky, cell_index::Int)
    children = Grids.refine!(sky.grid, cell_index)
    for (j, child_index) in enumerate(children)
        if j == 5 # midpoint will be the same as parent
            push!(sky.values, sky.values[cell_index])
        else
            child = sky.grid.cells[child_index]
            th = acos(child.pos[1])
            value = sky.calculate_value(th, child.pos[2])
            push!(sky.values, value)
        end
    end
end

"""
    trace_initial!(sky::AdaptiveSky; level = 3)

Initialise the sky to a given level of refinement.
"""
function trace_initial!(sky::AdaptiveSky; level = 3)
    for cell in sky.grid.cells
        th = acos(cell.pos[1])
        value = sky.calculate_value(th, cell.pos[2])
        push!(sky.values, value)
    end

    for level = 1:(level-1)
        N = lastindex(sky.grid.cells)
        for i = 1:N
            if !Grids.has_children(sky.grid, i)
                Grids.refine!(sky, i)
            end
        end
    end

    sky
end

"""
    find_need_refine(sky::AdaptiveSky, check_refine::Function)::Vector{Int}

Check all cells with `check_refine`, and return those cell IDs that would need
to refined.
"""
function find_need_refine(sky::AdaptiveSky, check_refine::Function)::Vector{Int}
    N = lastindex(sky.grid.cells)

    to_refine = Set{Int}()

    # TODO: multi-thread using a divide and conquer approach
    for cell_index = 1:N
        # skip those we've already refined
        if Grids.has_children(sky.grid, cell_index)
            continue
        end

        bordering = Grids.get_bordering(sky.grid, cell_index)

        for index in bordering
            if check_refine(sky, cell_index, index)
                push!(to_refine, index)
            end
        end
    end

    sort!(collect(to_refine))
end

"""
    function refine_and_trace!(
        sky::AdaptiveSky,
        cell_ids::Vector{Int};
        kwargs...,
    )

Given a list of cell IDs, refine them using [`refine!`](@ref), and then trace
the children to populate their values.

All kwargs are forwarded to [`trace_all!`](@ref).
"""
function refine_and_trace!(
    sky::AdaptiveSky{T,V},
    cell_ids::Vector{Int};
    kwargs...,
) where {T,V}
    N = lastindex(sky.grid.cells)

    to_trace = Int[]
    sizehint!(to_trace, length(cell_ids) * 8)

    # now that all have been reaped, apply refining
    for index in cell_ids
        @assert !Grids.has_children(sky.grid, index)
        _children = Grids.refine!(sky.grid, index)
        for (j, child_id) in enumerate(_children)
            if j == 5
                push!(sky.values, sky.values[index])
            else
                push!(to_trace, child_id)
                push!(sky.values, make_null(V))
            end
        end
    end

    if !isempty(to_trace)

        trace_all!(sky, to_trace; kwargs...)
    end

    sky
end

"""
    trace_step!(sky::AdaptiveSky; check_refine = sky.check_refine, verbose = false)

Apply the refinement metric across each cell boundary and refine the cells
where the metric is `true`.

A different refinemenet metric from the default can be used by passing the
`check_refine` kwarg, using the same function signature as documented in
[`AdaptiveSky`](@ref).

If `verbose` is true, a progress bar will be displayed during refinement.
"""
function trace_step!(sky::AdaptiveSky; check_refine = sky.check_refine, kwargs...)
    to_refine = find_need_refine(sky, check_refine)
    refine_and_trace!(sky, to_refine; kwargs...)
end

vector_average(distances, values) = sum(i -> i[1] * i[2], zip(distances, values))

# TODO: use a buffer
function vector_average(sky::AdaptiveSky{T,V}, point, cell_index, bordering) where {T,V}
    # computes the weighted average from a number of unstructured points
    w(pos) = sum(i -> i^2, point - pos)

    values = zeros(V, length(bordering) + 1)
    distances = zeros(Float64, length(bordering) + 1)

    values[1] = sky.values[cell_index]
    distances[1] = w(point - sky.grid.cells[cell_index].pos)

    for (i, b) in enumerate(bordering)
        vi = if isnan(sky.values[b].r)
            base
        else
            sky.values[b]
        end
        values[i+1] = vi
        distances[i+1] = w(sky.grid.cells[b].pos)
    end

    # normalise distance vector
    total_distances = sqrt(sum(i -> i^2, distances))
    @. distances = (distances / total_distances)^2

    @assert isapprox(sum(distances), 1) "$(sum(distances))"

    vector_average(distances, values)
end

function _refine_interpolate!(sky::AdaptiveSky, cell_index::Int, max_level::Int)
    if isnan(sky.values[cell_index].r)
        return
    end

    bordering = Grids.get_bordering(sky.grid, cell_index)
    children = Grids.refine!(sky.grid, cell_index)

    for (i, child) in enumerate(children)
        if i == 5 # middle inherits from parent
            push!(sky.values, sky.values[cell_index])
        else
            push!(
                sky.values,
                vector_average(sky, sky.grid.cells[child].pos, cell_index, bordering),
            )
        end
    end

    for child in children
        if sky.grid.cells[child].level < max_level
            _refine_interpolate!(sky, child, max_level)
        end
    end
end

function interpolate_refine!(sky::AdaptiveSky)
    # find the highest level
    # refine all cells below that level without children
    # interpolate the value from the neighbours
    max_level = maximum(i -> i.level, sky.grid.cells)

    N = lastindex(sky.grid.cells)
    for i = 1:N
        cell = sky.grid.cells[i]
        if (!Grids.has_children(sky.grid, i)) && (cell.level < max_level)
            _refine_interpolate!(sky, i, max_level)
        end
    end

end

"""
    unpack_sky(sky::AdaptiveSky)

Return three arrays, representing the `x`, `y` coordinates of each point, and
it's associated value `V`.
"""
function unpack_sky(sky::AdaptiveSky{T,V}) where {T,V}
    X = Float64[]
    Y = Float64[]
    Z = V[]
    for i in eachindex(sky.grid.cells)
        if !Grids.has_children(sky.grid, i)
            cell = sky.grid.cells[i]
            vals = sky.values[i]
            push!(Y, cell.pos[1])
            push!(X, cell.pos[2])
            push!(Z, vals)
        end
    end
    X, Y, Z
end

export AdaptiveSky, trace_initial!, unpack_sky, trace_step!
