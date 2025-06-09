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
vector_average(weight::Vector{<:Number}, values::Vector{V})::V where{V}
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

AdaptiveSky(V::Type, calc_v, check_refine; pole_offset = 1e-6) = AdaptiveSky(
    Grids.AdaptiveGrid(cos(pole_offset), cos(π - pole_offset), -π, π),
    V[],
    calc_v,
    check_refine,
)

"""
refine_all!(sky::AdaptiveSky, to_refine::Vector{Int})

Refine all cells in `to_refine` using a multi-threaded loop.
"""
function refine_all!(sky::AdaptiveSky, to_refine::Vector{Int})
    Threads.@threads for index in to_refine
        cell = sky.grid.cells[index]
        th = acos(cell.pos[1])
        value = sky.calculate_value(th, cell.pos[2])
        sky.values[index] = value
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
    trace_step!(sky::AdaptiveSky)

Apply the refinement metric across each cell boundary and refine the cells
where the metric is `true`.
"""
function trace_step!(sky::AdaptiveSky{T,V}) where {T,V}
    N = lastindex(sky.grid.cells)

    to_refine = Set{Int}()
    for cell_index = 1:N
        # skip those we've already refined
        if Grids.has_children(sky.grid, cell_index)
            continue
        end

        bordering = Grids.get_bordering(sky.grid, cell_index)

        r1 = sky.values[cell_index].r
        g1 = sky.values[cell_index].g
        J1 = sky.values[cell_index].J

        for index in bordering
            if sky.check_refine(sky, cell_index, index)
                push!(to_refine, index)
            end
        end
    end

    # now that all have been reaped, apply refining
    to_trace = map(collect(to_refine)) do index
        if Grids.has_children(sky.grid, index)
            Int[]
        else
            _children = Grids.refine!(sky.grid, index)
            for (j, _) in enumerate(_children)
                if j == 5
                    push!(sky.values, sky.values[index])
                else
                    push!(sky.values, make_null(V))
                end
            end
            # return all but the middle for tracing
            Int[_children[[1, 2, 3, 4, 6, 7, 8, 9]]...]
        end
    end

    if !isempty(to_trace)
        all_trace = reduce(vcat, to_trace)
        refine_all!(sky, all_trace)
    end

    sky
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
