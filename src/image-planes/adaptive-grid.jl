"""
A set of unit vectors that allows for iterating over all of the cardinal
directions.

The cardinal directions are
```
   | 1 |
---+---+---
 4 | x | 2
---+---+---
   | 3 |
```
"""
const CARDINAL_DIRECTIONS =
    SVector{2,Int}[SVector{2}(x, y) for y in (-1, 0, 1), x in (-1, 0, 1)]

# the indices facing in a given direction
const FACING_TOP = (1, 4, 7)
const FACING_RIGHT = (7, 8, 9)
const FACING_DOWN = (3, 6, 9)
const FACING_LEFT = (1, 2, 3)
const FACINGS = (FACING_TOP, FACING_RIGHT, FACING_DOWN, FACING_LEFT)
# maps which direction something is looking to which way the neighbouring cell
# ought to be facing
const FACING_LOOKUP = Int[3, 4, 1, 2]

"""
    AdaptiveCell{P}

A cell in the [`AdaptiveGrid`](@ref). `P` is the point type (conventionally an
SVector).
"""
struct AdaptiveCell{P}
    "Midpoint position of the cell"
    pos::P
    "The indices that this cell occupies for it's parent at different levels"
    locations::Vector{Int}
    # starts at 0 for the grid
    # 3^level gives the fraction of an exes that a cell covers
    "The refinement level of the grid"
    level::Int
end

"""
    AdaptiveGrid{P}(x1, x2, y1, y2; x_periodic = true)

Constructs an `AdaptiveGrid` with domain `x1` to `x2`, and `y1` to `y2`. If
`periodic = true` then the grid is periodic in the `x` direction.

An `AdaptiveGrid` is an implementation of a grid refinement structure, which
splits the domain into 3x3 cells which can be refined into 3x3 child cells and
so forth.

Each refinement is referred to as a `level`, with `level == 1` being the top
level (i.e. initial 9 cells). The indexing scheme for (child) cells is column
orientated:
```
 1 | 4 | 7
---+---+---
 2 | 5 | 8
---+---+---
 3 | 6 | 9
```
"""
struct AdaptiveGrid{P} <: Abstract2DGrid
    "All cells in the grid."
    cells::Vector{AdaptiveCell{P}}
    "The indexes of the children of the corresponding cell."
    children::Vector{Union{Nothing,NTuple{9,Int}}}
    "The indexes of the neighbours of the corresponding cell."
    neighbours::Vector{NTuple{4,Int}}
    "The limits of the grid domain."
    limits::NTuple{2,P}
    "An internal buffer used when gathering (boundary) cells indices."
    _buffer::Vector{Int}
end

function AdaptiveGrid(
    x_1::Number,
    x_2::Number,
    y_1::Number,
    y_2::Number;
    T = Float64,
    x_periodic = true,
)
    x1 = min(x_1, x_2)
    x2 = max(x_1, x_2)
    y1 = min(y_1, y_2)
    y2 = max(y_1, y_2)

    δ = SVector{2}(x2 - x1, y2 - y1)
    mid = δ/2 + SVector{2}(x1, y1)

    grid = AdaptiveGrid(
        AdaptiveCell{SVector{2,T}}[],
        Union{Nothing,NTuple{9,Int}}[],
        NTuple{4,Int}[],
        (SVector{2,T}(x1, x2), SVector{2,T}(y1, y2)),
        Int[],
    )

    # initialise the top level cells
    for i = 1:9
        pos = mid + (δ/3 .* CARDINAL_DIRECTIONS[i])
        _add_cell!(grid, AdaptiveCell(pos, [i], 1))
    end

    _set_neighbours!(grid, 1, top = 0, left = x_periodic ? 7 : 0, right = 4, bottom = 2)
    _set_neighbours!(grid, 2, top = 1, left = x_periodic ? 8 : 0, right = 5, bottom = 3)
    _set_neighbours!(grid, 3, top = 2, left = x_periodic ? 9 : 0, right = 6, bottom = 0)

    _set_neighbours!(grid, 4, top = 0, left = 1, right = 7, bottom = 5)
    _set_neighbours!(grid, 5, top = 4, left = 2, right = 8, bottom = 6)
    _set_neighbours!(grid, 6, top = 5, left = 3, right = 9, bottom = 0)

    _set_neighbours!(grid, 7, top = 0, left = 4, right = x_periodic ? 1 : 0, bottom = 8)
    _set_neighbours!(grid, 8, top = 7, left = 5, right = x_periodic ? 2 : 0, bottom = 9)
    _set_neighbours!(grid, 9, top = 8, left = 6, right = x_periodic ? 3 : 0, bottom = 0)

    grid
end

function _add_cell!(grid::AdaptiveGrid, cell::AdaptiveCell)
    push!(grid.cells, cell)
    push!(grid.children, nothing)
    push!(grid.neighbours, (0, 0, 0, 0))
    lastindex(grid.cells)
end

function _set_children!(grid::AdaptiveGrid, cell_index::Int, children::NTuple{9,Int})
    grid.children[cell_index] = children
end

function _set_neighbours!(grid::AdaptiveGrid, cell_index::Int, neighbours::NTuple{4,Int})
    grid.neighbours[cell_index] = neighbours
end

function _set_neighbours!(
    grid::AdaptiveGrid,
    cell_index::Int;
    top::Int,
    bottom::Int,
    right::Int,
    left::Int,
)
    _set_neighbours!(grid, cell_index, (top, right, bottom, left))
end

"""
    extent(grid::AdaptiveGrid)

Returns the width of the domain of the [`AdaptiveGrid`](@ref), i.e. `Δx`, `Δy`.
"""
function extent(grid::AdaptiveGrid)
    Δx = abs(grid.limits[1][2] - grid.limits[1][1])
    Δy = abs(grid.limits[2][2] - grid.limits[2][1])
    Δx, Δy
end

"""
    has_children(grid::AdaptiveGrid, cell_index::Int)

Returns `true` if the cell at `cell_index` is refined (i.e. has children).
"""
has_children(grid::AdaptiveGrid, cell_index::Int) = !isnothing(grid.children[cell_index])

"""
    refine!(grid::AdaptiveGrid, cell_index::Int)

Refine the cell at index `cell_index`. That is, split the cell into 9 child
cells. Returns a tuple of cell indices of the newly added children.
"""
function refine!(grid::AdaptiveGrid, cell_index::Int)
    has_children(grid, cell_index) && throw("$cell_index already refined")

    cell = grid.cells[cell_index]
    _neighbours = grid.neighbours[cell_index]

    fraction = 3^(cell.level + 1)
    Δx, Δy = extent(grid)
    δ = SVector{2}(Δx, Δy)/fraction

    children = map((1:9...,)) do i
        pos = cell.pos + (δ .* CARDINAL_DIRECTIONS[i])
        _add_cell!(grid, AdaptiveCell(pos, vcat(cell.locations, i), cell.level+1))
    end
    _set_children!(grid, cell_index, children)

    # top row of children
    _set_neighbours!(
        grid,
        children[1];
        top = _neighbours[1],
        right = children[4],
        left = _neighbours[4],
        bottom = children[2],
    )
    _set_neighbours!(
        grid,
        children[4];
        top = _neighbours[1],
        right = children[7],
        left = children[1],
        bottom = children[5],
    )
    _set_neighbours!(
        grid,
        children[7];
        top = _neighbours[1],
        right = _neighbours[2],
        left = children[4],
        bottom = children[8],
    )
    # middle row
    _set_neighbours!(
        grid,
        children[2];
        top = children[1],
        right = children[5],
        left = _neighbours[4],
        bottom = children[3],
    )
    _set_neighbours!(
        grid,
        children[5];
        top = children[4],
        right = children[8],
        left = children[2],
        bottom = children[6],
    )
    _set_neighbours!(
        grid,
        children[8];
        top = children[7],
        right = _neighbours[2],
        left = children[5],
        bottom = children[9],
    )
    # bottom row
    _set_neighbours!(
        grid,
        children[3];
        top = children[2],
        right = children[6],
        left = _neighbours[4],
        bottom = _neighbours[3],
    )
    _set_neighbours!(
        grid,
        children[6];
        top = children[5],
        right = children[9],
        left = children[3],
        bottom = _neighbours[3],
    )
    _set_neighbours!(
        grid,
        children[9];
        top = children[8],
        right = _neighbours[2],
        left = children[6],
        bottom = _neighbours[3],
    )

    map((1:9...,)) do i
        children[i]
    end
end

_add_facing!(buffer, index) = push!(buffer, index)

"""
    _get_position(location::Int, direction::Int)

Given the location (i.e. index between 1 and 9) and a facing direction, returns which position the cell is in.

For direction 2 or 4 (right or left)
- 1 is top, 2 is middle, 3 is bottom

For direction 1 or 3 (up or down)
- 1 is left, 2 is middle, 3 is right
"""
function _get_position(location::Int, direction::Int)
    if iseven(direction)
        rem(location - 1, 3) + 1
    else
        div(location - 1, 3) + 1
    end
end

_get_position(cell::AdaptiveCell, direction::Int) = _get_position(cell.location, direction)

function _get_facing!(
    buffer::AbstractVector{Int},
    grid::AdaptiveGrid,
    cell_index::Int,
    direction::Int,
    level::Int,
    locations::Vector{Int},
)
    children = grid.children[cell_index]

    if isnothing(children)
        _add_facing!(buffer, cell_index)
        return
    end

    self = grid.cells[cell_index]

    # if at the same level and same locations, add all children
    if self.level < level
        position = _get_position(locations[self.level+1], direction)
        index = FACINGS[direction][position]
        _get_facing!(buffer, grid, children[index], direction, level, locations)
    else
        for index in FACINGS[direction]
            _get_facing!(buffer, grid, children[index], direction, level, locations)
        end
    end
end

"""
    get_bordering(grid::AdaptiveGrid, cell_index::Int)

Get all of the bordering cells (i.e. the indices of all children surrounding
the given cell). Returns a vector of indices.

# Warning

This function uses an internal buffer to avoid excessive allocations. It is
consequently **not thread safe**, and if it must be called twice, the caller
should copy the memory between calls to avoid overwriting the previous values.

To make this threadsafe, you can pass the `buffer` keyword.
"""
function get_bordering(grid::AdaptiveGrid, cell_index::Int; buffer = grid._buffer)
    empty!(buffer)

    cell = grid.cells[cell_index]
    neighbours = grid.neighbours[cell_index]

    for i = 1:4
        if neighbours[i] != 0
            _get_facing!(
                buffer,
                grid,
                neighbours[i],
                FACING_LOOKUP[i],
                cell.level,
                cell.locations,
            )
        end
    end

    buffer
end

"""
    Base.contains(grid::AdaptiveGrid, cell_index::Int, x, y)

Returns `true` if the cell at `cell_index` contains the point `(x, y)`.
"""
function Base.contains(grid::AdaptiveGrid, cell_index::Int, x, y)
    cell = grid.cells[cell_index]

    Δx, Δy = get_cell_width(grid, cell_index) ./ 2

    x_good = ((cell.pos[1] - Δx) < x) && ((cell.pos[1] + Δx) > x)
    y_good = ((cell.pos[2] - Δy) < y) && ((cell.pos[2] + Δy) > y)
    x_good && y_good
end

function _relative_location(grid::AdaptiveGrid, cell_index::Int, x, y)
    cell = grid.cells[cell_index]
    if !contains(grid, cell_index, x, y)
        return nothing
    end

    # divide by 3 to get child size, then by 2 to get +/-
    Δx, Δy = get_cell_width(grid, cell_index) ./ 6

    x_offset = ((cell.pos[1]) - Δx > x) ? -1 : ((cell.pos[1]) + Δx < x) ? 1 : 0
    y_offset = ((cell.pos[2]) - Δy > y) ? -1 : ((cell.pos[2]) + Δy < y) ? 1 : 0

    ((1 + x_offset) * 3) + y_offset + 2
end

"""
    get_parent_index(grid::AdaptiveGrid, x, y)

Get the index of the most refined cell that contains the point `(x, y)`.
"""
function get_parent_index(grid::AdaptiveGrid, x, y)
    current_cell = 0

    for i = 1:9
        if contains(grid, i, x, y)
            current_cell = i
            break
        end
    end

    (current_cell == 0) && return 0

    while has_children(grid, current_cell)
        I = _relative_location(grid, current_cell, x, y)
        if isnothing(I)
            throw("Unreachable")
        end
        current_cell = grid.children[current_cell][I]
    end

    current_cell
end

function get_cell_width(grid::AdaptiveGrid, cell_index::Int)
    cell = grid.cells[cell_index]
    fraction = 3^cell.level
    Δx, Δy = extent(grid) ./ fraction
    Δx, Δy
end
