struct CunninghamTransferData{T}
    "``g^\\ast`` values."
    g✶::Vector{T}
    "Transfer function data."
    f::Vector{T}
    "Timing data."
    t::Vector{T}
    gmin::T
    gmax::T
    "Emission radius."
    rₑ::T
end

struct CunninghamTransferGrid{V<:AbstractVector,M<:AbstractMatrix}
    r_grid::V
    g✶_grid::V
    # the minimum g for each r
    g_min::V
    # the maximum g for each r
    g_max::V
    # the matrix is column major, and since we will be interpolating over g we
    # want to have each column be a different r
    lower_f::M
    upper_f::M
    lower_time::M
    upper_time::M
end

function Base.show(io::IO, ctg::CunninghamTransferGrid)
    T = eltype(ctg.r_grid)
    N = length(ctg.r_grid)
    print(
        io,
        "CunninghamTransferGrid{T=$T, N=$N, rmin=$(ctg.r_grid[1]), rmax=$(ctg.r_grid[end])}",
    )
end

function inner_radius(ctg::CunninghamTransferGrid)
    ctg.r_grid[1]
end

function outer_radius(ctg::CunninghamTransferGrid)
    ctg.r_grid[end]
end

# for whatever reason, `fieldnames` always seems to fail to type infer here
# so just unroll it by hand
function MultiLinearInterpolations.restructure(
    grid::CunninghamTransferGrid,
    vs::AbstractVector,
)
    @views begin
        start = 1

        stop = start + length(grid.r_grid) - 1
        r_grid = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.g✶_grid) - 1
        g✶_grid = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.g_min) - 1
        g_min = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.g_max) - 1
        g_max = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.lower_f) - 1
        lower_f = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.upper_f) - 1
        upper_f = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.lower_time) - 1
        lower_time = vs[start:stop]
        start = stop + 1

        stop = start + length(grid.upper_time) - 1
        upper_time = vs[start:end]

        CunninghamTransferGrid(
            r_grid,
            g✶_grid,
            g_min,
            g_max,
            reshape(lower_f, size(grid.lower_f)),
            reshape(upper_f, size(grid.upper_f)),
            reshape(lower_time, size(grid.lower_time)),
            reshape(upper_time, size(grid.upper_time)),
        )
    end
end

struct CunninghamTransferTable{N,T,I<:MultilinearInterpolator}
    params::NTuple{N,Vector{T}}
    grids::Array{CunninghamTransferGrid{Vector{T},Matrix{T}},N}
    cache::I
end

function CunninghamTransferTable(
    x::NTuple{N},
    grids::AbstractArray{<:CunninghamTransferGrid},
) where {N}
    _grids = reshape(grids, length.(x))
    cache = MultilinearInterpolator{N}(_grids)
    CunninghamTransferTable(x, _grids, cache)
end

(table::CunninghamTransferTable{1})(x::Number) = table((x,))
function (table::CunninghamTransferTable{N})(x::Tuple) where {N}
    @assert N == length(x)
    interpolate!(table.cache, table.params, table.grids, promote(x...))
end

struct TransferBranches{SameDomain,T,F}
    upper_f::F
    lower_f::F
    upper_t::F
    lower_t::F
    gmin::T
    gmax::T
    rₑ::T
    function TransferBranches{SameDomain}(
        upper_f::F,
        lower_f::F,
        upper_t::F,
        lower_t::F,
        gmin::T,
        gmax::T,
        rₑ::T,
    ) where {SameDomain,F,T}
        new{SameDomain,T,F}(upper_f, lower_f, upper_t, lower_t, gmin, gmax, rₑ)
    end
end

struct InterpolatingTransferBranches{T,F,SameDomain}
    branches::Vector{TransferBranches{SameDomain,T,F}}
    radii::Vector{T}
    gmin::Vector{T}
    gmax::Vector{T}
end

function outer_radius(itb::InterpolatingTransferBranches)
    itb.radii[end]
end

function inner_radius(itb::InterpolatingTransferBranches)
    itb.radii[1]
end

struct LagTransferFunction{T,X,E,P}
    max_t::T
    x::X
    image_plane_areas::Vector{T}
    coronal_geodesics::E
    observer_to_disc::Vector{P}
end

function Base.show(io::IO, ::MIME"text/plain", tf::LagTransferFunction)
    text = """LagTransferFunction for $(typeof(tf.coronal_geodesics.metric))
      . observer position
          $(tf.x)
      . model                         : $(typeof(tf.coronal_geodesics.model))
      . observer to disc photon count : $(length(tf.observer_to_disc))
      . source to disc photon count   : $(length(tf.coronal_geodesics.geodesic_points))
      Total memory: $(Base.format_bytes(Base.summarysize(tf)))
    """
    print(io, text)
end
