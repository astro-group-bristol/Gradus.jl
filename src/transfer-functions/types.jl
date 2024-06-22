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

struct CunninghamTransferGrid{T}
    r_grid::Vector{T}
    g✶_grid::Vector{T}
    # the minimum g for each r
    g_min::Vector{T}
    # the maximum g for each r
    g_max::Vector{T}
    # the matrix is column major, and since we will be interpolating over g we
    # want to have each column be a different r
    lower_f::Matrix{T}
    upper_f::Matrix{T}
    lower_time::Matrix{T}
    upper_time::Matrix{T}
end

# for whatever reason, `fieldnames` always seems to fail to type infer here
# so just unroll it by hand
function _set_value!(out::CunninghamTransferGrid, v::CunninghamTransferGrid)
    _set_value!(out.r_grid, v.r_grid)
    _set_value!(out.g✶_grid, v.g✶_grid)
    _set_value!(out.g_min, v.g_min)
    _set_value!(out.g_max, v.g_max)
    _set_value!(out.lower_f, v.lower_f)
    _set_value!(out.upper_f, v.upper_f)
    _set_value!(out.lower_time, v.lower_time)
    _set_value!(out.upper_time, v.upper_time)
end
function _linear_interpolate!(
    out::CunninghamTransferGrid,
    y1::CunninghamTransferGrid,
    y2::CunninghamTransferGrid,
    θ,
)
    _linear_interpolate!(out.r_grid, y1.r_grid, y2.r_grid, θ)
    _linear_interpolate!(out.g✶_grid, y1.g✶_grid, y2.g✶_grid, θ)
    _linear_interpolate!(out.g_min, y1.g_min, y2.g_min, θ)
    _linear_interpolate!(out.g_max, y1.g_max, y2.g_max, θ)
    _linear_interpolate!(out.lower_f, y1.lower_f, y2.lower_f, θ)
    _linear_interpolate!(out.upper_f, y1.upper_f, y2.upper_f, θ)
    _linear_interpolate!(out.lower_time, y1.lower_time, y2.lower_time, θ)
    _linear_interpolate!(out.upper_time, y1.upper_time, y2.upper_time, θ)
end

struct CunninghamTransferTable{N,T,CacheT}
    params::NTuple{N,Vector{T}}
    grids::Vector{CunninghamTransferGrid{T}}
    cache::CacheT
end

function CunninghamTransferTable(
    x::NTuple{N},
    grids::AbstractVector{<:CunninghamTransferGrid},
) where {N}
    _grids = reshape(grids, length.(x))
    cache = InterpolationCache{N}(_grids)
    CunninghamTransferTable(x, _grids, cache)
end

(table::CunninghamTransferTable{1})(x::Number) = table((x,))
function (table::CunninghamTransferTable{N})(x::NTuple{N}) where {N}
    interpolate!(table.cache, table.params, table.grids, x)
end

struct TransferBranches{T,F}
    upper_f::F
    lower_f::F
    upper_t::F
    lower_t::F
    gmin::T
    gmax::T
    rₑ::T
end

struct InterpolatingTransferBranches{T,F}
    branches::Vector{TransferBranches{T,F}}
    radii::Vector{T}
    gmin::Vector{T}
    gmax::Vector{T}
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
