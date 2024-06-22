struct NaNLinearInterpolator{X,Y}
    x::Vector{X}
    y::Vector{Y}
    default::Y
end

function _interpolate(interp::NaNLinearInterpolator{X}, x::X) where {X}
    idx = clamp(searchsortedlast(interp.x, x), 1, length(interp.x) - 1)
    x1, x2 = interp.x[idx], interp.x[idx+1]
    y1, y2 = interp.y[idx], interp.y[idx+1]

    w = (x - x1) / (x2 - x1)
    y = (1 - w) * y1 + w * y2

    if isnan(y)
        if (w < 0.5)
            isnan(y1) && return interp.default
            return y1
        else
            isnan(y2) && return interp.default
            return y2
        end
    end

    y
end

function (interp::NaNLinearInterpolator)(x)
    _interpolate(interp, x)
end

"""
    _make_interpolation(x, y)

Interpolate `y` over `x`

Utility method to wrap some interpolation library and provide the same interface
for our needs.
"""
function _make_interpolation(x, y)
    @assert size(x) == size(y) "size(x) = $(size(x)), size(y) = $(size(y))"
    @assert issorted(x) "x must be sorted!"
    # NaNLinearInterpolator(x, y, zero(eltype(y)))
    DataInterpolations.LinearInterpolation(y, x)
end

@inline function _enforce_interpolation_bounds(r::Number, r_min::Number, r_max::Number)
    if (r < r_min) || (r > r_max)
        @warn "Interpolation out of bounds $r ∉ [$(r_min),  $(r_max)]. Additional geodesic samples may be required."
    end
    clamp(r, r_min, r_max)
end

function _linear_interpolate(y1, y2, θ)
    @. (1 - θ) * y1 + θ * y2
end

function _linear_interpolate!(
    out::AbstractArray{<:Number},
    y1::Union{<:Number,AbstractArray{<:Number}},
    y2::Union{<:Number,AbstractArray{<:Number}},
    θ,
)
    @. out = (1 - θ) * y1 + θ * y2
end

function _linear_interpolate!(out::AbstractArray{T}, y1::T, y2::T, θ) where {T}
    _linear_interpolate!(out[1], y1, y2, θ)
end

function _linear_interpolate!(
    out::AbstractArray{T},
    y1::AbstractArray{T},
    y2::AbstractArray{T},
    θ,
) where {T}
    for i in eachindex(y1)
        _linear_interpolate!(out[i], y1[i], y2[i], θ)
    end
end


@inline function _linear_interpolate(arr::AbstractVector, idx, θ)
    _linear_interpolate(arr[idx], arr[idx+1], θ)
end

"""
    struct InterpolationCache{D,T,N}

A `D` dimensional interpolation cache.
"""
struct InterpolationCache{D,T,N}
    cache::Array{T,N}
    function InterpolationCache{D}(values::AbstractArray{T,1}) where {D,T}
        cache::Array{T,1} = [deepcopy(values[1])]
        new{D,T,1}(cache)
    end
    function InterpolationCache{D}(values::AbstractArray{T,N}) where {D,N,T}
        cache::Array{T,N - 1} = deepcopy(_get_dim_slice(values, Val{1}()))
        new{D,T,N - 1}(cache)
    end
end

@generated function _get_all_slice(values::AbstractArray{T,N}, i) where {T,N}
    rem = [:(:) for _ = 1:N-1]
    :(@views values[$(rem...), i])
end

@generated function _get_dim_slice(values::AbstractArray{T,N}, ::Val{M}) where {T,N,M}
    rem = [:(:) for _ = 1:(N-M)]
    inds = [:(1) for i = 1:M]
    :(@views values[$(rem...), $(inds...)])
end

function _make_cache_slices(cache::InterpolationCache{D}) where {D}
    itr = ((Val{i}() for i = 0:D-1)...,)
    map(itr) do i
        _get_dim_slice(cache.cache, i)
    end
end

function interpolate!(
    cache::InterpolationCache{D},
    grids::NTuple{D,<:AbstractArray},
    values::AbstractArray,
    x::NTuple{D},
) where {D}
    itr = (1:D...,)
    slices = _make_cache_slices(cache)
    vs = (values, slices...)
    _unroll_for(Val{D}(), (zip(slices, vs, itr)...,)) do K
        c, v, i = K
        _inplace_interpolate!(c, grids[i], v, x[i])
    end
    if D === 1
        first(cache.cache)
    else
        slices[D]
    end
end

function _set_value!(out::AbstractArray{<:Number}, value::AbstractArray{<:Number})
    @. out = value
end
function _set_value!(out::AbstractArray{<:Number}, v::Number)
    out[1] = v
end
function _set_value!(out::AbstractArray{<:T}, v::T) where {T}
    _set_value!(out[1], v)
end
function _set_value!(out::AbstractArray{<:T}, v::AbstractArray{<:T}) where {T}
    for (o, k) in zip(out, v)
        @views _set_value!(o, k)
    end
end

function _inplace_interpolate!(out, grid::AbstractArray, values::AbstractArray, x)
    i2 = searchsortedfirst(grid, x)
    if (i2 == 1)
        _set_value!(out, _get_all_slice(values, 1))
        return out
    end
    if i2 > lastindex(grid) || grid[i2] > grid[end]
        _set_value!(out, _get_all_slice(values, lastindex(grid)))
        return out
    end

    i1 = i2 - 1

    x1 = grid[i1]
    x2 = grid[i2]

    # interpolation weight
    θ = (x - x1) / (x2 - x1)
    y1 = _get_all_slice(values, i1)
    y2 = _get_all_slice(values, i2)
    _linear_interpolate!(out, y1, y2, θ)
    out
end
