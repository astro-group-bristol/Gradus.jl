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

@inline function _linear_interpolate(y1, y2, θ)
    @. (1 - θ) * y1 + θ * y2
end

@inline function _linear_interpolate!(out, y1, y2, θ)
    @. out = (1 - θ) * y1 + θ * y2
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
    function InterpolationCache{D}(values::AbstractArray{T,N}) where {D,N,T}
        cache::Array{T,N - 1} = zeros(T, size(values)[1:N-1])
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
    slices[D]
end

function _inplace_interpolate!(out, grid::AbstractArray, values::AbstractArray, x)
    i2 = searchsortedfirst(grid, x)
    if (i2 == 1)
        @. out = values[1]
        return out
    end
    if i2 > lastindex(grid) || grid[i2] > grid[end]
        @. out = values[end]
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
