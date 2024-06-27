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
        @warn "Interpolation out of bounds $r ∉ [$(r_min),  $(r_max)]. Additional geodesic samples may be required (will not log again)." maxlog =
            1
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

restructure(::Number, vals::AbstractVector) = first(vals)
restructure(::AbstractArray, vals::AbstractVector) = vals
restructure(::NTuple{N}, vals::AbstractVector) where {N} = ((vals[i] for i = 1:N)...,)

function _tuple_set(tuple::NTuple{N}, index, v)::NTuple{N} where {N}
    if index == 1
        (v, tuple[2:end]...)
    elseif index == N
        (tuple[1:end-1]..., v)
    else
        (tuple[1:index-1]..., v, tuple[index+1:end]...)
    end
end

function _dual_size(len::Int; chunk_size = ForwardDiff.pickchunksize(len))
    cs = prod((chunk_size * ones(Int, 1)) .+ 1)
    cs * len
end

struct MultilinearInterpolator{D,T}
    indices::Vector{NTuple{D,Int}}
    weights::Vector{T}
    output::Vector{T}
    output_len::Int
    function MultilinearInterpolator{D}(
        values::AbstractArray{C};
        T = Float64,
        kwargs...,
    ) where {D,C}
        size_points = 2^D
        indices = Vector{NTuple{D,Int}}(undef, size_points)
        weights = zeros(T, _dual_size(D; kwargs...))

        len = if C <: Number
            1
        elseif eltype(C) <: Number
            length(values[1])
        else
            sum(fieldnames(C)) do f
                length(getproperty(values[1], f))
            end
        end

        output = Vector{T}(undef, _dual_size(len + 1; kwargs...))
        new{D,T}(indices, weights, output, len)
    end
end

_reinterpret_dual(::Type, v::AbstractArray, n::Int) = view(v, 1:n)
function _reinterpret_dual(
    DualType::Type{<:ForwardDiff.Dual},
    v::AbstractArray{T},
    n::Int,
) where {T}
    n_elems = div(sizeof(DualType), sizeof(T)) * n
    if n_elems > length(v)
        @warn "Resizing..."
        resize!(v, n_elems)
    end
    reinterpret(DualType, view(v, 1:n_elems))
end


function update_indices!(
    cache::MultilinearInterpolator{D},
    axes::NTuple{D},
    x::NTuple{D,<:Number},
) where {D}
    its = ((1:D)...,)

    weights = _reinterpret_dual(typeof(first(x)), cache.weights, D)

    map(its) do I
        stride = 2^(D - I)
        ax = axes[I]
        x0 = x[I]

        i = min(lastindex(ax), searchsortedfirst(ax, x0))
        i1, i2 = if i == 1
            1, 2
        else
            i - 1, i
        end

        weights[I] = (x0 - ax[i1]) / (ax[i2] - ax[i1])

        for (q, j) in enumerate(range(1, lastindex(cache.indices), step = stride))
            for k = j:j+stride-1
                tup = cache.indices[k]
                if !iseven(q)
                    cache.indices[k] = _tuple_set(tup, I, i1)
                else
                    cache.indices[k] = _tuple_set(tup, I, i2)
                end
            end
        end

        nothing
    end
    cache.indices
end

function _build_multilinear_expression(D::Int, field_name)
    function _lerp(y1, y2, w)
        :($(y2) * $(w) + $(y1) * (1 - $(w)))
    end
    assignments = []
    _index(i) =
        if !isnothing(field_name)
            sym = Base.gensym()
            push!(
                assignments,
                :(
                    $(sym) = getproperty(
                        values[cache.indices[$(i)]...],
                        $(Meta.quot(field_name)),
                    )
                ),
            )
            sym
        else
            :(values[cache.indices[$(i)]...])
        end
    _weight(i) = :(weights[$(i)])

    weight_index = D

    # get the knots
    knots = map(1:2^(D-1)) do d
        i = d * 2
        _lerp(_index(i - 1), _index(i), _weight(weight_index))
    end

    while length(knots) > 1
        weight_index -= 1
        knots = map(range(1, lastindex(knots), step = 2)) do i
            _lerp(knots[i], knots[i+1], _weight(weight_index))
        end
    end

    assignments, first(knots)
end

@inline @generated function _interpolate_cache!(
    cache::MultilinearInterpolator{D},
    values::AbstractArray{T,D},
    p::NTuple,
) where {D,T}
    assignments = []
    exprs = if T <: Number || eltype(T) <: Number
        _, interp = _build_multilinear_expression(D, nothing)
        expr = quote
            @. output = $(interp)
        end
        [expr]
    else
        interps = (_build_multilinear_expression(D, i) for i in fieldnames(T))
        map(zip(interps, fieldnames(T))) do I
            (assign, interp), f = I
            append!(assignments, assign)
            sym = Base.gensym()
            quote
                start = stop + 1
                shape = size(getproperty(values[cache.indices[1]...], $(Meta.quot(f))))
                @views begin
                    $(sym) = if length(shape) > 0
                        stop = start + prod(shape) - 1
                        reshape(output[start:stop], shape)
                    else
                        stop = start
                        output[start:stop]
                    end
                    @. $(sym) = $(interp)
                end
            end
        end
    end
    quote
        begin
            weights = _reinterpret_dual(eltype(p), cache.weights, D)
            output = _reinterpret_dual(eltype(p), cache.output, cache.output_len)

            start::Int = 0
            stop::Int = 0
            $(assignments...)
            $(exprs...)

            output
        end
    end
end

function interpolate!(
    cache::MultilinearInterpolator{D},
    axes::NTuple{D},
    vals::AbstractArray{T,D},
    x::NTuple{D,<:Number},
) where {D,T}
    update_indices!(cache, axes, x)
    output = _interpolate_cache!(cache, vals, x)
    restructure(first(vals), output)
end

export MultilinearInterpolator, interpolate!
