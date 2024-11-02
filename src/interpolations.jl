struct NaNLinearInterpolator{V1,V2,Y}
    t::V1
    u::V2
    default::Y
end

function _interpolate(interp::NaNLinearInterpolator, x)
    idx = clamp(searchsortedlast(interp.t, x), 1, length(interp.t) - 1)
    x1, x2 = interp.t[idx], interp.t[idx+1]
    y1, y2 = interp.u[idx], interp.u[idx+1]

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
    # @assert issorted(x) "x must be sorted!"
    NaNLinearInterpolator(x, y, zero(eltype(y)))
    # DataInterpolations.LinearInterpolation(y, x)
end

@inline function _enforce_interpolation_bounds(r::Number, r_min::Number, r_max::Number)
    if (r < r_min) || (r > r_max)
        @warn "Interpolation out of bounds $r ∉ [$(r_min),  $(r_max)]. Additional geodesic samples may be required (will not log again)." maxlog =
            1
    end
    clamp(r, r_min, r_max)
end
