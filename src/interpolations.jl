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

function gaussian_kernel(kernel_size; σ = 1, domain = (-5, 5))
    function _gaussian(x, y)
        exp(-((x/σ)^2 + (y/σ)^2))
    end
    kernel = zeros(Float64, kernel_size)
    for (i, xi) in enumerate(range(domain..., kernel_size[1]))
        for (j, yj) in enumerate(range(domain..., kernel_size[2]))
            kernel[j, i] = _gaussian(xi, yj)
        end
    end
    kernel ./= sum(kernel)
    kernel
end

function constant_kernel(kernel_size)
    kernel = fill(1.0, kernel_size)
    kernel ./= sum(kernel)
    kernel
end

function kernel_interpolate!(
    data;
    kernel_size = (5, 5),
    kf = gaussian_kernel,
    verbose = false,
    kwargs...,
)
    output = deepcopy(data)
    kernel = kf(kernel_size; kwargs...)

    half_x = div(kernel_size[1], 2)
    half_y = div(kernel_size[2], 2)

    bar = Gradus.init_progress_bar("interpolation", size(output, 1), verbose)

    for i = (half_x+1):(size(output, 1)-half_x)
        for j = (half_y+1):(size(output, 2)-half_y)
            v = data[i, j]
            if isnan(v)
                slice = @views data[(i-half_x):(i+half_x), ((j-half_y):(j+half_y))]
                @assert size(kernel) == size(slice) "$(size(kernel)) != $(size(slice))"

                output[i, j] = 0
                weight = 0.0
                for (v, k) in zip(slice, kernel)
                    if !isnan(v)
                        output[i, j] += v * k
                        weight += k
                    end
                end

                if weight > 0
                    output[i, j] /= weight
                end
            end
        end
        Gradus.ProgressMeter.next!(bar)
    end

    Gradus.ProgressMeter.finish!(bar)

    data .= output
end
