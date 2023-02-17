"""
    $(TYPEDSIGNATURES)

Calculates initial four-velocity vectors from impact parameters ``\\alpha`` and
``\\beta``, with ``v^t = 0``, i.e. unconstrained, given some ``init_pos`` four-vector.

`α_generator` and `β_generator` should be any iterable of ``\\alpha`` and ``\\beta``
respectively, whereas ``\\beta`` is a single ``T`` value, such that a variety of geodesics
in the same plane may be easily traced.
"""
function calculate_velocities(m::AbstractMetricParams, init_pos, α_generator, β::Number)
    [map_impact_parameters(m, init_pos, α, β) for α in α_generator]
end

function calculate_velocities(m::AbstractMetricParams, init_pos, α_genetator, β_generator)
    [map_impact_parameters(m, init_pos, α, β) for α in α_genetator, β in β_generator]
end

"""
In-place specialisation, writing the four-velocities into `vs`.
"""
function calculate_velocities!(
    vs,
    m::AbstractMetricParams,
    init_pos,
    α_generator,
    β::Number,
)
    for (i, α) in enumerate(α_generator)
        vs[i] = map_impact_parameters(m, init_pos, α, β)
    end
end


"""
    $(TYPEDSIGNATURES)

Utility function for converting some `X` on an image plane into ``\\alpha``, given
the midpoint `x_mid` and field-of-view multiplier `fov`.
"""
# have to use a slight 0.001 offset to avoid integrating α=0.0 geodesics in first order methods
x_to_α(X, x_mid, fov) = (X + 1e-3 - x_mid) / fov

"""
    $(TYPEDSIGNATURES)

Utility function for converting some `Y` on an image plane into ``\\beta``, given
the midpoint `y_mid` and field-of-view multiplier `fov`.
"""
y_to_β(Y, y_mid, fov) = (Y - y_mid) / fov

function init_progress_bar(text, N, enabled)
    ProgressMeter.Progress(
        N;
        desc = text,
        dt = 0.5,
        barglyphs = BarGlyphs("[=> ]"),
        barlen = 40,
        color = :none,
        enabled = enabled,
        showspeed = true,
    )
end

function impact_axes(width, height, fov)
    x_mid = width ÷ 2
    y_mid = height ÷ 2
    α = [Gradus.x_to_α(X, x_mid, fov) for X = 1:width]
    β = [Gradus.y_to_β(Y, y_mid, fov) for Y = 1:height]
    α, β
end
