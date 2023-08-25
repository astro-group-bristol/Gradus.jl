"""
    calculate_velocities(m::AbstractMetric, init_pos, α_generator, β::Number)
    calculate_velocities(m::AbstractMetric, init_pos, α_genetator, β_generator)

Calculates initial four-velocity vectors from impact parameters ``\\alpha`` and
``\\beta``, with ``v^t = 0``, i.e. unconstrained, given some ``init_pos`` four-vector.

`α_generator` and `β_generator` should be any iterable of ``\\alpha`` and ``\\beta``
respectively, whereas ``\\beta`` is a single ``T`` value, such that a variety of geodesics
in the same plane may be easily traced.
"""
function calculate_velocities(m::AbstractMetric, init_pos, α_generator, β::Number)
    [map_impact_parameters(m, init_pos, α, β) for α in α_generator]
end
function calculate_velocities(m::AbstractMetric, init_pos, α_genetator, β_generator)
    [map_impact_parameters(m, init_pos, α, β) for α in α_genetator, β in β_generator]
end

"""
    calculate_velocities!(vs, m::AbstractMetric, init_pos, α_generator, β::Number)

In-place specialisation, writing the four-velocities into `vs`.
"""
function calculate_velocities!(vs, m::AbstractMetric, init_pos, α_generator, β::Number)
    for (i, α) in enumerate(α_generator)
        vs[i] = map_impact_parameters(m, init_pos, α, β)
    end
end

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

function impact_axes(width, height, αlims, βlims)
    α = [a for a in range(αlims..., width)]
    β = [b for b in range(βlims..., height)]
    α, β
end
