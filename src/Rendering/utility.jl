function calculate_velocities(
    m::AbstractMetricParams{T},
    init_pos,
    α_generator,
    β::T,
) where {T}
    [map_impact_parameters(m, init_pos, α, β) for α in α_generator]
end

function calculate_velocities(
    m::AbstractMetricParams{T},
    init_pos,
    α_genetator,
    β_generator,
) where {T}
    [map_impact_parameters(m, init_pos, α, β) for α in α_genetator, β in β_generator]
end

"""
In-place specialisation
"""
function calculate_velocities!(
    vs,
    m::AbstractMetricParams{T},
    init_pos,
    α_generator,
    β::T,
) where {T}
    for (i, α) in enumerate(α_generator)
        vs[i] = map_impact_parameters(m, init_pos, α, β)
    end
end

# have to use a slight 0.001 offset to avoid integrating α=0.0 geodesics
x_to_α(X, x_mid, fov_factor) = (X + 0.001 - x_mid) / fov_factor
y_to_β(Y, y_mid, fov_factor) = (Y - y_mid) / fov_factor
