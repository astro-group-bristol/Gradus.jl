function integrator_problem(
    m::AbstractMetricParams{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain
) where {S,T}
    SecondOrderODEProblem{false}(vel, pos, time_domain, m) do v, u, p, λ
        SVector(geodesic_eq(p, u, v)...)
    end
end

function integrator_problem(
    m::AbstractMetricParams{T},
    pos::AbstractVector{T},
    vel::AbstractVector{T},
    time_domain
) where {T}
    SecondOrderODEProblem{true}(vel, pos, time_domain, m) do dv, v, u, p, λ
        dv .= geodesic_eq(p, u, v)
    end
end
