function integrator_problem(
    m::AbstractMetricParams{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
) where {S,T}
    u_init = vcat(pos, vel)
    ODEProblem{false}(u_init, time_domain) do u, p, λ
        @inbounds let x = @view(u[1:4]), v = @view(u[5:8])
            dv = SVector{4}(geodesic_eq(m, x, v))
            SVector{8}(v[1], v[2], v[3], v[4], dv[1], dv[2], dv[3], dv[4])
        end
    end
end

function integrator_problem(
    m::AbstractMetricParams{T},
    pos::AbstractVector{T},
    vel::AbstractVector{T},
    time_domain,
) where {T}
    u_init = vcat(pos, vel)
    ODEProblem{true}(u_init, time_domain) do du, u, p, λ
        @inbounds let x = @view(u[1:4]), v = @view(u[5:8])
            dv = SVector{4}(geodesic_eq(m, x, v))

            du[1:4] = v
            du[5:8] = dv
        end
    end
end
