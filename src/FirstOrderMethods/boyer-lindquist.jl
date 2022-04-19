@with_kw struct BoyerLindquistFO{T} <: AbstractFirstOrderMetricParams{T}
    @deftype T
    M = 1.0
    a = 0.0
    E = 1.0
end

function four_velocity(u, E, M, a, p)
    let r = u[2], θ = u[3], L = p.L, Q = p.Q
        Σ₀ = BoyerLindquistFOCoords.Σ(r, a, θ)
        (
            BoyerLindquistFOCoords.Σδt_δλ(E, L, M, r, a, θ) / Σ₀,
            p.r * BoyerLindquistFOCoords.Σδr_δλ(E, L, M, Q, r, a) / Σ₀,
            p.θ * BoyerLindquistFOCoords.Σδθ_δλ(E, L, Q, a, θ) / Σ₀,
            BoyerLindquistFOCoords.Σδϕ_δλ(E, L, M, r, a, θ) / Σ₀,
        )
    end
end

calc_lq(m::BoyerLindquistFO{T}, pos, vel) where {T} =
    BoyerLindquistFOCoords.LQ(m.M, pos[2], m.a, pos[3], vel[3], vel[4])


inner_radius(m::BoyerLindquistFO{T}) where {T} = m.M + √(m.M^2 - m.a^2)
constrain(::BoyerLindquistFO{T}, u, v; μ = 0.0) where {T} = v[1]


function integrator_problem(
    m::BoyerLindquistFO{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
) where {S,T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{false}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do u, p, λ
        SVector(four_velocity(u, m.E, m.M, m.a, p)...)
    end
end

function integrator_problem(
    m::BoyerLindquistFO{T},
    pos::AbstractVector{T},
    vel::AbstractVector{T},
    time_domain,
) where {T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{true}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do du, u, p, λ
        du .= four_velocity(u, m.E, m.M, m.a, p)
    end
end

function alpha_beta_to_vel(m::BoyerLindquistFO{T}, u, α, β) where {T}
    sinΦ, sinΨ = BoyerLindquistFOCoords.sinΦsinΨ(m.M, u[2], m.a, u[3], α, β)
    (β < 0.0 ? -1.0 : 1.0, sinΦ, sinΨ)
end


function Vr(m::BoyerLindquistFO{T}, u, p) where {T}
    L, Q, _, _ = p
    BoyerLindquistFOCoords.Vr(m.E, L, m.M, Q, u[2], m.a)
end
function Vθ(m::BoyerLindquistFO{T}, u, p) where {T}
    L, Q, _, _ = p
    BoyerLindquistFOCoords.Vθ(m.E, L, Q, m.a, u[3])
end
