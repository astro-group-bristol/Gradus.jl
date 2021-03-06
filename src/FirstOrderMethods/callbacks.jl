function flip_radial_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.r = -p.r
    integrator.sol.prob.p.changes[1] = integrator.t[end]
end

function flip_angular_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.θ = -p.θ
    integrator.sol.prob.p.changes[2] = integrator.t[end]
end

function radial_negative_check(m::AbstractFirstOrderMetricParams{T}) where {T}
    (u, λ, integrator) -> Vr(m, u, integrator.p) < 0
end

function angular_negative_check(m::AbstractFirstOrderMetricParams{T}) where {T}
    (u, λ, integrator) -> Vθ(m, u, integrator.p) < 0
end
