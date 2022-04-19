function ensure_domain(m::AbstractMetricParams{T}) where {T}
    min_r = inner_radius(m)
    (u, λ, integrator) -> u[2] ≤ min_r * 1.1 || u[2] > 1200.0
end

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

function radial_negative_check(m::CarterMethodBL{T}) where {T}
    (u, λ, integrator) -> begin
        Vr(m, u, integrator.p) < 0
    end
end

function angular_negative_check(m::CarterMethodBL{T}) where {T}
    (u, λ, integrator) -> begin
        Vθ(m, u, integrator.p) < 0
    end
end
