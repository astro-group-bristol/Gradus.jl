
module __NoZMetric
using ..StaticArrays
using ..MuladdMacro

@muladd @fastmath begin

    epsilon(M, a, ϵ, y) = ϵ * M * a * y

    # the way this function must be defined is a little complex
    # but helps with type-stability
    function metric_components(M, a, ϵ, rθ)
        (r, θ) = rθ
        sinθ2 = sin(θ)^2

        y = cos(θ)
        _epsilon = epsilon(M, a, ϵ, y)

        tt =
            -1 + (
                (2M * r * (r^(2) + a^(2) * y^(2))) /
                ((r^(2) + a^(2) * y^(2))^(2) + (r^(2) - 2M * r + a^(2) * y^(2)) * _epsilon)
            )
        ϕϕ =
            (
                (1 - y^2) *
                (r^2 + a^2 * y^2 + _epsilon) *
                (
                    r^4 +
                    a^4 * y^2 +
                    r^2 * (a^2 + a^2 * y^2 + _epsilon) +
                    a^2 * _epsilon +
                    2M * r * (a^2 - a^2 * y^2 - _epsilon)
                )
            ) / ((r^2 + a^2 * y^2)^2 + (r^2 - 2M * r + a^2 * y^2) * _epsilon)
        rr = (r^(2) + a^(2) * y^(2) + _epsilon) / (r^(2) - 2M * r + a^(2))
        yy = (r^(2) + a^(2) * y^(2) + _epsilon) / (1 - y^(2))

        tϕ =
            -(2M * r * a * (1 - y^(2)) * (r^(2) + a^(2) * y^(2) + _epsilon)) /
            ((r^(2) + a^(2) * y^(2))^(2) + (r^(2) - 2M * r + a^(2) * y^(2)) * _epsilon)

        # include the coordinate transformation factor
        #   y = cos(θ)  -->  dy^2 = dθ^2 * sin(θ)^2
        @SVector [tt, rr, yy * sinθ2, ϕϕ, tϕ]
    end
end

end # module

# new structure for our spacetime
"""
    struct NoZMetric
"""
@with_kw struct NoZMetric{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "Deviation parameter"
    ϵ = 0.0
end

# implementation
metric_components(m::NoZMetric, rθ) = __NoZMetric.metric_components(m.M, m.a, m.ϵ, rθ)
inner_radius(m::NoZMetric) = m.M + √(m.M^2 - m.a^2)

function _solve_orbit_θ(m, r)
    function _objective(θ)
        rθ = SVector(r, θ)
        _, J = Gradus.metric_jacobian(m, rθ)
        ∂rg = J[:, 1]
        ∂θg = J[:, 2]
        Ωϕ = Gradus.CircularOrbits._Ω_analytic(∂rg, false)

        ∂θg[1] + 2 * ∂θg[5] * Ωϕ + ∂θg[4] * Ωϕ^2
    end

    Gradus.Roots.find_zero(_objective, π / 2)
end

function make_circular_velocity_function(
    m::NoZMetric{T};
    outer_radius = T(500.0),
    num_samples::Int = 200,
) where {T}
    isco = Gradus.isco(m)

    rs = collect(Grids._geometric_grid(isco, outer_radius, num_samples))
    θs = Gradus._solve_orbit_θ.(m, rs)
    interp = Gradus._make_interpolation(rs, θs)
    function _velocity_function(r)
        θ = interp(r)
        CircularOrbits.fourvelocity(m, SVector(r, θ))
    end
end

function isco(m::NoZMetric{T}) where {T}
    kerr_isco = isco(KerrMetric(M = m.M, a = m.a))
    rs = range(inner_radius(m), kerr_isco + 1.0, 100)
    thetas = map(rs) do r
        try
            _solve_orbit_θ(m, r)
        catch
            zero(T)
        end
    end

    # TODO: this is such a hack
    replace!(thetas, zero(T) => thetas[findfirst(!=(zero(T)), thetas)])

    interp = _make_interpolation(rs, thetas)

    dE(r) = ForwardDiff.derivative(x -> CircularOrbits.energy(m, SVector(x, interp(x))), r)
    # optimize from the Kerr equivalent metric
    Roots.find_zero(dE, kerr_isco)
end

export NoZMetric
