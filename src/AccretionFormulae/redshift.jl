module RedshiftFunctions
import ..AccretionFormulae.Gradus
import ..AccretionFormulae.Gradus: __BoyerLindquistFO
import ..AccretionFormulae: AbstractMetricParams, metric
using StaticArrays
using Tullio: @tullio

"""
    eⱽ(M, r, a, θ)
Modified from Cunningham et al. (1975) eq. (A2a):
```math
e^\\nu = \\sqrt{\\frac{\\Delta \\Sigma}{A}}.
```
"""
eⱽ(M, r, a, θ) = √(
    __BoyerLindquistFO.Σ(r, a, θ) * __BoyerLindquistFO.Δ(M, r, a) /
    __BoyerLindquistFO.A(M, r, a, θ),
)


"""
    eᶲ(M, r, a, θ) 
Modified from Cunningham et al. (1975) eq. (A2b):
```math
e^\\Phi = \\sin \\theta \\sqrt{\\frac{A}{\\Sigma}}.
```
"""
eᶲ(M, r, a, θ) =
    sin(θ) * √(__BoyerLindquistFO.A(M, r, a, θ) / __BoyerLindquistFO.Σ(r, a, θ))


"""
    ω(M, r, a, θ)
From Cunningham et al. (1975) eq. (A2c):
```math
\\omega = \\frac{2 a M r}{A}.
```
"""
ω(M, r, a, θ) = 2 * a * M * r / __BoyerLindquistFO.A(M, r, a, θ)


"""
    Ωₑ(M, r, a)
Coordinate angular velocity of an accreting gas. 
Taken from Cunningham et al. (1975) eq. (A7b):
```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} + r_e^{3/2}}.
```
# Notes
Fanton et al. (1997) use
```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} \\pm r_e^{3/2}},
```
where the sign is dependent on co- or contra-rotation. This function may be extended in the future to support this definition.
"""
Ωₑ(M, r, a) = √M / (r^1.5 + a * √M)


"""
    Vₑ(M, r, a, θ)  
Velocity of an accreting gas in a locally non-rotating reference frame (see Bardeen et al. 1973).
Taken from Cunningham et al. (1975) eq. (A7b):
```math
V_e = (\\Omega_e - \\omega) e^{\\Phi - \\nu}.
```
"""
Vₑ(M, r, a, θ) = (Ωₑ(M, r, a) - ω(M, r, a, θ)) * eᶲ(M, r, a, θ) / eⱽ(M, r, a, θ)


"""
    Lₑ(M, rms, a) 
Angular momentum of an accreting gas within ``r_ms``.
Taken from Cunningham et al. (1975) eq. (A11b):
```math
L_e = \\sqrt{M} \\frac{
        r_{\\text{ms}}^2 - 2 a \\sqrt{M r_{\\text{ms}}} + a^2
    }{
        r_{\\text{ms}}^{3/2} - 2 M \\sqrt{r_{\\text{ms}}} + a \\sqrt{M}
    }.
```
"""
Lₑ(M, rms, a) = √M * (rms^2 - 2 * a * √(M * rms) + a^2) / (rms^1.5 - 2 * M * √rms + a * √M)


"""
    H(M, rms, r, a)
Taken from Cunningham et al. (1975) eq. (A12e):
```math
H = \\frac{2 M r_e - a \\lambda_e}{\\Delta},
```
where we distinguing ``r_e`` as the position of the accreting gas. 
"""
H(M, rms, r, a) = (2 * M * r - a * Lₑ(M, rms, a)) / __BoyerLindquistFO.Δ(M, r, a)


"""
    γₑ(M, rms) 
Taken from Cunningham et al. (1975) eq. (A11c):
```math
\\gamma_e = \\sqrt{1 - \\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }}.
```
"""
γₑ(M, rms) = √(1 - (2 * M) / (3 * rms))


"""
    uʳ(M, rms, r)
Taken from Cunningham et al. (1975) eq. (A12b):
```math
u^r = - \\sqrt{\\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }} \\left(
        \\frac{ r_{\\text{ms}} }{r_e} - 1
    \\right)^{3/2}.
```
"""
uʳ(M, rms, r) = -√((2 * M) / (3 * rms)) * (rms / r - 1)^1.5


"""
    uᶲ(M, rms, r, a)
Taken from Cunningham et al. (1975) eq. (A12c):
```math
u^\\phi = \\frac{\\gamma_e}{r_e^2} \\left( 
        L_e + aH 
    \\right).
```
"""
uᶲ(M, rms, r, a) = γₑ(M, rms) / r^2 * (Lₑ(M, rms, a) + a * H(M, rms, r, a))


"""
    uᵗ(M, rms, r, a) 
Taken from Cunningham et al. (1975) eq. (A12b):
```math
u^t = \\gamma_e \\left(
        1 + \\frac{2 M (1 + H)}{r_e}
    \\right).
```
"""
uᵗ(M, rms, r, a) = γₑ(M, rms) * (1 + 2 * M * (1 + H(M, rms, r, a)) / r)


"""
Experimental API:
"""
function regular_pdotu_inv(L, M, r, a, θ)
    (eⱽ(M, r, a, θ) * √(1 - Vₑ(M, r, a, θ)^2)) / (1 - L * Ωₑ(M, r, a))
end

@inline function regular_pdotu_inv(m::AbstractMetricParams{T}, u, v) where {T}
    metric_matrix = metric(m, u)
    # reverse signs of the velocity vector
    # since we're integrating backwards
    p = @inbounds @SVector [-v[1], 0, 0, -v[4]]

    # TODO: this only works for Kerr
    disc_norm = (eⱽ(m.M, u[2], m.a, u[3]) * √(1 - Vₑ(m.M, u[2], m.a, u[3])^2))

    u_disc = @SVector [1 / disc_norm, 0, 0, Ωₑ(m.M, u[2], m.a) / disc_norm]

    # use Tullio to do the einsum
    @tullio g := metric_matrix[i, j] * (u_disc[i]) * p[j]
    1 / g
end

function plunging_p_dot_u(E, a, M, L, Q, rms, r, sign_r)
    inv(
        uᵗ(M, rms, r, a) - uᶲ(M, rms, r, a) * L -
        sign_r * uʳ(M, rms, r) * __BoyerLindquistFO.Σδr_δλ(E, L, M, Q, r, a) /
        __BoyerLindquistFO.Δ(M, r, a),
    )
end

function plunging_p_dot_u(m::AbstractMetricParams{T}, u, v, rms) where {T}
    metric_matrix = metric(m, u)

    # reverse signs of the velocity vector
    # since we're integrating backwards
    p = @inbounds @SVector [-v[1], v[2], 0, -v[4]]

    u_disc = @inbounds @SVector [
        uᵗ(m.M, rms, u[2], m.a),
        uʳ(m.M, rms, u[2]),
        0,
        uᶲ(m.M, rms, u[2], m.a),
    ]

    @tullio g := metric_matrix[i, j] * u_disc[i] * p[j]
    1 / g
end

@inline function redshift_function(m::Gradus.BoyerLindquistFO{T}, u, p, λ) where {T}
    isco = Gradus.isco(m)
    if u[2] > isco
        @inbounds regular_pdotu_inv(p.L, m.M, u[2], m.a, u[3])
    else
        # change sign if we're after the sign flip
        # TODO: this isn't robust to multiple sign changes 
        # TODO: i feel like there should be a better way than checking this with two conditions
        #       used to have λ > p.changes[1] to make sure we're ahead of the time flip (but when wouldn't we be???)
        #       now p.changes[1] > 0.0 to make sure there was a time flip at all
        sign_r = (p.changes[1] > 0.0 ? 1 : -1) * p.r
        @inbounds plunging_p_dot_u(m.E, m.a, m.M, p.L, p.Q, isco, u[2], sign_r)
    end
end

@inline function redshift_function(m::AbstractMetricParams{T}, u, v) where {T}
    isco = Gradus.isco(m)
    if u[2] > isco
        regular_pdotu_inv(m, u, v)
    else
        plunging_p_dot_u(m, u, v, isco)
    end
end

end # module

# point functions exports
function _redshift_guard(
    m::Gradus.BoyerLindquistFO{T},
    gp::FirstOrderGeodesicPoint{T},
    max_time,
) where {T}
    RedshiftFunctions.redshift_function(m, gp.u, gp.p, gp.t)
end
function _redshift_guard(m::AbstractMetricParams{T}, gp, max_time) where {T}
    RedshiftFunctions.redshift_function(m, gp.u, gp.v)
end
