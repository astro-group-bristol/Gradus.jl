# for use with static, axis-symmetric metrics
# create new abstract type for easy re-definition
"""
    AbstractAutoDiffMetricParams{T} <: AbstractMetricParams{T}

Abstract type for metrics using the 2nd-order integration method, with the automatic
differentiation backend.
"""
abstract type AbstractAutoDiffMetricParams{T} <: AbstractMetricParams{T} end

"""
    AbstractAutoDiffStaticAxisSymmetricParams{T} <: AbstractAutoDiffMetricParams{T}

Specialisation for static, axis-symmetric metrics. Here, the metric is of the form
```math
    g_{\\mu\\nu} =
    \\left( \\begin{matrix}
        g_{tt}     & 0      & 0                  & g_{t\\phi}     \\\\
        0          & g_{rr} & 0                  & 0              \\\\
        0          & 0      & g_{\\theta\\theta} & 0              \\\\
        g_{t\\phi} & 0      & 0                  & g_{\\phi\\phi}
    \\end{matrix} \\right),
```
where the only non-zero off axis elements are ``g_{t\\phi}``.

Required implementations:
- [`inner_radius`](@ref)
- [`metric_components`](@ref)
"""
abstract type AbstractAutoDiffStaticAxisSymmetricParams{T} <:
              AbstractAutoDiffMetricParams{T} end

"""
    $(TYPEDSIGNATURES)

Interface for [`AbstractAutoDiffStaticAxisSymmetricParams`](@ref). Should return
a vector or tuple with the elements
```math
\\left(
    g_{tt}, g_{rr}, g_{\\theta \\theta}, g_{\\phi \\phi}, g_{t\\phi}
\\right).
```
"""
metric_components(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, rθ) where {T} =
    error("Not implemented for $(typeof(m)).")

"""
    inverse_metric_components(g_comp)

Calculates ``g^{tt}``, ``g^{rr}``, ``g^{\\theta\\theta}``, ``g^{\\phi\\phi}``, ``g^{t\\phi}`` of a static,
axis-symmetric metric from ``g_{tt}``, ``g_{rr}``, ``g_{\\theta\\theta}``, ``g_{\\phi\\phi}``, ``g_{t\\phi}``
using a symbolically computed inverse matrix method.

# Notes

To recreate:

```julia
using Symbolics
@variables g[1:5] # non zero metric components
metric = [
    g[1] 0 0 g[5]
    0 g[2] 0 0
    0 0 g[3] 0
    g[5] 0 0 g[4]
]
inv(metric)
```
"""
@muladd @fastmath function inverse_metric_components(g_comp)
    @inbounds let g1 = g_comp[1],
        g2 = g_comp[2],
        g3 = g_comp[3],
        g4 = g_comp[4],
        g5 = g_comp[5]

        (
            (g2 * g3 * g4) / (g1 * g2 * g3 * g4 - (g5^2) * g2 * g3),
            (g1 * g3 * g4 - (g5^2) * g3) / (g1 * g2 * g3 * g4 - (g5^2) * g2 * g3),
            (g1 * g2 * g4 - (g5^2) * g2) / (g1 * g2 * g3 * g4 - (g5^2) * g2 * g3),
            (g1 * g2 * g3) / (g1 * g2 * g3 * g4 - (g5^2) * g2 * g3),
            (-g2 * g3 * g5) / (g1 * g2 * g3 * g4 - (g5^2) * g2 * g3),
        )
    end
end


"""
    $(TYPEDSIGNATURES)

Using the inverse metric `ginv`, the Jacobian of the metric for ``r`` and ``\\theta``,
`j1` and `j2` respectively, and velocity four-vector `v`, calculates the four-acceleration
via the geodesic equation.

Returns the components of ``\\frac{\\text{d}^2 u^\\mu}{\\text{d} \\lambda^2}`` via

```math
\\frac{\\text{d}^2 u^\\mu}{\\text{d} \\lambda^2}
    + \\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}
    \\frac{\\text{d}u^\\nu}{\\text{d} \\lambda}
    \\frac{\\text{d}u^\\sigma}{\\text{d} \\lambda}
= 0,
```

where ``x^\\mu`` is a position four-vector, ``\\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}``
are the Christoffel symbols of the second kind, and ``\\lambda`` the affine parameter
describing the curve.

The Christoffel symbols ``\\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}`` are defined

```math
\\Gamma^{\\mu}_{\\phantom{\\mu}\\nu\\sigma}
= \\frac{1}{2} g^{\\mu, \\rho} \\left(
    \\partial_{\\nu}g_{\\rho \\sigma}
    + \\partial_{\\sigma}g_{\\rho \\nu}
    - \\partial_{\\rho}g_{\\sigma \\nu}
\\right).
```

Limitations:
- currenly pre-supposes static, axis-symmetric metric.

# Notes
This function is symbolically pre-computed using the following code:

```julia
using Symbolics, Tullio
@variables ginv[1:5], j1[1:5], j2[1:5], v[1:4] # non zero metric components
inverse_metric = [
    ginv[1] 0 0 ginv[5]
    0 ginv[2] 0 0
    0 0 ginv[3] 0
    ginv[5] 0 0 ginv[4]
]
j1_mat = [
    j1[1] 0 0 j1[5]
    0 j1[2] 0 0
    0 0 j1[3] 0
    j1[5] 0 0 j1[4]
]
j2_mat = [
    j2[1] 0 0 j2[5]
    0 j2[2] 0 0
    0 0 j2[3] 0
    j2[5] 0 0 j2[4]
]
j0 = zeros(Float64, (4, 4))
jacobian = (j0, j1_mat, j2_mat, j0)
# christoffel symbols
@tullio Γ[i, k, l] :=
    1 / 2 *
    inverse_metric[i, m] *
    (jacobian[l][m, k] + jacobian[k][m, l] - jacobian[m][k, l])
# compute geodesic equation
@tullio δxδλ[i] := -v[j] * Γ[i, j, k] * v[k]
```
"""
@muladd @fastmath function compute_geodesic_equation(ginv, j1, j2, v)
    @inbounds let gi1 = ginv[1],
        gi2 = ginv[2],
        gi3 = ginv[3],
        gi4 = ginv[4],
        gi5 = ginv[5],
        j11 = j1[1],
        j12 = j1[2],
        j13 = j1[3],
        j14 = j1[4],
        j15 = j1[5],
        j21 = j2[1],
        j22 = j2[2],
        j23 = j2[3],
        j24 = j2[4],
        j25 = j2[5],
        v1 = v[1],
        v2 = v[2],
        v3 = v[3],
        v4 = v[4]

        (
            -2(0.5gi5 * j14 + 0.5gi1 * j15) * v2 * v4 -
            2(0.5gi1 * j11 + 0.5gi5 * j15) * v1 * v2 -
            2(0.5gi1 * j21 + 0.5gi5 * j25) * v1 * v3 -
            2(0.5gi5 * j24 + 0.5gi1 * j25) * v3 * v4,
            0.5(v1^2) * gi2 * j11 +
            0.5(v3^2) * gi2 * j13 +
            0.5(v4^2) * gi2 * j14 +
            gi2 * j15 * v1 * v4 - 0.5(v2^2) * gi2 * j12 - gi2 * j22 * v2 * v3,
            0.5(v1^2) * gi3 * j21 +
            0.5(v2^2) * gi3 * j22 +
            0.5(v4^2) * gi3 * j24 +
            gi3 * j25 * v1 * v4 - 0.5(v3^2) * gi3 * j23 - gi3 * j13 * v2 * v3,
            -2(0.5gi4 * j14 + 0.5gi5 * j15) * v2 * v4 -
            2(0.5gi4 * j24 + 0.5gi5 * j25) * v3 * v4 -
            2(0.5gi5 * j11 + 0.5gi4 * j15) * v1 * v2 -
            2(0.5gi4 * j25 + 0.5gi5 * j21) * v1 * v3,
        )
    end
end

"""
    $(FUNCTIONNAME)(g_comp, v, μ = 0.0, positive::Bool = true)

Constrains the time component of the four-velocity `v`, given metric components `g_comp` and
effective mass `μ`.

```math
g_{\\sigma\\nu} \\dot{u}^\\sigma \\dot{u}^\\nu = -\\mu^2,
```

for ``v^t``. The argument `positive` allows the sign of ``\\mu`` to be changed. `true`
corresponds to time-like geodesics, `false` to space-like.

This function should rarely be directly called, and instead is invoked by [`constrain`](@ref).

Limitations:
- currenly pre-supposes static, axis-symmetric metric.
"""
@muladd @fastmath function constrain_time(g_comp, v, μ = 0.0, positive::Bool = true)
    @inbounds begin
        discriminant = (
            -g_comp[1] * g_comp[2] * v[2]^2 - g_comp[1] * g_comp[3] * v[3]^2 -
            g_comp[1] * μ^2 - (g_comp[1] * g_comp[4] - g_comp[5]^2) * v[4]^2
        )
        if positive
            -(g_comp[5] * v[4] + √discriminant) / g_comp[1]
        else
            -(g_comp[5] * v[4] - √discriminant) / g_comp[1]
        end
    end
end

function constrain(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    u,
    v;
    μ::T = 0.0,
) where {T}
    rθ = (u[2], u[3])
    g_comps = metric_components(m, rθ)
    constrain_time(g_comps, v, μ)
end

"""
    metric_jacobian(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, rθ)

Calculate the value and Jacobian elements of the metric with respect to ``r`` and ``\\theta``.

Limitations:
- currenly pre-supposes static, axis-symmetric metric.

## Notes

Function body is equivalent to
```julia
f = x -> metric_components(m, x)
J = ForwardDiff.vector_mode_jacobian(f, rθ)
f(rθ), J
```
but non-allocating.
"""
function metric_jacobian(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ)
    f = x -> metric_components(m, x)
    T = typeof(ForwardDiff.Tag(f, eltype(rθ)))
    ydual = ForwardDiff.static_dual_eval(T, f, rθ)
    ForwardDiff.value.(T, ydual), ForwardDiff.extract_jacobian(T, ydual, rθ)
end

@inbounds function geodesic_eq(
    m::AbstractAutoDiffStaticAxisSymmetricParams,
    u::AbstractArray{T},
    v::AbstractArray{T},
) where {T}
    # get the only position components we need for this metric type
    rθ = SVector{2,T}(u[2], u[3])
    # calculate all non-zero components, and use AD to get derivatives
    g_comps, jacs = metric_jacobian(m, rθ)
    # calculate all non-zero inverse matric components
    ginv_comps = inverse_metric_components(g_comps)
    # calculate acceleration
    compute_geodesic_equation(ginv_comps, jacs[:, 1], jacs[:, 2], v)
end


function metric(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, u) where {T}
    rθ = (u[2], u[3])
    comps = metric_components(m, rθ)
    @SMatrix [
        comps[1] 0 0 comps[5]
        0 comps[2] 0 0
        0 0 comps[3] 0
        comps[5] 0 0 comps[4]
    ]
end

export AbstractAutoDiffMetricParams, AbstractAutoDiffStaticAxisSymmetricParams
