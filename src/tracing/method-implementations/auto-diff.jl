"""
    AbstractStaticAxisSymmetric{T}

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
AbstractStaticAxisSymmetric

"""
    $(TYPEDSIGNATURES)

Interface for [`AbstractStaticAxisSymmetric`](@ref). Should return
a vector or tuple with the elements
```math
\\left(
    g_{tt}, g_{rr}, g_{\\theta \\theta}, g_{\\phi \\phi}, g_{t\\phi}
\\right).
```
"""
metric_components(m::AbstractStaticAxisSymmetric, rθ) =
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

        Δ = inv(g1 * g2 * g3 * g4 - (g5^2) * g2 * g3)
        SVector{5}(
            (g2 * g3 * g4) * Δ,
            (g1 * g3 * g4 - (g5^2) * g3) * Δ,
            (g1 * g2 * g4 - (g5^2) * g2) * Δ,
            (g1 * g2 * g3) * Δ,
            (-g2 * g3 * g5) * Δ,
        )
    end
end
inverse_metric_components(m::AbstractStaticAxisSymmetric, rθ) =
    inverse_metric_components(metric_components(m, rθ))

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
= \\frac{1}{2} g^{\\mu\\rho} \\left(
    \\partial_{\\nu}g_{\\rho \\sigma}
    + \\partial_{\\sigma}g_{\\rho \\nu}
    - \\partial_{\\rho}g_{\\sigma \\nu}
\\right).
```

Limitations:
- currenly pre-supposes static, axis-symmetric metric.
"""
@generated function compute_geodesic_equation(
    _ginv::SVector,
    _j1::SVector,
    _j2::SVector,
    _v::SVector,
)
    Symbolics.@variables ginv[1:5], j1[1:5], j2[1:5], v[1:4] # non zero metric components
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
    j0 = zeros(Symbolics.Num, (4, 4))
    jacobian = (j0, j1_mat, j2_mat, j0)
    # christoffel symbols
    @tullio Γ[i, k, l] :=
        1 / 2 *
        inverse_metric[i, m] *
        (jacobian[l][m, k] + jacobian[k][m, l] - jacobian[m][k, l])
    quote
        @inbounds @muladd @fastmath let ginv = _ginv, j1 = _j1, j2 = _j2, v = _v
            Γ1 = SMatrix{4,4}($(Symbolics.toexpr.(Γ[1, :, :])...))
            Γ2 = SMatrix{4,4}($(Symbolics.toexpr.(Γ[2, :, :])...))
            Γ3 = SMatrix{4,4}($(Symbolics.toexpr.(Γ[3, :, :])...))
            Γ4 = SMatrix{4,4}($(Symbolics.toexpr.(Γ[4, :, :])...))

            -SVector{4}((Γ1 * v) ⋅ v, (Γ2 * v) ⋅ v, (Γ3 * v) ⋅ v, (Γ4 * v) ⋅ v)
        end
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

@inline function constrain(m::AbstractStaticAxisSymmetric{T}, u, v; μ::T = T(0.0)) where {T}
    rθ = (u[2], u[3])
    g_comps = metric_components(m, rθ)
    constrain_time(g_comps, v, μ)
end

const _static_dual_eval =
    Base.get_extension(ForwardDiff, :ForwardDiffStaticArraysExt).static_dual_eval
const _extract_jacobian =
    Base.get_extension(ForwardDiff, :ForwardDiffStaticArraysExt).extract_jacobian

"""
    metric_jacobian(m::AbstractStaticAxisSymmetric{T}, rθ)

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
function metric_jacobian(m::AbstractStaticAxisSymmetric, rθ)
    f = x -> metric_components(m, x)
    T = typeof(ForwardDiff.Tag(f, eltype(rθ)))
    ydual = _static_dual_eval(T, f, rθ)
    (ForwardDiff.value.(T, ydual), _extract_jacobian(T, ydual, rθ))
end

@inbounds function geodesic_equation(
    m::AbstractStaticAxisSymmetric,
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

function metric(m::AbstractStaticAxisSymmetric, u)
    rθ = (u[2], u[3])
    comps = metric_components(m, rθ)
    @SMatrix [
        comps[1] 0 0 comps[5]
        0 comps[2] 0 0
        0 0 comps[3] 0
        comps[5] 0 0 comps[4]
    ]
end

export AbstractStaticAxisSymmetric
