# for use with static, axis-symmetric metrics
# create new abstract type for easy re-definition
abstract type AbstractAutoDiffMetricParams{T} <: AbstractMetricParams{T} end
abstract type AbstractAutoDiffStaticAxisSymmetricParams{T} <:
              AbstractAutoDiffMetricParams{T} end

# interface
metric_components(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, rθ) where {T} =
    error("Not implemented for $(typeof(m)).")

"""
    inverse_metric_components(g_comp)

Calculates ``g^{tt}``, ``g^{rr}``, ``g^{\\theta\\theta}``, ``g^{\\phi\\phi}``, ``g^{t\\phi}`` of a static,
axis-symmetric metric from ``g_{tt}``, ``g_{rr}``, ``g_{\\theta\\theta}``, ``g_{\\phi\\phi}``, ``g_{t\\phi}``
using a symbolically computed inverse matrix method.

## Developer notes

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
@fastmath function inverse_metric_components(g_comp)
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
    calc_symb_geodesic_equation(g1inv, j1, j2, v)

Symbolically pre-computed geodesic equations for a static, axis-symmetric metric.

## Developer notes

To recreate:

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
@fastmath function compute_geodesic_equation(ginv, j1, j2, v)
    @inbounds let gi1 = ginv[1], gi2 = ginv[2], gi3 = ginv[3], gi4 = ginv[4], gi5 = ginv[5],
            j11 = j1[1], j12 = j1[2], j13 = j1[3], j14 = j1[4], j15 = j1[5],
            j21 = j2[1], j22 = j2[2], j23 = j2[3], j24 = j2[4], j25 = j2[5],
            v1 = v[1], v2 = v[2], v3 = v[3], v4 = v[4]
        (
            -2(0.5gi5 * j14 + 0.5gi1 * j15) * v2 * v4 -
            2(0.5gi1 * j11 + 0.5gi5 * j15) * v1 * v2 -
            2(0.5gi1 * j21 + 0.5gi5 * j25) * v1 * v3 -
            2(0.5gi5 * j24 + 0.5gi1 * j25) * v3 * v4,
            0.5(v1^2) * gi2 * j11 +
            0.5(v3^2) * gi2 * j13 +
            0.5(v4^2) * gi2 * j14 +
            gi2 * j15 * v1 * v4 - 0.5(v2^2) * gi2 * j12 -
            gi2 * j22 * v2 * v3,
            0.5(v1^2) * gi3 * j21 +
            0.5(v2^2) * gi3 * j22 +
            0.5(v4^2) * gi3 * j24 +
            gi3 * j25 * v1 * v4 - 0.5(v3^2) * gi3 * j23 -
            gi3 * j13 * v2 * v3,
            -2(0.5gi4 * j14 + 0.5gi5 * j15) * v2 * v4 -
            2(0.5gi4 * j24 + 0.5gi5 * j25) * v3 * v4 -
            2(0.5gi5 * j11 + 0.5gi4 * j15) * v1 * v2 -
            2(0.5gi4 * j25 + 0.5gi5 * j21) * v1 * v3,
        )
    end
end

@fastmath function constrain_time(g_comp, v, μ = 0.0, positive::Bool = true)
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

@inbounds function geodesic_eq(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    u,
    v,
) where {T}
    # get the only position components we need for this metric type
    rθ = @SVector [u[2], u[3]]
    # calculate all non-zero components
    g_comps = metric_components(m, rθ)
    # use AD to get derivatives
    jacs = ForwardDiff.vector_mode_jacobian(x -> metric_components(m, x), rθ)
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
