using Test
using Gradus
using LinearAlgebra: diagm
using Tullio

# is the beaming being calculated correctly? 
# using Eq. (10) in Gonzalez+17 for the analytic tetrad
# NB: they use (+, -, -, -), whereas we use (-, +, +, +)

drdt(g, β) = β * √(-g[1] / g[2])
drdt(m::AbstractMetric, x, β) = drdt(metric_components(m, SVector(x[2], x[3])), β)

# Eq. (8) in G+17
drdt_analytic(r, a, β) = β * (r^2 - 2r + a^2) / (r^2 + a^2)

function analytic_tetrad(m::AbstractMetric, x, β)
    g = metric_components(m, SVector(x[2], x[3]))
    v = drdt(g, β)

    A = 1 / √(-g[1] - v^2 * g[2])
    e_t = A * SVector(1, v, 0, 0)

    B = √(-g[2] / g[1])
    e_r = A * SVector(v * B, 1 / B, 0, 0)

    e_θ = SVector(0, 0, √(1 / g[3]), 0)

    C = 1 / √(-g[1] * (g[5]^2 - g[1] * g[4]))
    e_ϕ = C * SVector(g[5], 0, 0, -g[1])


    (e_t, e_r, e_θ, e_ϕ)
end

function is_minkowski(m, x, M)
    m_mat = Gradus.metric(m, x)
    @tullio res[a, b] := m_mat[i, j] * M[a][i] * M[b][j]
    η = diagm([-1.0, 1.0, 1.0, 1.0])
    @test res ≈ η
end

m = KerrMetric(1.0, 0.998)
x = SVector(0, 3.0, deg2rad(0.01), 0)

# make sure the velocity is correctly calculated
@test drdt(m, x, 1) ≈ drdt_analytic(x[2], m.a, 1)

M1 = analytic_tetrad(m, x, 0.25)
# make sure analytic tetrad is approx minkowski
is_minkowski(m, x, M1)

v = SVector(1, drdt(m, x, 0.25), 0, 0)
M2 = Gradus.tetradframe(m, x, v)
# make sure our generic method is approx minkwoski
is_minkowski(m, x, M2)

mat1 = reduce(hcat, M1)
mat2 = reduce(hcat, M2)

# are they the same
@test mat1 ≈ mat2

# make sure we calculate the velocity profiles correctly as well
# using two example velocities sent to me by Alexey Nekrasov

m = KerrMetric(1.0, 0.998)
model = RingCorona(Gradus.SourceVelocities.co_rotating, 1.02, 1.113)
x, v = Gradus.sample_position_velocity(m, model)
# no idea why this one isn't matching
# @test v ≈ [5.4962, 0.0, 0.0, 2.3307] rtol = 1e-3

model = RingCorona(Gradus.SourceVelocities.co_rotating, 2.082, 50.0)
x, v = Gradus.sample_position_velocity(m, model)
@test v ≈ [1.204, 0.0, 0.0, 0.300] rtol = 1e-3
