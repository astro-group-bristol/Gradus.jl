
@testset "tetradframe" begin
    all_metrics = (
        # can't do first order yet since no four velocity
        # BoyerLindquistFO(1.0, 0.998, 1.0),
        # BoyerLindquistFO(1.0, -0.998, 1.0),
        BoyerLindquistAD(1.0, 0.998),
        BoyerLindquistAD(1.0, -0.998),
        JohannsenAD(M = 1.0, a = 0.998, α13 = 1.0),
        JohannsenAD(M = 1.0, a = 0.998, α22 = 1.0),
        DilatonAxionAD(M = 1.0, a = 0.998, β = 0.2, b = 1.0),
    )
    radii = 5.0:0.8:10.0
    angles = 0.1:0.5:2π
    minkowski = @SMatrix [
        -1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    for m in all_metrics, r in radii, θ in angles
        u = @SVector([0.0, r, θ, 0.0])
        v = Gradus.constrain_all(m, u, CircularOrbits.fourvelocity(m, r), 1.0)

        # function that we are testing
        M = GradusBase.tetradframe(m, u, v)

        m_mat = Gradus.metric(m, u)
        @tullio res[a, b] := m_mat[i, j] * M[a][i] * M[b][j]
        # ensure it gives minkowski
        @test isapprox(res, minkowski, atol = 1e-13)
    end
end
