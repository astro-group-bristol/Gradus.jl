@testset "tetradframe" begin
    all_metrics = (
        # can't do first order yet since no four velocity
        # BoyerLindquistFO(1.0, 0.998, 1.0),
        # BoyerLindquistFO(1.0, -0.998, 1.0),
        BoyerLindquistAD(1.0, 0.998),
        BoyerLindquistAD(1.0, -0.998),
        JohannsenAD(M = 1.0, a = 0.998, α13 = 1.0),
        JohannsenAD(M = 1.0, a = 0.998, α22 = 1.0),
        # DilatonAxionAD(M = 1.0, a = 0.998, β = 0.2, b = 1.0),
    )
    radii = 5.0:0.8:10.0
    angles = 0.1:0.5:π
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
        M = Gradus.GradusBase.tetradframe(m, u, v)

        m_mat = Gradus.metric(m, u)
        @tullio res[a, b] := m_mat[i, j] * M[a][i] * M[b][j]
        # ensure it gives minkowski
        @test isapprox(res, minkowski, atol = 1e-13)
    end

    @testset "lnrf" begin
        for m in all_metrics, r in radii, θ in angles
            u = @SVector([0.0, r, θ, 0.0])
    
            # function that we are testing
            M = Gradus.GradusBase.lnrframe(m, u)
    
            m_mat = Gradus.metric(m, u)
            @tullio res[a, b] := m_mat[i, j] * M[a][i] * M[b][j]
            # ensure it gives minkowski
            @test isapprox(res, minkowski, atol = 1e-13)
        end
    end

    @testset "kerr-lnrf" begin
        # these test explicity check the LNRF calculations for the known
        # theory of the Kerr metric

        function kerr_lnrframe(m, u)
            A = Gradus.__BoyerLindquistFO.A(m.M, u[2], m.a, u[3])
            Σ = Gradus.__BoyerLindquistFO.Σ(u[2], m.a, u[3])
            Δ = Gradus.__BoyerLindquistFO.Δ(m.M, u[2], m.a)
            ω = 2* m.M *m.a * u[2] / A
        
            et = √(A / (Σ * Δ)) * @SVector [1.0, 0.0, 0.0, ω]
            er = √(Δ / Σ) * @SVector [0.0, 1.0, 0.0, 0.0]
            eθ = √(1 / Σ) * @SVector [0.0, 0.0, 1.0, 0.0]
            eϕ = √(Σ / A) * (1/sin(u[3])) * @SVector [0.0, 0.0, 0.0, 1.0]

            vecs = (et, er, eθ, eϕ)
            reduce(hcat, vecs)
        end
        
        function numerical_lnrframe(m, u)
            vecs = Gradus.GradusBase.lnrframe(m, u)
            reduce(hcat, vecs)
        end

        for M in 0.2:0.8:2.0, a in -M:0.5:M
            m = BoyerLindquistAD(M, a)
            r = inner_radius(m) + 4.2
            for θ in angles
                u = @SVector [0.0, r, θ, 0.0] 
                expected = kerr_lnrframe(m, u)
                calculated = numerical_lnrframe(m, u)

                @test isapprox(calculated, expected, atol=1e-13)
            end
        end
    end

    @testset "kerr-lnrf-basis" begin
        function kerr_lnrbasis(m, u)
            A = Gradus.__BoyerLindquistFO.A(m.M, u[2], m.a, u[3])
            Σ = Gradus.__BoyerLindquistFO.Σ(u[2], m.a, u[3])
            Δ = Gradus.__BoyerLindquistFO.Δ(m.M, u[2], m.a)
            ω = 2* m.M *m.a * u[2] / A
        
            et = √(Σ * Δ / A) * @SVector [1.0, 0.0, 0.0, 0.0]
            er = √(Σ / Δ) * @SVector [0.0, 1.0, 0.0, 0.0]
            eθ = √Σ * @SVector [0.0, 0.0, 1.0, 0.0]
            eϕ = √(A / Σ) * sin(u[3]) * @SVector [-ω, 0.0, 0.0, 1.0]

            vecs = (et, er, eθ, eϕ)
            reduce(hcat, vecs)
        end
        
        function numerical_lnrbasis(m, u)
            vecs = Gradus.GradusBase.lnrbasis(m, u)
            reduce(hcat, vecs)
        end

        for M in 0.2:0.8:2.0, a in -M:0.5:M
            m = BoyerLindquistAD(M, a)
            r = inner_radius(m) + 0.3
            for θ in angles
                u = @SVector [0.0, r, θ, 0.1] 
                expected = kerr_lnrbasis(m, u)
                calculated = numerical_lnrbasis(m, u)

                @test isapprox(calculated, expected, atol=1e-10) 
            end
        end
    end
end
