using RecipesBase

@userplot Plot_Paths
@recipe function f(p::Plot_Paths; range = (0, 20), n_points = 400, projection = :none)
    sol = p.args[1]
    projection := projection
    range = (-range[2], range[2])
    if projection == :polar
        xlims --> range
        ylims --> range
    else
        xlims --> range
        ylims --> range
        aspect_ratio --> 1
    end

    itr = if !(typeof(sol) <: SciMLBase.EnsembleSolution)
        (; u = (sol,))
    else
        sol
    end

    for s in itr.u
        k = length(s.u) ÷ 2
        λ_range =
            Base.range(max(s.t[k] - 100.0, s.t[1]), min(s.t[k] + 100.0, s.t[end]), n_points)
        r = [s(t)[2] for t in λ_range]
        ϕ = [s(t)[4] for t in λ_range]
        coords = if projection == :polar
            ϕ, r
        else
            x = @. r * cos(ϕ)
            y = @. r * sin(ϕ)
            x, y
        end
        @series begin
            coords
        end
    end
end

@userplot Plot_Horizon
@recipe function f(p::Plot_Horizon; projection = :none, n_points = 100)
    projection := projection
    R = inner_radius(p.args[1])

    ϕ = collect(range(0.0, 2π, n_points))
    r = fill(R, size(ϕ))

    if projection == :polar
        ϕ, r
    else
        x = @. r * cos(ϕ)
        y = @. r * sin(ϕ)
        x, y
    end
end

@recipe function f(p::PolarPlane, u = SVector(0.0, 1000.0, π / 2, 0.0))
    a, b = image_plane(p, u)
    seriestype --> :scatter
    color --> :black
    legend --> false
    xlabel --> "α"
    ylabel --> "β"
    (a, b)
end
