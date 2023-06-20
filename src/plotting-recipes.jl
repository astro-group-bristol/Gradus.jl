using RecipesBase

function _extract_path(sol, n_points, projection = :none)
    mid_i = max(1, length(sol.u) ÷ 2)

    start_t = max(sol.t[mid_i] - 100.0, sol.t[1])
    end_t = min(sol.t[mid_i] + 100.0, sol.t[end])

    t_range = range(start_t, end_t, n_points)

    r = [sol(t)[2] for t in t_range]
    θ = [sol(t)[3] for t in t_range]
    ϕ = [sol(t)[4] for t in t_range]

    if projection == :polar
        r, θ, ϕ
    else
        x = @. r * cos(ϕ) * sin(θ)
        y = @. r * sin(ϕ) * sin(θ)
        z = @. r * cos(θ) 
        x, y, z
    end
end

@userplot Plot_Paths_3D
@recipe function f(p::Plot_Paths_3D; range = 20, n_points = 400)
    sol = p.args[1]
    range = (-range, range)
    xlims --> range
    ylims --> range
    zlims --> range
    aspect_ratio --> 1

    itr = if !(typeof(sol) <: SciMLBase.EnsembleSolution)
        (; u = (sol,))
    else
        sol
    end

    for s in itr.u
        path = _extract_path(s, n_points, :none)
        @series begin
            path
        end
    end
end

@userplot Plot_Paths
@recipe function f(p::Plot_Paths; range = (0, 20), n_points = 400, projection = :none)
    sol = p.args[1]
    projection := projection
    range = (-range, range)
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
        path = _extract_path(s, n_points, projection)
        coords = if projection == :polar
            path[3], path[1]
        else
            path[1], path[2]
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
