using RecipesBase

function _extract_path(sol, n_points; projection = :none, t_span = 30.0)
    status = Gradus.get_status_code(sol.prob.p)
    mid_i =
        if (status == StatusCodes.IntersectedWithGeometry) ||
           (status == StatusCodes.WithinInnerBoundary)
            lastindex(sol.t)
        else
            max(1, length(sol.u) ÷ 2)
        end

    start_t = max(sol.t[mid_i] - t_span, sol.t[1])
    end_t = min(sol.t[mid_i] + t_span, sol.t[end])

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
@recipe function f(p::Plot_Paths_3D; extent = 20, n_points = 400, t_span = 100.0)
    sol = p.args[1]
    _range = (-extent, extent)
    xlims --> _range
    ylims --> _range
    zlims --> _range
    aspect_ratio --> 1

    itr = if typeof(sol) <: SciMLBase.EnsembleSolution
        sol.u
    else
        sol
    end

    for s in itr
        path = _extract_path(s, n_points, projection = :none, t_span = t_span)
        @series begin
            path
        end
    end
end

@userplot Plot_Paths
@recipe function f(
    p::Plot_Paths;
    extent = 20,
    n_points = 400,
    projection = :none,
    t_span = 100.0,
    plane = :XY,
)
    sol = p.args[1]
    projection := projection
    _range = (-extent, extent)
    if projection == :polar
        xlims --> _range
        ylims --> _range
    else
        xlims --> _range
        ylims --> _range
        aspect_ratio --> 1
    end

    itr = if !(typeof(sol) <: SciMLBase.EnsembleSolution)
        (; u = (sol,))
    else
        sol
    end

    for s in itr.u
        path = _extract_path(s, n_points, projection = projection, t_span = t_span)
        coords = if projection == :polar
            path[3], path[1]
        else
            if plane == :XY
                path[1], path[2]
            elseif plane == :XZ
                path[1], path[3]
            elseif plane == :YZ
                path[2], path[3]
            else
                error("Unknown plane: $(plane). Only XY, XZ, or YZ are implemented.")
            end
        end
        @series begin
            coords
        end
    end
end

@userplot Plot_Horizon_3D
@recipe function f(p::Plot_Horizon_3D; projection = :none, n_points = 32)
    projection := projection
    R = inner_radius(p.args[1])

    u = range(0, stop = 2π, length = n_points)
    v = range(0, stop = π, length = n_points)
    x = R * @. cos(u) * sin(v)'
    y = R * @. sin(u) * sin(v)'
    z = R .* repeat(cos.(v)', outer = [n_points, 1])

    seriestype := :surface

    x, y, z
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

@recipe function f(p::AbstractImagePlane, u = SVector(0.0, 1000.0, π / 2, 0.0))
    a, b = image_plane(p, u)
    seriestype --> :scatter
    color --> :black
    legend --> false
    xlabel --> "α"
    ylabel --> "β"
    (a, b)
end

@recipe function f(p::RadialDiscProfile; normalize = nothing)
    xlabel --> "r (rg)"
    ylabel --> "ε (arb.)"
    xscale --> :log10
    yscale --> :log10

    y = p.ε[2:(end-1)]

    if !isnothing(normalize)
        y = normalize(y)
    end
    p.radii[2:(end-1)], y
end

@recipe function f(ctf::CunninghamTransferData; h = 1e-4)
    xlabel --> "g✶"
    ylabel --> "f"
    mask = @. (ctf.g✶ > h) & (ctf.g✶ < 1 - h)
    ctf.g✶[mask], ctf.f[mask]
end

@userplot Plot_Emissivity_Index
@recipe function f(p::Plot_Emissivity_Index)
    prof = p.args[1]
    if !(prof isa RadialDiscProfile)
        error("Plotting argument for emissivity index must be a radial disc profile.")
    end

    xscale --> :log10
    xlabel --> "r"
    ylabel --> "α"
    title --> "ε ∝ r⁻ᵅ"

    y = prof.ε
    x = prof.radii

    dydx = diff(y) ./ diff(x)
    x = x[2:end]
    y = y[2:end]
    em_index = @. dydx * x / y

    x[1:(end-1)], -em_index[1:(end-1)]
end
