abstract type AbstractImagePlane{G} end
image_plane(plane::AbstractImagePlane, u) = error("Not implemented for $plane")
trajectory_count(plane::AbstractImagePlane) = error("Not implemented for $plane")
unnormalized_areas(plane::AbstractImagePlane) = error("Not implemented for $plane")

function impact_parameters(plane::AbstractImagePlane, u)
    αs, βs = image_plane(plane, u)
    vec(αs), vec(βs)
end

struct PolarPlane{G,T} <: AbstractImagePlane{G}
    grid::G
    Nr::Int
    Nθ::Int
    # domain
    r_min::T
    r_max::T
    θ_min::T
    θ_max::T
end

function PolarPlane(
    grid::AbstractImpactParameterGrid;
    Nr = 400,
    Nθ = 100,
    r_min = 1.0,
    r_max = 250.0,
    θ_min = 0.0,
    θ_max = 2π,
)
    PolarPlane(grid, Nr, Nθ, r_min, r_max, θ_min, θ_max)
end

function trajectory_count(plane::PolarPlane)
    plane.Nr * plane.Nθ
end

function image_plane(plane::PolarPlane, u)
    rs = plane.grid(plane.r_min, plane.r_max, plane.Nr)
    δθ = (plane.θ_max - plane.θ_min) / (plane.Nθ)
    θs = range(plane.θ_min, plane.θ_max - δθ, plane.Nθ)

    αs = [r * cos(θ) for r in rs, θ in θs]
    βs = [r * sin(θ) * sin(u[3]) for r in rs, θ in θs]

    αs, βs
end

function unnormalized_areas(plane::PolarPlane)
    rs = collect(plane.grid(plane.r_min, plane.r_max, plane.Nr))
    # todo: do we need dr here?
    # dr = diff(rs) 
    # push!(dr, 0.0)
    A = @. rs^2 # * abs(dr)
    repeat(A, inner = (1, plane.Nθ))
end

struct CartesianPlane{G,T} <: AbstractImagePlane{G}
    grid::G
    Nx::Int
    Ny::Int
    x_min::T
    x_max::T
    y_min::T
    y_max::T
end

function CartesianPlane(
    grid::AbstractImpactParameterGrid;
    Nx = 150,
    Ny = 150,
    x_min = 0.0,
    x_max = 150.0,
    y_min = 0.0,
    y_max = 150.0,
)
    CartesianPlane(grid, Nx, Ny, x_min, x_max, y_min, y_max)
end

function trajectory_count(plane::CartesianPlane)
    (2 * (plane.Ny ÷ 2) - 1) * (2 * (plane.Nx ÷ 2) - 1)
end

function image_plane(plane::CartesianPlane, u)
    xs = collect(plane.grid(plane.x_min, plane.x_max, plane.Nx ÷ 2))
    ys = collect(plane.grid(plane.y_min, plane.y_max, plane.Ny ÷ 2))

    X_size = 2 * (plane.Ny ÷ 2) - 1
    Y_size = 2 * (plane.Nx ÷ 2) - 1

    X = @views repeat(xs[2:end]', inner = (X_size, 1))
    Y = @views repeat(ys[2:end], inner = (1, Y_size))
    αs = hcat(-reverse(X, dims = 2), fill(xs[1], X_size), X)
    βs = vcat(-reverse(Y, dims = 1), fill(ys[1], Y_size)', Y)

    # αs[:, 1:size(αs,2)÷2] = -reverse(αs[:, 1:size(αs,2)÷2], dims=2)
    # βs[1:size(βs,1)÷2, :] = -reverse(βs[1:size(βs,1)÷2, :], dims=1)

    αs, βs
end

function tracegeodesics(
    m::AbstractMetricParams,
    observer_position,
    plane::AbstractImagePlane,
    args...;
    kwargs...,
)
    αs, βs = impact_parameters(plane, observer_position)
    velfunc(i) = map_impact_parameters(m, observer_position, αs[i], βs[i])
    tracegeodesics(
        m,
        observer_position,
        velfunc,
        args...;
        trajectories = trajectory_count(plane),
        kwargs...,
    )
end

export PolarPlane, CartesianPlane, image_plane
