abstract type AbstractImagePlane{G} end
image_plane(plane::AbstractImagePlane, u) = error("Not implemented for $plane")
trajectory_count(plane::AbstractImagePlane) = error("Not implemented for $plane")

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
function trajectory_count(plane::PolarPlane)
    plane.Nr * plane.Nθ
end
function image_plane(plane::PolarPlane, u)
    rs = plane.grid(plane.r_min, plane.r_max, plane.Nr)
    θs = range(plane.θ_min, plane.θ_max, plane.Nθ)

    αs = [r * cos(θ) for r in rs, θ in θs]
    βs = [r * sin(θ) * sin(u[3]) for r in rs, θ in θs]

    αs, βs
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

function tracegeodesics(
    m::AbstractMetricParams{T},
    observer_position,
    plane::AbstractImagePlane,
    time_domain::Tuple{T,T};
    kwargs...,
) where {T}
    αs, βs = impact_parameters(plane, observer_position)
    velfunc(i) = map_impact_parameters(m, observer_position, αs[i], βs[i])
    tracegeodesics(
        m,
        observer_position,
        velfunc,
        time_domain;
        trajectories = trajectory_count(plane),
        kwargs...,
    )
end

export PolarPlane, image_plane
