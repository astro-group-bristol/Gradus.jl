"""
    AbstractImagePlane

An plane abstraction used to represent the observer's image plane. These are
particularly relevant for [`BinningMethod`](@ref) computations.

Concerete implementations include:
- [`PolarPlane`](@ref)
- [`CartesianPlane`](@ref)

An `AbstractImagePlane` must implement the following functions:
- [`image_plane`](@ref)
- [`trajectory_count`](@ref)
- [`unnormalized_areas`](@ref)
"""
abstract type AbstractImagePlane{G} end

"""
    image_plane(plane::AbstractImagePlane)
    image_plane(plane::AbstractImagePlane, x::SVector{4})

Return two vectors, representing the ``\\alpha`` and ``\\beta`` impact
parameters that parameterise the image plane. Each pair of ``\\alpha`` and
``\\beta`` must correspond to an [`unnormalized_areas`](@ref).

"""
image_plane(plane::AbstractImagePlane) = image_plane(plane, SVector(0, 10000, π / 2, 0))
# todo: remove this, as i don't think it will be needed
image_plane(plane::AbstractImagePlane, x) = error("Not implemented for $plane")

"""
    trajectory_count(plane::AbstractImagePlane)

Return an integer that counts how many unique geodesics need to be calculated to
map the plane. This should be equal to `length(first(image_plane(plane)))`, and
is used to pre-allocate buffers.
"""
trajectory_count(plane::AbstractImagePlane) = error("Not implemented for $plane")

"""
    unnormalized_areas(plane::AbstractImagePlane)

Return a vector where each element is the area (number) of a given geodesic
element on the image plane. For a pixel image plane, each area will be a
constant `1`. These are used to weight the contributions of each region when
calculating observational results.
"""
unnormalized_areas(plane::AbstractImagePlane) = error("Not implemented for $plane")

function impact_parameters(plane::AbstractImagePlane, x)
    αs, βs = image_plane(plane, x)
    vec(αs), vec(βs)
end


"""
    function PolarPlane(
        grid::Abstract2DGrid;
        Nr = 400,
        Nθ = 100,
        r_min = 1.0,
        r_max = 250.0,
        θ_min = 0.0,
        θ_max = 2π,
    )

Divide the image plane into a polar grid centered at ``\\alpha = 0`` and
``\\beta = 0``.
"""
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
    grid::Abstract2DGrid;
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

function image_plane(plane::PolarPlane, x)
    rs = plane.grid(plane.r_min, plane.r_max, plane.Nr)
    δθ = (plane.θ_max - plane.θ_min) / (plane.Nθ)
    θs = range(plane.θ_min, plane.θ_max - δθ, plane.Nθ)

    αs = [r * cos(θ) for r in rs, θ in θs]
    βs = [r * sin(θ) for r in rs, θ in θs]

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

"""
    function CartesianPlane(
        grid::Abstract2DGrid;
        Nx = 150,
        Ny = 150,
        x_min = 0.0,
        x_max = 150.0,
        y_min = 0.0,
        y_max = 150.0,
    )

Represent the image plane as equi-rectangular regions.
"""
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
    grid::Abstract2DGrid;
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

function image_plane(plane::CartesianPlane, x)
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

function unnormalized_areas(plane::CartesianPlane{<:LinearGrid,T}) where {T}
    X_size = 2 * (plane.Ny ÷ 2) - 1
    Y_size = 2 * (plane.Nx ÷ 2) - 1
    ones(T, (Y_size, X_size))
end

function promote_velfunc(m::AbstractMetric, position, plane::AbstractImagePlane, _unused)
    αs, βs = impact_parameters(plane, position)
    velfunc(i) = map_impact_parameters(m, position, αs[i], βs[i])
    velfunc, trajectory_count(plane)
end

export PolarPlane, CartesianPlane, image_plane
