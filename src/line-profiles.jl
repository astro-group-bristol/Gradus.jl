abstract type AbstractLineProfileAlgorithm end
struct CunninghamLineProfile <: AbstractLineProfileAlgorithm end
struct BinnedLineProfile <: AbstractLineProfileAlgorithm end

@inline function lineprofile(
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry;
    bins = collect(range(0.1, 1.5, 180)),
    kwargs...,
)
    lineprofile(bins, r -> r^-3, m, u, d; kwargs...)
end

@inline function lineprofile(
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry,
    ep::AbstractDiscProfile;
    bins = collect(range(0.1, 1.5, 180)),
    algorithm = BinnedLineProfile(),
    kwargs...,
)
    lineprofile(bins, r -> ep.f.ε(r), m, u, d; algorithm = algorithm, kwargs...)
end

@inline function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry;
    algorithm::AbstractLineProfileAlgorithm = CunninghamLineProfile(),
    kwargs...,
)
    lineprofile(bins, ε, m, u, d, algorithm; kwargs...)
end

function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry,
    ::CunninghamLineProfile,
    ;
    minrₑ = isco(m) + 1e-2, # delta to avoid numerical instabilities
    maxrₑ = 50,
    numrₑ = 100,
    verbose = false,
    h = 2e-8,
    Nr = 1000,
    kwargs...,
)
    radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ)
    itfs = interpolated_transfer_branches(m, u, d, radii; verbose = verbose, kwargs...)
    flux = integrate_lineprofile(ε, itfs, radii, bins; h = h, Nr = Nr)
    bins, flux
end

function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry,
    ::BinnedLineProfile;
    plane = PolarPlane(GeometricGrid(); Nr = 450, Nθ = 1300, r_max = 50.0),
    λ_max = 2 * u[2],
    redshift_pf = ConstPointFunctions.redshift(m, u),
    verbose = false,
    minrₑ = isco(m),
    maxrₑ = T(50),
    # this 5maxrₑ is a heuristic that is insulting
    # todo: make it friendly
    plane = PolarPlane(GeometricGrid(); Nr = 450, Nθ = 1300, r_max = 5maxrₑ),
    solver_args...,
)
    progress_bar = init_progress_bar("Lineprofile: ", trajectory_count(plane), verbose)

    gps = tracegeodesics(
        m,
        u,
        plane,
        d,
        (0.0, λ_max);
        save_on = false,
        verbose = verbose,
        progress_bar = progress_bar,
        callback = domain_upper_hemisphere(),
        ensemble = EnsembleEndpointThreads(),
        solver_args...,
    )

    I = intersected_with_geometry(gps, x -> (minrₑ <= _equatorial_project(x) <= maxrₑ))
    points = @views gps[I]
    areas = unnormalized_areas(plane)[I]

    # calculate physical flux
    g = redshift_pf.(m, points, λ_max)
    r = map(p -> _equatorial_project(p.x), points)

    f = @. ε(r) * g^3 * areas
    # bin
    flux = bucket(Simple(), g, f, bins)
    bins, flux ./ sum(flux)
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
