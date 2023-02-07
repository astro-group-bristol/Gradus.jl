abstract type AbstractLineProfileAlgorithm end
struct CunninghamLineProfile <: AbstractLineProfileAlgorithm end
struct BinnedLineProfile <: AbstractLineProfileAlgorithm end

@inline function lineprofile(m::AbstractMetricParams, u, d::AbstractAccretionGeometry; kwargs...)
    bins = collect(range(0.1, 1.3, 120))
    lineprofile(bins, r -> r^-3, m, u, d; kwargs...)
end

@inline function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetricParams,
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
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    ::CunninghamLineProfile,
    ;
    minrₑ = isco(m) + 1e-2, # delta to avoid numerical instabilities
    maxrₑ = 50,
    numrₑ = 100,
    verbose = false,
    Ng✶ = 67,
    h = 2e-8,
    kwargs...,
)
    radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ)
    itfs = interpolated_transfer_branches(m, u, d, radii; verbose = verbose, kwargs...)
    flux = integrate_drdg✶(ε, itfs, radii, bins; h = h, Ng✶ = Ng✶)
    bins, flux
end

function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    ::BinnedLineProfile;
    plane = PolarPlane(GeometricGrid(); Nr = 250, Nθ = 1300, r_max = 50.0),
    λ_max = 2 * u[2],
    redshift_pf = ConstPointFunctions.redshift(m, u),
    kwargs...,
)
    simsols = tracegeodesics(m, u, plane, d, (0.0, λ_max); save_on = false)
    gps = process_solution.(m, simsols.u)

    I = [i.status == StatusCodes.IntersectedWithGeometry for i in gps]
    points = @views gps[I]
    areas = unnormalized_areas(plane)[I]

    # calculate physical flux
    g = redshift_pf.(m, points, λ_max)
    r = [p.u2[2] for p in points]

    f = @. ε(r) * g^3 * areas
    # bin
    F = bucket(Simple(), g, f, bins)
    bins, F ./ sum(F)
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
