_warn_disc_integration_limits(::AbstractAccretionGeometry, _, _) = nothing
function _warn_disc_integration_limits(d::ThinDisc, minrₑ, maxrₑ)
    if !(d.inner_radius ≈ minrₑ) && (d.inner_radius != 0)
        @warn """
        `ThinDisc` inner radius does not equal `minrₑ`. Integration occurs only over 
        `minrₑ`, and inner radius of the disc is ignored. Set the inner radius to `0` to
        supress this message.
        """
    end
    if !(d.inner_radius ≈ maxrₑ) && (d.outer_radius != Inf)
        @warn """
        `ThinDisc` outer radius does not equal `maxrₑ`. Integration occurs only over 
        `maxrₑ`, and outer radius of the disc is ignored. Set the outer radius to `Inf` to
        supress this message.
        """
    end
end

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
    profile::AbstractDiscProfile;
    bins = collect(range(0.1, 1.5, 180)),
    method = BinningMethod(),
    kwargs...,
)
    lineprofile(bins, r -> emissivity_at(profile, r), m, u, d; method = method, kwargs...)
end

@inline function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry;
    method::AbstractComputationalMethod = TransferFunctionMethod(),
    kwargs...,
)
    lineprofile(bins, ε, m, u, d, method; kwargs...)
end

function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric,
    u,
    d::AbstractAccretionGeometry,
    ::TransferFunctionMethod,
    ;
    minrₑ = isco(m) + 1e-2, # delta to avoid numerical instabilities
    maxrₑ = 50,
    numrₑ = 100,
    verbose = false,
    h = 2e-8,
    Nr = 1000,
    kwargs...,
)
    _warn_disc_integration_limits(d, minrₑ, maxrₑ)
    radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ)
    itfs = interpolated_transfer_branches(m, u, d, radii; verbose = verbose, kwargs...)
    flux = integrate_lineprofile(ε, itfs, radii, bins; h = h, Nr = Nr)
    bins, flux
end

function lineprofile(
    bins,
    ε::Function,
    m::AbstractMetric{T},
    u,
    d::AbstractAccretionGeometry,
    ::BinningMethod;
    λ_max = 2 * u[2],
    redshift_pf = ConstPointFunctions.redshift(m, u),
    verbose = false,
    minrₑ = isco(m),
    maxrₑ = T(50),
    # this 5maxrₑ is a heuristic that is insulting
    # todo: make it friendly
    plane = PolarPlane(GeometricGrid(); Nr = 450, Nθ = 1300, r_max = 5maxrₑ),
    solver_args...,
) where {T}
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

export AbstractComputationalMethod, BinningMethod, TransferFunctionMethod, lineprofile
