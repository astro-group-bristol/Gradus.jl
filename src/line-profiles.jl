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

"""
    lineprofile(m::AbstractMetric, x::SVector, d::AbstractAccretionGeometry; kwargs...)
    lineprofile(
        m::AbstractMetric, 
        x::SVector, 
        d::AbstractAccretionDisc, 
        profile::AbstractDiscProfile; 
        kwargs...
    )

Compute a line profile for a given set of model components. Returns both the
grid and flux.

The dispatch that includes a `profile` can be used with an
[`AbstractDiscProfile`](@ref), conventionally calculated by a coronal model.

If no profile is specified, the emissivity is assumed to be a power-law
``\\varepsilon(r) = r^{-3}``.

This function accepts a number of keyword arguments depending on method used to
compute the line profile. Common arguments with their defaults are:

- `bins = collect(range(0.1, 1.5, 180))`: the (energy) bins used as the domain of the line profile.
- `method = TransferFunctionMethod()`: used to select which method to use when
   computing the line profile. Alternatives include [`BinningMethod`](@ref).

The [`TransferFunctionMethod`](@ref) dispatch additionally accepts the following keyword arguments:

- `minrₑ = isco(m)`: the innermost transfer function radius to compute / inner
   radius of line profile integration.
- `maxrₑ = 50`: the outermost transfer function radius to compute / outer radius
   of line profile integration.
- `numrₑ = 100`: the number of transfer functions to calculate.
- `verbose = false`: show a progress bar.
- `h = 2e-8`: an integration padding value to avoid numerical instabilities. See
   Dauser et al., 2010 for details.
- `Nr = 1000`: the number of radial steps used to interpolate between transfer
   function branches when integrating.

The [`BinningMethod`](@ref) dispatch additionally accepts the following keyword arguments:

- `λ_max = 2 * x[2]`: the maximum integration time (affine parameter).
- `redshift_pf = ConstPointFunctions.redshift(m, x)`: the function used to
   compute redshift values.
- `verbose = false`: show a progress bar.
- `minrₑ = isco(m)`:` the inner radius of the line profile.
- `maxrₑ = T(50)`:` the outer radius of the line profile.
- `plane = PolarPlane(GeometricGrid(); Nr = 450, Nθ = 1300, r_max = 5maxrₑ)`:
   the image plane used in the calculation.

All other keyword arguments are passed to [`tracegeodesics`](@ref).

There is an additional dispatch that does not accept `bins` as a keyword argument:

    lineprofile(
        bins,
        ε::Function,
        m::AbstractMetric,
        u,
        d::AbstractAccretionGeometry;
        kwargs...
    )

This dispatch is special as it be used to pass any arbitrary function to act as
the emissivity profile of the disc.
"""
lineprofile

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
    n_radii = 1000,
    kwargs...,
)
    _warn_disc_integration_limits(d, minrₑ, maxrₑ)
    tfs = transferfunctions(
        m,
        u,
        d;
        numrₑ = numrₑ,
        verbose = verbose,
        minrₑ = minrₑ,
        maxrₑ = maxrₑ,
        kwargs...,
    )
    flux = integrate_lineprofile(ε, tfs, bins; h = h, n_radii = n_radii)
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
    callback = domain_upper_hemisphere(),
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
        callback = callback,
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
