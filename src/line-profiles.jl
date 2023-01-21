abstract type AbstractLineProfileAlgorithm end
struct CunninghamLineProfile <: AbstractLineProfileAlgorithm end
struct BinnedLineProfile <: AbstractLineProfileAlgorithm end

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
    m::AbstractMetricParams{T},
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
) where {T}
    radii = Grids._inverse_grid(minrₑ, maxrₑ, numrₑ)
    itfs = interpolated_transfer_branches(m, u, d, radii; verbose = verbose, kwargs...)
    flux = integrate_drdg✶(ε, itfs, radii, bins; h = h, Ng✶ = Ng✶)
    bins, flux
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
