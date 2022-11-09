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
    lineprofile(bins, algorithm, m, u, d, ε; kwargs...)
end

function lineprofile(
    bins,
    ::CunninghamLineProfile,
    m::AbstractMetricParams{T},
    u,
    d::AbstractAccretionGeometry,
    ε;
    minrₑ = isco(m) + 1e-2, # delta to avoid numerical instabilities
    maxrₑ = 50,
    numrₑ = 50,
    verbose = false,
    Ng✶ = 10,
    h = 2e-8,
    kwargs...,
) where {T}
    radii = weighted_rₑ_grid(minrₑ, maxrₑ, numrₑ)
    itfs = interpolated_transfer_branches(m, u, d, radii; verbose = verbose, kwargs...)
    flux = integrate_drdg✶(ε, itfs, radii, bins; h = h, Ng✶=Ng✶)
    bins, flux
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
