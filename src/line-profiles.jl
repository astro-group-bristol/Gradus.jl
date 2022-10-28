abstract type AbstractLineProfileAlgorithm end
struct CunninghamLineProfile <: AbstractLineProfileAlgorithm end
struct BinnedLineProfile <: AbstractLineProfileAlgorithm end

@inline function lineprofile(
    ε::Function,
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry;
    algorithm::AbstractLineProfileAlgorithm = CunninghamLineProfile(),
    kwargs...,
)
    lineprofile(algorithm, m, u, d, ε; kwargs...)
end

function _change_interval(f, a, b)
    α = (b - a) / 2
    β = (a + b) / 2
    (α, x -> f(α * x + β))
end

function lineprofile(
    ::CunninghamLineProfile,
    m::AbstractMetricParams{T},
    u,
    d::AbstractAccretionGeometry,
    ε;
    num_points = 100,
    min_re = isco(m) + 1e-2, # delta to avoid numerical instabilities
    max_re = 20,
    num_re = 100,
    bins = range(0.0, 1.5, 100),
    verbose = false,
    offset = 1e-7,
    kwargs...,
) where {T}
    # this is just a placeholder: desire a distribution that favours
    # small radii over large radii, and this one does that quite well
    # radii here are emission radii rₑ
    radii = exp.(range(log(1), log(1000), num_re))
    _max_radii = maximum(radii)
    # rescale inplace
    @. radii = (radii / _max_radii) * (max_re - min_re) * min_re

    ictbs = _calculate_interpolated_transfer_branches(
        m,
        u,
        d,
        radii;
        num_points = num_points,
        verbose = verbose,
        offset = offset,
        kwargs...,
    )

    integrate_transfer_function(ε, ictbs, bins)
end


export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
