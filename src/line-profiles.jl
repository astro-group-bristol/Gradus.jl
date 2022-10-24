abstract type AbstractLineProfileAlgorithm end
struct CunninghamLineProfile <: AbstractLineProfileAlgorithm end
struct BinnedLineProfile <: AbstractLineProfileAlgorithm end

@inline function lineprofile(
    intensity_func::Function,
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry;
    algorithm::AbstractLineProfileAlgorithm = CunninghamLineProfile(),
    kwargs...,
)
    lineprofile(algorithm, m, u, d, intensity_func; kwargs...)
end

function lineprofile(
    ::CunninghamLineProfile,
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    intensity;
    num_points = 100,
    radii = range(isco(m), 30.0, 100),
    bins = range(0.0, 1.5, 100),
    kwargs...,
)
    ictbs = map(radii) do re
        ctf = cunningham_transfer_function(
            m,
            u,
            d,
            re,
            2000.0;
            num_points = num_points,
            offset_max = 1.5 * re + 1.0,
            kwargs...,
        )
        _interpolate_branches(ctf)
    end
    y = map(bins) do i
        _integrate_tranfer_function_branches(ictbs, i, intensity)
    end

    (bins, y)
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
