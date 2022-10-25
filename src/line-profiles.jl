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
    m::AbstractMetricParams,
    u,
    d::AbstractAccretionGeometry,
    ε;
    num_points = 100,
    min_re = isco(m),
    max_re = 20,
    num_re = 100,
    bins = range(0.0, 1.5, 100),
    kwargs...,
)
    radii, w = gausslobatto(num_re)
    f = rₑ -> begin
        ctf = cunningham_transfer_function(
            m,
            u,
            d,
            rₑ,
            2000.0;
            num_points = num_points,
            offset_max = rₑ + 20.0,
            kwargs...,
        )
        Gradus._interpolate_branches(ctf)
    end

    α, 𝔉 = Gradus._change_interval(f, min_re, max_re)

    bin_extrema = extrema(bins)
    y = map(bins) do i
        α * sum(
            (i, r) -> begin
                ictb = 𝔉(r)
                w[i] * ε(ictb.radius) * _integrate_tranfer_function_branches(ictb, i, bin_extrema...)
            end,
            enumerate(radii)
        )
    end

    (bins, y)
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
