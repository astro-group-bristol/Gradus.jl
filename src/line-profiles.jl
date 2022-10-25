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
    min_re = isco(m),
    max_re = 20,
    num_re = 100,
    bins = range(0.0, 1.5, 100),
    kwargs...,
) where {T}
    # Gauss-Lobatto to include [-1, 1]
    radii, w = gausslobatto(num_re)

    αs = zeros(T, num_points)
    βs = zeros(T, num_points)
    Js = zeros(T, num_points)
    f =
        rₑ -> begin
            # this feels like such a bad practice
            # calling GC on young objects to clear the temporarily allocated memory
            # but we get a significant speedup
            GC.gc(false)
            ctf = cunningham_transfer_function!(
                αs,
                βs,
                Js,
                m,
                u,
                d,
                rₑ,
                2000.0;
                offset_max = rₑ + 20.0,
                kwargs...,
            )
            Gradus._interpolate_branches(ctf)
        end

    α, 𝔉 = _change_interval(f, min_re, max_re)
    ictbs = map(𝔉, radii)

    bin_extrema = extrema(bins)
    y = map(bins) do g
        α * sum(
            ((i, ictb),) -> begin
                w[i] *
                ε(ictb.radius) *
                _integrate_tranfer_function_branches(ictb, g, bin_extrema...)
            end,
            enumerate(ictbs),
        )
    end

    (bins, y)
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
