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

function _cunningham_line_profile_integrand(ictb)
    gs -> begin
        if !Gradus._gs_in_limits(ictb, gs)
            0.0
        else
            gmin, gmax = ictb.g_extrema
            g = gstar_to_g(gs, gmin, gmax)
            w = g^3 / √(gs * (1 - gs))
            f = ictb.lower(gs) + ictb.upper(gs)
            w * f / (gmax - gmin)
        end
    end
end

function _calculate_interpolated_transfer_branches(
    m,
    u,
    d,
    radii;
    num_points = 100,
    verbose = false,
    offset = 1e-7,
    kwargs...,
)
    # pre-allocate arrays
    αs = zeros(T, num_points)
    βs = zeros(T, num_points)
    Js = zeros(T, num_points)

    progress_bar = init_progress_bar("Transfer functions:", length(radii), verbose)

    # IILF for calculating the interpolated branches
    𝔉 =
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
                offset_max = 2rₑ + 20.0,
                kwargs...,
            )
            ProgressMeter.next!(progress_bar)
            _interpolate_branches(ctf; offset = offset)
        end

    # calculate interpolated transfer functions for each emission radius
    map(𝔉, radii)
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

    integrate_line_profile(ε, ictbs, bins)
end

function integrate_line_profile(
    ε,
    ictbs::Vector{<:InterpolatedCunninghamTransferBranches{T}},
    bins,
) where {T}
    # global min and max
    ggmin = maximum(i -> first(i.g_limits), ictbs)
    ggmax = minimum(i -> last(i.g_limits), ictbs)

    radii = map(i -> i.radius, ictbs)

    r_low, r_high = extrema(radii)

    y = map(eachindex(@view(bins[1:end-1]))) do index
        bin_low = bins[index]
        bin_high = bins[index+1]
        integrated_bin_for_radii =
            _areas_under_transfer_functions(ε, ictbs, bin_low, bin_high, ggmin, ggmax)
        intp = DataInterpolations.LinearInterpolation(integrated_bin_for_radii, radii)
        res, _ = quadgk(intp, r_low, r_high)
        res
    end

    (bins, y)
end

function _areas_under_transfer_functions(ε, ictbs, bin_low, bin_high, ggmin, ggmax)
    map(ictbs) do ictb
        𝒮 = _cunningham_line_profile_integrand(ictb)
        Sg = g -> begin
            gs = gstar(g, ictb.g_extrema...)
            if ggmin < gs < ggmax
                𝒮(gs)
            else
                0.0
            end
        end
        res, _ = quadgk(Sg, bin_low, bin_high)
        r = ictb.radius
        res * r * ε(r)
    end
end

export AbstractLineProfileAlgorithm, BinnedLineProfile, CunninghamLineProfile, lineprofile
