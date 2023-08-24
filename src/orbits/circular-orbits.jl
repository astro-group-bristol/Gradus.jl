module CircularOrbits
import ..StaticArrays: SVector, @SVector
import ..Gradus:
    AbstractStaticAxisSymmetric,
    metric_components,
    metric_jacobian,
    inverse_metric_components,
    MuladdMacro

MuladdMacro.@muladd begin
    @inline function _Ω_analytic(∂rg, contra_rotating)
        Δ = √(∂rg[5]^2 - ∂rg[1] * ∂rg[4])
        if contra_rotating
            -(∂rg[5] + Δ) / ∂rg[4]
        else
            -(∂rg[5] - Δ) / ∂rg[4]
        end
    end

    function Ω(m::AbstractStaticAxisSymmetric, rθ; contra_rotating = false)
        _, jacs = metric_jacobian(m, rθ)
        ∂rg = jacs[:, 1]
        _Ω_analytic(∂rg, contra_rotating)
    end

    function ut_uϕ(𝛺::Number, ginv)
        𝒜 = -(𝛺 * ginv[1] - ginv[5])
        ℬ = (𝛺 * ginv[5] - ginv[4])

        denom = ℬ^2 * ginv[1] + 2 * 𝒜 * ℬ * ginv[5] + 𝒜^2 * ginv[4]
        d = -sign(denom) * sqrt(inv(abs(denom)))
        ut = ℬ * d
        uϕ = 𝒜 * d

        SVector(ut, uϕ)
    end

    function ut_uϕ(
        m::AbstractStaticAxisSymmetric,
        rθ,
        ginv = inverse_metric_components(m, rθ);
        kwargs...,
    )
        ut_uϕ(Ω(m, rθ; kwargs...), ginv)
    end

    # these 4 functions can be overwritten for a specific
    # metric, e.g. Kerr-Newman
    function energy(::AbstractStaticAxisSymmetric, rθ, utuϕ)
        -utuϕ[1]
    end
    function angmom(::AbstractStaticAxisSymmetric, rθ, utuϕ)
        utuϕ[2]
    end
    # this component doesn't actually seem to correctly constrain the geodesic
    # to being light- / null-like, or timelike. maybe revist? else call constrain before returning
    vt(::AbstractStaticAxisSymmetric, rθ, ginv, utuϕ) =
        ginv[1] * utuϕ[1] + ginv[5] * utuϕ[2]
    vϕ(::AbstractStaticAxisSymmetric, rθ, ginv, utuϕ) =
        ginv[5] * utuϕ[1] + ginv[4] * utuϕ[2]

    # dispatch helpers
    energy(
        m::AbstractStaticAxisSymmetric,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    ) = energy(m, rθ, ut_uϕ(m, rθ; contra_rotating = contra_rotating, kwargs...); kwargs...)
    angmom(
        m::AbstractStaticAxisSymmetric,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    ) = angmom(m, rθ, ut_uϕ(m, rθ; contra_rotating = contra_rotating, kwargs...); kwargs...)
    energy(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        energy(m, SVector(r, π / 2); kwargs...)
    angmom(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        angmom(m, SVector(r, π / 2); kwargs...)

    function energy_angmom(m::AbstractStaticAxisSymmetric, rθ::SVector; kwargs...)
        utuϕ = ut_uϕ(m, rθ; kwargs...)
        energy(m, rθ, utuϕ; kwargs...), angmom(m, rθ, utuϕ; kwargs...)
    end
    energy_angmom(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        energy_angmom(m, SVector(r, π / 2); kwargs...)

    function vt(
        m::AbstractStaticAxisSymmetric,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    )
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; contra_rotating = contra_rotating, kwargs...)
        vt(m, ginv, rθ, utuϕ)
    end
    vt(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        vt(m, SVector(r, π / 2); kwargs...)

    function vϕ(
        m::AbstractStaticAxisSymmetric,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    )
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; contra_rotating = contra_rotating, kwargs...)
        vϕ(m, rθ, ginv, utuϕ)
    end
    vϕ(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        vϕ(m, SVector(r, π / 2); kwargs...)

    function fourvelocity(m::AbstractStaticAxisSymmetric, rθ; kwargs...)
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; kwargs...)

        SVector(vt(m, rθ, ginv, utuϕ), 0, 0, vϕ(m, rθ, ginv, utuϕ))
    end
    fourvelocity(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        fourvelocity(m, SVector(r, π / 2); kwargs...)
    fourvelocity(m::AbstractStaticAxisSymmetric, x::SVector{4}; kwargs...) =
        fourvelocity(m, SVector(x[2], x[3]); kwargs...)

    function plunging_fourvelocity(
        m::AbstractStaticAxisSymmetric,
        rθ;
        contra_rotating = false,
        kwargs...,
    )
        g = metric_components(m, rθ)
        ginv = inverse_metric_components(g)
        utuϕ = ut_uϕ(m, rθ, ginv; contra_rotating = contra_rotating, kwargs...)
        E = energy(m, rθ, utuϕ; kwargs...)
        L = angmom(m, rθ, utuϕ; kwargs...)
        _vt = vt(m, rθ, ginv, utuϕ; kwargs...)
        _vϕ = vϕ(m, rθ, ginv, utuϕ; kwargs...)

        nom = ginv[1] * E^2 - 2ginv[5] * E * L + ginv[4] * L^2 + 1
        denom = -g[2]

        @SVector[_vt, -sqrt(abs(nom / denom)), 0.0, _vϕ]
    end
    plunging_fourvelocity(m::AbstractStaticAxisSymmetric, r::Number; kwargs...) =
        plunging_fourvelocity(m, @SVector([r, π / 2]); kwargs...)
    plunging_fourvelocity(m::AbstractStaticAxisSymmetric, x::SVector{4}; kwargs...) =
        fourvelocity(m, SVector(x[2], x[3]); kwargs...)

end # mulladd macro
end # module

_keplerian_velocity_projector(m::AbstractMetric, ::AbstractAccretionGeometry) =
    _keplerian_velocity_projector(m)
function _keplerian_velocity_projector(m::AbstractMetric)
    r_isco = isco(m)
    interp = interpolate_plunging_velocities(m)

    function _keplerian_project(x::SVector{4})
        r = _equatorial_project(x)
        if r < r_isco
            vtemp = interp(r)
            SVector(vtemp[1], -vtemp[2], vtemp[3], vtemp[4])
        else
            CircularOrbits.fourvelocity(m, r)
        end
    end
end

export CircularOrbits
