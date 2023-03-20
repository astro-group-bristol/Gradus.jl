module CircularOrbits
import ..StaticArrays: SVector, @SVector
import ..Gradus:
    AbstractAutoDiffStaticAxisSymmetricParams,
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

    function Ω(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; contra_rotating = false)
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
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ,
        ginv = inverse_metric_components(m, rθ);
        kwargs...,
    )
        ut_uϕ(Ω(m, rθ; kwargs...), ginv)
    end

    # these 4 functions can be overwritten for a specific
    # metric, e.g. Kerr-Newman
    function energy(::AbstractAutoDiffStaticAxisSymmetricParams, rθ, utuϕ)
        -utuϕ[1]
    end
    function angmom(::AbstractAutoDiffStaticAxisSymmetricParams, rθ, utuϕ)
        utuϕ[2]
    end
    # this component doesn't actually seem to correctly constrain the geodesic
    # to being light- / null-like, or timelike. maybe revist? else call constrain before returning
    vt(::AbstractAutoDiffStaticAxisSymmetricParams, rθ, ginv, utuϕ) =
        ginv[1] * utuϕ[1] + ginv[5] * utuϕ[2]
    vϕ(::AbstractAutoDiffStaticAxisSymmetricParams, rθ, ginv, utuϕ) =
        ginv[5] * utuϕ[1] + ginv[4] * utuϕ[2]

    # dispatch helpers
    energy(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    ) = energy(m, rθ, ut_uϕ(m, rθ; contra_rotating = contra_rotating, kwargs...); kwargs...)
    angmom(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    ) = angmom(m, rθ, ut_uϕ(m, rθ; contra_rotating = contra_rotating, kwargs...); kwargs...)
    energy(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        energy(m, SVector(r, π / 2); kwargs...)
    angmom(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        angmom(m, SVector(r, π / 2); kwargs...)

    function energy_angmom(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        kwargs...,
    )
        utuϕ = ut_uϕ(m, rθ; kwargs...)
        energy(m, rθ, utuϕ; kwargs...), angmom(m, rθ, utuϕ; kwargs...)
    end
    energy_angmom(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        energy_angmom(m, SVector(r, π / 2); kwargs...)

    function vt(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    )
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; contra_rotating = contra_rotating, kwargs...)
        vt(m, ginv, rθ, utuϕ)
    end
    vt(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        vt(m, SVector(r, π / 2); kwargs...)

    function vϕ(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        contra_rotating = false,
        kwargs...,
    )
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; contra_rotating = contra_rotating, kwargs...)
        vϕ(m, rθ, ginv, utuϕ)
    end
    vϕ(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        vϕ(m, SVector(r, π / 2); kwargs...)

    function fourvelocity(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        rθ::SVector;
        kwargs...,
    )
        ginv = inverse_metric_components(m, rθ)
        utuϕ = ut_uϕ(m, rθ, ginv; kwargs...)

        SVector(vt(m, rθ, ginv, utuϕ), 0, 0, vϕ(m, rθ, ginv, utuϕ))
    end
    fourvelocity(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
        fourvelocity(m, SVector(r, π / 2); kwargs...)

    function plunging_fourvelocity(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
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
    plunging_fourvelocity(
        m::AbstractAutoDiffStaticAxisSymmetricParams,
        r::Number;
        kwargs...,
    ) = plunging_fourvelocity(m, @SVector([r, π / 2]); kwargs...)

end # mulladd macro
end # module

export CircularOrbits
