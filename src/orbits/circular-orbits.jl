module CircularOrbits
import ..StaticArrays: @SVector
import ..Gradus:
    AbstractAutoDiffStaticAxisSymmetricParams,
    metric_components,
    metric_jacobian,
    inverse_metric_components

function Ω(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; contra_rotating = false)
    _, jacs = metric_jacobian(m, rθ)
    ∂rg = jacs[:, 1]

    Δ = √(∂rg[5]^2 - ∂rg[1] * ∂rg[4])
    if contra_rotating
        -(∂rg[5] + Δ) / ∂rg[4]
    else
        -(∂rg[5] - Δ) / ∂rg[4]
    end
end

function __energy(g, Ωϕ)
    @inbounds -(g[1] + g[5] * Ωϕ) / √(-g[1] - 2g[5] * Ωϕ - g[4] * Ωϕ^2)
end

function energy(m, rθ; contra_rotating = false)
    Ωϕ = Ω(m, rθ; contra_rotating = contra_rotating)
    __energy(metric_components(m, rθ), Ωϕ)
end
energy(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    energy(m, @SVector([r, π / 2]); kwargs...)

function __angmom(g, Ωϕ, prograde)
    @inbounds res = (g[5] + g[4] * Ωϕ) / √(-g[1] - 2g[5] * Ωϕ - g[4] * Ωϕ^2)
    if prograde
        res
    else
        -res
    end
end

function angmom(m, rθ; contra_rotating = false, prograde = true)
    Ωϕ = Ω(m, rθ; contra_rotating = contra_rotating)
    __angmom(metric_components(m, rθ), Ωϕ, prograde)
end
angmom(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    angmom(m, @SVector([r, π / 2]); kwargs...)

function g_ginv_energy_angmom(
    m::AbstractAutoDiffStaticAxisSymmetricParams,
    rθ;
    contra_rotating = false,
    prograde = true,
)
    g = metric_components(m, rθ)
    ginv = inverse_metric_components(g)

    Ωϕ = Ω(m, rθ; contra_rotating = contra_rotating)
    E = __energy(g, Ωϕ)
    L = __angmom(g, Ωϕ, prograde)

    (g, ginv, E, L)
end

function __vϕ(ginv, E, L)
    -E * ginv[5] + L * ginv[4]
end

function vϕ(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; kwargs...)
    _, ginv, E, L = g_ginv_energy_angmom(m, rθ; kwargs...)
    __vϕ(ginv, E, L)
end
vϕ(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    vϕ(m, @SVector([r, π / 2]); kwargs...)

# this component doesn't actually seem to correctly constrain the geodesic
# to being light- / null-like, or timelike. maybe revist? else call constrain before returning
function __vt(ginv, E, L)
    -E * ginv[1] + L * ginv[5]
end

function vt(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; kwargs...)
    _, ginv, E, L = g_ginv_energy_angmom(m, rθ; kwargs...)
    __vt(ginv, E, L)
end
vt(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    vt(m, @SVector([r, π / 2]); kwargs...)

function fourvelocity(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; kwargs...)
    _, ginv, E, L = g_ginv_energy_angmom(m, rθ; kwargs...)

    vt = __vt(ginv, E, L)
    vϕ = __vϕ(ginv, E, L)

    @SVector [vt, 0.0, 0.0, vϕ]
end
fourvelocity(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    fourvelocity(m, @SVector([r, π / 2]); kwargs...)

function plunging_fourvelocity(m::AbstractAutoDiffStaticAxisSymmetricParams, rθ; kwargs...)
    g, ginv, E, L = g_ginv_energy_angmom(m, rθ; kwargs...)
    vt = __vt(ginv, E, L)
    vϕ = __vϕ(ginv, E, L)

    nom = ginv[1] * E^2 - 2ginv[5] * E * L + ginv[4] * L^2 + 1
    denom = -g[2]

    @SVector[vt, -sqrt(abs(nom / denom)), 0.0, vϕ]
end
plunging_fourvelocity(m::AbstractAutoDiffStaticAxisSymmetricParams, r::Number; kwargs...) =
    plunging_fourvelocity(m, @SVector([r, π / 2]); kwargs...)

end # module

export CircularOrbits
