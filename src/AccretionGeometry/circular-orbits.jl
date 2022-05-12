module CircularOrbits
import ..StaticArrays: @SVector
import ..AccretionGeometry:
    AbstractAutoDiffStaticAxisSymmetricParams,
    metric_components,
    metric_jacobian,
    inverse_metric_components

function Ω(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    rθ::AbstractVector{T},
    pos,
) where {T}
    jacs = metric_jacobian(m, rθ)
    ∂rg = jacs[:, 1]

    Δ = √(∂rg[5]^2 - ∂rg[1] * ∂rg[4])
    if pos
        -(∂rg[5] + Δ) / ∂rg[4]
    else
        -(∂rg[5] - Δ) / ∂rg[4]
    end
end

function __energy(g, Ωϕ) where {T}
    -(g[1] + g[5] * Ωϕ) / √(-g[1] - 2g[5] * Ωϕ - g[4] * Ωϕ^2)
end

energy(m, r; kwargs...) =
    let rθ = @SVector([r, π / 2])
        Ωϕ = Ω(m, rθ, contra_rotating)
        __energy(metric_components(m), Ωϕ; kwargs...)
    end

function __angmom(g, Ωϕ, prograde) where {T}
    res = (g[5] + g[4] * Ωϕ) / √(-g[1] - 2g[5] * Ωϕ - g[4] * Ωϕ^2)
    if prograde
        res
    else
        -res
    end
end

angmom(m, r; contra_rotating = false, prograde = true) =
    let rθ = @SVector([r, π / 2])
        Ωϕ = Ω(m, rθ, contra_rotating)
        __angmom(metric_components(m), Ωϕ, prograde)
    end

function g_ginv_energy_angmom(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    r;
    contra_rotating = false,
    prograde = true,
) where {T}
    let rθ = @SVector([r, π / 2])

        g = metric_components(m, rθ)
        ginv = inverse_metric_components(g)

        Ωϕ = Ω(m, rθ, contra_rotating)
        E = __energy(g, Ωϕ)
        L = __angmom(g, Ωϕ, prograde)

        (g, ginv, E, L)
    end
end

function __vϕ(ginv, E, L)
    -E * ginv[5] + L * ginv[4]
end

function vϕ(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, r; kwargs...) where {T}
    _, ginv, E, L = g_ginv_energy_angmom(m, r; kwargs...)
    __vϕ(ginv, E, L)
end

function __vt(ginv, E, L)
    -E * ginv[1] + L * ginv[5]
end

function vt(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    r;
    contra_rotating = false,
    prograde = true,
) where {T}
    _, ginv, E, L = g_ginv_energy_angmom(m, r; kwargs...)
    __vt(ginv, E, L)
end

function fourvelocity(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    r;
    kwargs...,
) where {T}
    _, ginv, E, L = g_ginv_energy_angmom(m, r; kwargs...)

    vt = __vt(ginv, E, L)
    vϕ = __vϕ(ginv, E, L)

    @SVector [vt, 0.0, 0.0, vϕ]
end

function plunging_fourvelocity(
    m::AbstractAutoDiffStaticAxisSymmetricParams{T},
    r;
    kwargs...,
) where {T}
    g, ginv, E, L = g_ginv_energy_angmom(m, r; kwargs...)
    vt = __vt(ginv, E, L)
    vϕ = __vϕ(ginv, E, L)

    nom = ginv[1] * E^2 - 2ginv[5] * E * L + ginv[4] * L^2 + 1
    denom = -g[2]

    @SVector[vt, -sqrt(abs(nom / denom)), 0.0, vϕ]
end

end # module
