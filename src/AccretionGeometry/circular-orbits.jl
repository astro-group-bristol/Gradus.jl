module CircularOrbits
import ..StaticArrays: @SVector
import ..AccretionGeometry: AbstractAutoDiffStaticAxisSymmetricParams, metric_components, metric_jacobian, inverse_metric_components

function Ω(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, rθ::AbstractVector{T}, pos) where {T}
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
    -(g[1] + g[5] * Ωϕ) / √(-g[1] - 2g[5]*Ωϕ - g[4]*Ωϕ^2)
end

energy(m, r; kwargs...) = let rθ = @SVector([r, π/2])
    Ωϕ = Ω(m, rθ, contra_rotating)
    __energy(metric_components(m), Ωϕ; kwargs...)
end

function __angmom(g, Ωϕ, prograde) where {T}
    res = (g[5] + g[4] * Ωϕ) / √(-g[1] - 2g[5]*Ωϕ - g[4]*Ωϕ^2)
    if prograde
        res
    else
        -res
    end
end

angmom(m, r; contra_rotating=false, prograde=true) = let rθ = @SVector([r, π/2])
    Ωϕ = Ω(m, rθ, contra_rotating)
    __angmom(metric_components(m), Ωϕ, prograde)
end

function vϕ(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, r; contra_rotating=false, prograde=true) where {T}
    let rθ = @SVector([r, π/2])

        g = metric_components(m, rθ)
        ginv = inverse_metric_components(g)

        Ωϕ = Ω(m, rθ, contra_rotating)
        E = __energy(g, Ωϕ)
        L = __angmom(g, Ωϕ, prograde)

        -E * ginv[5] + L * ginv[4]
    end
end

end # module