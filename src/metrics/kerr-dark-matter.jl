module __KerrDarkMatter
using ..StaticArrays
using ..MuladdMacro
using ..Gradus: _smooth_interpolate

@muladd @fastmath begin
    Σ(r, a, θ) = r^2 + a^2 * cos(θ)^2
    Δ(r, R, a) = r^2 + a^2 - R * r

    function G(r, Δr, rₛ)
        dr = (r - rₛ) / Δr
        (3 - 2 * dr) * dr^2
    end

    function dark_matter_mass(r, Δr, rₛ, M)
        if r < rₛ
            zero(r)
        elseif r < rₛ + Δr
            M * G(r, Δr, rₛ)
        else
            M
        end
    end

    # the way this function must be defined is a little complex
    # but helps with type-stability
    function metric_components(M_bh, a, M_dm, Δr, rₛ, rθ)
        (r, θ) = rθ

        M = M_bh + dark_matter_mass(r, Δr, rₛ, M_dm)

        R = 2M
        sinθ2 = sin(θ)^2
        cosθ2 = (1 - sinθ2)
        # slightly faster, especially when considering AD evals
        Σ₀ = r^2 + a^2 * cosθ2

        tt = -(1 - (R * r) / Σ₀)
        rr = Σ₀ / Δ(r, R, a)
        θθ = Σ₀
        ϕϕ = sinθ2 * (r^2 + a^2 + (sinθ2 * R * r * a^2) / Σ₀)

        tϕ = (-R * r * a * sinθ2) / Σ₀
        @SVector [tt, rr, θθ, ϕϕ, tϕ]
    end
end

end # module

"""
    struct KerrDarkMatter

https://arxiv.org/pdf/2003.06829.pdf
"""
@with_kw struct KerrDarkMatter{T} <: AbstractStaticAxisSymmetric{T}
    @deftype T
    "Black hole mass."
    M = 1.0
    "Black hole spin."
    a = 0.0
    "Dark matter mass."
    M_dark_matter = 2.0
    Δr = 20.0
    rₛ = 10.0
end

# implementation
metric_components(m::KerrDarkMatter{T}, rθ) where {T} =
    __KerrDarkMatter.metric_components(m.M, m.a, m.M_dark_matter, m.Δr, m.rₛ, rθ)
inner_radius(m::KerrDarkMatter{T}) where {T} = m.M + √(m.M^2 - m.a^2)

export KerrDarkMatter
