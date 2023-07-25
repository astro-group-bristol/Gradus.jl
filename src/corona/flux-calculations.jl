"""
    local_velocity(m::AbstractMetric, x, v, component::Int)

Compute the component of the velocity `v` in the LNRF given by `component`,
where in Boyer-Lindquist `(1, 2, 3, 4)` corresponds to `(r, t, θ, ϕ)`.

Evaluates Bardeen+73 (3.9),

```math
\\mathscr{V}^{(i)} = \\frac{v^\\mu e_\\mu^{(i)}}{v^\\nu e_\\nu^{(t)}}.
```
"""
function local_velocity(m::AbstractMetric, x, v, component::Int)
    es = GradusBase.lnrbasis(m, x)
    𝒱t = _fast_dot(es[1], v)
    𝒱i = _fast_dot(es[component], v)
    𝒱i / 𝒱t
end

"""
    lorentz_factor(m::AbstractMetric, g::AbstractAccretionGeometry, x)

Calculate the Lorentz factor of a patch of the accretion geometry `g` at
position `x` via

```math
\\left( 1 - v^{\\phi}^2 \\right)^{\\frac{-1}{2}}
```
"""
lorentz_factor(m::AbstractMetric, ::AbstractAccretionGeometry, x; kwargs...) =
    lorentz_factor(m, x, CircularOrbits.fourvelocity(m, x[2]); kwargs...)
function lorentz_factor(m::AbstractMetric, x, v; component = 4)
    𝒱ϕ = local_velocity(m, x, v, component)
    inv(√(1 - 𝒱ϕ^2))
end

function flux_source_to_disc(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    vdp::AbstractDiscProfile;
    kwargs...,
)
    error(
        "Not implemented for metric $(typeof(m)) with model $(typeof(model)) and disc profile $(typeof(vdp)).",
    )
end

function flux_source_to_disc(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    points,
    areas::AbstractVector;
    α = 1.0,
)
    v_source = source_velocity(m, model)

    total_area = sum(areas)

    isco_r = isco(m)
    intp = interpolate_plunging_velocities(m)

    disc_velocity(r) =
        if r < isco_r
            vtemp = intp(r)
            SVector(vtemp[1], -vtemp[2], vtemp[3], vtemp[4])
        else
            CircularOrbits.fourvelocity(m, r)
        end

    flux = args -> begin
        (i, gp) = args
        g_1 = metric(m, gp.x_init)
        g_2 = metric(m, gp.x)

        # energy at source
        @tullio E_s := -g_1[i, j] * gp.v_init[i] * v_source[j]

        # energy at disc
        v_disc = disc_velocity(gp.x[2])
        @tullio E_d := -g_2[i, j] * gp.v[i] * v_disc[j]

        # relative redshift source to disc
        g_sd = E_d / E_s

        # area element
        dA = inv(√(g_2[2, 2] * g_2[4, 4]))

        γ = lorentz_factor(g_2, gp.x, v_disc)
        f_sd = inv(areas[i] / total_area)
        # total reflected flux 
        g_sd^(1 + α) * E_d^(-α) * dA * f_sd / γ
    end

    map(flux, points)
end

function energy_ratio(m, gp, v_src)
    g_src = metric(m, gp.x_init)
    @tullio e_src := g_src[i, j] * gp.v_init[i] * v_src[j]
    # at the disc
    g_disc = metric(m, gp.x)
    v_disc = CircularOrbits.fourvelocity(m, SVector(gp.x[2], gp.x[3]))
    @tullio e_disc := g_disc[i, j] * gp.v[i] * v_disc[j]
    # ratio g = E_source / E_disc
    e_src / e_disc
end

function flux_source_to_disc(
    m::AbstractMetric,
    model::AbstractCoronaModel,
    vdp::VoronoiDiscProfile;
    kwargs...,
)
    areas = getareas(vdp)
    flux_source_to_disc(m, model, vdp.geodesic_points, areas; kwargs...)
end

export flux_source_to_disc
