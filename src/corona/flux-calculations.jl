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
    es = lnrbasis(m, x)
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
    lorentz_factor(m, x, CircularOrbits.fourvelocity(m, _equatorial_project(x)); kwargs...)
function lorentz_factor(m::AbstractMetric, x, v; component = 4)
    𝒱ϕ = local_velocity(m, x, v, component)
    inv(√(1 - 𝒱ϕ^2))
end

# todo: these are currently all unused
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
        E_s = dotproduct(g_1, gp.v_init, v_source)

        # energy at disc
        v_disc = disc_velocity(_equatorial_project(gp.x))
        E_d = dotproduct(g_2, gp.v, v_disc)

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

energy_ratio(m::AbstractMetric, gp::GeodesicPoint, v_src) =
    energy_ratio(m, gp, v_src, CircularOrbits.fourvelocity(m, _equatorial_project(gp.x)))
function energy_ratio(m::AbstractMetric, gp::GeodesicPoint, v_src, v_disc)
    # at the source
    g_src = metric(m, gp.x_init)
    e_src = dotproduct(g_src, gp.v_init, v_src)
    # at the disc
    g_disc = metric(m, gp.x)
    e_disc = dotproduct(g_disc, gp.v, v_disc)

    # TODO: this should be the otherway round to be consistent with g = E_end /
    # E_start but changing this will require making sure all 1/g are replaced
    # with g and vice versa in the emissivity calculations.
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
