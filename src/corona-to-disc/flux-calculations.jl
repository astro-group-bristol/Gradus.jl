"""
    lorentz_factor(g::AbstractMatrix, isco, u, v) 

Calculate Lorentz factor in LNRF of `u`.
"""
function lorentz_factor(g::AbstractMatrix, isco_r, u, v)
    frame = Gradus.GradusBase.lnrbasis(g)
    B = reduce(hcat, frame)
    denom = B[:, 1] ⋅ v

    𝒱ϕ = (B[:, 4] ⋅ v) / denom

    if u[2] < isco_r
        𝒱r = (B[:, 2] ⋅ v) / denom
        inv(√(1 - 𝒱r^2 - 𝒱ϕ^2))
    else
        inv(√(1 - 𝒱ϕ^2))
    end
end

function flux_source_to_disc(
    m::AbstractMetricParameters,
    model::AbstractCoronaModel,
    vdp::AbstractDiscProfile;
    kwargs...,
)
    error(
        "Not implemented for metric $(typeof(m)) with model $(typeof(model)) and disc profile $(typeof(vdp)).",
    )
end

function flux_source_to_disc(
    m::AbstractMetricParameters,
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
            @SVector [vtemp[1], -vtemp[2], vtemp[3], vtemp[4]]
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

        γ = lorentz_factor(g_2, isco_r, gp.x, v_disc)
        f_sd = inv(areas[i] / total_area)
        # total reflected flux 
        g_sd^(1 + α) * E_d^(-α) * dA * f_sd / γ
    end

    map(flux, points)
end

function flux_source_to_disc(
    m::AbstractMetricParameters,
    model::AbstractCoronaModel,
    vdp::VoronoiDiscProfile;
    kwargs...,
)
    areas = getareas(vdp)
    flux_source_to_disc(m, model, vdp.geodesic_points, areas; kwargs...)
end

export flux_source_to_disc
