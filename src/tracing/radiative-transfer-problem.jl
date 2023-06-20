function radiative_transfer(m::AbstractMetric, x, k, geometry, I, ν, invν3, r_isco)
    a_ν, j_ν, u = covariant_absorption_emission_velocity(m, x, ν, geometry, r_isco)
    g = metric(m, x)
    # cache inv(ν)^3 to avoid costly division
    k_rev = SVector(k[1], k[2], k[3], -k[4])
    -(g * k_rev) ⋅ u * (-a_ν * I + j_ν * invν3)
end

function covariant_absorption_emission_velocity(
    m::AbstractMetric,
    x,
    ν,
    d::AbstractAccretionGeometry,
    r_isco,
)
    u = if x[2] > r_isco
        CircularOrbits.fourvelocity(m, x[2])
    else
        CircularOrbits.plunging_fourvelocity(m, x[2])
    end

    (absorption_coefficient(m, d, x, ν), emissivity_coefficient(m, d, x, ν), u)
end

absorption_coefficient(m::AbstractMetric, d::AbstractAccretionGeometry, x, ν) = 0.0
emissivity_coefficient(m::AbstractMetric, d::AbstractAccretionGeometry, x, ν) = 0.0

mutable struct RadiativeTransferIntegrationParameters <: AbstractIntegrationParameters
    status::StatusCodes.T
    within_geometry::Bool
end

function update_integration_parameters!(
    p::RadiativeTransferIntegrationParameters,
    new::RadiativeTransferIntegrationParameters,
)
    p.status = new.status
    p.within_geometry = new.within_geometry
    p
end

function _radiative_goemtry_collision_callback(
    g::AbstractAccretionDisc;
    gtol,
    interp_points = 8,
)
    function _radiative_transfer_geometry_callback(u, λ, integrator)
        dist = distance_to_disc(g, u; gtol = gtol)
        integrator.p.within_geometry ? -dist : dist
    end
    function _radiate_transfer_geometry_effect!(integrator)
        # flip
        integrator.p.within_geometry = !integrator.p.within_geometry
    end

    ContinuousCallback(
        _radiative_transfer_geometry_callback,
        _radiate_transfer_geometry_effect!,
        interp_points = interp_points,
        save_positions = (false, false),
    )
end

function geometry_collision_callback(
    g::AbstractAccretionDisc,
    ::TraceRadiativeTransfer;
    gtol,
    interp_points = 8,
)
    # callbacks should only be defined on optically thick geometries
    if is_optically_thin(g)
        return _thin_geometry_collision_callback(
            g;
            gtol = gtol,
            interp_points = interp_points,
        )
    else
        return _radiative_goemtry_collision_callback(g; gtol = gtol, interp_points = interp_points)
    end
end

function radiative_transfer_ode_problem(
    trace::TraceRadiativeTransfer,
    m::AbstractMetric,
    pos,
    vel,
    time_domain,
    callback,
    geometry,
)
    invν3 = inv(trace.ν)^3
    r_isco = Gradus.isco(m)

    function f(u::SVector{9,T}, p, λ) where {T}
        @inbounds let x = SVector{4,T}(u[1:4]), k = SVector{4,T}(u[5:8]), I = u[9]

            dk = SVector{4,T}(geodesic_equation(m, x, k))
            dI = if p.within_geometry
                radiative_transfer(m, x, k, geometry, I, trace.ν, invν3, r_isco)
            else
                0.0
            end
            vcat(k, dk, SVector(dI))
        end
    end

    # add our new quantity that we want to trace
    u_init = vcat(pos, vel, SVector(trace.I₀))
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        RadiativeTransferIntegrationParameters(StatusCodes.NoStatus, false);
        callback = callback,
    )
end

function assemble_tracing_problem(
    trace::TraceRadiativeTransfer,
    config::TracingConfiguration,
)
    if isnothing(config.geometry)
        error(
            "Cannot calculate radiative transfer without geometry (geometry is `nothing`).",
        )
    end
    # create the callback set for the problem
    cbs = create_callback_set(config.metric, config.callback, config.chart)
    # wrap a function that can build the ODE problems
    function _problem_builder(x, v)
        radiative_transfer_ode_problem(
            trace,
            config.metric,
            x,
            v,
            config.λ_domain,
            cbs,
            config.geometry,
        )
    end
    wrap_arguments(config, _problem_builder, trace.μ)
end
