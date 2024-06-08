function radiative_transfer(m::AbstractMetric, x, k, geometry, I, ν₀, r_isco, λ)
    # fluid velocity in at the current point
    u = fluid_velocity(m, geometry, x, r_isco, λ)
    # prefactor ds / dλ
    dsdλ = -dotproduct(m, x, k, u)
    # radiative transfer coefficients
    ν = ν₀ * dsdλ
    aν, jν = fluid_absorption_emission(m, geometry, x, ν, u)
    # covariant radiative transfer
    dsdλ * (-aν * I + jν / ν^3)
end

function fluid_velocity(m::AbstractMetric, ::AbstractAccretionGeometry, x, r_isco, λ)
    r = _equatorial_project(x)
    if r > r_isco
        CircularOrbits.fourvelocity(m, r)
    else
        CircularOrbits.plunging_fourvelocity(m, r)
    end
end

function fluid_absorption_emission(m::AbstractMetric, d::AbstractAccretionGeometry, x, ν, u)
    absorption_coefficient(m, d, x, ν), emissivity_coefficient(m, d, x, ν)
end

absorption_coefficient(m::AbstractMetric, d::AbstractAccretionGeometry, x, ν) = 0.0
emissivity_coefficient(m::AbstractMetric, d::AbstractAccretionGeometry, x, ν) = 0.0

struct RadiativeTransferIntegrationParameters{M,V} <: AbstractIntegrationParameters{M}
    metric::M
    status::MutStatusCode
    within_geometry::V
    RadiativeTransferIntegrationParameters(
        metric::M,
        status,
        within_geometry::V,
    ) where {M,V} = new{M,V}(metric, MutStatusCode(status), within_geometry)
end

function _radiative_transfer_integration_parameters(
    metric::AbstractMetric,
    status::StatusCodes.T,
    geometry,
)
    within_geometry = map(!is_finite_disc, geometry)
    RadiativeTransferIntegrationParameters(metric, status, within_geometry)
end

set_status_code!(params::RadiativeTransferIntegrationParameters, status::StatusCodes.T) =
    params.status[1] = status
get_status_code(params::RadiativeTransferIntegrationParameters) = params.status[1]
get_metric(params::RadiativeTransferIntegrationParameters) = params.metric


function update_integration_parameters!(
    p::RadiativeTransferIntegrationParameters,
    new::RadiativeTransferIntegrationParameters,
)
    set_status_code!(p, get_status_code(new))
    p.within_geometry .= new.within_geometry
    p
end

function _intensity_delta(
    m,
    x,
    k::AbstractArray{T},
    geometry::CompositeGeometry{K,D},
    within,
    I,
    ν₀,
    r_isco,
    λ,
) where {T,K,D}
    indexed_geom = (enumerate(geometry.geometry)...,)
    sum(indexed_geom) do args
        i, geom = args
        _intensity_delta(
            m,
            x,
            k,
            geom,
            within,
            I,
            ν₀,
            r_isco,
            λ;
            index = i,
        )
    end
end

function _intensity_delta(
    m,
    x,
    k::AbstractArray{T},
    geometry::AbstractAccretionGeometry,
    within,
    I,
    ν₀,
    r_isco,
    λ;
    index = 1,
) where {T}
    if within[index]
        radiative_transfer(m, x, k, geometry, I, ν₀, r_isco, λ)
    else
        zero(T)
    end
end

function intersection_callbacks(
    g::AbstractAccretionGeometry,
    ::TraceRadiativeTransfer,
    i::Int;
    gtol,
)
    # callbacks should only be defined on optically thick geometries
    if is_optically_thin(g)
        return _intersection_condition(g; gtol = gtol), _intersection_affect(g)
    end

    function _radiative_transfer_condition(u, λ, integrator)
        dist = distance_to_disc(g, u; gtol = gtol)
        integrator.p.within_geometry[i] ? -dist : dist
    end
    function _radiative_transfer_affect!(integrator)
        # flip
        integrator.p.within_geometry[i] = !integrator.p.within_geometry[i]
    end

    _radiative_transfer_condition, _radiative_transfer_affect!
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
    r_isco = Gradus.isco(m)

    function f(u::SVector{9,T}, p, λ) where {T}
        @inbounds @views let x = SVector{4,T}(u[1:4]), k = SVector{4,T}(u[5:8]), I = u[9]
            _m = get_metric(p)
            dk = SVector{4,T}(geodesic_equation(_m, x, k))
            dI = _intensity_delta(
                _m,
                x,
                k,
                geometry,
                p.within_geometry,
                I,
                trace.ν,
                r_isco,
                λ,
            )
            vcat(k, dk, SVector(dI))
        end
    end

    # add our new quantity that we want to trace
    u_init = vcat(pos, vel, SVector(trace.I₀))
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        _radiative_transfer_integration_parameters(m, StatusCodes.NoStatus, geometry);
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
