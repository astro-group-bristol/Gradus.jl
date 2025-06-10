struct CoronaGridValues{T}
    "Time"
    t::T
    "Radius on the disc"
    r::T
    "Radial angle on the disc"
    ϕ::T
    "Redshift"
    g::T
    "|∂(theta, phi) / ∂(r, ϕ)| * sin(theta)"
    J::T
end

_to_vector(v::CoronaGridValues) = SVector{5}(v.t, v.r, v.ϕ, v.g, v.J)
_from_vector(v::AbstractVector) = CoronaGridValues(v[1], v[2], v[3], v[4], v[5])

make_null(::Type{T}) where {T<:CoronaGridValues} = T(NaN, NaN, NaN, NaN, NaN)

vector_average(
    weights::AbstractVector{<:Number},
    values::AbstractVector{<:CoronaGridValues},
) = _from_vector(sum(i -> i[1] * _to_vector(i[2]), zip(weights, values)))

function _make_emissivity_tracer(
    m::AbstractMetric,
    corona::AbstractCoronaModel,
    d::AbstractAccretionDisc;
    t_max = 1_000_000.0,
    kwargs...,
)
    x_src, v_src = sample_position_velocity(m, corona)

    function trace_angles(th, ph)
        v = sky_angles_to_velocity(m, x_src, v_src, th, ph)
        # TODO: make this a reusable integrator
        sol = tracegeodesics(
            m,
            x_src,
            v,
            d,
            t_max;
            save_on = false,
            callback = domain_upper_hemisphere(),
            chart = chart_for_metric(m, 2*t_max),
            kwargs...,
        )
        unpack_solution(sol)
    end

    # function for obtaining keplerian velocities
    _disc_velocity = _keplerian_velocity_projector(m, d)

    function trace_jacobian(th::T, ph::T) where {T}
        function _f(x_)
            gp = trace_angles(x_[1], x_[2])
            if gp.status != StatusCodes.IntersectedWithGeometry
                SVector{4,eltype(x_)}(NaN, NaN, NaN, NaN)
            else
                v_disc = _disc_velocity(gp.x)
                _redshift = energy_ratio(m, gp, v_src, v_disc)
                r = _equatorial_project(gp.x)
                SVector{4,eltype(x_)}(r, gp.x[4], gp.x[1], _redshift)
            end
        end

        x0 = SVector(th, ph)

        _Tag = typeof(ForwardDiff.Tag(_f, T))
        ydual = _static_dual_eval(_Tag, _f, x0)

        res = ForwardDiff.value.(_Tag, ydual)
        jac = _extract_jacobian(_Tag, SVector{2}(ydual[1], ydual[2]), x0)

        CoronaGridValues{T}(res[3], res[1], res[2], res[4], abs(inv(det(jac))) * sin(th))
    end
end

function check_refine(
    sky::AdaptiveSky{T,<:CoronaGridValues},
    i1::Int,
    i2::Int;
    percentage = 2,
) where {T}
    v1 = sky.values[i1]
    v2 = sky.values[i2]

    if isnan(v1.r) && isnan(v2.r)
        return false
    end

    g_too_coarse = !isapprox(v1.g, v2.g, atol = percentage / 100)
    J_too_coarse = !isapprox(v1.J, v2.J, atol = percentage / 100)

    g_too_coarse || J_too_coarse
end

function AdaptiveSky(
    m::AbstractMetric{T},
    corona::AbstractCoronaModel,
    d::AbstractAccretionDisc;
    kwargs...,
) where {T}
    AdaptiveSky(
        CoronaGridValues{T},
        _make_emissivity_tracer(m, corona, d; kwargs...),
        check_refine,
    )
end

