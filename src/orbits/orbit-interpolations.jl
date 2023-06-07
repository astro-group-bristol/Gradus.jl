struct PlungingInterpolation{M,_interp_type}
    m::M
    t::_interp_type
    r::_interp_type
    ϕ::_interp_type

    function PlungingInterpolation(m::M, sol) where {M<:AbstractMetric{T}} where {T}
        # sort relevant features
        I = sortperm(@view(sol[2, :]))[2:end]

        r = sol[2, :][I]

        vt = sol[5, :][I]
        vr = sol[6, :][I]
        vϕ = sol[8, :][I]

        rinterp = DataInterpolations.LinearInterpolation(vt, r)

        new{M,typeof(rinterp)}(
            m,
            rinterp,
            DataInterpolations.LinearInterpolation(vr, r),
            DataInterpolations.LinearInterpolation(vϕ, r),
        )
    end
end

function Base.show(io::IO, ::PlungingInterpolation{M}) where {M}
    write(io, "PlungingInterpolation for $M")
end

function (pintrp::PlungingInterpolation)(r)
    vt = pintrp.t(r)
    vr = pintrp.r(r)
    vϕ = pintrp.ϕ(r)
    SVector(vt, vr, 0, vϕ)
end

function interpolate_plunging_velocities(
    m::AbstractMetric{T};
    max_time = 50_000,
    contra_rotating = false,
    reltol = 1e-9,
    kwargs...,
) where {T}
    isco = Gradus.isco(m)

    # rule of thumb to achieve desired error
    δr = (reltol * 10)
    u = @SVector([0.0, isco - δr, deg2rad(90), 0.0])
    v = Gradus.CircularOrbits.plunging_fourvelocity(
        m,
        isco;
        contra_rotating = contra_rotating,
    )
    sol = Gradus.tracegeodesics(
        m,
        u,
        v,
        (0.0, T(max_time));
        μ = 1.0,
        reltol = reltol,
        # ensure we gets sufficiently close to the event horizon
        chart = chart_for_metric(m; closest_approach = 1.000001),
        kwargs...,
    )

    PlungingInterpolation(m, sol)
end

export PlungingInterpolation, interpolate_plunging_velocities
