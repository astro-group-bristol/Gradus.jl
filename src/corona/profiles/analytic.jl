
struct AnalyticRadialDiscProfile{E,T} <: AbstractDiscProfile
    ε::E
    t::T

    function AnalyticRadialDiscProfile(emissivity, time)
        new{typeof(emissivity),typeof(time)}(emissivity, time)
    end
end

function AnalyticRadialDiscProfile(emissivity, cg::CoronaGeodesics)
    J = sortperm(cg.geodesic_points; by = i -> _equatorial_project(i.x))
    radii = @views [_equatorial_project(i.x) for i in cg.geodesic_points[J]]
    times = @views [i.x[1] for i in cg.geodesic_points[J]]
    t = _make_interpolation(radii, times)
    AnalyticRadialDiscProfile(emissivity, t)
end

emissivity_at(prof::AnalyticRadialDiscProfile, r::Number) = prof.ε(r)
emissivity_at(prof::AnalyticRadialDiscProfile, gp::AbstractGeodesicPoint) =
    emissivity_at(prof, _equatorial_project(gp.x))

coordtime_at(prof::AnalyticRadialDiscProfile, r::Number) = prof.t(r)
coordtime_at(prof::AnalyticRadialDiscProfile, gp::AbstractGeodesicPoint) =
    coordtime_at(prof, _equatorial_project(gp.x)) + gp.x[1]
