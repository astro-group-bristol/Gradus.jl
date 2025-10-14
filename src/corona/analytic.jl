
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

function emissivity_at(prof::AnalyticRadialDiscProfile, r::Number)
    r_bounded = _enforce_interpolation_bounds(r, prof)
    prof.ε(r)
end
emissivity_at(prof::AnalyticRadialDiscProfile, gp::AbstractGeodesicPoint) =
    emissivity_at(prof, _equatorial_project(gp.x))

function coordtime_at(prof::AnalyticRadialDiscProfile, r::Number)
    r_bounded = _enforce_interpolation_bounds(r, prof)
    prof.t(r_bounded)
end
coordtime_at(prof::AnalyticRadialDiscProfile, gp::AbstractGeodesicPoint) =
    coordtime_at(prof, _equatorial_project(gp.x)) + gp.x[1]

function _enforce_interpolation_bounds(r::Number, prof::AnalyticRadialDiscProfile)
    r_min = first(prof.t.t)
    r_max = last(prof.t.t)
    _enforce_interpolation_bounds(r, r_min, r_max)
end
