export VoronoiDiscProfile,
    getareas, getproperarea, getbarycenter, RadialDiscProfile, get_emissivity

# exported interface
function emitted_flux(profile::AbstractDiscProfile, gps)
    error("Not implemented for $(typeof(profile))")
end
function delay(profile::AbstractDiscProfile, gps)
    error("Not implemented for $(typeof(profile))")
end

# tuple so the calculation may be combined if desired
function delay_flux(profile::AbstractDiscProfile, gps)
    (delay(profile, gps), emitted_flux(profile, gps))
end

include("profiles/radial.jl")
include("profiles/voronoi.jl")
