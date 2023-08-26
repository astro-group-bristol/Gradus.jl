export emissivity_at, coordtime_at

emissivity_at(prof::AbstractDiscProfile, ::AbstractGeodesicPoint) =
    error("Not implemented for $(typeof(prof))")

coordtime_at(prof::AbstractDiscProfile, ::AbstractGeodesicPoint) =
    error("Not implemented for $(typeof(prof))")

emissivity_at(prof::AbstractDiscProfile, v::AbstractArray) =
    map(i -> emissivity_at(prof, i), v)
coordtime_at(prof::AbstractDiscProfile, v::AbstractArray) =
    map(i -> coordtime_at(prof, i), v)

include("profiles/radial.jl")
include("profiles/voronoi.jl")
include("profiles/analytic.jl")
