module GradusBase

import Parameters: @with_kw
import SciMLBase
import Base
import StaticArrays: SVector

include("metric-params.jl")
include("physical-quantities.jl")
include("geodesic-solutions.jl")
include("geometry.jl")

export AbstractMetricParams, metric_params, metric, getgeodesicpoint
GeodesicPoint,
vector_to_local_sky,
AbstractMetricParams,
geodesic_eq,
geodesic_eq!,
constrain,
inner_radius,
metric_type

end # module
