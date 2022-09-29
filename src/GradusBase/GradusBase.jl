module GradusBase

import Parameters: @with_kw
import SciMLBase
import Base
import StaticArrays: SVector, MMatrix, SMatrix, @SVector
import Tullio: @tullio

# for doc bindings
import ..Gradus

using DocStringExtensions

include("metric-params.jl")
include("physical-quantities.jl")
include("geodesic-solutions.jl")
include("geometry.jl")

export AbstractMetricParams, metric_params, metric, getgeodesicpoint
GeodesicPoint,
AbstractGeodesicPoint,
vector_to_local_sky,
AbstractMetricParams,
geodesic_eq,
geodesic_eq!,
constrain,
inner_radius,
metric_type,
metric_components,
inverse_metric_components

end # module
