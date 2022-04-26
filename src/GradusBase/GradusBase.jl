module GradusBase

import Parameters: @with_kw
import SciMLBase
import Base

include("metric-params.jl")
include("physical-quantities.jl")
include("geodesic-solutions.jl")
include("geometry.jl")

export AbstractMetricParams,
    metric_params,
    metric,
    get_endpoint,
    get_startpoint,
    get_point,
    GeodesicPoint,
    vector_to_local_sky,
    AbstractMetricParams,
    geodesic_eq,
    geodesic_eq!,
    constrain,
    on_chart,
    inner_radius,
    metric_type

end # module
