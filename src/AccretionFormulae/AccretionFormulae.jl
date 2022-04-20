module AccretionFormulae

import ..Gradus
import ..GradusBase: AbstractMetricParams, metric

using ..FirstOrderMethods: FirstOrderGeodesicPoint
using ..Rendering: PointFunction

include("redshift.jl")

const redshift = PointFunction(_redshift_guard)

export redshift

end # module
