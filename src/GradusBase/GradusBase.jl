module GradusBase

import Base
import SciMLBase

using Parameters: @with_kw
using StaticArrays: SVector, MMatrix, SMatrix, @SVector
using Tullio: @tullio
using LinearAlgebra: norm, inv
using EnumX

# for doc bindings
import ..Gradus

using DocStringExtensions

@enumx StatusCodes begin
    OutOfDomain
    WithinInnerBoundary
    IntersectedWithGeometry
    NoStatus
end

include("abstract-metric.jl")
include("geodesic-solutions.jl")
include("geometry.jl")
include("physical-quantities.jl")

export AbstractMetric,
    AbstractCoordinates,
    AbstractTrace,
    BoyerLindquist,
    AbstractStaticAxisSymmetric,
    metric,
    unpack_solution,
    unpack_solution_full,
    GeodesicPoint,
    AbstractGeodesicPoint,
    vector_to_local_sky,
    AbstractMetric,
    geodesic_equation,
    constrain,
    inner_radius,
    metric_type,
    metric_components,
    inverse_metric_components,
    dotproduct,
    propernorm,
    tetradframe,
    lnrbasis,
    lnrframe,
    lowerindices,
    raiseindices,
    StatusCodes,
    AbstractIntegrationParameters,
    IntegrationParameters,
    update_integration_parameters!,
    restrict_ensemble

end # module
