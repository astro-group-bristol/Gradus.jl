module Rendering

import Base.Threads: @threads

using Parameters

import SciMLBase

using ..GradusBase
using ..GeodesicTracer

include("utility.jl")
include("render-cache.jl")
include("point-functions.jl")
include("render.jl")


function rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    vf = ConstValueFunctions.shadow,
    kwargs...,
) where {T}
    __rendergeodesics(
        m,
        init_pos;
        image_width = 350,
        image_height = 250,
        fov_factor = 3.0,
        max_time = convert(T, max_time),
        vf = vf,
        kwargs...,
    )
end

function prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    cache = DefaultCache(),
    kwargs...,
) where {T}
    __prerendergeodesics(
        m,
        init_pos,
        cache;
        image_width = 350,
        image_height = 250,
        fov_factor = 3.0,
        max_time = convert(T, max_time),
        kwargs...,
    )
end


export (
    rendergeodesics, prerendergeodesics,
    PointFunction, FilterPointFunction, ConstPointFunctions
)

end # module
