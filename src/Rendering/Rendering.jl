module Rendering

import Base.Threads: @threads
import ThreadsX

using Parameters
using ProgressMeter

import SciMLBase

using ..GradusBase
using ..GeodesicTracer

# for doc bindings
import ..Gradus

using DocStringExtensions

include("utility.jl")
include("render-cache.jl")
include("point-functions.jl")
include("render.jl")


function rendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    pf = PointFunction((m, gp, mt) -> gp.t2) ∘
         FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t2 < max_time, NaN),
    image_width = 350,
    image_height = 250,
    fov_factor = 3.0,
    kwargs...,
) where {T}
    image = zeros(T, (image_height, image_width))
    render_into_image!(
        image,
        m,
        init_pos,
        max_time;
        pf = pf,
        image_width = image_width,
        image_height = image_height,
        fov_factor = fov_factor,
        kwargs...,
    )
    image
end

function prerendergeodesics(
    m::AbstractMetricParams{T},
    init_pos,
    max_time;
    cache = EndpointCache(),
    image_width = 350,
    image_height = 250,
    fov_factor = 3.0,
    kwargs...,
) where {T}
    __prerendergeodesics(
        m,
        init_pos,
        max_time,
        cache;
        image_width = image_width,
        image_height = image_height,
        fov_factor = fov_factor,
        kwargs...,
    )
end


export rendergeodesics,
    prerendergeodesics, AbstractPointFunction, PointFunction, FilterPointFunction, apply


end # module
