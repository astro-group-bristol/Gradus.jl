abstract type AbstractCacheStrategy end

struct SolutionCache <: AbstractCacheStrategy end
struct EndpointCache <: AbstractCacheStrategy end

abstract type AbstractRenderCache{M,T} end

struct SolutionRenderCache{M,T,G} <: AbstractRenderCache{M,T}
    # metric
    m::M

    # max time
    max_time::T

    # size information
    height::Int
    width::Int

    # geodesics themselves in 2d array
    geodesics::AbstractMatrix{G}

    function SolutionRenderCache(
        m::AbstractMetricParams{T},
        max_time::T,
        height,
        width,
        cache::AbstractVector{SciMLBase.EnsembleSolution{T,N,Vector{O}}},
    ) where {T,N,O}

        geodesics = Matrix{O}(undef, (height, width))

        # populate store
        for (col, simsol) in enumerate(cache)
            for (row, sol) in enumerate(simsol)
                geodesics[col, row] = sol
            end
        end

        # return instance 
        new{typeof(m),T,O}(m, max_time, height, width, geodesics)
    end
end

struct EndpointRenderCache{M,T,P} <: AbstractRenderCache{M,T}
    # metric
    m::M

    # max time
    max_time::T

    # size information
    height::Int
    width::Int

    points::Matrix{P}
end

function Base.show(io::IO, cache::AbstractRenderCache{M}) where {M}
    repr = "$(Base.typename(typeof(cache)).name){$M} (dimensions $(cache.height)x$(cache.width))"
    write(io, repr)
end
