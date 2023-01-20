module Grids
import ..Gradus
export AbstractImpactParameterGrid, GeometricGrid, InverseGrid, LinearGrid

abstract type AbstractImpactParameterGrid end

struct GeometricGrid <: AbstractImpactParameterGrid end

function (grid::GeometricGrid)(min::T, max::T, N) where {T}
    _geometric_grid(min, max, N)
end

function _geometric_grid(min, max, N)
    K = (max / min)^(1 / N)
    (min * K^(i - 1) for i = 1:N)
end

struct InverseGrid <: AbstractImpactParameterGrid end

function (grid::InverseGrid)(min::T, max::T, N) where {T}
    _inverse_grid(min, max, N)
end

function _inverse_grid(min, max, N)
    Iterators.reverse(inv(x) for x in range(1 / max, 1 / min, N))
end

struct LinearGrid <: AbstractImpactParameterGrid end

function (grid::LinearGrid)(min::T, max::T, N) where {T}
    range(min, max, N)
end

end # module

using .Grids

export GeometricGrid, InverseGrid, LinearGrid
