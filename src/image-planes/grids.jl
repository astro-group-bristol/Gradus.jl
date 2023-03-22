module Grids
import ..Gradus
export AbstractImpactParameterGrid, GeometricGrid, InverseGrid, LinearGrid

abstract type AbstractImpactParameterGrid end

struct GeometricGrid <: AbstractImpactParameterGrid end

function (grid::GeometricGrid)(min, max, N)
    _geometric_grid(min, max, N)
end

function _geometric_grid(min, max, N)
    K = (max / min)^(1 / N)
    (min * K^(i - 1) for i = 1:N)
end

struct InverseGrid <: AbstractImpactParameterGrid end

function (grid::InverseGrid)(min, max, N)
    _inverse_grid(min, max, N)
end

function _inverse_grid(min, max, N)
    Iterators.reverse(inv(x) for x in range(1 / max, 1 / min, N))
end

struct LinearGrid <: AbstractImpactParameterGrid end

function (grid::LinearGrid)(min, max, N)
    range(min, max, N)
end

struct SinGrid <: AbstractImpactParameterGrid end

function (grid::SinGrid)(min, max, N)
    _sin_grid(min, max, N)
end

function _sin_grid(min, max, N)
    (((sin(p) + 1) / 2) * (max - min) + min for p in range(-π / 2, π / 2, N))
end

struct LogisticGrid{T} <: AbstractImpactParameterGrid
    k::T
end

function (grid::LogisticGrid)(min, max, N)
    _logistic_grid(min, max, N; k = grid.k)
end

function _logistic_grid(min, max, N; k = 0.5)
    f(x) = (max - min) * inv(1 + exp(-k * x)) + min
    (f(y) for y in range(-10, 10, N))
end

end # module

using .Grids

export GeometricGrid, InverseGrid, LinearGrid
