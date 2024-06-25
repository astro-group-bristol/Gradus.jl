using Test
using Gradus

v = zeros(Float64, (3, 2, 2))

@inferred MultilinearInterpolator{2}(v)

X1 = [0.0, 1.0]
X2 = [1.0, 2.0]
vals = reshape(
    NTuple{3,Float64}[
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, -1.0),
        (2.0, 2.0, 2.0),
        (-2.0, -2.0, -2.0),
    ],
    (2, 2),
)

cache = MultilinearInterpolator{2}(vals)
interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

@test Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.0)) == (1.0, 1.0, 1.0)
@test Gradus.interpolate!(cache, (X1, X2), vals, (1.0, 1.5)) == (-1.5, -1.5, -1.5)
@test Gradus.interpolate!(cache, (X1, X2), vals, (0.5, 1.5)) == (0.0, 0.0, 0.0)

@inferred Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

@allocated Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

ff(x, y) = 3x^2 + x * y - sin(y)
ff(x) = SVector(ff(x[1], x[2]))
X1 = collect(range(0, 1, 1000))
X2 = collect(range(0, 1, 1000))
vals = reshape([ff(x, y) for y in X1, x in X2], (length(X1), length(X2)))
cache = MultilinearInterpolator{2}(vals)

# check that dual cache works too
function _interpolate_wrapper(cache, X1, X2, vals)
    function _f(x)
        x2, x1 = x
        SVector(interpolate!(cache, (X1, X2), vals, (x1, x2)))
    end
end
interp_f = _interpolate_wrapper(cache, X1, X2, vals)
x0 = SVector(X1[3], X2[5])
@test Gradus.ForwardDiff.jacobian(interp_f, x0) ≈ Gradus.ForwardDiff.jacobian(ff, x0) atol =
    1e-2

# now single dimension edge case
X1 = [0.0, 1.0]
vals = [-1.0, 0.0]
cache = MultilinearInterpolator{1}(vals)

@test Gradus.interpolate!(cache, (X1,), vals, (0.0,)) == -1.0
@test Gradus.interpolate!(cache, (X1,), vals, (0.6,)) == -0.4

# check that dual cache works too
function _interpolate_wrapper(cache, X1, vals)
    function _f(x)
        interpolate!(cache, (X1,), vals, (x,))
    end
end
interp_f = _interpolate_wrapper(cache, X1, vals)
@test Gradus.ForwardDiff.derivative(interp_f, 0.6) == 1
@inferred Gradus.ForwardDiff.derivative(interp_f, 0.6)

# let's try it with arbitrary data structures
struct Thing{V<:AbstractVector,M<:AbstractMatrix}
    a::V
    b::M
end
Thing(a::Number, b::Number) = Thing([a], [b 0; 0 b])

function Gradus.restructure(data::Thing, vals::AbstractVector)
    @views Thing(vals[1:length(data.a)], reshape(vals[length(data.a)+1:end], size(data.b)))
end

X1 = [0.0, 1.0]
X2 = [1.0, 2.0]
vals = reshape(
    [Thing(1.0, 1.0), Thing(-1.0, -1.0), Thing(2.0, 2.0), Thing(-2.0, -2.0)],
    (2, 2),
)

cache = MultilinearInterpolator{2}(vals)
intp = Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))
# do it twice to make sure no side-effects
intp = Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))
@test intp.a == [1.5]
@test intp.b == [1.5 0; 0 1.5]
