using Test
using Gradus

v = zeros(Float64, (3, 2, 2))

@test size(Gradus._get_dim_slice(v, Val{2}())) == (3,)

cache = Gradus.InterpolationCache{2}(v)
@test Gradus._make_cache_slices(cache) ==
      @views (cache.cache[:, :, 1], cache.cache[:, 1, 1])
@inferred Gradus._make_cache_slices(cache)

X1 = [0.0, 1.0]
X2 = [1.0, 2.0]
vals = [[[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]], [[-1.0, -1.0, -1.0], [-2.0, -2.0, -2.0]]]
vals = reshape(reduce(vcat, reduce(vcat, vals)), (3, 2, 2))

cache = Gradus.InterpolationCache{2}(vals)
Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

@test Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.0)) == [1.0, 1.0, 1.0]
@test Gradus.interpolate!(cache, (X1, X2), vals, (1.0, 1.5)) == [-1.5, -1.5, -1.5]
@test Gradus.interpolate!(cache, (X1, X2), vals, (0.5, 1.5)) == [0.0, 0.0, 0.0]

@inferred Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

@allocated Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))

# let's try it with arbitrary data structures

struct Thing
    a::Vector{Float64}
    b::Vector{Float64}
end
Thing(a::Number, b::Number) = Thing([a], [b])

function Gradus._set_value!(out::Thing, v::Thing)
    out.a[1] = v.a[1]
    out.b[1] = v.b[1]
end
function Gradus._linear_interpolate!(out::Thing, y1::Thing, y2::Thing, θ)
    Gradus._linear_interpolate!(out.a, y1.a, y2.a, θ)
    Gradus._linear_interpolate!(out.b, y1.b, y2.b, θ)
end

vals = reshape(
    [Thing(1.0, 1.0), Thing(2.0, 2.0), Thing(-1.0, -1.0), Thing(-2.0, -2.0)],
    (2, 2),
)

cache = Gradus.InterpolationCache{2}(vals)
intp = Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))
# do it twice to make sure no side-effects
intp = Gradus.interpolate!(cache, (X1, X2), vals, (0.0, 1.5))
@test intp.a == [1.5]
