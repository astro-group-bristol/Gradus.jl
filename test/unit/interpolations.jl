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
