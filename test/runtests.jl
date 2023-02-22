using Test
using Gradus

using Aqua

@testset "smoke-tests" verbose = true begin
    include("smoke-tests/rendergeodesics.jl")
    include("smoke-tests/tracegeodesics.jl")
    include("smoke-tests/pointfunctions.jl")
    include("smoke-tests/circular-orbits.jl")
    include("smoke-tests/disc-profiles.jl")
    include("smoke-tests/special-radii.jl")
    include("smoke-tests/cunningham-transfer-functions.jl")
end

@testset "metric-geometry" verbose = true begin
    include("unit/gradusbase.geometry.jl")
    include("test-special-radii.jl")
end

@testset "integration" verbose = true begin
    include("integration/test-inference.jl")
    include("integration/test-charts.jl")
end

@testset "transfer-functions" verbose = true begin
    include("transfer-functions/test-2d.jl")
end

@testset "image-planes" verbose = true begin
    include("image-planes/test-polar-grids.jl")
    include("image-planes/test-cartesian-grids.jl")
end

@testset "line-profiles" verbose = true begin
    include("line-profiles/test-cunningham.jl")
    include("line-profiles/test-binning.jl")
end

# little bit of aqua
Aqua.test_undefined_exports(Gradus)
Aqua.test_unbound_args(Gradus)
Aqua.test_stale_deps(Gradus)
