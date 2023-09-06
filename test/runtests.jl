using Test
using Gradus

using Aqua

@time @testset "smoke-tracing" verbose = true begin
    include("smoke-tests/tracegeodesics.jl")
    include("smoke-tests/rendergeodesics.jl")
    include("smoke-tests/prerendergeodesics.jl")
end

@time @testset "smoke-utility" verbose = true begin
    include("smoke-tests/pointfunctions.jl")
    include("smoke-tests/circular-orbits.jl")
    include("smoke-tests/disc-profiles.jl")
    include("smoke-tests/special-radii.jl")
    include("smoke-tests/cunningham-transfer-functions.jl")
end

@time @testset "smoke-reverberation" verbose = true begin
    include("smoke-tests/reverberation.jl")
end

@time @testset "metric-geometry" verbose = true begin
    include("unit/orthonormalization.jl")
    include("unit/metrics.kerr-newman.jl")
    include("test-special-radii.jl")
end

@time @testset "corona" verbose = true begin
    include("unit/coronal-beaming.jl")
    include("smoke-tests/coronal-spectra.jl")
    include("unit/emissivity.jl")
    include("disc-profiles/test-beamedpointsource.jl")
end

@time @testset "integration" verbose = true begin
    include("integration/test-inference.jl")
    include("integration/test-charts.jl")
    include("integration/test-precision.jl")
end

@time @testset "transfer-functions" verbose = true begin
    include("transfer-functions/test-2d.jl")
    include("transfer-functions/test-thick-disc.jl")
end

@time @testset "image-planes" verbose = true begin
    include("image-planes/test-polar-grids.jl")
    include("image-planes/test-cartesian-grids.jl")
end

@time @testset "line-profiles" verbose = true begin
    include("line-profiles/test-cunningham.jl")
    include("line-profiles/test-binning.jl")
end

@time @testset "geometry" verbose = true begin
    include("discs/test-polish-doughnut.jl")
end

# little bit of aqua
Aqua.test_undefined_exports(Gradus)
Aqua.test_unbound_args(Gradus)
Aqua.test_stale_deps(Gradus)
