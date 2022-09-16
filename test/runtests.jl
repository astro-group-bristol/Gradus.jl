using Gradus, StaticArrays
using Test
using Aqua

include("smoke-tests/rendergeodesics.jl")
include("smoke-tests/tracegeodesics.jl")
include("smoke-tests/pointfunctions.jl")
include("smoke-tests/circular-orbits.jl")
include("smoke-tests/disc-profiles.jl")
include("smoke-tests/special-radii.jl")
include("smoke-tests/cunningham-transfer-functions.jl")

# little bit of aqua
Aqua.test_undefined_exports(Gradus)
Aqua.test_unbound_args(Gradus)
Aqua.test_stale_deps(Gradus)
