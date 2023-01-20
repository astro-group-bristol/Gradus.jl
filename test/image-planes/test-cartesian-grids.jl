using Test
using Gradus
using StaticArrays

include("../utils.jl")

m = BoyerLindquistAD()
u = @SVector [1.0, 1e3, π / 2, 0.0]

# check the trace bootstrap works for each grid type
plane = CartesianPlane(LinearGrid(), x_min = 0.1, y_min = 0.1, Nx = 12, Ny = 12)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 1

plane = CartesianPlane(GeometricGrid(), x_min = 0.1, y_min = 0.1, Nx = 12, Ny = 12)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 45

plane = CartesianPlane(InverseGrid(), x_min = 0.1, y_min = 0.1, Nx = 12, Ny = 12)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 81
