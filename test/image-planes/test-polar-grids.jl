using Test
using Gradus
using StaticArrays

m = BoyerLindquistAD()
u = @SVector [1.0, 1e3, π / 2, 0.0]

function count_inner_boundary(m, simsols)
    points = getgeodesicpoint.(m, simsols.u)
    count(i -> i.status == StatusCodes.WithinInnerBoundary, points)
end

# check the trace bootstrap works for each grid type
plane = PolarPlane(LinearGrid(), Nr = 10, Nθ = 10)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 10

plane = PolarPlane(GeometricGrid(), Nr = 10, Nθ = 10)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 30

plane = PolarPlane(InverseGrid(), Nr = 10, Nθ = 10)
simsols = tracegeodesics(m, u, plane, (0.0, 2000.0))
@test count_inner_boundary(m, simsols) == 80
