abstract type AbstractAccretionGeometry{T} end
abstract type AbstractAccretionDisc{T} <: AbstractAccretionGeometry{T} end

function to_cartesian(vec::AbstractVector{T}) where {T}
    SVector{3,T}(to_cartesian(vec[2], vec[3], vec[4]))
end

function to_cartesian(r, ϕ, θ)
    sinϕ = sin(ϕ)
    (r * sinϕ * cos(θ), r * sinϕ * sin(θ), r * cos(ϕ))
end

"""
    to_cartesian(gp::GradusBase.AbstractGeodesicPoint{T})

Return the position of `gp` in Cartesian coordinates.
"""
function to_cartesian(gp::GradusBase.AbstractGeodesicPoint{T}) where {T}
    let r = gp.u[2], ϕ = gp.u[4]
        x = r * cos(ϕ)
        y = r * sin(ϕ)
        GeometryBasics.Point2{T}(x, y)
    end
end

function cartesian_line_element(u, integrator)
    (to_cartesian(integrator.uprev), to_cartesian(u))
end

function line_element(u, integrator)
    @inbounds (@view(integrator.uprev[2:4]), @view(u[2:4]))
end

# most of these functions are utility methods for GeometryBasics
# which may or may not be included. GeometryBasics seems to be about 
# to undergo major breaking changes, which would refine the package
# so may wait before using another libraries version of these:

@fastmath function getorientation(line::GeometryBasics.Line{2,T}, p = GeometryBasics.Point2{T}(0.0, 0.0)) where {T}
    let lp1 = line.points[1], lp2 = line.points[2]
        o = p - lp1
        b = lp1 - lp2
        t = (b[2] * o[1]) - (b[1] * o[2])
        # branchless : t < 0 ? 1 : -1
        (t < 0) * 2 - 1
    end
end

function getcycliclines(poly::GeometryBasics.Polygon{2, T, GeometryBasics.Point2{T}, L, V}) where {T,L,V}
    lines = poly.exterior.points
    vcat(lines, GeometryBasics.Line(lines[end].points[2], lines[1].points[1]))
end

function getpoints(poly::GeometryBasics.Polygon{2, T, GeometryBasics.Point2{T}, L, V}) where {T,L,V}
    lines = poly.exterior.points
    ps = first.(lines)
    push!(ps, lines[end][2])
    ps
end

function getcyclicpoints(poly::GeometryBasics.Polygon{2, T, GeometryBasics.Point2{T}, L, V}) where {T,L,V}
    points = getpoints(poly)
    push!(points, points[1])
    points
end

function inpolygon(poly::GeometryBasics.Polygon{2, T, GeometryBasics.Point2{T}, L, V}, p::GeometryBasics.Point2{T}) where {T,L,V}
    lines = getcycliclines(poly)
    side = getorientation(lines[1], p)
    for i = 2:length(lines)
        if getorientation(lines[i], p) != side
            return false
        end
    end
    true
end

# shoelace formula
function getarea(poly::GeometryBasics.Polygon{2, T, GeometryBasics.Point2{T}, L, V}) where {T,L,V}
    a = 0.0
    ppoints = getcyclicpoints(poly)
    @inbounds for i = 2:length(ppoints)
        p1 = ppoints[i-1]
        p2 = ppoints[i]
        # (x1 * y2) - (x2 * y1)
        a += (p1[1] * p2[2]) - (p2[1] * p1[2])
    end
    abs(a / 2)
end
