function to_cartesian(vec::SVector{8,T}) where {T}
    SVector{3,T}(to_cartesian(vec[2], vec[3], vec[4]))
end

function to_cartesian(vec::SVector{4,T}) where {T}
    SVector{3,T}(to_cartesian(vec[2], vec[3], vec[4]))
end

function to_cartesian(vec::SVector{3,T}) where {T}
    SVector{3,T}(to_cartesian(vec[1], vec[2], vec[3]))
end

function to_cartesian(r, ϕ, θ)
    sinϕ = sin(ϕ)
    (r * sinϕ * cos(θ), r * sinϕ * sin(θ), r * cos(ϕ))
end

"""
    to_cartesian(gp::GradusBase.AbstractGeodesicPoint{T})

Return the end position of `gp` in Cartesian coordinates.
"""
function to_cartesian(gp::GradusBase.AbstractGeodesicPoint{T}) where {T}
    @inbounds let r = gp.u2[2], ϕ = gp.u2[4]
        x = r * cos(ϕ)
        y = r * sin(ϕ)
        SVector{2,T}(x, y)
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

getorientation(line::GeometryBasics.AbstractPolygon, p) = getorientation(line.points, p)

@fastmath function getorientation(line, p)
    let p1 = line[1], p2 = line[2]
        o = p - p1
        b = p1 - p2
        t = (b[2] * o[1]) - (b[1] * o[2])
        # branchless : t < 0 ? 1 : -1
        (t < 0) * 2 - 1
    end
end

function getcycliclines(
    poly::GeometryBasics.Polygon{2,T,GeometryBasics.Point2{T},L,V},
) where {T,L,V}
    lines = poly.exterior.points
    vcat(lines, GeometryBasics.Line(lines[end].points[2], lines[1].points[1]))
end

function getcycliclines(poly)
    (i == 1 ? (poly[end], poly[i]) : (poly[i-1], poly[i]) for i = 1:length(poly))
end

function getpoints(poly::GeometryBasics.Polygon)
    lines = poly.exterior.points
    ps = first.(lines)
    push!(ps, lines[end][2])
    ps
end

getpoints(poly) = (i for i in poly)

function getcyclicpoints(poly)
    pts = getpoints(poly)
    (i for i in Iterators.flatten((pts, (first(pts),))))
end

function inpolygon(poly, p)
    lines = getcycliclines(poly)
    side = getorientation(first(lines), p)
    for line in lines
        if getorientation(line, p) != side
            return false
        end
    end
    true
end

# shoelace formula
function getarea(poly)
    a = 0.0
    p1, ppoints = Iterators.peel(getcyclicpoints(poly))
    @inbounds for p2 in ppoints
        # (x1 * y2) - (x2 * y1)
        a += (p1[1] * p2[2]) - (p2[1] * p1[2])
        p1 = p2
    end
    abs(a / 2)
end

function getbarycenter(poly::GeometryBasics.Polygon)
    pts = getpoints(poly)
    n_points = length(pts)

    x_sum = 0.0
    y_sum = 0.0
    for pt in pts
        x_sum += pt[1]
        y_sum += pt[2]
    end

    @SVector [x_sum / n_points, y_sum / n_points]
end

export getpoints, getarea
