"""
    MeshAccretionGeometry(mesh)

$(FIELDS)
"""
struct MeshAccretionGeometry{T} <: AbstractAccretionGeometry{T}
    mesh::Vector{Tuple{SVector{3,T},SVector{3,T},SVector{3,T}}}
    x_extent::Tuple{T,T}
    y_extent::Tuple{T,T}
    z_extent::Tuple{T,T}
end

function MeshAccretionGeometry(mesh)
    static_mesh = map(mesh) do triangle
        Tuple(SVector(p[1], p[2], p[3]) for p in triangle)
    end
    MeshAccretionGeometry(static_mesh, bounding_box(mesh)...)
end

# naive implementation
function bounding_box(mesh::Union{GeometryBasics.Mesh{3,T}}) where {T}
    bounding_box(T, mesh)
end

function bounding_box(mesh::Vector{Tuple{SVector{3,T},SVector{3,T},SVector{3,T}}}) where {T}
    bounding_box(T, mesh)
end

function bounding_box(T, mesh)
    xmin = typemax(T)
    xmax = -typemax(T)
    ymin = typemax(T)
    ymax = -typemax(T)
    zmin = typemax(T)
    zmax = -typemax(T)
    for t in mesh
        for p in t
            let x = p[1], y = p[2], z = p[3]
                (x > xmax) && (xmax = x)
                (x < xmin) && (xmin = x)
                (y > ymax) && (ymax = y)
                (y < ymin) && (ymin = y)
                (z > zmax) && (zmax = z)
                (z < zmin) && (zmin = z)
            end
        end
    end
    (xmin, xmax), (ymin, ymax), (zmin, zmax)
end

function in_nearby_region(m::MeshAccretionGeometry{T}, line_element) where {T}
    p = line_element[2]
    @inbounds m.x_extent[1] < p[1] < m.x_extent[2] &&
              m.y_extent[1] < p[2] < m.y_extent[2] &&
              m.z_extent[1] < p[3] < m.z_extent[2]
end

function has_intersect(m::MeshAccretionGeometry{T}, line_element) where {T}
    for triangle in m.mesh
        dist_sq = sum(@.((triangle[1] - line_element[2])^2))
        # assume line element and mesh triangle are small; check if we're within a
        # certain distance before running jsr
        if dist_sq < 3.0 && jsf_algorithm(triangle..., line_element...)
            return true
        end
    end
    false
end

function collision_callback(m::MeshAccretionGeometry{T}) where {T}
    (u, λ, integrator) -> intersects_geometry(m, cartesian_line_element(u, integrator))
end
