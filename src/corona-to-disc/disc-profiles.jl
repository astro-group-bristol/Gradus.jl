function tracegeodesics(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    time_domain::Tuple{T,T};
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, time_domain; kwargs...)
end

function tracegeodesics(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T},
    time_domain::Tuple{T,T},
    ;
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, d, time_domain; kwargs...)
end

function renderprofile(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T},
    max_time::T;
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
) where {T}
    __renderprofile(m, model, d, n_samples, (0.0, max_time); kwargs...)
end

"""
    AbstractDiscProfile

Abstract type for binning structures over discs (e.g., radial bins, voronoi).
"""
abstract type AbstractDiscProfile end

# evaluate(aap::AbstractAccretionProfile, p) =
#     error("Not implemented for `$(typeof(aap))` yet.")

struct VoronoiDiscProfile{D,V} <: AbstractDiscProfile
    disc::D
    polys::Vector{Vector{V}}
    generators::Vector{V}

    function VoronoiDiscProfile(
        d::D,
        polys::Vector{Vector{V}},
        gen::Vector{V},
    ) where {V<:AbstractArray{T}} where {D,T}
        if !isapprox(d.inclination, π / 2)
            return error(
                "Currently only supported for discs in the equitorial plane (θ=π/2).",
            )
        end
        new{D,V}(d, polys, gen)
    end
end

function Base.show(io::IO, vdp::VoronoiDiscProfile{D,T}) where {D,T}
    write(io, "VoronoiDiscProfile for $D with $(length(vdp.generators)) generators")
end

function VoronoiDiscProfile(
    m::AbstractMetricParams{T},
    d::AbstractAccretionDisc{T},
    endpoints::AbstractVector{GeodesicPoint{T,V}};
    padding = 1,
) where {T,V}
    dim = d.outer_radius + padding
    rect = VoronoiCells.Rectangle(
        GeometryBasics.Point2{T}(-dim, -dim),
        GeometryBasics.Point2{T}(dim, dim),
    )

    generators = to_cartesian.(endpoints)

    generators_points = convert.(GeometryBasics.Point2{T}, generators)

    tess = VoronoiCells.voronoicells(generators_points, rect)
    polys = GeometryBasics.Polygon.(tess.Cells)

    # cut the polygons down to shape
    cutchart!(polys, m, d)

    polys_vecs = unpack_polys(polys)

    VoronoiDiscProfile(d, polys_vecs, generators)
end

function VoronoiDiscProfile(
    m::AbstractMetricParams{T},
    d::AbstractAccretionDisc{T},
    sols::AbstractArray{S},
) where {T,S<:SciMLBase.AbstractODESolution}
    VoronoiDiscProfile(m, d, map(sol -> getgeodesicpoint(m, sol), sols))
end

function VoronoiDiscProfile(
    m::AbstractMetricParams{T},
    d::AbstractAccretionDisc{T},
    simsols::SciMLBase.EnsembleSolution{T},
) where {T}
    VoronoiDiscProfile(m, d, filter(i -> i.retcode == :Intersected, simsols.u))
end

@noinline function findindex(vdp::VoronoiDiscProfile{D,V}, p) where {D,V}
    for (i, poly) in enumerate(vdp.polys)
        # check if we're at all close
        let p1 = poly[1]
            d = @inbounds (p1[1] - p[1])^2 + (p1[2] - p[2])^2

            if inpolygon(poly, p)
                return i
            end
        end
    end
    -1
end

@inline function findindex(
    vdp::VoronoiDiscProfile{D,V},
    gp::GeodesicPoint{T,P},
) where {D,V,T,P}
    p = to_cartesian(gp)
    findindex(vdp, p)
end

function findindex(
    vdp::VoronoiDiscProfile{D,V},
    gps::AbstractArray{GeodesicPoint{T,P}};
    kwargs...,
) where {D,V,T,P}
    ret = fill(-1, size(gps))
    Threads.@threads for i = 1:length(gps)
        gp = gps[i]
        ret[i] = findindex(vdp, gp)
    end
    ret
end

getareas(vdp::VoronoiDiscProfile{D,T}) where {D,T} = getarea.(vdp.polys)

function getproperarea(
    poly::P,
    m::AbstractMetricParams{T},
) where {P<:AbstractArray{V}} where {V,T}
    A = getarea(poly)
    c = getbarycenter(poly)
    m_params = metric_components(m, (sqrt(c[1]^2 + c[2]^2), π / 2))
    sqrt(m_params[2] * m_params[3]) * A
end

getproperarea(vdp::VoronoiDiscProfile{D,K}, m::AbstractMetricParams{T}) where {D,K,T} =
    map(p -> getproperarea(p, m), vdp.polys)

function unpack_polys(
    polys::AbstractVector{GeometryBasics.Polygon{2,T,GeometryBasics.Point2{T},L,V}},
) where {T,L,V}
    map(polys) do poly
        map(SVector{2,T}, getpoints(poly))
    end
end

"""
    getbarycenter(poly)

Calculate the barycentric point of `poly`. Given that `poly` is an array of 2D points, this
method calculates

```math
\\vec{c} = \\frac{1}{N} \\sum_i^N \\vec{p}_i,
```

where ``\\vec{p}_i`` are the vectors to point ``i`` of the polygon.
"""
function getbarycenter(poly)
    xs = sum(i -> first(i), poly)
    ys = sum(i -> last(i), poly)
    n = length(poly)
    (xs / n, ys / n)
end


"""
    cutchart!(polys, m::AbstractMetricParams{T}, d::AbstractAccretionDisc{T})

Trims each polygon in `polys` to match the chart of the disc. Walks the edges of each
polygon until some intersection ``A`` is found. Continues walking until a second intersection
``B`` is found, and then interpolates an arc between ``A`` and ``B``.
"""
function cutchart!(polys, m::AbstractMetricParams{T}, d::AbstractAccretionDisc{T}) where {T}
    rs = inner_radius(m)
    inside_radius = rs < d.inner_radius ? d.inner_radius : rs
    outside_radius = d.outer_radius
    for i in eachindex(polys)
        poly = polys[i]
        poly = _cut_polygon(inside_radius, poly; inner = true)
        poly = _cut_polygon(outside_radius, poly; inner = false)
        polys[i] = poly
    end
    # filter!(!isnothing, polys)
end

function _cut_polygon(radius, poly; inner = true)
    clines = getcycliclines(poly)
    (i1, t1) = _circ_path_intersect(radius, clines)
    if t1 > 0
        (i2, t2) = _circ_path_intersect(radius, @view(clines[i1+1:end]))
        if inner
            _circ_cut(radius, clines, i1, t1, i2 + i1, t2; outer = false)
        else
            _circ_cut(radius, clines, i1, t1, i2 + i1, t2)
        end
    else
        poly
    end
end

function _circ_cut(radius, clines, i1, t1, i2, t2; outer = true)
    # get intersection points
    A = clines[i1][1] .+ t1 .* (clines[i1][2] - clines[i1][1])
    B = clines[i2][1] .+ t2 .* (clines[i2][2] - clines[i2][1])

    @assert length(clines) > 2

    # check if we started in or outside of the circle
    dist = √sum(clines[1][1] .^ 2)
    if (outer && (dist ≤ radius)) || (!outer && (dist ≥ radius))
        # inside circle
        points = first.(@view(clines[1:i1]))
        push!(points, A)
        push!(points, B)
        points = vcat(points, first.(@view(clines[i2+1:end])))

        return GeometryBasics.Polygon(points)
    else
        # outside circle
        points = first.(@view(clines[i1+1:i2]))
        push!(points, B)
        push!(points, A)

        return GeometryBasics.Polygon(points)
    end
end


function _circ_path_intersect(
    radius,
    lines::AbstractArray{L},
) where {L<:GeometryBasics.Line}
    for (i, line) in enumerate(lines)
        t = _circ_line_intersect(radius, line)
        if t > 0
            return i, t
        end
    end
    -1, -1
end

@fastmath function _circ_line_intersect(radius, line::GeometryBasics.Line)
    let A = line.points[1], B = line.points[2]
        d = B .- A
        d2 = sum(d .* d)
        A2 = sum(A .* A)
        da = sum(d .* A)
        Δ = da^2 - d2 * (A2 - radius^2)

        if Δ > 0
            tp = (-da + √Δ) / d2
            tn = (-da - √Δ) / d2

            if (1 ≥ tp ≥ 0)
                return tp
            end
            if (1 ≥ tn ≥ 0)
                return tn
            end
        end
    end
    -1
end


function __renderprofile(
    m::AbstractMetricParams{T},
    model::AbstractCoronaModel{T},
    d::AbstractAccretionGeometry{T},
    time_domain;
    sampler,
    kwargs...,
) where {T}
    # TODO
    error("Not implemented for $(typeof(d)).")
end
