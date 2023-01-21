function tracegeodesics(
    m::AbstractMetricParams,
    model::AbstractCoronaModel,
    time_domain::NTuple{2},
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
)
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, time_domain; kwargs...)
end

function tracegeodesics(
    m::AbstractMetricParams,
    model::AbstractCoronaModel,
    d::AbstractAccretionGeometry,
    time_domain::NTuple{2},
    ;
    n_samples = 1024,
    sampler = WeierstrassSampler(res = 100.0),
    kwargs...,
)
    us = sample_position(m, model, n_samples)
    vs = sample_velocity(m, model, sampler, us, n_samples)
    tracegeodesics(m, us, vs, d, time_domain; kwargs...)
end

struct VoronoiDiscProfile{D,V,G} <: AbstractDiscProfile
    disc::D
    polys::Vector{Vector{V}}
    generators::Vector{V}
    geodesic_points::Vector{G}

    function VoronoiDiscProfile(
        d::D,
        polys::Vector{Vector{V}},
        gen::Vector{V},
        gps::Vector{G},
    ) where {D<:AbstractAccretionDisc,V<:AbstractArray,G}
        if !isapprox(d.inclination, π / 2)
            return error(
                "Currently only supported for discs in the equitorial plane (θ=π/2).",
            )
        end
        new{D,V,G}(d, polys, gen, gps)
    end
end

function Base.show(io::IO, vdp::VoronoiDiscProfile{D}) where {D}
    write(io, "VoronoiDiscProfile for $D with $(length(vdp.generators)) generators")
end

function VoronoiDiscProfile(
    m::AbstractMetricParams{T},
    d::AbstractAccretionDisc{T},
    endpoints::AbstractVector{<:GeodesicPoint{T}};
    padding = 1,
) where {T}
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

    VoronoiDiscProfile(d, polys_vecs, generators, endpoints)
end

function VoronoiDiscProfile(
    m::AbstractMetricParams,
    d::AbstractAccretionDisc,
    sols::AbstractArray{S},
) where {S<:SciMLBase.AbstractODESolution}
    VoronoiDiscProfile(m, d, map(sol -> process_solution(m, sol), sols))
end

function VoronoiDiscProfile(
    m::AbstractMetricParams,
    d::AbstractAccretionDisc,
    simsols::SciMLBase.EnsembleSolution,
)
    VoronoiDiscProfile(
        m,
        d,
        filter(i -> i.prob.p.status == StatusCodes.IntersectedWithGeometry, simsols.u),
    )
end

@noinline function findindex(vdp::VoronoiDiscProfile, p)
    for (i, poly) in enumerate(vdp.polys)
        # check if we're at all close
        let p1 = poly[1]
            # todo: what's going on here?
            d = @inbounds (p1[1] - p[1])^2 + (p1[2] - p[2])^2

            if inpolygon(poly, p)
                return i
            end
        end
    end
    -1
end

@inline function findindex(vdp::VoronoiDiscProfile, gp::GeodesicPoint)
    p = to_cartesian(gp)
    findindex(vdp, p)
end

function findindex(vdp::VoronoiDiscProfile, gps::AbstractArray{<:GeodesicPoint}; kwargs...)
    ret = fill(-1, size(gps))
    Threads.@threads for i = 1:length(gps)
        gp = gps[i]
        ret[i] = findindex(vdp, gp)
    end
    ret
end

getareas(vdp::VoronoiDiscProfile) = getarea.(vdp.polys)

function getproperarea(poly::AbstractArray, m::AbstractMetricParams)
    A = getarea(poly)
    c = getbarycenter(poly)
    # get value of metric at the radius of the polygon's barycenter, and in the equitorial plane
    m_params = metric_components(m, (sqrt(c[1]^2 + c[2]^2), π / 2))
    # need radial and azimuthal components of the metric
    sqrt(m_params[2] * m_params[4]) * A
end

getproperarea(vdp::VoronoiDiscProfile, m::AbstractMetricParams) =
    map(p -> getproperarea(p, m), vdp.polys)

function unpack_polys(
    polys::AbstractVector{GeometryBasics.Polygon{2,T,GeometryBasics.Point2{T},L,V}},
) where {T,L,V}
    map(polys) do poly
        map(SVector{2,T}, getpoints(poly))
    end
end

struct SurrogateDiscProfile{D,S}
    disc::D
    surrogate::S
end

function SurrogateDiscProfile(
    disc::D,
    x::AbstractArray{<:SVector{N}},
    y;
    S = Surrogates.RadialBasis,
) where {D,N}
    bounds = map(1:N) do i
        extrema(el -> el[i], x)
    end
    lower_bounds = first.(bounds)
    upper_bounds = last.(bounds)
    surrogate = S(x, y, lower_bounds, upper_bounds)
    SurrogateDiscProfile{D,typeof(surrogate)}(disc, surrogate)
end

function evaluate(sdp::SurrogateDiscProfile, x::AbstractArray{<:AbstractArray})
    map(i -> evaluate(sdp, i), x)
end

function evaluate(sdp::SurrogateDiscProfile, x)
    sdp.surrogate(x)
end

function Base.show(io::IO, ::MIME"text/plain", sdp::SurrogateDiscProfile{D,S}) where {D,S}
    print(
        io,
        "SurrogateDiscProfile\n - disc      : $D\n - surrogate : $(Base.typename(S).name) ($(length(sdp.surrogate.x)) points)",
    )
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
    SVector{2}(xs / n, ys / n)
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

function _circ_path_intersect(radius, lines::AbstractArray{<:GeometryBasics.Line})
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

export VoronoiDiscProfile,
    getareas, getproperarea, getbarycenter, findindex, SurrogateDiscProfile, evaluate
