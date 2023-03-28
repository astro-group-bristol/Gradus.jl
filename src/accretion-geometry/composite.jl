
struct CompositeGeometry{T,G} <: AbstractAccretionGeometry{T}
    geometry::G
    function CompositeGeometry(g::Vararg{<:AbstractAccretionGeometry{T}}) where {T}
        new{T,typeof(g)}(g)
    end
end

function Base.show(io::IO, ::MIME"text/plain", cg::CompositeGeometry)
    buf = IOBuffer()
    println(buf, "CompositeGeometry:")
    for g in cg.geometry
        println(buf, "  - ", typeof(g).name.name)
    end
    print(io, String(take!(buf)))
end
