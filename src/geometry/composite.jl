function CompositeGeometry()
    error("Must provide at least one disc as argument to constructor")
end

function CompositeGeometry(g::Vararg{<:AbstractAccretionGeometry{T}}) where {T}
    CompositeGeometry{T,typeof(g)}(g)
end

Base.length(cg::CompositeGeometry) = length(cg.geometry)

function Base.show(io::IO, ::MIME"text/plain", cg::CompositeGeometry)
    buf = IOBuffer()
    println(buf, "CompositeGeometry:")
    for g in cg.geometry
        println(buf, "  - ", typeof(g).name.name)
    end
    print(io, String(take!(buf)))
end

Base.:∘(d1::AbstractAccretionGeometry, d2::AbstractAccretionGeometry) =
    CompositeGeometry(d1, d2)
