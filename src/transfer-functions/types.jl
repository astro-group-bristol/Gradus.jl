struct CunninghamTransferFunction{T}
    "``g^\\ast`` values."
    g✶::Vector{T}
    "Transfer function data."
    f::Vector{T}
    gmin::T
    gmax::T
    "Emission radius."
    rₑ::T
end

struct InterpolatedCunninghamTransferFunction{T,U,L}
    upper_f::U
    lower_f::L
    gmin::T
    gmax::T
    rₑ::T
end

struct LagTransferFunction{T,M,U,P}
    max_t::T
    metric::M
    u::U
    image_plane_areas::Vector{T}
    source_to_disc::Vector{P}
    observer_to_disc::Vector{P}
end

function Base.show(io::IO, ::MIME"text/plain", tf::LagTransferFunction)
    text = """LagTransferFunction for $(typeof(tf.metric)) 
      . observer position      
          $(tf.u)
      . observer to disc photon count : $(length(tf.observer_to_disc))
      . source to disc photon count   : $(length(tf.source_to_disc))
      Total memory: $(Base.format_bytes(Base.summarysize(tf)))
    """
    print(io, text)
end
