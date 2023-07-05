struct CunninghamTransferFunction{T}
    "``g^\\ast`` values."
    g✶::Vector{T}
    "Transfer function data."
    f::Vector{T}
    "Timing data."
    t::Vector{T}
    gmin::T
    gmax::T
    "Emission radius."
    rₑ::T
end

struct InterpolatedCunninghamTransferFunction{T,U,L}
    upper_f::U
    lower_f::L
    upper_t::L
    lower_t::U
    gmin::T
    gmax::T
    rₑ::T
end

struct LagTransferFunction{T,X,E,P}
    max_t::T
    x::X
    image_plane_areas::Vector{T}
    emissivity_profile::E
    observer_to_disc::Vector{P}
end

function Base.show(io::IO, ::MIME"text/plain", tf::LagTransferFunction)
    text = """LagTransferFunction for $(typeof(tf.emissivity_profile.metric)) 
      . observer position      
          $(tf.x)
      . model                         : $(typeof(tf.emissivity_profile.model))
      . observer to disc photon count : $(length(tf.observer_to_disc))
      . source to disc photon count   : $(length(tf.emissivity_profile.geodesic_points))
      Total memory: $(Base.format_bytes(Base.summarysize(tf)))
    """
    print(io, text)
end
