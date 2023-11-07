struct CunninghamTransferData{T}
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

struct TransferBranches{T,F}
    upper_f::F
    lower_f::F
    upper_t::F
    lower_t::F
    gmin::T
    gmax::T
    rₑ::T
end

struct InterpolatingTransferBranches{T,F}
    branches::Vector{TransferBranches{T,F}}
    radii::Vector{T}
    gmin::Vector{T}
    gmax::Vector{T}
end

struct LagTransferFunction{T,X,E,P}
    max_t::T
    x::X
    image_plane_areas::Vector{T}
    coronal_geodesics::E
    observer_to_disc::Vector{P}
end

function Base.show(io::IO, ::MIME"text/plain", tf::LagTransferFunction)
    text = """LagTransferFunction for $(typeof(tf.coronal_geodesics.metric)) 
      . observer position      
          $(tf.x)
      . model                         : $(typeof(tf.coronal_geodesics.model))
      . observer to disc photon count : $(length(tf.observer_to_disc))
      . source to disc photon count   : $(length(tf.coronal_geodesics.geodesic_points))
      Total memory: $(Base.format_bytes(Base.summarysize(tf)))
    """
    print(io, text)
end
