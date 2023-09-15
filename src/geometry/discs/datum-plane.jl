struct DatumPlane{T} <: AbstractAccretionDisc{T}
    height::T
end
optical_property(::Type{<:DatumPlane}) = OpticallyThin()

function distance_to_disc(d::DatumPlane{T}, x4; gtol) where {T}
    # no abs on h since datum planes have no underside
    h = _spinaxis_project(x4, signed = true)
    h - d.height
end

inner_radius(d::DatumPlane) = 0

function datumplane(disc::AbstractThickAccretionDisc, ρ::T) where {T}
    h = cross_section(disc, ρ)
    DatumPlane(h)
end
datumplane(disc::AbstractThickAccretionDisc, x::SVector{4}) =
    datumplane(disc, _equatorial_project(x))

export DatumPlane
