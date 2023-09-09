struct DatumPlane{T} <: AbstractAccretionDisc{T}
    height::T
end
optical_property(::Type{<:DatumPlane}) = OpticallyThin()

function distance_to_disc(d::DatumPlane{T}, x4; gtol) where {T}
    h = x4[2] * cos(x4[3])
    abs(h) - d.height - (gtol * x4[2])
end

function datumplane(disc::AbstractThickAccretionDisc, ρ::T) where {T}
    h = cross_section(disc, ρ)
    DatumPlane(h)
end
datumplane(disc::AbstractThickAccretionDisc, x::SVector{4}) =
    datumplane(disc, _equatorial_project(x))

export DatumPlane
