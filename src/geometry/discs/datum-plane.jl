struct DatumPlane{T} <: AbstractAccretionDisc{T}
    height::T
end
optical_property(::Type{<:DatumPlane}) = OpticallyThin()

function distance_to_disc(d::DatumPlane{T}, x4; gtol) where {T}
    h = @inbounds let r = x4[2], θ = x4[3], ϕ = x4[4]
        r * cos(θ)
    end
    abs(h) - d.height - (gtol * x4[2])
end

function datumplane(disc::AbstractThickAccretionDisc, r::T) where {T}
    h = r_cross_section(disc, r)
    DatumPlane(h)
end
function datumplane(disc::AbstractThickAccretionDisc, x::SVector{4})
    h = cross_section(disc, x)
    DatumPlane(h)
end

export DatumPlane
