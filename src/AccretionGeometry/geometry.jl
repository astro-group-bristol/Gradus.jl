abstract type AbstractAccretionGeometry{T} end
abstract type AbstractAccretionDisc{T} <: AbstractAccretionGeometry{T} end

function to_cartesian(vec::AbstractVector{T}) where {T}
    SVector{3,T}(to_cartesian(vec[2], vec[3], vec[4])...)
end

function to_cartesian(r, ϕ, θ)
    sinϕ = sin(ϕ)
    (r * sinϕ * cos(θ), r * sinϕ * sin(θ), r * cos(ϕ))
end

function cartesian_line_element(u::ArrayPartition{F,T}, integrator) where {F,T}
    (to_cartesian(integrator.uprev.x[2]), to_cartesian(u.x[2]))
end

function cartesian_line_element(u, integrator)
    (to_cartesian(integrator.uprev), to_cartesian(u))
end

function line_element(u::ArrayPartition{F,T}, integrator) where {F,T}
    @inbounds (@view(integrator.uprev.x[2][2:4]), @view(u.x[2][2:4]))
end

function line_element(u, integrator)
    @inbounds (@view(integrator.uprev[2:4]), @view(u[2:4]))
end
