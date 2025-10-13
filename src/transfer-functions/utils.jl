mutable struct _TransferDataAccumulator{T}
    data::Matrix{T}
    cutoff::Int
    mask::BitVector
end

function _TransferDataAccumulator(T::Type, M::Int, cutoff::Int; dims = 4)
    mat = zeros(T, (dims, M))
    mask = BitVector(undef, M)
    _TransferDataAccumulator(mat, cutoff, mask)
end

function Base.getproperty(t::_TransferDataAccumulator, f::Symbol)
    if f === :θs
        @views getfield(t, :data)[1, :]
    elseif f === :gs
        @views getfield(t, :data)[2, :]
    elseif f === :Js
        @views getfield(t, :data)[3, :]
    elseif f === :ts
        @views getfield(t, :data)[4, :]
    else
        getfield(t, f)
    end
end

Base.eachindex(data::_TransferDataAccumulator) = eachindex(data.θs)
Base.lastindex(data::_TransferDataAccumulator) = size(data.data, 2)

function Base.sort!(data::_TransferDataAccumulator)
    I = sortperm(data.θs)
    Base.permutecols!!(data.data, I)
end

function _search_visibility_range(vec)
    i1::Int = 0
    i2::Int = 0
    if all(==(true), vec)
        return 0, 0
    end
    s = first(vec)
    for (i, v) in enumerate(vec)
        if v != s
            if i1 == 0
                i1 = i
            elseif i2 == 0
                i2 = i - 1
                break
            end
            s = !s
        end
    end
    if first(vec) == 0
        if i2 == 0
            i2 = lastindex(vec)
        else
            i1 -= 1
        end
    else
        i1 -= 1
        i2 += 1
    end
    i1, i2
end

function insert_data!(data::_TransferDataAccumulator, i, θ, vals::NTuple{3})
    data.mask[i] = true
    data.data[:, i] .= (θ, vals...)
    data
end
function insert_data!(data::_TransferDataAccumulator, i, θ, vals::NTuple{4})
    data.mask[i] = true
    data.data[:, i] .= (θ, vals[1:3]...)
    if vals[4] == 0
        data.Js[i] = NaN
    end
    data
end

function _check_gmin_gmax(_gmin, _gmax, rₑ, gs)
    gmin = _gmin
    gmax = _gmax
    # sometimes the interpolation seems to fail if there are duplicate knots
    if isnan(gmin)
        @warn "gmin is NaN for rₑ = $rₑ (gmin = $gmin, extrema(gs) = $(extrema(gs)). Using extrema."
        gmin = minimum(gs)
    end
    if isnan(gmax)
        @warn "gmax is NaN for rₑ = $rₑ (gmax = $gmax, extrema(gs) = $(extrema(gs)). Using extrema."
        gmax = maximum(gs)
    end
    if gmin == gmax
        @warn (
            "gmin == gmax for rₑ = $rₑ (gmin = gmax = $gmin, extrema(gs) = $(extrema(gs)). Using extrema."
        )
        gmin, gmax = extrema(gs)
        if gmin == gmax
            error("Cannot use extrema")
        end
    end

    # check if they are better than the estimates from the array
    if gmin > minimum(gs)
        gmin = minimum(gs)
    end

    if gmax < maximum(gs)
        gmax = maximum(gs)
    end

    return gmin, gmax
end

function _rθ_to_αβ(r, θ; α₀ = 0, β₀ = 0)
    sinθ, cosθ = sincos(θ)
    α = (r * cosθ) + α₀
    β = (r * sinθ) + β₀
    (α, β)
end

function _normalize!(flux::AbstractVector{T}, grid) where {T}
    Σflux = zero(T)
    @inbounds for i = 1:(length(grid)-1)
        ḡ = (grid[i+1] + grid[i])
        flux[i] = flux[i] / ḡ
        Σflux += flux[i]
    end
    if Σflux > 0
        @. flux = flux / Σflux
    end
    flux
end

function _normalize!(flux::AbstractMatrix{T}, grid) where {T}
    Σflux = zero(T)
    @views @inbounds for i = 1:(length(grid)-1)
        ḡ = (grid[i+1] + grid[i])
        @. flux[i, :] = flux[i, :] / ḡ
        Σflux += sum(flux[i, :])
    end
    if Σflux > 0
        @. flux = flux / Σflux
    end
    # normalize so that the maximal row sums to 1
    flux = flux ./ maximum(sum(flux, dims = 2))
    flux
end
