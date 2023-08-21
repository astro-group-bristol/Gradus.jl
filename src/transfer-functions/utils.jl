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

function remove_unused_elements!(data::_TransferDataAccumulator)
    data.data = data.data[:, data.mask]
    data
end

function Base.sort!(data::_TransferDataAccumulator)
    I = sortperm(data.θs)
    Base.permutecols!!(data.data, I)
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
    if gmin > minimum(gs)
        @warn (
            "Inferred minima > array minimum rₑ = $rₑ (gmin = $gmin, extrema(gs) = $(extrema(gs)). Using minimum."
        )
        gmin = minimum(gs)
    end
    if gmax < maximum(gs)
        @warn (
            "Inferred maximum < array maximum rₑ = $rₑ (gmax = $gmax, extrema(gs) = $(extrema(gs)). Using maximum."
        )
        gmax = maximum(gs)
    end
    return gmin, gmax
end
