struct _TransferDataAccumulator{T}
    θs::Vector{T}
    gs::Vector{T}
    Js::Vector{T}
    ts::Vector{T}
    cutoff::Int
end

function _TransferDataAccumulator(T::Type, M::Int, cutoff::Int)
    _TransferDataAccumulator(zeros(T, M), zeros(T, M), zeros(T, M), zeros(T, M), cutoff)
end

Base.eachindex(data) = eachindex(data.θs)
Base.lastindex(data) = lastindex(data.θs)
function remove_unused_elements!(data)
    I = findall(==(0), data.θs)
    deleteat!(data.θs, I)
    deleteat!(data.gs, I)
    deleteat!(data.Js, I)
    deleteat!(data.ts, I)
end
function Base.sort!(data)
    I = sortperm(data.θs)
    @. data.θs = data.θs[I]
    @. data.gs = data.gs[I]
    @. data.Js = data.Js[I]
    @. data.ts = data.ts[I]
end

function insert_data!(data::_TransferDataAccumulator, i, θ, g, J, t)
    @assert i <= lastindex(data.θs) "$i > $(lastindex(data.θs))"
    data.θs[i] = θ
    data.gs[i] = g
    data.Js[i] = J
    data.ts[i] = t
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
