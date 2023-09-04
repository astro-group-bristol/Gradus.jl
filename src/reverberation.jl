# first some utility functions

"""
Assumes x is equally spaced. Will not raise an error if not, but results
will be invalid.
"""
function extend_domain_with_zeros(x, y, x_max)
    Δx = x[2] - x[1]

    x̄ = range(minimum(x), x_max, step = Δx) |> collect
    ȳ = zeros(eltype(y), size(x̄))
    ȳ[firstindex(y):lastindex(y)] .= y

    x̄, ȳ
end

sum_impulse_response(f::AbstractMatrix) =
    sum(i -> isnan(i) ? zero(eltype(f)) : i, f, dims = 1) |> vec

phase_difference(𝔉ψ::AbstractArray{<:Complex}) = @. atan(imag(𝔉ψ) / (1 + real(𝔉ψ)))

function time_lag!(out, freq, ϕ)
    for i in eachindex(out)
        out[i] = ϕ[i] / (2π * freq[i])
    end
    out
end

function lag_frequency(t, ψ; R = 1)
    freq = FFTW.fftfreq(length(t), 1 / (t[2] - t[1]))
    𝔉ψ = R .* FFTW.fft(ψ)
    # pick only positive frequencies
    I = 1:(length(freq)÷2)
    ϕ = @views phase_difference(𝔉ψ[I])
    # store the time lag in ϕ directly
    time_lag!(ϕ, freq, ϕ)
    freq[I], -ϕ
end

function lag_frequency(t, f::AbstractMatrix; flo = 5e-5, kwargs...)
    ψ = sum_impulse_response(f)
    t_extended, ψ_extended = extend_domain_with_zeros(t, ψ, 1 / flo)

    lag_frequency(t_extended, ψ_extended; kwargs...)
end

export lag_frequency
