abstract type AbstractSkyDomain end

struct LowerHemisphere <: AbstractSkyDomain end
struct BothHemispheres <: AbstractSkyDomain end

abstract type AbstractDirectionSampler{AbstractSkyDomain} end

struct RandomSampler{D} <: AbstractDirectionSampler{D}
    RandomSampler(domain::AbstractSkyDomain) = new{domain}()
    RandomSampler(; domain::AbstractSkyDomain = LowerHemisphere()) = new{typeof(domain)}()
end
struct EvenSampler{D} <: AbstractDirectionSampler{D}
    EvenSampler(domain::AbstractSkyDomain) = new{domain}()
    EvenSampler(; domain::AbstractSkyDomain = LowerHemisphere()) = new{typeof(domain)}()
end
struct WeierstrassSampler{D} <: AbstractDirectionSampler{D}
    resolution::Float64
    WeierstrassSampler(res, domain::AbstractSkyDomain) = new{typeof(domain)}(res)
    WeierstrassSampler(; res = 100.0, domain::AbstractSkyDomain = LowerHemisphere()) =
        WeierstrassSampler(res, domain)
end

sample_azimuth(sm::AbstractDirectionSampler{D}, i) where {D} =
    error("Not implemented for $(typeof(sm)).")
sample_elevation(sm::AbstractDirectionSampler{D}, i) where {D} =
    error("Not implemented for $(typeof(sm)).")
sample_angles(sm::AbstractDirectionSampler{D}, i, N) where {D} =
    (sample_elevation(sm, i / N), sample_azimuth(sm, i))

# even
sample_azimuth(::RandomSampler{D}, i) where {D} = rand(Float64) * 2π
sample_elevation(::RandomSampler{LowerHemisphere}, i) = acos(1 - rand(Float64))
sample_elevation(::RandomSampler{BothHemispheres}, i) = acos(1 - 2 * rand(Float64))

# even
sample_azimuth(::EvenSampler{D}, i) where {D} = π * (1 + √5) * i
sample_elevation(::EvenSampler{LowerHemisphere}, i) = acos(1 - i)
sample_elevation(::EvenSampler{BothHemispheres}, i) = acos(1 - 2i)

# uniform
sample_azimuth(::WeierstrassSampler{D}, i) where {D} = π * (1 + √5) * i
sample_elevation(sm::WeierstrassSampler{LowerHemisphere}, i) =
    π - 2atan(√(sm.resolution / i))
function sample_elevation(sm::WeierstrassSampler{BothHemispheres}, i)
    ϕ = 2atan(√(sm.resolution / i))
    if iseven(i)
        ϕ
    else
        π - ϕ
    end
end
sample_angles(sm::WeierstrassSampler{D}, i, N) where {D} =
    (sample_elevation(sm, i), sample_azimuth(sm, i))

export LowerHemisphere, BothHemispheres, RandomSampler, EvenSampler, WeierstrassSampler
