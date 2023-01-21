
struct LowerHemisphere <: AbstractSkyDomain end
struct BothHemispheres <: AbstractSkyDomain end

struct RandomGenerator <: AbstractGenerator end
struct GoldenSpiralGenerator <: AbstractGenerator end
struct EvenGenerator <: AbstractGenerator end

struct EvenSampler{D,G} <: AbstractDirectionSampler{D,G}
    EvenSampler(domain::AbstractSkyDomain, generator::AbstractGenerator) =
        new{typeof(domain),typeof(generator)}()
    EvenSampler(;
        domain::AbstractSkyDomain = LowerHemisphere(),
        generator::AbstractGenerator = GoldenSpiralGenerator(),
    ) = EvenSampler(domain, generator)
end
struct WeierstrassSampler{D,G} <: AbstractDirectionSampler{D,G}
    resolution::Float64
    WeierstrassSampler(res, domain::AbstractSkyDomain, generator::AbstractGenerator) =
        new{typeof(domain),typeof(generator)}(res)
    WeierstrassSampler(;
        res = 100.0,
        domain::AbstractSkyDomain = LowerHemisphere(),
        generator::AbstractGenerator = GoldenSpiralGenerator(),
    ) = WeierstrassSampler(res, domain, generator)
end

@inline geti(::AbstractDirectionSampler{D,EvenGenerator}, i, N) where {D} = i
@inline geti(::AbstractDirectionSampler{D,GoldenSpiralGenerator}, i, N) where {D} = i
@inline geti(::AbstractDirectionSampler{D,RandomGenerator}, i, N) where {D} =
    rand(Float64) * N

@inline getangle(::AbstractDirectionSampler{D,EvenGenerator}, i) where {D} = 2π * i
@inline getangle(::AbstractDirectionSampler{D,GoldenSpiralGenerator}, i) where {D} =
    π * (1 + √5) * i
@inline getangle(::AbstractDirectionSampler{D,RandomGenerator}, i) where {D} = 2π * i

sample_angles(sm::AbstractDirectionSampler{D,G}, i, N) where {D,G} =
    (sample_elevation(sm, i / N), mod2pi(getangle(sm, i)))
@inline sample_angles(sm::WeierstrassSampler{D}, i, N) where {D} =
    (sample_elevation(sm, i), mod2pi(getangle(sm, i)))


sample_elevation(sm::AbstractDirectionSampler, i) =
    error("Not implemented for $(typeof(sm)).")
@inline sample_elevation(::EvenSampler{LowerHemisphere}, i) = acos(1 - i)
@inline sample_elevation(::EvenSampler{BothHemispheres}, i) = acos(1 - 2i)
@inline sample_elevation(sm::WeierstrassSampler{LowerHemisphere}, i) =
    2atan(√(sm.resolution / i))
@inline function sample_elevation(sm::WeierstrassSampler{BothHemispheres}, i)
    ϕ = 2atan(√(sm.resolution / i))
    if iseven(i)
        ϕ
    else
        π - ϕ
    end
end

function _cart_to_spher_jacobian(θ, ϕ)
    @SMatrix [
        sin(θ)*cos(ϕ) sin(θ)*sin(ϕ) cos(θ)
        cos(θ)*cos(ϕ) cos(θ)*sin(ϕ) -sin(θ)
        -sin(ϕ) cos(ϕ) 0
    ]
end

function _cart_local_direction(θ, ϕ)
    @SVector [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]
end

export LowerHemisphere,
    BothHemispheres, RandomGenerator, GoldenSpiralGenerator, EvenSampler, WeierstrassSampler
