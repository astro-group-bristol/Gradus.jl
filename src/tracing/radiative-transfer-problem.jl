
covariant_absorption_emission_velocity(args...) = (0, 0, SVector(0.0, 0.0, 0.0, 0.0))

function radiative_transfer(m::AbstractMetricParameters, x, k, geometry, I, ν, invν3)
    a_ν, j_ν, u = covariant_absorption_emission_velocity(m, x, ν, geometry)
    g = metric(m, x)
    # cache inv(ν)^3 to avoid costly division
    -(g * k) ⋅ u * (-a_ν * I + j_ν * invν3)
end

function radiative_transfer_ode_problem(
    m::AbstractMetricParameters,
    pos::StaticVector,
    vel::StaticVector,
    geometry,
    time_domain;
    ν = 1.0,
    kwargs...,
)
    invν3 = inv(ν)^3
    println("YOU MADE IT HERE!")
    function f(u::SVector{9}, p, λ)
        @inbounds let x = SVector{4,T}(@view(u[1:4])),
            k = SVector{4,T}(@view(u[5:8])),
            I = u[9]

            dv = SVector{4,T}(geodesic_eq(m, x, k))
            dI = radiative_transfer(m, x, k, geometry, I, ν, invν3)
            vcat(v, dv, SVector(dI))
        end
    end

    u_init = vcat(pos, vel, SVector(0.0))
    ODEProblem{false}(
        f,
        u_init,
        time_domain,
        IntegrationParameters(StatusCodes.NoStatus);
        kwargs...,
    )
end
