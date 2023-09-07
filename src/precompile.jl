Base.precompile(
    Tuple{
        typeof(_second_order_ode_f),
        SVector{8,Float64},
        IntegrationParameters{KerrMetric{Float64}},
        Float64,
    },
)   # time: 5.720331
Base.precompile(
    Tuple{typeof(lineprofile),KerrMetric{Float64},SVector{4,Float64},ThinDisc{Float64}},
)   # time: 2.5695212
Base.precompile(
    Tuple{
        typeof(tracegeodesics),
        KerrMetric{Float64},
        SVector{4,Float64},
        SVector{4,Float64},
        Vararg{Any},
    },
)   # time: 1.4936496
let fbody = try
        __lookup_kwbody__(
            which(tracegeodesics, (KerrMetric{Float64}, SVector{4,Float64}, Vararg{Any})),
        )
    catch missing
    end
    if !ismissing(fbody)
        precompile(
            fbody,
            (
                Float64,
                Float64,
                TraceGeodesic{Float64},
                Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}},
                typeof(tracegeodesics),
                KerrMetric{Float64},
                SVector{4,Float64},
                Vararg{Any},
            ),
        )
    end
end   # time: 1.0984381
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:n_samples,),Tuple{Int64}},
        typeof(emissivity_profile),
        KerrMetric{Float64},
        ThinDisc{Float64},
        LampPostModel{Float64},
    },
)   # time: 0.9623503
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:n_samples,),Tuple{Int64}},
        typeof(tracecorona),
        KerrMetric{Float64},
        ThinDisc{Float64},
        LampPostModel{Float64},
    },
)   # time: 0.4676826
let fbody = try
        __lookup_kwbody__(
            which(tracegeodesics, (KerrMetric{Float64}, SVector{4,Float64}, Vararg{Any})),
        )
    catch missing
    end
    if !ismissing(fbody)
        precompile(
            fbody,
            (
                Float64,
                Float64,
                TraceRadiativeTransfer{Float64},
                Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}},
                typeof(tracegeodesics),
                KerrMetric{Float64},
                SVector{4,Float64},
                Vararg{Any},
            ),
        )
    end
end   # time: 0.3491159
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:n_samples,),Tuple{Int64}},
        typeof(tracegeodesics),
        KerrMetric{Float64},
        LampPostModel{Float64},
        ThinDisc{Float64},
        Vararg{Any},
    },
)   # time: 0.3443357
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:image_width, :image_height),Tuple{Int64,Int64}},
        typeof(rendergeodesics),
        KerrMetric{Float64},
        SVector{4,Float64},
        ThinDisc{Float64},
        Vararg{Any},
    },
)   # time: 0.19957814
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:trace,),Tuple{TraceRadiativeTransfer{Float64}}},
        typeof(tracegeodesics),
        KerrMetric{Float64},
        SVector{4,Float64},
        SVector{4,Float64},
        Vararg{Any},
    },
)   # time: 0.13094269
Base.precompile(Tuple{Type{RadialDiscProfile},Vector{Float64},Vector{Float64},Vector{Any}})   # time: 0.08506817
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{
            (:trace, :image_width, :image_height),
            Tuple{TraceRadiativeTransfer{Float64},Int64,Int64},
        },
        typeof(rendergeodesics),
        KerrMetric{Float64},
        SVector{4,Float64},
        ThinDisc{Float64},
        Vararg{Any},
    },
)   # time: 0.0841068
let fbody = try
        __lookup_kwbody__(
            which(
                tracegeodesics,
                (
                    TraceGeodesic{Float64},
                    KerrMetric{Float64},
                    SVector{4,Float64},
                    Vararg{Any},
                ),
            ),
        )
    catch missing
    end
    if !ismissing(fbody)
        precompile(
            fbody,
            (
                Base.Pairs{
                    Symbol,
                    Any,
                    Tuple{Symbol,Symbol},
                    NamedTuple{(:save_on, :chart),Tuple{Bool,PolarChart{Float64}}},
                },
                typeof(tracegeodesics),
                TraceGeodesic{Float64},
                KerrMetric{Float64},
                SVector{4,Float64},
                Vararg{Any},
            ),
        )
    end
end   # time: 0.036170956
Base.precompile(
    Tuple{
        typeof(tracegeodesics),
        KerrMetric{Float64},
        Vector{SVector{4,Float64}},
        Vector{SVector{4,Float64}},
        Vararg{Any},
    },
)   # time: 0.033662505
Base.precompile(
    Tuple{
        typeof(tracing_configuration),
        TraceGeodesic{Float64},
        KerrMetric{Float64},
        SVector{4,Float64},
        SVector{4,Float64},
        ThinDisc{Float64},
        Float64,
    },
)   # time: 0.030030003
Base.precompile(
    Tuple{
        typeof(tracing_configuration),
        TraceRadiativeTransfer{Float64},
        KerrMetric{Float64},
        SVector{4,Float64},
        SVector{4,Float64},
        ThinDisc{Float64},
        Float64,
    },
)   # time: 0.01967225
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{
            (:trajectories, :save_on, :ensemble),
            Tuple{Int64,Bool,EnsembleEndpointThreads},
        },
        typeof(tracing_configuration),
        TraceGeodesic{Float64},
        KerrMetric{Float64},
        SVector{4,Float64},
        Function,
        ThinDisc{Float64},
        Float64,
    },
)   # time: 0.010621044
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{
            (:trajectories, :save_on, :ensemble),
            Tuple{Int64,Bool,EnsembleEndpointThreads},
        },
        typeof(tracing_configuration),
        TraceRadiativeTransfer{Float64},
        KerrMetric{Float64},
        SVector{4,Float64},
        Function,
        ThinDisc{Float64},
        Float64,
    },
)   # time: 0.008629832
Base.precompile(
    Tuple{
        typeof(tracing_configuration),
        TraceGeodesic{Float64},
        KerrMetric{Float64},
        Vector{SVector{4,Float64}},
        Vector{SVector{4,Float64}},
        ThinDisc{Float64},
        Float64,
    },
)   # time: 0.007926913
Base.precompile(
    Tuple{
        typeof(update_integration_parameters!),
        RadiativeTransferIntegrationParameters{KerrMetric{Float64},Vector{Bool}},
        RadiativeTransferIntegrationParameters{KerrMetric{Float64},Vector{Bool}},
    },
)   # time: 0.007825662
Base.precompile(
    Tuple{typeof(point_source_equatorial_disc_emissivity),Float64,Any,Float64,Float64},
)   # time: 0.006162999
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:gtol,),Tuple{Float64}},
        typeof(distance_to_disc),
        ThinDisc{Float64},
        SVector{9,Float64},
    },
)   # time: 0.003949623
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:gtol,),Tuple{Float64}},
        typeof(distance_to_disc),
        ThinDisc{Float64},
        SVector{8,Float64},
    },
)   # time: 0.001979128
Base.precompile(
    Tuple{
        typeof(set_status_code!),
        IntegrationParameters{KerrMetric{Float64}},
        Gradus.StatusCodes.T,
    },
)   # time: 0.001029539
