Base.precompile(
    Tuple{
        typeof(lineprofile),
        KerrMetric{Float64},
        SVector{4,Float64},
        GeometricThinDisc{Float64},
    },
)   # time: 4.8591228
Base.precompile(
    Tuple{
        typeof(tracegeodesics),
        KerrMetric{Float64},
        SVector{4,Float64},
        SVector{4,Float64},
        Vararg{Any},
    },
)   # time: 1.1855857
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
end   # time: 0.9053782
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
        GeometricThinDisc{Float64},
        Vararg{Any},
    },
)   # time: 0.11555993
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
)   # time: 0.04722405
Base.precompile(
    Tuple{
        typeof(tracegeodesics),
        KerrMetric{Float64},
        Vector{SVector{4,Float64}},
        Vector{SVector{4,Float64}},
        Vararg{Any},
    },
)   # time: 0.034065828
Base.precompile(
    Tuple{
        typeof(Core.kwcall),
        NamedTuple{(:n_samples,),Tuple{Int64}},
        typeof(tracegeodesics),
        KerrMetric{Float64},
        LampPostModel{Float64},
        GeometricThinDisc{Float64},
        Vararg{Any},
    },
)   # time: 0.025394669
Base.precompile(
    Tuple{
        typeof(update_integration_parameters!),
        RadiativeTransferIntegrationParameters{Vector{Bool}},
        RadiativeTransferIntegrationParameters{Vector{Bool}},
    },
)   # time: 0.007881916
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
        GeometricThinDisc{Float64},
        Float64,
    },
)   # time: 0.007585501
Base.precompile(
    Tuple{
        var"#200#threadsfor_fun#119"{
            var"#200#threadsfor_fun#118#120"{
                KerrMetric{Float64},
                Matrix{Float64},
                Vector{GeodesicPoint{Float64,SVector{1,Float64}}},
                PointFunction{
                    var"#154#156"{
                        var"#154#155#157"{
                            FilterPointFunction{var"#110#116"{var"#110#111#117"},Float64},
                            var"#109#115",
                            var"#110#116"{var"#110#111#117"},
                        },
                    },
                },
                Float64,
                Base.OneTo{Int64},
            },
        },
        Int64,
    },
)   # time: 0.001631875
