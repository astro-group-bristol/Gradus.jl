var documenterSearchIndex = {"docs":
[{"location":"overview/metrics/#Implemented-Metrics","page":"Available metrics","title":"Implemented Metrics","text":"","category":"section"},{"location":"overview/metrics/","page":"Available metrics","title":"Available metrics","text":"CurrentModule = Gradus","category":"page"},{"location":"overview/metrics/","page":"Available metrics","title":"Available metrics","text":"Gradus.jl implements a library of metrics ready to use for integrations and rendering.","category":"page"},{"location":"overview/metrics/","page":"Available metrics","title":"Available metrics","text":"note: Note\nTo implement your own custom metrics, please see Implementing a new metric.","category":"page"},{"location":"overview/metrics/#First-Order","page":"Available metrics","title":"First-Order","text":"","category":"section"},{"location":"overview/metrics/","page":"Available metrics","title":"Available metrics","text":"BoyerLindquistFO","category":"page"},{"location":"overview/metrics/#Gradus.BoyerLindquistFO","page":"Available metrics","title":"Gradus.BoyerLindquistFO","text":"A first-order implementation of BoyerLindquistAD.\n\nM\nBlack hole mass. Default: 1.0\na\nBlack hole spin. Default: 0.0\nE\nGeodesic energy (a consant of motion). Default: 1.0\n\n\n\n\n\n","category":"type"},{"location":"overview/metrics/#Second-Order","page":"Available metrics","title":"Second-Order","text":"","category":"section"},{"location":"overview/metrics/","page":"Available metrics","title":"Available metrics","text":"BoyerLindquistAD\nJohannsenAD\nKerrRefractiveAD\nDilatonAxionAD\nMorrisThorneAD","category":"page"},{"location":"overview/metrics/#Gradus.BoyerLindquistAD","page":"Available metrics","title":"Gradus.BoyerLindquistAD","text":"struct BoyerLindquistAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n\nThe Kerr metric in Boyer-Lindquist coordinates, describing a black hole with mass M and angular spin a:\n\nbeginalign*\n    textds^2 =\n    - left( 1 - frac2 M rSigma right)textdt^2\n    - frac2M r a sin^2(theta)Sigma textdt textdphi\n    \n    + fracSigmaDelta textdr^2\n    + Sigma textdtheta^2\n    + left(r^2 + a^2 + frac2 M r a^2 sin^2(theta)Sigma right) sin^2(theta) textdphi^2\nendalign*\n\nwhere\n\nSigma = r^2 + a^2 cos^2 (theta)\nquad textand quad\nDelta = r^2 - 2Mr + a^2\n\nParameters\n\nM\nBlack hole mass. Default: 1.0\na\nBlack hole spin. Default: 0.0\n\n\n\n\n\n","category":"type"},{"location":"overview/metrics/#Gradus.JohannsenAD","page":"Available metrics","title":"Gradus.JohannsenAD","text":"struct JohannsenAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n\nThe Johannsen (20xx) metric.\n\nM\nBlack hole mass. Default: 1.0\na\nBlack hole spin. Default: 0.0\nα13\nalpha_13 deviation parameter. Default: 0.0\nα22\nalpha_22 deviation parameter. Default: 0.0\nα52\nalpha_52 deviation parameter. Default: 0.0\nϵ3\nepsilon_3 deviation parameter. Default: 0.0\n\n\n\n\n\n","category":"type"},{"location":"overview/metrics/#Gradus.KerrRefractiveAD","page":"Available metrics","title":"Gradus.KerrRefractiveAD","text":"struct KerrRefractiveAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n\nKerr metric in Boyer-Lindquist coordintes with a path-length ansatz, equivalent to a refractive index n, within the coronal radius corona_radius.\n\nM\nBlack hole mass. Default: 1.0\na\nBlack hole spin. Default: 0.0\nn\nRefractive index within the corona. Default: 1.0\ncorona_radius\nRadius of the corona. Default: 20.0\n\n\n\n\n\n","category":"type"},{"location":"overview/metrics/#Gradus.DilatonAxionAD","page":"Available metrics","title":"Gradus.DilatonAxionAD","text":"struct DilatonAxionAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n\nEinstein-Maxwell-Dilaton-Axion metric.\n\nM\nSingularity mass. Default: 1.0\na\nSingularity spin. Default: 0.0\nβ\nDilaton coupling strength. Default: 0.0\nb\nAxion coupling strength. Default: 1.0\n\n\n\n\n\n","category":"type"},{"location":"overview/metrics/#Gradus.MorrisThorneAD","page":"Available metrics","title":"Gradus.MorrisThorneAD","text":"struct MorrisThorneAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n\nMorris-Thorne wormhole metric.\n\nb\nThroat size. Default: 1.0\n\n\n\n\n\n","category":"type"},{"location":"overview/geodesic-integration/#Geodesic-integration-strategies","page":"Geodesic integration","title":"Geodesic integration strategies","text":"","category":"section"},{"location":"overview/geodesic-integration/#Second-Order","page":"Geodesic integration","title":"Second-Order","text":"","category":"section"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"The motivation behind the second-order methods is to permit the computation of geodesics in generic spacetimes, via the geodesic equation","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"fractextd^2 x^mutextd lambda^2 \n    + Gamma^mu_phantommunusigma\n    fractextdx^nutextd lambda\n    fractextdx^sigmatextd lambda\n= 0","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"where x^mu is a position four-vector, Gamma^mu_phantommunusigma are the Christoffel symbols of the second kind, and lambda the affine parameter describing the curve.","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"The above can be solved as a second-order ODE, subject to an initial position and initial velocity","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"u^mu = left(t r theta phi right)\nquad\ntextand\nquad\ndotu^mu  \n    = left( dott dotr dottheta dotphi right)","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"where the dot refers to the derivative with respect to lambda. In general, the spatial components of the initial velocity are known a priori, and the time-component is determined via the constraint:","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"g_sigmanu dotu^sigma dotu^nu = -mu^2","category":"page"},{"location":"overview/geodesic-integration/","page":"Geodesic integration","title":"Geodesic integration","text":"with the metric tensor g_munu, and where mu is related to the effective mass associated with the geodesic.","category":"page"},{"location":"api-documentation/Gradus/#Gradus-API","page":"Gradus","title":"Gradus API","text":"","category":"section"},{"location":"api-documentation/Gradus/","page":"Gradus","title":"Gradus","text":"CurrentModule = Gradus","category":"page"},{"location":"api-documentation/Gradus/#Special-radii","page":"Gradus","title":"Special radii","text":"","category":"section"},{"location":"api-documentation/Gradus/","page":"Gradus","title":"Gradus","text":"isco\nr_ph\nr_mb\nr_s","category":"page"},{"location":"api-documentation/Gradus/#Gradus.isco","page":"Gradus","title":"Gradus.isco","text":"Innermost stable circular orbit (ISCO), defined by\n\n    fractextdtextdr left( fracEmu right) = 0\n\nUses analytic solutions if known for that metric, else uses a root finder to calculate the radius at which the defining condition is met.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/Gradus/#Gradus.r_ph","page":"Gradus","title":"Gradus.r_ph","text":"Photon orbit radius, defined as the radius for which\n\n    fracEmu rightarrow infty \n\n\n\n\n\n","category":"function"},{"location":"api-documentation/Gradus/#Gradus.r_mb","page":"Gradus","title":"Gradus.r_mb","text":"Marginally bound orbit\n\n    fracEmu = 1\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/Gradus/#Gradus.r_s","page":"Gradus","title":"Gradus.r_s","text":"Event horizon radius, often equivalent to GradusBase.inner_radius, however remains distinct, such that the latter may still be an arbitrary chart cutoff.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#GradusBase-API","page":"GradusBase","title":"GradusBase API","text":"","category":"section"},{"location":"api-documentation/GradusBase/","page":"GradusBase","title":"GradusBase","text":"CurrentModule = Gradus.GradusBase","category":"page"},{"location":"api-documentation/GradusBase/","page":"GradusBase","title":"GradusBase","text":"AbstractMetricParams\ngeodesic_eq\ngeodesic_eq!\nconstrain\ninner_radius\nmetric\nvector_to_local_sky","category":"page"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.AbstractMetricParams","page":"GradusBase","title":"Gradus.GradusBase.AbstractMetricParams","text":"abstract type AbstractMetricParams{T} end\n\nAbstract type used to dispatch different geodesic problems.\n\n\n\n\n\n","category":"type"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.geodesic_eq","page":"GradusBase","title":"Gradus.GradusBase.geodesic_eq","text":"geodesic_eq(m::AbstractMetricParams{T}, u, v)\ngeodesic_eq!(m::AbstractMetricParams{T}, u, v)\n\nCalculate the acceleration components of the geodesic equation given a position u, a velocity v, and a metric m.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.constrain","page":"GradusBase","title":"Gradus.GradusBase.constrain","text":"constrain(m::AbstractMetricParams{T}, u, v; μ::T=0.0)\n\nCalculate time component v^t which would constrain a velocity vector v at position x as a geodesic with mass μ.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.inner_radius","page":"GradusBase","title":"Gradus.GradusBase.inner_radius","text":"inner_radius(m::AbstractMetricParams{T})\n\nReturn the innermost valid coordinate relative to the origin, for use in geodesic tracing.\n\nThis usually represents some property of the metric, e.g. event horizon radius in Kerr/Schwarzschild metrics, or throat diameter in worm hole metrics.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.metric","page":"GradusBase","title":"Gradus.GradusBase.metric","text":"metric(m::AbstractMetricParams{T}, u)\n\nNumerically evaluate the corresponding metric for AbstractMetricParams, given parameter values m and some point u.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Physical-Quantities","page":"GradusBase","title":"Physical Quantities","text":"","category":"section"},{"location":"api-documentation/GradusBase/","page":"GradusBase","title":"GradusBase","text":"E\nLz","category":"page"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.E","page":"GradusBase","title":"Gradus.GradusBase.E","text":"E(m::AbstractMatrix{T}, v)\nE(m::AbstractMetricParams{T}, u, v)\n\nCompute the energy for a numerically evaluated metric, and some velocity four vector v,\n\nE = - p_t = - g_tnu p^nu\n\nFor null geodesics, the velocity is the momentum v^nu = p^nu. For massive geodesics, the mass mu needs to be known to compute mu v^nu = p^nu.\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Gradus.GradusBase.Lz","page":"GradusBase","title":"Gradus.GradusBase.Lz","text":"Lz(m::AbstractMatrix{T}, v)\nLz(m::AbstractMetricParams{T}, u, v)\n\nCompute the angular momentum for a numerically evaluated metric, and some velocity four vector v.\n\nL_z = p_phi = - g_phinu p^nu\n\n\n\n\n\n","category":"function"},{"location":"api-documentation/GradusBase/#Geodesic-Points","page":"GradusBase","title":"Geodesic Points","text":"","category":"section"},{"location":"api-documentation/GradusBase/","page":"GradusBase","title":"GradusBase","text":"GeodesicPoint\ngetgeodesicpoint\nunpack_solution","category":"page"},{"location":"overview/point-functions/#Point-functions","page":"Point Functions","title":"Point functions","text":"","category":"section"},{"location":"internals/custom-metrics/#Implementing-a-new-metric","page":"Implementing new metrics","title":"Implementing a new metric","text":"","category":"section"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"CurrentModule = Gradus","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"Gradus.jl is able to integrate any 3+1 dimensional metric. A new metric may be defined by implementing one of the abstract types with a concrete type, and defining a number of methods. Depending on what you want to be able to do with a metric, different functions need to be implemented. ","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"Gradus also provides a few derivative abstract types to implement to ensure the most efficient code is executed for a given metric (see Metric parameter types below).","category":"page"},{"location":"internals/custom-metrics/#Example:-Schwarzschild","page":"Implementing new metrics","title":"Example: Schwarzschild","text":"","category":"section"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"As a minimal example, here is how the Schwarzschild metric may be implemented. First, we must define what the metric parameters for this metric are. These are effectively constants of the spacetime, representing physical quantities that appear in the metric expression. For the Schwarzschild metric, this is only the black hole mass M, but e.g. the Kerr metric also has the black hole spin a.","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"We can choose the integration strategy by sub-typing an abstract type representing different classes of spacetimes. For the Schwarzschild metric, we will use the static, axis-symmetric class, with the automatic differentiation (AD) backend. With AD, we only need to specify the non-zero components of the metric as Julia functions, and the rest is done for us.","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"For ease, we choose the Eddington-Finkelstein coordinates of the Schwarzschild solution, which may be written","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"textds^2 =\n    - left( 1 - frac2 Mr right) textdt^2\n    + left( 1 - frac2 Mr right)^-1 textdr^2\n    + r^2 textdtheta^2\n    + r^2 sin^2(theta) textdphi^2","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"Here is a possible implementation for Gradus.jl:","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"using Gradus\n\n@with_kw struct EddingtonFinkelsteinAD{T} <: AbstractAutoDiffStaticAxisSymmetricParams{T}\n    M = 1.0\nend\n\nfunction GeodesicTracer.metric_components(m::EddingtonFinkelsteinAD{T}, rθ) where {T}\n    (r, θ) = rθ\n    M = m.M\n\n    tt = 1 - (2M / r)\n    rr = inv(tt)\n    θθ = r^2\n    ϕϕ = r^2 * sin(θ)^2\n\n    (tt, rr, θθ, ϕϕ, T(0.0))\nend\n\nGradusBase.inner_radius(m::BoyerLindquistAD{T}) where {T} = 2 * m.M","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"A few notes:","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"We use @with_kw from Parameters.jl to define various utility constructors for us.\nGeodesicTracer.metric_components must return five elements for AbstractAutoDiffStaticAxisSymmetricParams, where the last element is the off-axis g_t phi matrix element, which in this case is always 0.\nThe GradusBase.inner_radius function defines the inner-radius of integration chart. This defines where the integration should terminate to avoid running indefinitely, and is, in this case, set to the event-horizon of our metric.","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"That's all we need! This metric is now ready to be traced in the usual way.","category":"page"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"note: Note\nFor more examples of how to implement different metrics, click on the \"source\" button of a metric in Implemented Metrics. Alternatively, view the source code directly here.","category":"page"},{"location":"internals/custom-metrics/#Metric-parameter-types","page":"Implementing new metrics","title":"Metric parameter types","text":"","category":"section"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"The following types may be implemented to add new metrics. Each type has different requirements for its interface.","category":"page"},{"location":"internals/custom-metrics/#First-Order","page":"Implementing new metrics","title":"First-Order","text":"","category":"section"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"FirstOrderMethods.AbstractFirstOrderMetricParams\nFirstOrderMethods.four_velocity\nFirstOrderMethods.calc_lq\nFirstOrderMethods.Vr\nFirstOrderMethods.Vθ","category":"page"},{"location":"internals/custom-metrics/#Gradus.FirstOrderMethods.AbstractFirstOrderMetricParams","page":"Implementing new metrics","title":"Gradus.FirstOrderMethods.AbstractFirstOrderMetricParams","text":"AbstractFirstOrderMetricParams{T} <: AbstractMetricParams{T}\n\nAbstract type for metrics using the 1st-order integration method. The 1st-order methods reuse the velocity vector as a parameter vector, where only element vel[2] and vel[3] are used, and are local observer ratios sin Theta and sin Phi respectively.\n\nRequire implementation of\n\nGradusBase.inner_radius\nGradusBase.constrain\nFirstOrderMethods.four_velocity\nFirstOrderMethods.calc_lq\nFirstOrderMethods.Vr\nFirstOrderMethods.Vθ\nGeodesicTracer.alpha_beta_to_vel\n\n\n\n\n\n","category":"type"},{"location":"internals/custom-metrics/#Gradus.FirstOrderMethods.four_velocity","page":"Implementing new metrics","title":"Gradus.FirstOrderMethods.four_velocity","text":"four_velocity(u, m::AbstractFirstOrderMetricParams{T}, p) -> NTuple{4, Any}\n\n\nCalculate the four-velocity at a point u, given a set of metric parameters and the constants of motion in p.\n\n\n\n\n\n","category":"function"},{"location":"internals/custom-metrics/#Gradus.FirstOrderMethods.calc_lq","page":"Implementing new metrics","title":"Gradus.FirstOrderMethods.calc_lq","text":"Calculate constants of motion L and Q, given a set of metric parameters, the geodesic position, and the param vector.\n\n\n\n\n\n","category":"function"},{"location":"internals/custom-metrics/#Gradus.FirstOrderMethods.Vr","page":"Implementing new metrics","title":"Gradus.FirstOrderMethods.Vr","text":"Effective potential in the radial direction. Used only to track sign changes.\n\n\n\n\n\n","category":"function"},{"location":"internals/custom-metrics/#Gradus.FirstOrderMethods.Vθ","page":"Implementing new metrics","title":"Gradus.FirstOrderMethods.Vθ","text":"Vθ(m::AbstractFirstOrderMetricParams{T}, u, p) -> Any\n\n\nEffective potential in the angular direction. Used only to track sign changes.\n\n\n\n\n\n","category":"function"},{"location":"internals/custom-metrics/#Second-Order","page":"Implementing new metrics","title":"Second-Order","text":"","category":"section"},{"location":"internals/custom-metrics/","page":"Implementing new metrics","title":"Implementing new metrics","text":"GeodesicTracer.AbstractAutoDiffMetricParams\nGeodesicTracer.metric_components\nGeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams","category":"page"},{"location":"internals/custom-metrics/#Gradus.GeodesicTracer.AbstractAutoDiffMetricParams","page":"Implementing new metrics","title":"Gradus.GeodesicTracer.AbstractAutoDiffMetricParams","text":"AbstractAutoDiffMetricParams{T} <: AbstractMetricParams{T}\n\nAbstract type for metrics using the 2nd-order integration method, with the automatic differentiation backend.\n\n\n\n\n\n","category":"type"},{"location":"internals/custom-metrics/#Gradus.GeodesicTracer.metric_components","page":"Implementing new metrics","title":"Gradus.GeodesicTracer.metric_components","text":"metric_components(m::AbstractAutoDiffStaticAxisSymmetricParams{T}, rθ) -> Any\n\n\nInterface for GeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams. Should return a vector or tuple with the elements\n\nleft(\n    g_tt g_rr g_theta theta g_phi phi g_tphi\nright)\n\n\n\n\n\n","category":"function"},{"location":"internals/custom-metrics/#Gradus.GeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams","page":"Implementing new metrics","title":"Gradus.GeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams","text":"AbstractAutoDiffStaticAxisSymmetricParams{T} <: AbstractAutoDiffMetricParams{T}\n\nSpecialisation for static, axis-symmetric metrics. Here, the metric is of the form\n\n    g_munu =\n    left( beginmatrix\n        g_tt      0       0                   g_tphi     \n        0           g_rr  0                   0              \n        0           0       g_thetatheta  0              \n        g_tphi  0       0                   g_phiphi\n    endmatrix right)\n\nwhere the only non-zero off axis elements are g_tphi.\n\nRequired implementations:\n\nGradusBase.inner_radius\nGeodesicTracer.metric_components\n\n\n\n\n\n","category":"type"},{"location":"#Gradus.jl-Documentation","page":"Home","title":"Gradus.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A pure Julia geodesic integration system built on DifferentialEquations.jl using automatic differentiation (AD) and computer algebra systems (CAS) to efficiently compute the geodesic equation. This package requires only a specification of the non-zero metric components in order to solve the 2nd order geodesic system. Alternatively, an implementation of the four velocity components may be specified to integrate a regular 1st order system.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The motivation behind this package began with an interest in studying reverberation lags around accreting black holes, however the scope has since expanded to facilitate the exploration of generic metrics through time-like, space-like, and null geodesics. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Our aim is to make testing modified Kerr metrics and alternative gravity theories fast.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<p align=\"center\" pa=\"0\" ma=\"0\">\n<img width=\"30%\" src=\"assets/uob-astro-grey.png\">\n</p>","category":"page"},{"location":"","page":"Home","title":"Home","text":"Gradus.jl allows for drastically different relativistic simulations to be computed with a composable and reusable API, permitting an end user to simply and expressively calculate physical formulae, create observational signatures, and interface with other popular astrophysics tools. Gradus.jl implements a number of high level abstractions, on the path towards a fully parallelized, high performance numerical relativity ecosystem, scalable from personal computers to super computers.","category":"page"},{"location":"#About","page":"Home","title":"About","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Gradus.jl is a research tool for calculating geodesic paths in arbitrary space-times. It is currently work-in-progress, and breaking changes are frequent, as the interface is redesigned to match changing use-cases.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is part of a larger developing eco-system of Strong Gravity Codes, created by members of the University of Bristol Astrophysics Group","category":"page"},{"location":"","page":"Home","title":"Home","text":"Fergus Baker (PhD Student)\nDr. Andrew Young (Associate Professor)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For more University of Bristol Astrophysics Group codes, see our GitHub organisation.","category":"page"}]
}
