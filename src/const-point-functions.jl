module ConstPointFunctions
import ..Rendering: PointFunction, FilterPointFunction
import ..AccretionFormulae: _redshift_guard

const filter_early_term =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t < max_time, NaN)

const filter_intersected =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.retcode == :Intersected, NaN)

const affine_time = PointFunction((m, gp, max_time; kwargs...) -> gp.t)

const shadow = affine_time ∘ filter_early_term

const redshift = PointFunction(_redshift_guard)

end # module
