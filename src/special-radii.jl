"""
$(TYPEDSIGNATURES)

Innermost stable circular orbit.

From Bardeen et al. (1972) eq. (2.21):

```math
r_\\text{ms} = M \\left\\{ 3 + Z_2 \\pm \\sqrt{(3 - Z_1)(3 + Z_1 + 2 Z_2)} \\right\\}.
```

The choice of ``\\pm`` is chosen by the sign of ``a``.
"""
isco(M, a, ±) = M * (3 + Z₂(M, a) ± √((3 - Z₁(M, a)) * (3 + Z₁(M, a) + 2 * Z₂(M, a))))
isco(M, a) = a > 0.0 ? isco(M, a, -) : isco(M, a, +)
