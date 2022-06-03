# Geodesic integration strategies


## Second-Order

The motivation behind the second-order methods is to permit the computation of geodesics in generic spacetimes, via the geodesic equation

```math
\frac{\text{d}^2 x^\mu}{\text{d} \lambda^2} 
    + \Gamma^{\mu}_{\phantom{\mu}\nu\sigma}
    \frac{\text{d}x^\nu}{\text{d} \lambda}
    \frac{\text{d}x^\sigma}{\text{d} \lambda}
= 0,
```

where ``x^\mu`` is a position four-vector, ``\Gamma^{\mu}_{\phantom{\mu}\nu\sigma}`` are the Christoffel symbols of the second kind, and ``\lambda`` the affine parameter describing the curve.

The above can be solved as a second-order ODE, subject to an initial position and initial velocity

```math
u^\mu = \left(t, r, \theta, \phi \right),
\quad
\text{and}
\quad
\dot{u}^\mu  
    = \left( \dot{t}, \dot{r}, \dot{\theta}, \dot{\phi} \right),
```
where the dot refers to the derivative with respect to ``\lambda``. In general, the spatial components of the initial velocity are known _a priori_, and the time-component is determined via the constraint:

```math
g_{\sigma\nu} \dot{u}^\sigma \dot{u}^\nu = -\mu^2,
```

with the metric tensor ``g_{\mu\nu}``, and where ``\mu`` is related to the effective mass associated with the geodesic.