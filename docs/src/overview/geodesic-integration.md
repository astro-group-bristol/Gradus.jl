# Geodesic integration strategies

```@meta
CurrentModule = Gradus
```

## Second-Order

The motivation behind the second-order methods is to permit the computation of geodesics in generic spacetimes, via the geodesic equation:

```@docs
Gradus.compute_geodesic_equation
```

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

```@docs
Gradus.constrain_time
```

## Using callbacks