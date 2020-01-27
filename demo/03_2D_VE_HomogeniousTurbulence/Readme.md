## Two dimensional Vorticity Euation

* Governing equation

  Two-dimensional vorticity equation

* Boundary condition

  Doubly-periodic

* Spacial discretization 

  Fourier-Galerkin

* Time discretization

  nonlinear term: Classical 4 stage 4th order Runge-Kutta method

  linear term: Integrating factor technique

## Problem setting

Two-dimensional decaying turbulence[^san]

# License

This software is released under the MIT License, see [LICENSE](https://github.com/toya42/psm4cfd/blob/master/LICENSE).

## References

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Fundamentals in Single Domain*, Springer-Verlag (2006).

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer-Verlag (2007).

[^san]: San,O.,Staples,A.E.:High-order methods for decaying two-dimensional homogeneous isotropic turbulence. Comput. Fluids 63, 105–127 (2012)