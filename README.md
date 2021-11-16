pde1d
=====

Qt program to compare solution methods for 1-dimensional PDEs.

Has several solvers for linear advection equation and some solvers for linear viscous diffusion and for viscous and inviscid Burger's Equation.

The plots are dynamic so solutions can be visualized in progress.  The instability growth can be seen.

The various solvers can be compared by running at the same time.  Tabs with solver specific values can be moved to independent windows.

An Error table shows various error and comparison terms.

Current dependencies include gsl, gslcblas, Qt4, qtw and blas.  The double part of superLU is included in the latest repository.

