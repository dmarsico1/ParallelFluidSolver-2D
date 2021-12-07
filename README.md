# ParallelFluidSolver-2D
Numerically solve the Boussinesq equations in parallel on a two dimensional domain.  The domain is decomposed into rectangular strips that extend the height of the domain, and the elliptic pressure equation is solved with a parallel conjugate gradient method.

To compile: make fluid_par

To run: mpirun -n (# of processors) fluid_par
