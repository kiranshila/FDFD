# A Finite Difference Frequency Domain Solver in Julia
This code follows the procedure of the FDFD method outlined in Dr. Raymond Rumpf's lecture notes from the EE 5337 course at UTEP.
This code was written as a final project for EEL 6935 Linear Algebra Course at USF.

## Curent Status
As of the latest commit, the solver is working for a 24 GHz binary diffraction grating using a slight modification to how the incident wave vector is being calculated. It is currently off by a scale of 0.0012.

## Prerequisites
* Julia 1.0
* IterativeSolvers
* SparseArrays
* LinearAlgebra
* Makie (for plotting)

## Running
Simply run the script
```sh
$ julia 2DFDFD.jl
```
## Future Work
I intend to discover why the k_inc section of the code is incorrect. Once I do that, I will add support for 3D solving and importing step file geometry. Included is an STL rasterizer in Utilities.jl. Additionally, I will add iterative convergance, a better iterative solver using Incomplete LU preconditioning, and calculations for transmitted and reflected power.
