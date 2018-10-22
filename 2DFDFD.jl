# 2D FDFD Solver
using LinearAlgebra
using GeometryTypes
include("Yee.jl")
include("Utilities.jl")

# Physical Constants
c_0 = 299792458
μ_0 = 1.256637061E-6
ϵ_0 = 1/(c_0^2*μ_0)
η_0 = sqrt(μ_0/ϵ_0)

# Solution Setup
Nx = 100
Ny = 100
NGRID = (Nx,Ny) # This is the entire solution space
NPML = (0,0,20,20) # This sets up the PML boundary, xlow, xhigh, ylow, yhigh
kinc = (2.2214,4.4429)
NGRID = (3,5)
RES = (0.5, 0.4)
thisBC = (Dirichlet,Periodic)

calculate_PML_2D(NGRID,NPML)
DEX,DEY = yee_grid_derivative(NGRID,RES,thisBC,kinc)
# Using copy to eagerly evaluate the transpose
DHX = copy(-transpose(DEX))
DHY = copy(-transpose(DEY))
