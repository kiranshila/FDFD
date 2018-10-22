# 2D FDFD Solver
using LinearAlgebra
using GeometryTypes
using Makie
include("Yee.jl")
include("Utilities.jl")

# Physical Constants
c_0 = 299792458
μ_0 = 1.256637061E-6
ϵ_0 = 1/(c_0^2*μ_0)
η_0 = sqrt(μ_0/ϵ_0)

# Solution Setup
freq = 38e9 # Source frequency
λ_0 = c_0/freq
Nx = 100
Ny = 100
NGRID = (Nx,Ny) # This is the entire solution space in number of boxes in the "2x" grid
RES = ((5*λ_0)/Nx,(5*λ_0)/Ny) # This is the grid resolution, test setup here is 5λ
NPML = (10,10,10,10) # This sets up the PML boundary, xlow, xhigh, ylow, yhigh - size in "1x" grid
θ = 20 # Angle of incidence in degrees
Polarization = E # E or H - Mode

# Setup calculations
kinc = (2.2214,4.4429)
Nx2 = 2*Nx
Ny2 = 2*Ny
dx2 = RES[1]/2
dy2 = RES[2]/2

# Test 1 - Air
ϵ_r_2X_Grid = ones(ComplexF64,Nx2,Ny2)
μ_r_2X_Grid = ones(ComplexF64,Nx2,Ny2)


# Step 2 - Calculate PML parameters for 2X grid
sx,sy = calculate_PML_2D(NGRID,NPML)
# Step 3 - Incorporate PML into the 2X grid
μ_r_x = μ_r_2X_Grid ./ sx .* sy
μ_r_y = μ_r_2X_Grid .* sx ./ sy
ϵ_r_x = ϵ_r_2X_Grid ./ sx .* sy
ϵ_r_y = ϵ_r_2X_Grid .* sx ./ sy
# Step 4 - Overlay materials onto 1X grid
μ_r_x = μ_r_x[1:2:Nx2,2:2:Ny2]


thisBC = (Dirichlet,Dirichlet)


DEX,DEY,DHX,DHY = yee_grid_derivative(NGRID,RES,thisBC,kinc)

# Visualize Data
heatmap(map(x->abs(x),ϵ_r_2X_Grid))
