# 2D FDFD Solver
using LinearAlgebra
using GeometryTypes
using Makie
include("Yee.jl")
include("Utilities.jl")

scene = Scene(resolution = (2048, 2048)); # Setup scene for visualization

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
NGRID = (Nx,Ny) # This is the entire solution space in number of boxes in the "1x" grid, Nx is rows, Ny is columns
RES = [(5*λ_0)/Nx,(5*λ_0)/Ny] # This is the grid resolution, test setup here is 5λ
NPML = (10,10,10,10) # This sets up the PML boundary, xlow, xhigh, ylow, yhigh - size in "1x" grid
θ = 20 # Angle of incidence in degrees
Polarization = E # E or H - H is TM, E is TE
n_reflected = 100 # FIXME what are these?
n_transmitted = 100
thisBC = (Dirichlet,Dirichlet) # Boundary conditions for x,y

# Setup calculations
Nx2 = 2*Nx
Ny2 = 2*Ny
dx2 = RES[1]/2
dy2 = RES[2]/2

### Test 1 - Air ###

# Initialize 2X Grid
ϵ_r_2X_Grid = ones(ComplexF64,Nx2,Ny2)
μ_r_2X_Grid = ones(ComplexF64,Nx2,Ny2)

# Step 2 - Calculate PML parameters for 2X grid
sx,sy = calculate_PML_2D(2 .* NGRID,2 .* NPML)

# Step 3 - Incorporate PML into the 2X grid
μ_r_x = μ_r_2X_Grid ./ sx .* sy
ϵ_r_x = ϵ_r_2X_Grid ./ sx .* sy
μ_r_y = μ_r_2X_Grid .* sx ./ sy
ϵ_r_y = ϵ_r_2X_Grid .* sx ./ sy
μ_r_z = μ_r_2X_Grid .* sx .* sy
ϵ_r_z = ϵ_r_2X_Grid .* sx .* sy

# Step 4 - Overlay materials onto 1X grid
μ_r_x = μ_r_x[1:2:Nx2,2:2:Ny2]
μ_r_y = μ_r_y[2:2:Nx2,1:2:Ny2]
ϵ_r_x = ϵ_r_x[2:2:Nx2,1:2:Ny2]
ϵ_r_y = ϵ_r_y[1:2:Nx2,2:2:Ny2]
ϵ_r_z = ϵ_r_z[1:2:Nx2,1:2:Ny2]
μ_r_z = μ_r_z[2:2:Nx2,2:2:Ny2]

# Step 5 - Compute wave vector terms
k₀ = (2*pi)/λ_0
k_inc = k₀ .* RES .* [cosd(θ);sind(θ)]
m = [-floor(Nx/2):floor(Nx/2)]
# FIXME

# Step 6 - Diagonalize material matrices
ϵ_r_x = spdiagm(0 => ϵ_r_x[:])
ϵ_r_y = spdiagm(0 => ϵ_r_y[:])
ϵ_r_z = spdiagm(0 => ϵ_r_z[:])
μ_r_z = spdiagm(0 => μ_r_z[:])
μ_r_x = spdiagm(0 => μ_r_x[:])
μ_r_y = spdiagm(0 => μ_r_y[:])

# Step 7 - Construct derivative matricies
DEX,DEY,DHX,DHY = yee_grid_derivative(NGRID,k₀*RES, thisBC, k_inc/k₀) # Normalize k terms

# Step 8 - Compute wave vector matrix A # Hard
A_E = DHX/μ_r_y*DEX + DHY/μ_r_x*DEY + ϵ_r_z
A_H = DEX/ϵ_r_y*DHX + DEY/ϵ_r_x*DHY + μ_r_z

# Step 9 - Compute source field
F_Src = [exp(1.0im*(k_inc[1]*i + k_inc[2]*j)) for i = 1:Nx, j = 1:Ny]

# Step 10 - Compute scattered-field masking matrix
 # FIXME, simple mask
Q = zeros(Ny,Nx)
for i = 1:div(Ny,2)
    for j = 1:Nx
        Q[i,j] = 1
    end
end

Q = spdiagm(0 => Q[:])

# Step 11 - Compute source vector b
if Polarization == H
    b = (Q*A_E - A_E*Q)*F_Src[:]
elseif Polarization == E
    b = (Q*A_H - A_H*Q)*F_Src[:]
end

# Step 12 - Solve
f = Array(A_E)^-1 * b

f = reshape(f,(Nx,Ny))
# Step 13 - Post Process

# Visualize Data
heatmap(map(x->real(x),f))
