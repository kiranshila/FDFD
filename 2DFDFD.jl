# 2D FDFD Solver
using LinearAlgebra
using GeometryTypes
using Makie
include("Yee.jl")
include("Utilities.jl")
include("BDG.jl")

println("Setting up physical constants")
# Physical Constants
c_0 = 299792458
μ_0 = 1.256637061E-6
ϵ_0 = 1/(c_0^2*μ_0)
η_0 = sqrt(μ_0/ϵ_0)

# Solution Setup
freq = 24e9
λ_0 = c_0/freq
Resolution = 40
θ = 15 # Angle of incidence in degrees
Polarization = E # E or H - H is TM, E is TE
thisBC = (Dirichlet,Periodic) # Boundary conditions for x,y

println("Building Binary Diffraction Grating")
# Initialize 2X Grid with Binary Diffraction Grating
ϵ_r_2X_Grid,μ_r_2X_Grid,NGRID,NPML,Nx2,Ny2,RES,Q_Limit = BDG(Resolution)

println("Calculating PML Parameters")
# Step 2 - Calculate PML parameters for 2X grid
sx,sy = calculate_PML_2D(2 .* NGRID,2 .* NPML)

println("Incorporating PML into 2X grid")
# Step 3 - Incorporate PML into the 2X grid
μ_r_x = μ_r_2X_Grid ./ sx .* sy
ϵ_r_x = ϵ_r_2X_Grid ./ sx .* sy
μ_r_y = μ_r_2X_Grid .* sx ./ sy
ϵ_r_y = ϵ_r_2X_Grid .* sx ./ sy
μ_r_z = μ_r_2X_Grid .* sx .* sy
ϵ_r_z = ϵ_r_2X_Grid .* sx .* sy

println("Overlaying materials onto 1X grid")
# Step 4 - Overlay materials onto 1X grid
μ_r_x = μ_r_x[1:2:Nx2,2:2:Ny2]
μ_r_y = μ_r_y[2:2:Nx2,1:2:Ny2]
ϵ_r_x = ϵ_r_x[2:2:Nx2,1:2:Ny2]
ϵ_r_y = ϵ_r_y[1:2:Nx2,2:2:Ny2]
ϵ_r_z = ϵ_r_z[1:2:Nx2,1:2:Ny2]
μ_r_z = μ_r_z[2:2:Nx2,2:2:Ny2]

println("Computing wave vector terms")
# Step 5 - Compute wave vector terms
k₀ = (2*pi)/λ_0
k_inc = k₀ .* RES .* [cosd(θ);sind(θ)]
# FIXME

println("Diagonalizing material matricies")
# Step 6 - Diagonalize material matrices
ϵ_r_x = spdiagm(0 => ϵ_r_x[:])
ϵ_r_y = spdiagm(0 => ϵ_r_y[:])
ϵ_r_z = spdiagm(0 => ϵ_r_z[:])
μ_r_z = spdiagm(0 => μ_r_z[:])
μ_r_x = spdiagm(0 => μ_r_x[:])
μ_r_y = spdiagm(0 => μ_r_y[:])

println("Constructing derivative matricies")
# Step 7 - Construct derivative matricies
DEX,DEY,DHX,DHY = yee_grid_derivative(NGRID,k₀*RES, thisBC, k_inc/k₀) # Normalize k terms

println("Computing wave vector matricies - this may take a while")
# Step 8 - Compute wave vector matrix A # Hard
A_E = DHX/μ_r_y*DEX + DHY/μ_r_x*DEY + ϵ_r_z
A_H = DEX/ϵ_r_y*DHX + DEY/ϵ_r_x*DHY + μ_r_z

println("Computing source field")
# Step 9 - Compute source field
F_Src = [exp(1.0im*(k_inc[1]*i + k_inc[2]*j)) for i = 1:NGRID[1], j = 1:NGRID[2]]

println("Computing source field mask")
# Step 10 - Compute scattered-field masking matrix
# Mask is top of solution space, past PML, a few cells into the free space
Q = zeros(NGRID[1],NGRID[2])
for i = 1:Q_Limit
    for j = 1:NGRID[2]
        Q[i,j] = 1
    end
end

Q = spdiagm(0 => Q[:])

println("Computing source vector b")
# Step 11 - Compute source vector b
if Polarization == H
    b = (Q*A_E - A_E*Q)*F_Src[:]
elseif Polarization == E
    b = (Q*A_H - A_H*Q)*F_Src[:]
end

println("Solving FDFD problem - this may take a while")
# Step 12 - Solve
f = Array(A_E)^-1 * b

f = reshape(f,(NGRID[1],NGRID[2]))
# Step 13 - Post Process

println("Visualizing")
# Visualize Data
ϵr_vis = heatmap(ϵ_r_2X_Grid,scale_plot = false)
fields = heatmap(map(x->real(x),f),scale_plot = false, interpolate = true)
scene = AbstractPlotting.vbox(ϵr_vis, fields)
display(scene)

println("Simulation complete")
