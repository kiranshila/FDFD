# 2D FDFD Solver
using LinearAlgebra

# Physical Constants
c_0 = 299792458
μ_0 = 1.256637061E-6
ϵ_0 = 1/(c_0^2*μ_0)
η_0 = sqrt(μ_0/ϵ_0)

# Solution Setup
Nx = 7
Ny = 4
NGRID = (Nx,Ny) # This is the entire solution space

NPML = (2,3,1,2) # This sets up the PML boundary, xlow, xhigh, ylow, yhigh

function calculate_PML(grid_size, PML_size)
    a_max = 3
    p = 3
    σ_prime_max = 1

    # There are two cases, 2D and 3D - assume 2D #FIXME
    sx = fill(1,grid_size) # Fill soulution space with 1s
    sy = fill(1,grid_size)
    
    # Add xlow PML
    for nx = 1:PML_size[1]
        ax = 1 + a_max * (nx/PML_size[1])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[1])))^2
        sx[PML_size[1]-nx+1,:] = fill(ax * (1 + 1.0im * η_0 * σ_prime_x),(grid_size[1],1))
    end
    # Add xhigh PML
    for nx = 1:PML_size[2]
        ax = 1 + a_max * (nx/PML_size[2])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[2])))^2
        sx[grid_size[1]-PML_size[2]+nx,:] = fill(ax * (1 + 1.0im * η_0 * σ_prime_x),(grid_size[1],1))
    end
    # Add ylow PML
    #=
    for ny = 1:PML_size[3]
        ay = 1 + a_max * (ny/PML_size[3])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[3])))^2
        sy[:,PML_size[3]-ny+1] = fill(ay * (1 + 1.0im * η_0 * σ_prime_y),(1,grid_size[2]))
    end
    # Add yhigh PML
    for ny = 1:PML_size[4]
        ay = 1 + a_max * (ny/PML_size[4])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[4])))^2
        sy[:,grid_size[2]-PML_size[4]+ny] = fill(ay * (1 + 1.0im * η_0 * σ_prime_y),(1,grid_size[2]))
    end
    =#
    return [sx sy]
end
