@enum BC Periodic=1 Dirichlet=2

function yee_grid_derivative(grid_size,grid_resolution,boundary_condition::Tuple{BC,BC},k_inc = [0 0])
    # Generates the Yee Grid Derivative on the 2D Grid #FIXME for 3D


    # To create a sparse array, we have to make I a vector of row idicies,
    # J a vector of column indicies, and V a vector of values
    # We do this to create our large derivative array for the two boundary condition cases

    # For the case when we supply an incident wave vector but the boundary is Dirichlet,
    # we need to set the k_inc back to zero
    if boundary_condition[1] == Dirichlet
        k_inc[1] = 0
    end
    if boundary_condition[2] == Dirichlet
        k_inc[2] = 0
    end

    # Storing results in a sparse array of row index I, column index J, and value V
    _I = Int64[] # I is reserved for the UniformScaling constant for Identity
    J = Int64[]
    V = Complex{Float64}[]
    DEX = nothing
    DEY = nothing

    # Check if Nx = 1
    if grid_size[1] == 1
        DEX = 1.0im*k_inc[1]*sparse(I,grid_size[1]*grid_size[2],grid_size[1]*grid_size[2])
    else
        # Every element is multiplied by 1/Δx
        entry_x = 1/grid_resolution[1]
        for i = 1:grid_size[1]*grid_size[2]
            # Main Diagonal, all -1*1/ΔX
            push!(_I,i)
            push!(J,i)
            push!(V,-1*entry_x)
        end
        for i = 1:grid_size[1]*grid_size[2]-1
            # Secondary Diagonal, above first diagonal, all 1*1/ΔX
            push!(_I,i)
            push!(J,i+1)
            push!(V,entry_x)
        end
        # Create sparse matrix for DEX
        DEX = sparse(_I,J,V)
        if boundary_condition[1] == Dirichlet
            # Every Nx element on secondary diagonal is 0
            # This makes sense as the value outside the solution space is forced
            # to zero
            for i = grid_size[1]:grid_size[1]:grid_size[1]*grid_size[2]-1
                DEX[i,i+1] = 0
            end
        end
        if boundary_condition[1] == Periodic
            # Calculate entry from Bloch's theorum
            Λ_x = grid_size[1] * grid_resolution[1]
            entry_periodic_x = exp(1.0im*k_inc[1]*Λ_x)

            # Place entry at boundary
            for i = grid_size[1]:grid_size[1]:grid_size[1]*grid_size[2]
                DEX[i,i-grid_size[1]+1] = entry_periodic_x*entry_x
                if i+1<= grid_size[1]*grid_size[2]
                    DEX[i,i+1] = 0
                end
            end
        end
    end

    # Clear Variables
    _I = Int64[]
    J = Int64[]
    V = Complex{Float64}[]

    # Check if Ny = 1
    if grid_size[2] == 1
        DEY = 1.0im*k_inc[2]*sparse(I,grid_size[1]*grid_size[2],grid_size[1]*grid_size[2])
    else
        # Every element is multiplied by 1/Δy
        entry_y = 1/grid_resolution[2]
        for i = 1:grid_size[1]*grid_size[2]
            # Main Diagonal
            push!(_I,i)
            push!(J,i)
            push!(V,-entry_y)
        end
        for i = 1:grid_size[1]*grid_size[2]-grid_size[1]
            # Secondary Diagonal
            push!(_I,i)
            push!(J,i+grid_size[1]) # The diagonal is located Nx past the main diagonal
            push!(V,entry_y)
        end
        DEY = sparse(_I,J,V)
        # There is no manifestation of the Dirichlet bondary condition in the Dy matrix
        if boundary_condition[2] == Periodic
            # For the periodic boundary condition in y, we only place elements when the boundary exists
            # This forms another diagonal, starting in column 1 after the roll over from the second diagonal
            # This roll over occurs when row Nx*Ny-Nx hits the edge of the matrix starting from Nx
            # This will only happen when the number of rows is greater than Ny
            # Or in other words, when Nx > 1
            # There will only be Nx entries then

            if grid_size[1] > 1
                # Bloch theorem entry
                Λ_y = grid_size[2] * grid_resolution[2]
                entry_periodic_y = exp(1.0im*k_inc[2]*Λ_y)

                # Place at boundary
                for i = 1:grid_size[1]
                    DEY[i+grid_size[1]*grid_size[2]-grid_size[1],i] = entry_periodic_y*entry_y
                end
            end
        end
    end
    # Using copy to eagerly evaluate the transpose
    return DEX,DEY,copy(-transpose(DEX)),copy(-transpose(DEY))
end

function calculate_PML_2D(grid_size, PML_size)
    # Creates a Uniaxial PML layer that matches in size of the solution space grid_size
    # at the boundaries defined by PML_size
    # For lossy maxwell's equations ∇XE = k[μ_r][S]H

    # These variables setup the behavior of the PML layer
    a_max = 3
    p = 3 # Drop off rate
    σ_prime_max = 1 # Maximum conductivity
    η = 376.73031333108594
    # There are two cases, 2D and 3D - assume 2D #FIXME
    # Nx is number of columns, Ny is number of rows
    sx = fill(1.0+0.0im,grid_size) # Fill soulution space with 1s
    sy = fill(1.0+0.0im,grid_size)
    # Add xlow PML
    for nx = 1:PML_size[1]
        ax = 1 + a_max * (nx/PML_size[1])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[1])))^2
        sx[PML_size[1]-nx+1,:] = fill(ax * (1 + 1.0im * η * σ_prime_x),(grid_size[2],1))
    end
    # Add xhigh PML
    for nx = 1:PML_size[2]
        ax = 1 + a_max * (nx/PML_size[2])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[2])))^2
        sx[grid_size[1]-PML_size[2]+nx,:] = fill(ax * (1 + 1.0im * η * σ_prime_x),(grid_size[2],1))
    end
    # Add ylow PML
    for ny = 1:PML_size[3]
        ay = 1 + a_max * (ny/PML_size[3])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[3])))^2
        sy[:,PML_size[3]-ny+1] = fill(ay * (1 + 1.0im * η * σ_prime_y),(1,grid_size[1]))
    end
    # Add yhigh PML
    for ny = 1:PML_size[4]
        ay = 1 + a_max * (ny/PML_size[4])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[4])))^2
        sy[:,grid_size[2]-PML_size[4]+ny] = fill(ay * (1 + 1.0im * η * σ_prime_y),(1,grid_size[1]))
    end
    return sx,sy
end
