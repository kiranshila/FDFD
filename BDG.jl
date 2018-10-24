## Build Binary Defraction Grating Device ##
function BDG(resolution)

    # Solution Setup
    λ_0 = c_0/8e9
    RES = [λ_0/resolution,λ_0/resolution] # This is the grid resolution on 2X

    # Some parameters to set up device in physical dimension
    L = 0.6755*λ_0
    t = 0.0510*λ_0
    d = 0.2405*λ_0
    totalH = t+d
    x1 = 0.1040*λ_0
    x2 = 0.0175*λ_0
    x3 = 0.1080*λ_0

    # Device material
    ϵ_r = 10.0
    μ_r = 1.0

    # Padding Setup
    y_pad = λ_0
    y_PML = 6 # PML space to add on top and bottom

    # Material Space Setup
    Nx = ceil(Int64,L/RES[1]) # Just the length of the device
    Ny = ceil(Int64,(totalH + 2*y_pad)/RES[2]) + y_PML*2 # Total height of device, plus 1 lambda spacer region, + PML
    NGRID = (Nx,Ny)
    # This is the entire solution space in number of boxes in the "1x" grid, Nx is rows, Ny is columns

    # Setup calculations
    Nx2 = 2*Nx
    Ny2 = 2*Ny
    dx2 = RES[1]/2
    dy2 = RES[2]/2

    # Initialize 2X Grid
    ϵ_r_2X_Grid = ones(Float64,Nx2,Ny2)
    μ_r_2X_Grid = ones(Float64,Nx2,Ny2)

    # Build substrate
    ny2 = ceil(Int64,Ny2/2 - t/(2*dy2))
    ny3 = ceil(Int64,Ny2/2 + t/(2*dy2))
    ϵ_r_2X_Grid[:,ny2:ny3] = fill(ϵ_r,Nx2,ny3-ny2+1)

    # Build fins
    ny1 = ceil(Int64,Ny2/2 - t/(2*dy2) - d/(dy2))
    nx1 = ceil(Int64,Nx2/2 - x2/(2*dx2) - x1/dx2)
    nx2 = ceil(Int64,Nx2/2 - x2/(2*dx2))
    nx3 = ceil(Int64,Nx2/2 + x2/(2*dx2))
    nx4 = ceil(Int64,Nx2/2 + x2/(2*dx2) + x3/dx2)

    ϵ_r_2X_Grid[nx1:nx2,ny1:ny2] = fill(ϵ_r,nx2-nx1+1,ny2-ny1+1)
    ϵ_r_2X_Grid[nx3:nx4,ny1:ny2] = fill(ϵ_r,nx4-nx3+1,ny2-ny1+1)

    # Decide where to build Q to
    Q = ceil(Int64,y_PML + y_pad/(4*dy2)) # Build the Q mask up to halfway into the padding regions

    # Just for visualization # FIXME FIXME REMOVE
    ϵ_r_2X_Grid[:,1:y_PML] = fill(ϵ_r,Nx2,y_PML)
    ϵ_r_2X_Grid[:,end-y_PML+1:end] = fill(ϵ_r,Nx2,y_PML)

    return ϵ_r_2X_Grid',μ_r_2X_Grid',reverse(NGRID),(20,20,0,0),Ny2,Nx2,reverse(RES),Q
end
