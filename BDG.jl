## Build Binary Defraction Grating Device ##
function BDG(resolution,NPML,highF)

    # Device Setup
    c_0 = 299792458
    λ = c_0/8e9
    λ_smallest = c_0/highF

    # Some parameters to set up device in physical dimension
    L = 0.6755*λ
    t = 0.0510*λ
    d = 0.2405*λ
    totalH = t+d
    x1 = 0.1040*λ
    x2 = 0.0175*λ
    x3 = 0.1080*λ

    # Device material
    ϵ_r = 1.0
    μ_r = 1.0

    # Padding Setup
    y_pad = λ
    y_PML = NPML # PML space to add on top and bottom on 2x grid

    # Physical size of device on the 1X grid
    X_Size = L
    Y_Size = totalH

    ## Optimize 1x grid size ##
    # How many cells does the smallest feature get? I'm picking 1
    # Pick smaller size, smallest wavelength or smallest feature
    initalRes = min(λ_smallest/(sqrt(ϵ_r)*Resolution),x2)
    ## Resolving Critical Dimensions ##
    N = ceil(Int64,L/initalRes)
    initalRes = L/N # So now L/N is an integer

    ## Calculate 1X Grid Size
    Nx = floor(Int64,L/initalRes)
    Ny = ceil(Int64,totalH/initalRes) + NPML + ceil(Int64,y_pad/initalRes)

    # Material Space Setup - 2X Grid
    dx2 = initalRes / 2
    dy2 = initalRes / 2
    Nx2 = Nx*2
    Ny2 = Ny*2

    # Initialize 2X Grid
    ϵ_r_2X_Grid = ones(Float64,Nx2,Ny2)
    μ_r_2X_Grid = ones(Float64,Nx2,Ny2)

    # Build substrate
    ny1 = ceil(Int64,Ny2/2 - t/(2*dy2) - d/(dy2) + totalH/(2*dy2))
    ny2 = ceil(Int64,Ny2/2 - t/(2*dy2) + totalH/(2*dy2))
    ny3 = ceil(Int64,Ny2/2 + t/(2*dy2) + totalH/(2*dy2))
    ϵ_r_2X_Grid[:,ny2:ny3] = fill(ϵ_r,Nx2,ny3-ny2+1)

    # Build fins
    nx1 = ceil(Int64,Nx2/2 - x2/(2*dx2) - x1/dx2)
    nx2 = ceil(Int64,Nx2/2 - x2/(2*dx2))
    nx3 = ceil(Int64,Nx2/2 + x2/(2*dx2))
    nx4 = ceil(Int64,Nx2/2 + x2/(2*dx2) + x3/dx2)

    ϵ_r_2X_Grid[nx1:nx2,ny1:ny2] = fill(ϵ_r,nx2-nx1+1,ny2-ny1+1)
    ϵ_r_2X_Grid[nx3:nx4,ny1:ny2] = fill(ϵ_r,nx4-nx3+1,ny2-ny1+1)

    # Decide where to build Q to on 1x grid
    Q = ceil(Int64,y_PML+3) # Build the Q mask 3 blocks from PML

    # Just for visualization # FIXME FIXME REMOVE
    #ϵ_r_2X_Grid[:,1:y_PML] = fill(ϵ_r,Nx2,y_PML)
    #ϵ_r_2X_Grid[:,end-y_PML+1:end] = fill(ϵ_r,Nx2,y_PML)

    return ϵ_r_2X_Grid',μ_r_2X_Grid',(y_PML,y_PML,y_PML,y_PML),Ny2,Nx2,[initalRes initalRes],Q
end
