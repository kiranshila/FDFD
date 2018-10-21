import Base.floor
using SparseArrays

function boundingBox(mesh::HomogenousMesh)
    matrix_verts = hcat(convert(Array{Array{Float32}},mesh.vertices)...)
    [[minimum(matrix_verts[1,:]),maximum(matrix_verts[1,:])],
     [minimum(matrix_verts[2,:]),maximum(matrix_verts[2,:])],
     [minimum(matrix_verts[3,:]),maximum(matrix_verts[3,:])]]
end

function boundingBox(triangle::Array{Float64,2})
    [[minimum(triangle[:,1]) maximum(triangle[:,1])];
     [minimum(triangle[:,2]) maximum(triangle[:,2])];
     [minimum(triangle[:,3]) maximum(triangle[:,3])]]
end

function xyzSize(limits)
    return (abs(limits[1][2] - limits[1][1]),
    abs(limits[2][2] - limits[2][1]),
    abs(limits[3][2] - limits[3][1]))
end

struct tomasCube
    C::Vec{3,Float64}
    halfSize::Float64
end

struct planeNormDist
    normal::Array{Float64,1}
    distanceFromOrigin::Float64
end

function AABBPlaneIntersection(AABB::tomasCube,plane::planeNormDist)
    # Compute the projection interval radius of b onto L(t) = box_center + t * plane_normal
    R = AABB.halfSize*(abs(plane.normal[1]) + abs(plane.normal[2]) + abs(plane.normal[3]))
    # Intersection occurs when distance s falls within [-r,+r] interval
    return abs(plane.distanceFromOrigin) <= R
end

function triangleIntersectsCube(triangle::Array{Float64,2},cube::tomasCube)
    ## Tomas Akenine-Moller Algorithm ##
    # 0 - Move box to test and triangle such that the box is centered at the origin
    x = cube.C[1]
    y = cube.C[2]
    z = cube.C[3]
    triangle = triangle - [x y z;x y z;x y z]; # Translation Matrix assuming triangle is [x y z; x y z; x y z]
    F = [triangle[2,:]-triangle[1,:],triangle[3,:]-triangle[2,:],triangle[1,:]-triangle[3,:]]

    # So now, box extents are
    E = [cube.halfSize 0.0 0.0;0.0 cube.halfSize 0.0; 0.0 0.0 cube.halfSize]

    # 1 - Test overlaping bounding box in x,y,z with min and max of the triangle vertices
    for i = 1:3 # For each axis
        max = findmax(triangle[:,i])[1] # Grab value of max for this axis
        min = findmin(triangle[:,i])[1] # Grab value of min for this axis
        if min > cube.halfSize || max < -cube.halfSize
            return false # The triangle's bounding box isn't within the dimensions of the test cube
        end
    end

    # 2 - Test if the box intersects the plane of the trianle
    tNorm = normalize(cross(F[1],F[2]))
    planeD = dot(triangle[1,:],tNorm)
    trianglePlane = planeNormDist(tNorm,planeD)
    if !AABBPlaneIntersection(cube,trianglePlane)
        return false
    end

    # Find the normal of the triangle

    # 3 - Auxilary Tests
    for i = 1:3
        for j = 1:3
            # a is the vector normal to one axis of the cube and an edge of the triangle
            a = cross(E[i,:],F[j])
            # Project the triangle's verticies onto a
            p = [dot(a,triangle[1,:]),dot(a,triangle[2,:]),dot(a,triangle[3,:])]
            # Find a box "radius"
            r = cube.halfSize*(abs(a[1])+abs(a[2])+abs(a[3]))
            min = findmin(p)[1]
            max = findmax(p)[1]
            if min > r || max < -r
                return false
            end
        end
    end
    # If all of these passed
    return true
end

@inline function cielBB(input::Array{Float64,2})
    output = zeros(Int64, 3, 2)
    for i = 1:3
        for j = 1:2
            output[i,j] = ceil(Int,input[i,j])
        end
    end
    output
end

function meshDevice!(cubes::Array{YeeSubCube,3},myMesh::HomogenousMesh,cubeSize::Float64,material::Float64)
    # 1 - Center device in solution space - Assuming the mesh was centered about (0,0,0)
    dims = size(cubes)
    numCubesPlotted = 0
    numCubesTested = 0
    x,y,z = dims .* cubeSize .*0.5
    z = 0.0 # These are the x,y,z positions of the center of the model
    # 2 - For every triangle
    for i = 1:3:length(myMesh.vertices)
        triangle = hcat(convert(Array{Float64,1},myMesh.vertices[i]),
            convert(Array{Float64,1},myMesh.vertices[i+1]),
            convert(Array{Float64,1},myMesh.vertices[i+2]))'
        triangle = triangle + [x y z;x y z;x y z]; # Translation Matrix assuming triangle is [x y z; x y z; x y z]
        #   - For every cube that has the potential to be intersecting with the triangle
        #FIXME
        triangleBBIndex = map(x->ceil(Int64,x+eps(Float64)),boundingBox(triangle) ./ cubeSize)
        for i = triangleBBIndex[1,1]:triangleBBIndex[1,2]
            for j = triangleBBIndex[2,1]:triangleBBIndex[2,2]
                for k = triangleBBIndex[3,1]:triangleBBIndex[3,2]
                    #   - Check if cube intersects triangle, if it does, apply material property
                    # Construct tomas cube from YeeCube
                    thisCube = tomasCube(cubeSize .* (i, j, k) ./ 2.0, cubeSize/2.0)
                    numCubesTested = numCubesTested + 1
                    #if triangleIntersectsCube(triangle,thisCube)
                        cubes[i,j,k] = YeeSubCube(E,material,1.0) # Modify given array
                        numCubesPlotted = numCubesPlotted + 1
                    #end
                end
            end
        end
    end
end

function calculate_PML_2D(grid_size, PML_size)
    # Creates a Uniaxial PML layer that matches in size of the solution space grid_size
    # at the boundaries defined by PML_size
    # For lossy maxwell's equations ∇XE = k[μ_r][S]H

    # These variables setup the behavior of the PML layer
    a_max = 3
    p = 3 # Drop off rate
    σ_prime_max = 1 # Maximum conductivity

    # There are two cases, 2D and 3D - assume 2D #FIXME
    sx = fill(1.0+0.0im,grid_size) # Fill soulution space with 1s
    sy = fill(1.0+0.0im,grid_size)
    # Add xlow PML
    for nx = 1:PML_size[1]
        ax = 1 + a_max * (nx/PML_size[1])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[1])))^2
        sx[PML_size[1]-nx+1,:] = fill(ax * (1 + 1.0im * η_0 * σ_prime_x),(grid_size[2],1))
    end
    # Add xhigh PML
    for nx = 1:PML_size[2]
        ax = 1 + a_max * (nx/PML_size[2])^p
        σ_prime_x = σ_prime_max * (sin((pi*nx)/(2*PML_size[2])))^2
        sx[grid_size[1]-PML_size[2]+nx,:] = fill(ax * (1 + 1.0im * η_0 * σ_prime_x),(grid_size[2],1))
    end
    # Add ylow PML
    for ny = 1:PML_size[3]
        ay = 1 + a_max * (ny/PML_size[3])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[3])))^2
        sy[:,PML_size[3]-ny+1] = fill(ay * (1 + 1.0im * η_0 * σ_prime_y),(1,grid_size[1]))
    end
    # Add yhigh PML
    for ny = 1:PML_size[4]
        ay = 1 + a_max * (ny/PML_size[4])^p
        σ_prime_y = σ_prime_max * (sin((pi*ny)/(2*PML_size[4])))^2
        sy[:,grid_size[2]-PML_size[4]+ny] = fill(ay * (1 + 1.0im * η_0 * σ_prime_y),(1,grid_size[1]))
    end
    return sy
end

@enum BC Periodic=1 Dirichlet=2

function yee_grid_derivative(grid_size,grid_resolution,boundary_condition::Tuple{BC,BC},k_inc = [0 0])
    # Generates the Yee Grid Derivative on the 2D Grid #FIXME for 3D


    # To create a sparse array, we have to make I a vector of row idicies,
    # J a vector of column indicies, and V a vector of values
    # We do this to create our large derivative array for the two boundary condition cases

    # Storing results in a sparse array of row index I, column index J, and value V
    I = Int64[]
    J = Int64[]
    V = Complex{Float64}[]
    DEX = nothing
    DEY = nothing

    # Check if Nx = 1
    if grid_size[1] == 1
        DEX = 1.0im*k_inc[1]*sparse(I,grid_size[1]*grid_size[2])
    else
        # Every element is multiplied by 1/Δx
        entry_x = 1/grid_resolution[1]
        for i = 1:grid_size[1]*grid_size[2]
            # Main Diagonal
            push!(I,i)
            push!(J,i)
            push!(V,-1)
        end
        for i = 1:grid_size[1]*grid_size[2]-1
            # Secondary Diagonal
            push!(I,i)
            push!(J,i+1)
            push!(V,1)
        end
        DEX = sparse(I,J,V)
        if boundary_condition[1] == Periodic
            Λ_x = grid_size[1] * grid_resolution[1]
            entry_periodic_x = exp(1.0im*k_inc[1]*Λ_x)
            for i = grid_size[1]:grid_size[1]:grid_size[1]*grid_size[2]-1
                DEX[i,Int(i/grid_size[1])] = entry_periodic_x
                DEX[i,i+1] = 0
            end
        end
    end

    # Clear Variables
    I = Int64[]
    J = Int64[]
    V = Float64[]

    # Check if Ny = 1
    if grid_size[1] == 1
        DEY = 1.0im*k_inc[2]*sparse(I,grid_size[1]*grid_size[2])
    else
        # Every element is multiplied by 1/Δx
        entry_y = 1/grid_resolution[2]
        for i = 1:grid_size[1]*grid_size[2]
            # Main Diagonal
            push!(I,i)
            push!(J,i)
            push!(V,-1)
        end
        for i = 1:grid_size[1]*grid_size[2]-1
            # Secondary Diagonal
            push!(I,i)
            push!(J,i+1)
            push!(V,1)
        end
        if boundary_condition[2] == Periodic
            # FIXME
        end
        DEY = sparse(I,J,V)
    end
    return DEX,DEY
end
