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
