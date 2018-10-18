using FileIO, GeometryTypes
using Makie
using LinearAlgebra
include("Yee.jl")
include("Utilities.jl")

# Set constant scene
#scene = Scene(resolution = (2048, 2048))

# Load STL File of Geometry
geometryMesh = load("Lens.stl")
#geometryMesh = Makie.loadasset("cat.obj")

# Set units of STL
stl_units = 1 # milimeters
# Extract vertices from mesh data
verts = geometryMesh.vertices

## Solver setup ##
# Voxel Resolution
resolution = 50
# Solution Frequency
frequency = 38e9
# Solution Space Scale
scale = 3
# Dielectric of model
ϵ_r = 2.2;

# Calculations of voxel dimensions and space setup
c_0 = 299792458.0
λ_0 = c_0/frequency
Δ_length = λ_0/resolution
geometry_size = boundingBox(geometryMesh)
solution_size = scale .* geometry_size .* stl_units
solution_indices = map(x -> floor(Int, x),xyzSize(solution_size) ./ Δ_length)

# Fill space with vacuum
solution_space = fill(YeeSubCube(E,1.0,1.0),solution_indices)

# Visualize Mesh
mesh(geometryMesh)

# Rasterize voxel grid to setup device in solution space
@time meshDevice!(solution_space,geometryMesh,Δ_length,ϵ_r)

# Visualize cubes
plotYee(solution_space,Δ_length)
