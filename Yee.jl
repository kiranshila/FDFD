@enum Field E=1 H=2

struct YeeSubCube
    ThisField::Field
    ϵ_r::Float64
    μ_r::Float64
end

struct YeeCube
    Ex::YeeSubCube
    Ey::YeeSubCube
    Ez::YeeSubCube
    Hx::YeeSubCube
    Hy::YeeSubCube
    Hz::YeeSubCube
end

function plotYee(cubes::Array{YeeSubCube,3},sizeOfCube)
    l,w,h = size(cubes)
    positions = [sizeOfCube*Vec3f0(i, j, k) for i=1:l, j=1:w, k=1:h if cubes[i,j,k].ϵ_r > 1]
    println("$(length(positions)) cubes plotted.")
    # Iterate through YeeCubes, add to position vector if ϵ_r > 1
    scatter(positions,
        color = RGBAf0(1.0, 0.0, 0.0, 0.5),
        marker = FRect3D(Vec3f0(0),
        Vec3f0(Δ_length)),
        show_axis = true,
        markersize=1,
        limits = FRect3D(Vec3f0(0),
        size(cubes) .* Δ_length))
end
