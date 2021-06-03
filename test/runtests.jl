using SafeTestsets

@safetestset "Parametric entity" begin include("parametricentity_test.jl") end

@safetestset "Simple shapes" begin include("simpleshapes_test.jl") end

@safetestset "Parametric elements" begin include("parametricelement_test.jl") end

@safetestset "Mesh" begin include("mesh_test.jl") end
