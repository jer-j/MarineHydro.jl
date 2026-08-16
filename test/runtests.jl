using Test

@testset "MarineHydro.jl test suite" begin
    include("./benchmark_meshes.jl")
    include("./maneuvering.jl")
    include("./boundary_layers.jl")
    include("./three_dimensional_boundary_layer.jl")
    include("./maneuvering_validation.jl")
    include("./consistency_with_Capytaine.jl")
    include("./greens_function.jl")
    include("./greens_function_differentiation.jl")
    include("./matrix_assembly.jl")
    include("./matrix_assembly_differentiation.jl")
    include("./consistency_with_analytical_solutions.jl")
    include("./outputs_differentiation.jl")
    include("./forward_speed_tests.jl")
end
