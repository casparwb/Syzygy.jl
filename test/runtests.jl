using Syzygy
using Test
# using LinearAlgebra
@testset "Syzygy.jl" begin
    include("initialization_tests.jl")
    include("initial_condition_tests.jl")
    include("solver_tests.jl")
    include("potential_tests.jl")
    
end
