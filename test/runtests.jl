using Syzygy
using Test

@testset "Syzygy.jl" begin
    include("initialization_tests.jl")
    include("initial_condition_tests.jl")
    include("solver_tests.jl")

    
end
