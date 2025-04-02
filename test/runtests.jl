using Syzygy
using Test
# using LinearAlgebra
@testset "Syzygy.jl" begin
    include("initialization_tests.jl")
    include("initial_condition_tests.jl")
    include("solver_tests.jl")
    include("potential_pair_tests.jl")
    include("potential_value_tests.jl")

end
