
module ODEAlgorithms

export ODESolvers

struct SyzygyODESolvers{T}
    solvers::T
end

function Base.getproperty(s::SyzygyODESolvers, p::Symbol)
    solvers = Base.getfield(s, :solvers)
    if haskey(solvers, p)
        return solvers[p]
    else 
        Base.getfield(s, p)
    end
end


import OrdinaryDiffEqRKN, OrdinaryDiffEqSymplecticRK
const core = OrdinaryDiffEqRKN.OrdinaryDiffEqCore

RKN_solvers = Symbol[]
for p in propertynames(OrdinaryDiffEqRKN)
    prop = getproperty(OrdinaryDiffEqRKN, p)
    prop_type = typeof(prop)

    try
        if prop <: core.OrdinaryDiffEqAlgorithm
            push!(RKN_solvers, p)
        end
    catch e
        continue
    end
end

symplectic_RK_solvers = Symbol[]
for p in propertynames(OrdinaryDiffEqSymplecticRK)
    prop = getproperty(OrdinaryDiffEqSymplecticRK, p)
    prop_type = typeof(prop)

    try
        if prop <: core.OrdinaryDiffEqAlgorithm
            push!(symplectic_RK_solvers, p)
        end
    catch e
        continue
    end
end


RKN_solvers = Dict(s => getproperty(OrdinaryDiffEqRKN, s)() for s in RKN_solvers)
symplectic_RK_solvers = Dict(s => getproperty(OrdinaryDiffEqSymplecticRK, s)() for s in symplectic_RK_solvers)

solvers = merge(RKN_solvers, symplectic_RK_solvers)

const ODESolvers = SyzygyODESolvers(solvers)

end