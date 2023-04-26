# __precompile__()

module FewBodySimulator

    include("units.jl")
    include("constants.jl")
    
    include("setup/potentials.jl")
    include("setup/framework.jl")
    
    include("physics/physics.jl")
    include("physics/orbital_elements.jl")
    include("physics/hierarchy_setup.jl")
    include("physics/kepler.jl")
    
    include("setup/initialization.jl")
    
    include("simulation/solver.jl")
    
    include("analysis/simulation_analysis.jl")
    include("analysis/visualization.jl")

    # EquilibriumTidalPotential(system::MultiBodySystem; G, τ) = EquilibriumTidalPotential([u"Msun"(p.mass).val for p in system.particles], G=G, τ=τ)

    export multibodysystem, getparticle, getbinary, getparent, getbinaries
    export PureGravitationalPotential, DynamicalTidalPotential, EquilibriumTidalPotential
    export simulation, simulate
    export analyse_simulation
    export 𝒢, pI, bI, ParticleIndex, BinaryIndex

end
