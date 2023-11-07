__precompile__(false)

module Syzygy
    using PrecompileTools

    function __init__()
        Unitful.register(Syzygy)
    end

    include("units.jl")
    include("constants.jl")
    
    include("setup/stellar_types.jl")
    
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


    export multibodysystem
    export PureGravitationalPotential, DynamicalTidalPotential, EquilibriumTidalPotential, StaticEquilibriumTidalPotential
    export simulation, simulate
    export analyse_simulation
    export 𝒢, pI, bI, ParticleIndex, BinaryIndex

    @compile_workload begin
        triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.4, 0.2])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision"])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision", "escape"])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision", "escape", "rlof", "tidal_disruption"])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision", "escape", "rlof", "tidal_disruption"],
                 potential=[PureGravitationalPotential(), DynamicalTidalPotential(G=𝒢.val, n=4, γ=[1.5, 1.5, 1.5])])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision", "escape", "rlof", "tidal_disruption"],
                 potential=[PureGravitationalPotential(), EquilibriumTidalPotential(𝒢.val)])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=["collision", "escape", "rlof", "tidal_disruption"],
                 potential=[PureGravitationalPotential(), StaticEquilibriumTidalPotential(triple)])

        sol = analyse_simulation(res)
    end

end
