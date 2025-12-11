__precompile__(false)

module Syzygy
    using PrecompileTools, Reexport, Preferences
    @reexport using DynamicQuantities

    if !(:Rsun in DynamicQuantities.UNIT_SYMBOLS) 
        @register_unit Rsun 6.975e8u"m" 
    end
    if !(:Msun in DynamicQuantities.UNIT_SYMBOLS) 
        @register_unit Msun 1.9884754153381438e30u"kg" 
    end
    if !(:Lsun in DynamicQuantities.UNIT_SYMBOLS) 
        @register_unit Lsun 3.828e26u"W" 
    end
    if !(:AU in DynamicQuantities.UNIT_SYMBOLS)
        @register_unit AU 149597870700u"m" 
    end

    const default_unit_length, default_unit_mass, default_unit_time = u"m", u"kg", u"s"

    include("constants.jl")


    include("ode_solvers.jl")
    
    include("setup/stellar_types.jl")
    include("setup/potentials.jl")
    include("setup/framework.jl")
    
    include("physics/physics.jl")
    include("physics/orbital_elements.jl")
    include("physics/hierarchy_setup.jl")
    include("physics/kepler.jl")
    
    include("setup/initialization.jl")
    
    include("simulation/solver.jl")
    
    include("analysis/postprocessing.jl")
    include("analysis/visualization.jl")
    include("io.jl")
    

    public centre_of_mass, center_of_mass, centre_of_mass_velocity, center_of_mass_velocity, potential_energy, 
           kinetic_energy, specific_orbital_energy, reduced_mass, gravitational_radius,
           schwarzschild_radius, roche_radius, octupole_parameter, is_unstable,
           stability_criterion_ma01, eccentricity

    using .ODEAlgorithms
    export ODESolvers

    export multibodysystem
    export PureGravitationalPotential, DynamicalTidalPotential, 
           TimeDependentEquilibriumTidalPotential, EquilibriumTidalPotential,
           PN1Potential, PN2Potential, PN2p5Potential, PNPotential
    export TidalSimulationParams

    export simulation, simulate
    export to_solution
    export GRAVCONST
    
    export CollisionCB, EscapeCB, RocheLobeOverflowCB, CPUTimeCB, SavingCB,
           CentreOfMassCB, HubbleTimeCB, DemocraticCheckCB, IonizationCB

    const postprocess = to_solution
    export postprocess
   
    export Units

    # @compile_workload begin
    #     triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.4, 0.2])
    #     res = simulate(triple, t_sim=1, save_everystep=false)
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
    #     #          callbacks=[CollisionCB()])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100)])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100)], potential=[PureGravitationalPotential(), PNPotential()])
    #     # res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
    #     #          potential=[PureGravitationalPotential(), DynamicalTidalPotential(4, [1.5, 1.5, 1.5])])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
    #     #          potential=[PureGravitationalPotential(), EquilibriumTidalPotential(triple)])

    #     sol = to_solution(res)
    # end

end
