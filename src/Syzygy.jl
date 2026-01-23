# __precompile__(false)

module Syzygy
    using PrecompileTools, Reexport
    @reexport using DynamicQuantities

    @register_unit Rsun 6.975e8u"m" 

    @register_unit Msun 1.9884754153381438e30u"kg" 

    @register_unit Lsun 3.828e26u"W" 

    @register_unit AU 149597870700u"m" 

    # if !(:Rsun in DynamicQuantities.UNIT_SYMBOLS) 
    #     @register_unit Rsun 6.975e8u"m" 
    # end
    # if !(:Msun in DynamicQuantities.UNIT_SYMBOLS) 
    #     @register_unit Msun 1.9884754153381438e30u"kg" 
    # end
    # if !(:Lsun in DynamicQuantities.UNIT_SYMBOLS) 
    #     @register_unit Lsun 3.828e26u"W" 
    # end
    # if !(:AU in DynamicQuantities.UNIT_SYMBOLS)
    #     @register_unit AU 149597870700u"m" 
    # end
    
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

    @compile_workload begin
        triple = multibodysystem([1.0, 1.0, 1.0]Msun, a=[0.1, 1.0]AU, e=[0.4, 0.2])
        res = simulate(triple, t_sim=0.1, save_everystep=false)

        threebody = multibodysystem(rand(3)u"kg", rand(3)u"m", rand(3)u"m/s", nbody_units=true)
        simulate(threebody, t_sim=0.1threebody.units.u_time)

        postprocess(res)
    end

end
