# __precompile__(false)

module Syzygy
    using PrecompileTools, Reexport, Preferences
    @reexport using Unitful, UnitfulAstro
    include("units.jl")

    function __init__()
        Unitful.register(Syzygy)
        merge!(Unitful.promotion, localpromotion)
    end


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
    
#     include("simulation/callbacks.jl")
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
           EquilibriumTidalPotential, StaticEquilibriumTidalPotential,
           PN1Potential, PN2Potential, PN2p5Potential, PNPotential,
           PN1SpinPrecessionPotential, PN1p5SpinPrecessionPotential, PN2SpinPrecessionPotential,
           SpinPrecessionPotential

    export simulation, simulate
    export to_solution
    export GRAVCONST
    
    export CollisionCB, EscapeCB, RocheLobeOverflowCB, CPUTimeCB, SavingCB,
           CentreOfMassCB, HubbleTimeCB, DemocraticCheckCB, IonizationCB
   

    @compile_workload begin
        triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.4, 0.2])
        res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
                 callbacks=[CollisionCB()])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
                 callbacks=[CollisionCB(), EscapeCB(100, 100)])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, max_cpu_time=1,
                 callbacks=[CollisionCB(), EscapeCB(100, 100)], potential=[PureGravitationalPotential(), PNPotential()])
        # res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)])
        # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
        #          potential=[PureGravitationalPotential(), DynamicalTidalPotential(4, [1.5, 1.5, 1.5])])
        # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
        #          potential=[PureGravitationalPotential(), StaticEquilibriumTidalPotential(triple)])

        sol = to_solution(res)
    end

end
