# __precompile__(false)

module Syzygy
    using PrecompileTools, Reexport, Preferences
    @reexport using Unitful, UnitfulAstro

    function set_units(units)
        allowed_units = ["SI", "Solar", "CGS"]

        if !(units in allowed_units)
            throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
        end

        @set_preferences!("units" => units)
        @info "New default units set; restart your Julia session for this change to take effect!"
    end

    function units_from_unit_system(unit_system)
        allowed_units = ["SI", "Solar", "SI"]
        units = if unit_system == "SI"
            u"kg, m, s"
        elseif unit_system == "CGS"
            u"g, cm, s"
        elseif unit_system == "Solar"
            u"Rsun, Msun, yr"
        else
            throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
        end

        return units
    end

    const unit_system = @load_preference("units", "Solar")
    const units = units_from_unit_system(unit_system)

    Unitful.preferunits(units...)
    const localpromotion = copy(Unitful.promotion)
    function __init__()
        Unitful.register(Syzygy)
        merge!(Unitful.promotion, localpromotion)
    end


    # include("units.jl")
    include("constants.jl")
    
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
    

    public centre_of_mass, center_of_mass, centre_of_mass_velocity, potential_energy, 
           kinetic_energy, specific_orbital_energy, reduced_mass, gravitational_radius,
           schwarzschild_radius, roche_radius, octupole_parameter, is_unstable,
           stability_criterion_ma01, eccentricity

    export multibodysystem
    export PureGravitationalPotential, DynamicalTidalPotential, 
           EquilibriumTidalPotential, StaticEquilibriumTidalPotential,
           PN1Potential, PN2Potential, PN2p5Potential, PNPotential,
           PN1SpinPrecessionPotential, PN1p5SpinPrecessionPotential, PN2SpinPrecessionPotential,
           SpinPrecessionPotential

    export simulation, simulate
    export to_solution
    export GRAVCONST
    
    export CollisionCB, EscapeCB, RocheLobeOverflowCB, CPUTimeCB, 
           CentreOfMassCB, HubbleTimeCB, DemocraticCheckCB, IonizationCB
   

    @compile_workload begin
        triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.4, 0.2])
        res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[CollisionCB()])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[CollisionCB(), EscapeCB(100, 100)])
        simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
                 callbacks=[CollisionCB(), EscapeCB(100, 100)], potential=[PureGravitationalPotential(), PNPotential()])
        # res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)])
        # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
        #          potential=[PureGravitationalPotential(), DynamicalTidalPotential(G=GRAVCONST.val, n=4, γ=[1.5, 1.5, 1.5])])
        # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
        #          potential=[PureGravitationalPotential(), EquilibriumTidalPotential(GRAVCONST.val)])
        # res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
        #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
        #          potential=[PureGravitationalPotential(), StaticEquilibriumTidalPotential(triple)])

        sol = to_solution(res)
    end

end
