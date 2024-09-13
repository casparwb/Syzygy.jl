__precompile__(false)

module Syzygy
    using PrecompileTools, Reexport 
    @reexport using Unitful, UnitfulAstro

    Unitful.preferunits(u"Rsun, Msun, yr"...)

    # function __init__()
    #     Unitful.register(Syzygy)
    # end

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

    export multibodysystem
    export PureGravitationalPotential, DynamicalTidalPotential, 
           EquilibriumTidalPotential, StaticEquilibriumTidalPotential
    export simulation, simulate
    export to_solution
    export GRAVCONST
    
    export CollisionCB, EscapeCB, RocheLobeOverflowCB, CPUTimeCB, 
           CentreOfMassCB, HubbleTimeCB, DemocraticCheckCB, IonizationCB
   

    # @compile_workload begin
    #     triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.4, 0.2])
    #     simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #              callbacks=[])
    #     simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #              callbacks=[CollisionCB()])
    #     simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #              callbacks=[CollisionCB(), EscapeCB(100, 100)])
    #     res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #              callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
    #     #          potential=[PureGravitationalPotential(), DynamicalTidalPotential(G=GRAVCONST.val, n=4, γ=[1.5, 1.5, 1.5])])
    #     # simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
    #     #          potential=[PureGravitationalPotential(), EquilibriumTidalPotential(GRAVCONST.val)])
    #     # res = simulate(triple, t_sim=10, save_everystep=false, showprogress=false, 
    #     #          callbacks=[CollisionCB(), EscapeCB(100, 100), RocheLobeOverflowCB(100)],
    #     #          potential=[PureGravitationalPotential(), StaticEquilibriumTidalPotential(triple)])

    #     sol = to_solution(res)
    # end

end
