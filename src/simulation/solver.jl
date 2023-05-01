

include("./setup.jl")
include("./callbacks.jl")

using OrdinaryDiffEq
using Unitful, UnitfulAstro, StaticArrays
using ProgressMeter


"""
    simulate(simulation::FewBodySimulation)

Simulate the given simulation setup. 
"""
function simulate(simulation::FewBodySimulation)

    args        = simulation.args
    diffeq_args = simulation.diffeq_args

    # retcode = Set()
    retcodes = Dict{Symbol, Any}(:Success => false)
    
    cbs = setup_callbacks(args[:callbacks], simulation.ic, simulation.params, 
                          retcodes, simulation.system.potential[:PureGravitationalPotential].G, args)
    callbacks = isnothing(cbs) ? nothing : CallbackSet(cbs...)
    
    ode_problem = SecondOrderODEProblem(simulation)
    integrator = OrdinaryDiffEq.init(ode_problem, args[:alg], saveat=args[:saveat], 
                                     callback=callbacks, maxiters=args[:maxiters], 
                                     abstol=args[:abstol], reltol=args[:reltol], dt=args[:dt]; 
                                     diffeq_args...)

    
    prog = ProgressUnknown("Evolving system:", showspeed=true, spinner=true, enabled=args[:showprogress])


    maxtime = simulation.tspan[end]
    runtime = time()
    try
        if args[:showprogress]
            # println("hello")
            for i in integrator
                next!(prog; showvalues=[(Symbol("System time"), u"kyr"(integrator.t * u"s")),
                                        (Symbol("System %"), integrator.t/maxtime*100)], 
                                        spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
            end
            solve!(integrator)
        else
            solve!(integrator)
        end
        retcodes[:Success] = true
        retcodes[:DiffEq] = integrator.sol.retcode
    catch e
        if e isa InterruptException
            @info "Stopped at t = $(u"kyr"(integrator.t * u"s"))"
            
        else
            throw(e)
            retcodes[:DiffEq] = integrator.sol.retcode
            # retcode = []
        end
    end

    ProgressMeter.finish!(prog)
    runtime =  (time() - runtime) * u"s"

    if args[:verbose]#retcode[1][1] != :Success &&
        if retcodes[:Success]
            @info "Simulation successful."
        else
            outcome = collect(keys(retcodes))
            t = u"kyr"(integrator.t * u"s")
            @info "$outcome at t = " t
        end
    end

    result = SimulationResult(integrator.sol, simulation, 
                              retcodes, runtime, 
                              merge(args, diffeq_args))

    if args[:verbose]
        println()

        E0 = total_energy(result, simulation.tspan[1])
        E1 = total_energy(result, simulation.tspan[2])

        @info "Energy loss: " 1 - E1/E0

        runtime = runtime
        @info "Total runtime: " runtime
        println()
    end

    return result
end

"""
    simulate(system::FewBodyInitialConditions; args...)

Allows call to `simulate` directly without setting up simulation first.
"""
function simulate(system::FewBodyInitialConditions; args...)
    sim = simulation(system; args...)
    if get(args, :verbose, :false)
        show(sim)
    end

    simulate(simulation(system; args...))
end

function get_masses(system::FewBodySystem)
    return [b.mass for b in system.bodies]
end

function total_energy(result::SimulationResult, time)
    masses = get_masses(result.simulation.system)
    idx = argmin(abs.(result.solution.t .- time))

    total_energy([result.solution.u[idx].x[2][:,i] for i ∈ eachindex(masses)], 
                 [result.solution.u[idx].x[1][:,i] for i ∈ eachindex(masses)],
                 masses)  
end
