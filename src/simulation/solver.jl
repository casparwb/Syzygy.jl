

include("./setup.jl")
include("./callbacks.jl")

using OrdinaryDiffEq
using Unitful, UnitfulAstro, StaticArrays
using ProgressMeter


"""
    simulate(simulation::MultiBodySimulation)

Simulate the given simulation setup. 
"""
function simulate(simulation::MultiBodySimulation)

    args        = simulation.args
    diffeq_args = simulation.diffeq_args

    retcodes = Dict{Symbol, Any}(:Success => false)
    start_time = time()
    if !isinf(args[:max_cpu_time])
        push!(args[:callbacks], "max_cpu_time")
    end
    
    cbs = setup_callbacks(args[:callbacks], simulation.ic, simulation.params, 
                          retcodes, simulation.system.potential[:PureGravitationalPotential].G, args, start_time=start_time)
    callbacks = isnothing(cbs) ? nothing : CallbackSet(cbs...)


    ##############################################################################################################
    #              This block allows the full simulation to run without allocations. Don't know why.             #
    ##############################################################################################################
    let
        ode_prob_static = sodeprob_static(simulation)
        integrator_static = OrdinaryDiffEq.init(ode_prob_static, args[:alg], saveat=args[:saveat], 
                                                callback=callbacks, maxiters=args[:maxiters], 
                                                abstol=args[:abstol], reltol=args[:reltol], dt=args[:dt]; 
                                                diffeq_args...)
        step!(integrator_static)
        terminate!(integrator_static)

    end
    ##############################################################################################################
    
    acc_funcs = gather_accelerations_for_potentials(simulation)
    ode_problem = SecondOrderODEProblem(simulation, acc_funcs)

    integrator = OrdinaryDiffEq.init(ode_problem, args[:alg], saveat=args[:saveat], 
                                     callback=callbacks, maxiters=args[:maxiters], 
                                     abstol=args[:abstol], reltol=args[:reltol], dt=args[:dt]; 
                                     diffeq_args...)
    prog = ProgressUnknown("Evolving system:", showspeed=true, spinner=true, enabled=args[:showprogress])
    maxtime = simulation.tspan[end]
    try
        if args[:showprogress]
            for i in integrator
                next!(prog; showvalues=[(Symbol("System time"), u"kyr"(integrator.t * u"s")),
                                        (Symbol("System %"), (integrator.t - simulation.tspan[1])/maxtime*100)], 
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
    runtime =  (time() - start_time) * u"s"

    if args[:verbose]
        if retcodes[:Success]
            @info "Simulation successful."
        else
            outcome = collect(keys(retcodes))
            t = u"kyr"(integrator.t * u"s")
            @info "$outcome at t = " t
        end
    end

    result = SimulationResult(integrator.sol, simulation, 
                              retcodes, runtime, integrator.p,
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
    simulate(system::MultiBodyInitialConditions; args...)

Allows call to `simulate` directly without setting up simulation first.
"""
function simulate(system::MultiBodyInitialConditions; args...)
    sim = simulation(system; args...)
    if get(args, :verbose, :false)
        show(sim)
    end

    simulate(simulation(system; args...))
end

function get_masses(system::MultiBodyODESystem)
    return [b.mass for b in system.bodies]
end

function total_energy(result::SimulationResult, time)
    masses = get_masses(result.simulation.system)
    idx = argmin(abs.(result.solution.t .- time))

    total_energy([result.solution.u[idx].x[2][:,i] for i ∈ eachindex(masses)], 
                 [result.solution.u[idx].x[1][:,i] for i ∈ eachindex(masses)],
                 masses)  
end
