include("./callbacks.jl")
include("./setup.jl")

using DiffEqBase, StaticArrays, ProgressMeter

function initialize_integrator(simulation)
    acc_funcs = gather_accelerations_for_potentials(simulation)

    args = simulation.args
    diffeq_args = simulation.diffeq_args
    ode_problem = SecondOrderODEProblem(simulation, acc_funcs, args[:dtype])

    retcodes = Dict{Symbol, Any}()

    if !isinf(args[:max_cpu_time]) && !(CPUTimeCB in typeof.(args[:callbacks]))
        push!(args[:callbacks], CPUTimeCB())
    end

    cbs = setup_callbacks(
        args[:callbacks],
        simulation.ic,
        simulation.params,
        retcodes, args
    )

    callbacks = isnothing(cbs) ? nothing : CallbackSet(cbs...)

    integrator = DiffEqBase.init(
        ode_problem, args[:alg], saveat = args[:saveat],
        callback = callbacks, maxiters = args[:maxiters],
        abstol = args[:abstol], reltol = args[:reltol], dt = args[:dt];
        diffeq_args...
    )

    return integrator, retcodes
end


"""
    simulate(simulation::MultiBodySimulation)

Simulate the given simulation setup. 
"""
function simulate(simulation::MultiBodySimulation)

    unit_length, unit_mass, unit_time = simulation.ic.units.u_length, simulation.ic.units.u_mass, simulation.ic.units.u_time

    args = simulation.args
    diffeq_args = simulation.diffeq_args

    # ##############################################################################################################
    # #              This block allows the full simulation to run without allocations. Don't know why.             #
    # ##############################################################################################################
    if simulation.ic.n < 10
        let
            ode_prob_static = sodeprob_static(simulation, args[:dtype])

            integrator_static = DiffEqBase.init(
                ode_prob_static, args[:alg], saveat = args[:saveat], maxiters = args[:maxiters],
                abstol = args[:abstol], reltol = args[:reltol], dt = args[:dt];
                diffeq_args...
            )
            try
                step!(integrator_static)
            catch err
                nothing
            finally
                terminate!(integrator_static)
            end
        end
    end
    # ##############################################################################################################

    integrator, retcodes = initialize_integrator(simulation)

    start_time = time()

    prog = ProgressUnknown(desc = "Evolving system:", showspeed = true, spinner = true, enabled = args[:showprogress])
    maxtime = simulation.tspan[end]
    try
        if args[:showprogress]
            for i in integrator
                next!(
                    prog; showvalues = [
                        (Symbol("System time"), integrator.t * unit_time),
                        (Symbol("System %"), (integrator.t - simulation.tspan[1]) / maxtime * 100),
                    ],
                    spinner = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
                )
            end
            solve!(integrator)
        else
            solve!(integrator)
        end

        retcodes[:DiffEq] = Symbol("$(integrator.sol.retcode)")
    catch e
        if e isa InterruptException
            @info "Stopped at t = $(integrator.t * unit_time)"
            savevalues!(integrator, true)
            retcodes[:DiffEq] = :Interrupted
        else
            throw(e)
        end
    end


    ProgressMeter.finish!(prog)
    runtime = (time() - start_time) * u"s"

    if args[:verbose]
        if retcodes[:DiffEq] == :Success
            @info "Simulation successful."
        else
            outcome = retcodes
            t = integrator.t * unit_time
            @info "Outcome: " t outcome
        end
    end

    result = SimulationResult(
        integrator.sol, simulation,
        retcodes, runtime, integrator.p,
        merge(args, diffeq_args)
    )

    if args[:verbose]
        println()

        E0 = total_energy(result, result.solution.t[1])
        E1 = total_energy(result, result.solution.t[end])

        @info "Energy loss: " E1 / E0 - 1

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

    return simulate(simulation(system; args...))
end

function simulate(res::SimulationResult, time; args...)
    t_idx = findmin(x -> abs(x - time), res.solution.t)[2]
    t = res.solution.t[t_idx] * unit_time
    u = res.solution[t_idx]
    rs = u.x[2] .* unit_length
    vs = u.x[1] .* unit_length / unit_time

    ic = res.simulation.ic
    masses = ic.particles.mass
    R = ic.particles.R
    # stellar_types = ic.particles.stellar_type

    structure_args = propertynames(ic.particles[1].structure)
    stellar_structure_params = Dict(sa => getproperty(ic.particles, sa) for sa in structure_args)

    pop!(stellar_structure_params, :mass)
    stellar_structure_params[:core_masses] = pop!(stellar_structure_params, :core_mass)
    stellar_structure_params[:envelope_masses] = pop!(stellar_structure_params, :envelope_mass)
    stellar_structure_params[:spins] = pop!(stellar_structure_params, :spin)
    stellar_structure_params[:radii] = pop!(stellar_structure_params, :radius)
    stellar_structure_params[:stellar_types] = getfield.(pop!(stellar_structure_params, :stellar_type), :number)
    stellar_structure_params[:core_radii] = pop!(stellar_structure_params, :core_radius)
    stellar_structure_params[:envelope_radii] = pop!(stellar_structure_params, :envelope_radius)
    stellar_structure_params[:L] = pop!(stellar_structure_params, :luminosity)

    sys = multibodysystem(masses, eachcol(rs), eachcol(vs), time = t; stellar_structure_params...)

    return simulate(sys; args...)
end


function total_energy(result::SimulationResult, time)
    masses = result.ode_params.masses
    idx = findmin(x -> abs(x - time), result.solution.t)[2]


    G = get_G_in_system_units(result.simulation.ic)
    return total_energy(
        [result.solution.u[idx].x[2][1:3, i] for i in eachindex(masses)],
        [result.solution.u[idx].x[1][1:3, i] for i in eachindex(masses)],
        masses, G
    )
end
