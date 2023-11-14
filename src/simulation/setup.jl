using Printf, Unitful, UnitfulAstro


function bodies(nbody::T where T <: FewBodyInitialConditions)

    positions = [zeros(3) for i = 1:nbody.n]
    velocities = [zeros(3) for i = 1:nbody.n]
    masses = zeros(nbody.n)

    for (key, particle) in nbody.particles
        positions[key] .= upreferred.(particle.position) |> ustrip
        velocities[key] .= upreferred.(particle.velocity) |> ustrip
        masses[key] = upreferred(particle.mass) |> ustrip
    end

    SA[[MassBody(SA[r...], SA[v...], m) for (r, v, m) in zip(positions, velocities, masses)]...]
end


function FewBodySystem(bodies, params, potential::FewBodyPotential)
    pot_dict = Dict{Symbol, FewBodyPotential}(nameof(typeof(potential)) => potential)
    system = FewBodySystem(bodies, pot_dict)
    return system
end

function FewBodySystem(bodies, params, potential::Vector)
    pot_dict = Dict{Symbol, FewBodyPotential}()
    for pot in potential
        pot_dict[nameof(typeof(pot))] = pot
    end
    system = FewBodySystem(bodies, pot_dict)
    return system
end


function FewBodySimulation(system::T where T <: FewBodyInitialConditions, tspan,
                          potential_params,
                          potential,
                          ode_params, args, diffeq_args)
    massbodies = bodies(system)
    ode_system = FewBodySystem(massbodies, potential_params, potential)
    FewBodySimulation(system, ode_system, tspan, 
                      ode_params, args, diffeq_args)
end


function get_kwarg!(kwargs::Dict, arg::Symbol, default_value)
    haskey(kwargs, arg) ? pop!(kwargs, arg) : default_value
end

function parse_arguments!(kwargs::Dict)

    default_args = Dict(
                        :t0 => nothing, :dt => 1/10, :t_sim => 1.0,
                        :alg => DPRKN8(), :saveat => [], :npoints => 0,
                        :maxiters => Inf,
                        :abstol => 1.0e-10, :reltol => 1.0e-10,
                        :potential => PureGravitationalPotential(), :potential_params => [𝒢.val],
                        :callbacks => ["collision"], :showprogress => false,
                        :verbose => false, :max_cpu_time => Inf
                        )

    args = copy(default_args)
    for (k, v) in kwargs
        if haskey(args, k)
            args[k] = pop!(kwargs, k)
        end
    end

    # args = merge(args, kwargs)
    return args
end


"""
    simulation(system::MultiBodySystem; <simulation kwargs>)

Setup a simulation with a given system and simulation arguments. Returns a 
[`FewBodySimulation`](@ref) object.

...
# Arguments
- `alg = DPRKN8()`: ODE solver to use for simulation. See [DiffEq website](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) for list of all solvers.
- `t_sim = 1.0`: simulation time. Can either be a ``Quantity`` of time, or a ``Number``, in which
           case ``t_sim`` is a multiple of the outermost binary in the system.
- `npoints = 0`: number of datapoints to save. The temporal locations of the snapshots
                 are then ``range(tspan..., length=npoints)``.
- `dt = 1/10`: time step for the ODE solver in multiples of the innermost binary. Only
               relevant if solver is symplectic.
- `potential = PureGravitationalPotential()`: potential to use in simulation. Can either be a single object
                                             of abstract type ``FewBodyPotential``, or a ``Vector{FewBodyPotential}``, in which case
                                             the total acceleration will be the sum of all acceleration functions for each potential. For 
                                             all potentials see [`potentials.jl`](@ref).
- `callbacks::Vector = ["collision"]`: callbacks to use in the integration. Can be used to define stopping conditions or other checks.
                                        Value should be an array containing a mix of strings, referencing a pre-defined callback,
                                         or a custom callback from the `DifferentialEquations.jl` ecosystem. See [`Event Handling and Callback Functions`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) for more.
- `showprogress = false`: whether to display the progress of the simulation.
- `max_cpu_time = Inf`: maximum time (in seconds) allowed for the simulation to run. Will terminate if it runs for longer than `max_cpu_time`.  
...

The function also accepts all keyword arguments supported by the `CommonSolve.solve` interface from
`DifferentialEquations.jl`. See ['Common Solver Options'](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options) for more information.


# Example
```jldoctest
julia> triple = multibodysystem([1.0, 1.0, 1.0]u"Msun")
julia> sim = simulation(triple, t_sim=10) # simulate for 10 outer orbits
julia> sim = simulation(triple, t_sim=500u"kyr") # simulate for 500 kyr
julia> sim = simulation(triple, t_sim=1, npoints=1000) # save 1000 datapoints
julia> sim = simulation(triple, t_sim=1, alg=Syzygy.McAte5(), dt=1.0) # use the McAte5 symplectic integrator with a timestep of 1*innermost period.
```
"""
function simulation(system::MultiBodySystem; kwargs...)


    kwargs = Dict{Symbol, Any}(kwargs)
    args = parse_arguments!(kwargs)

    particles = system.particles

    # Setup time step (only used if using symplectic integrator)
    periods = [bin.elements.P |> upreferred for bin in values(system.binaries)]
    P_in, P_out = extrema(periods)
    args[:dt] *= P_in.val # time step is multiple of inner period

    # Setup time span
    t0 = args[:t0]
    t0 = isnothing(t0) ? ustrip(upreferred(u"s"), system.time) : ustrip(upreferred(u"s"), t0)
    args[:t0] = t0
    t_sim = args[:t_sim]
    tspan = setup_timespan(t0, t_sim, P_out)
    args[:tspan] = tspan

    # Setup optional saving points
    if !iszero(args[:npoints])
        args[:saveat] = range(tspan..., length=args[:npoints])
    end

    # Setup parameters
    ode_params = setup_params(system.time, system.binaries, particles)
   
    simulation = FewBodySimulation(system, args[:tspan], 
                                          args[:potential_params], 
                                          args[:potential], ode_params,
                                          args, kwargs)

    return simulation
end


"""
    simulation(masses, positions, velocities; multiple_system_args=(;), kwargs...)

Setup a simulation with just `masses`, `positions`, and `velocities`, and optionally any other argument
accepted by [`multibodysystem`](@ref). The state vectors should be arrays of length `n`, where `n` is the number
of bodie in the system, in which each element is the state of each component.

See [`simulation`](@ref) for a complete overview of simulation-specific arguments.

# Example
```jldoctest
julia> positions = [[1.5, 0.09, 0.0], [1.5, -0.09, 0.0], [-1.5, 0.00, 0.0]]u"AU"
julia> velocities = [[-22.2, 17.2, 0.0], [22.2, 17.2, 0.0], [0.0, -17.2, 0.0]]u"km/s"
julia> masses = [1.0, 1.0, 2.0]u"Msun"
julia> sim = simulation(masses, positions, velocities, t_sim=500u"kyr") # simulate for 500 kyr
julia> sim = simulation(masses, positions, velocities, multibodysystem_args = (;R = [1.0, 1.0, 5.0]u"Rsun")) # include stellar structure arguments

```
"""
function simulation(masses, positions, velocities; multibodysystem_args=(;), kwargs...)

    kwargs = Dict{Symbol, Any}(kwargs)
    args = parse_arguments!(kwargs)

    n = length(masses)
    multiple_system = multibodysystem(masses; multibodysystem_args...)

    pos_body = [zeros(3) for i ∈ 1:n]
    vel_body = [zeros(3) for i ∈ 1:n]
    mass_body = ustrip(upreferred.(masses))

    for i ∈ 1:n
        pos_body[i] .= upreferred.(positions[i]) |> ustrip
        vel_body[i] .= upreferred.(velocities[i]) |> ustrip
    end

    periods = [bin.elements.P |> upreferred for bin in values(multiple_system.binaries)]
    P_in, P_out = extrema(periods)

    if !(args[:dt] isa Quantity)
        args[:dt] *= P_in.val
    else
        args[:dt] = upreferred(args[:dt])
    end

    # Setup time span
    t0 = args[:t0]
    t0 = isnothing(t0) ? ustrip(upreferred(u"s"), multiple_system.time) : ustrip(upreferred(u"s"), t0)
    args[:t0] = t0
    t_sim = args[:t_sim]

    tspan = setup_timespan(t0, t_sim, P_out)
    args[:tspan] = tspan

    bodies  = SA[[MassBody(SA[r...], SA[v...], m) for (r, v, m) in zip(pos_body, vel_body, mass_body)]...]
    
    ode_params = setup_params(multiple_system.time, multiple_system.binaries, 
                              multiple_system.particles)

    system = FewBodySystem(bodies, args[:potential_params], args[:potential])
    sim = FewBodySimulation(multiple_system, system, args[:tspan], ode_params, args, kwargs)

end

function setup_timespan(t0, t_sim, P_out)
    if t_sim isa Quantity
        t_sim *= 1.0
        t_sim = ustrip(upreferred(u"s"), t_sim) + t0
        tspan = (t0, t_sim)
    else
        tspan = (t0, t0 + t_sim*P_out.val)
        t_sim = t_sim*P_out.val
    end

    return tspan
end


function setup_params(time, binaries, particles)
    semi_major_axes = typeof(upreferred(1.0u"m"))[]
    masses = typeof(upreferred(1.0u"Msun"))[]
    luminosities = typeof(upreferred(1.0u"Lsun"))[]
    radii = typeof(upreferred(1.0u"Rsun"))[]
    spins = typeof(upreferred(1.0u"1/yr"))[]
    types = typeof(1.0u"stp")[]
    core_masses = typeof(upreferred(1.0u"Msun"))[]
    core_radii = typeof(upreferred(1.0u"Rsun"))[]
    ages = typeof(upreferred(1.0u"yr"))[]

    particle_keys = keys(particles) |> collect |> sort
    for i in particle_keys
        parent = binaries[particles[i].parent.i]
        sma = parent.elements.a |> upreferred
        push!(semi_major_axes, sma)
    end

    for i in particle_keys
        p = particles[i]
        mass = p.structure.m |> upreferred 
        luminosity = p.structure.L |> upreferred
        radius = p.structure.R |> upreferred 
        spin = p.structure.S |> upreferred 
        stellar_type = p.structure.type.index * u"stp"
        core_mass = p.structure.m_core |> upreferred
        core_radius = p.structure.R_core |> upreferred

        push!(core_masses, core_mass)
        push!(core_radii, core_radius)
        push!(ages, upreferred(time))
        push!(masses, mass)
        push!(luminosities, luminosity)
        push!(radii, radius)
        push!(spins, spin)
        push!(types, stellar_type)
    end

    semi_major_axes = MVector(semi_major_axes...)
    masses = MVector(masses...)
    luminosities = MVector(luminosities...)
    radii = MVector(radii...)
    spins = MVector(spins...)
    types = MVector(types...)
    core_masses = MVector(core_masses...)
    core_radii = MVector(core_radii...)
    ages = MVector(ages...)

    all_params = Dict(:a => semi_major_axes, :R => radii, 
                      :M => masses, :S => spins,
                      :L => luminosities,
                      :stellar_type => types,
                      :core_masses => core_masses,
                      :core_radii => core_radii,
                      :ages => ages)

    ode_params = DefaultSimulationParams(all_params[:a], all_params[:R], 
                                         all_params[:M], all_params[:L], 
                                         all_params[:S], all_params[:stellar_type],
                                         all_params[:core_masses], all_params[:core_radii], 
                                         all_params[:ages])

    return ode_params
end

function Base.show(io::IO, sim::FewBodySimulation)

    println(io, "\nSimulation setup\n-------------------------------")

    print(io, "\nTime span: ")
    timespan = u"kyr".(sim.tspan .* u"s")
    print(io, timespan)

    println("\n")

    println(io, "Arguments ")
    for (arg, val) in sim.args
        arg == :potential && continue 
        @printf(io, "   %-16s %s", "$arg", "$val")
        # @printf(io, "%7s ", "$val")
        println(io)
    end

    println(io)

    print(io, "DiffEq arguments: \n")
    for (arg, val) in sim.diffeq_args
        @printf(io, "   %-16s %s", "$arg", "$val")
        println(io)
    end

    print(io, "\nPotential(s): \n")
    if sim.args[:potential] isa Vector
        for pot in sim.args[:potential]
            pot = nameof(typeof(pot))
            @printf(io, "   %-16s", "$pot")
            println(io)
        end
    else
        pot = nameof(typeof(sim.args[:potential]))
        @printf(io, "   %-16s", "$pot")
    end

    println(io)
end

function Base.show(io::IO, sim::SimulationResult)

    println(io, "Simulation result:\n====================================\n")
    println(io, "Retcodes: ")
    for (k, v) in sim.retcode
        println(io, "   $k", " ", v)
    end

    println(io, "\nRuntime: $(sim.runtime)\n")

    println(io, "Number of datapoints: $(length(sim.solution.t))\n")
    println(io, "ODE Retcode: $(sim.solution.retcode)\n")
    println(io, "$(sim.solution.destats)")
end

