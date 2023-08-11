using Printf, LabelledArrays, Unitful, UnitfulAstro


function bodies(nbody::T where T <: FewBodyInitialConditions, frame="com")

    positions = [zeros(3) for i = 1:nbody.n]
    velocities = [zeros(3) for i = 1:nbody.n]
    masses = zeros(nbody.n)

    for (key, particle) in nbody.particles
        positions[key] .= upreferred.(particle.position) |> ustrip
        velocities[key] .= upreferred.(particle.velocity) |> ustrip
        masses[key] = upreferred(particle.mass) |> ustrip
    end
    # @show norm.(positions)
    SA[[MassBody(SA[r...], SA[v...], m) for (r, v, m) in zip(positions, velocities, masses)]...]
end


function FewBodySystem(bodies, params, potential::FewBodyPotential)#::T) where T <: FewBodyPotential
    # pot = potential(params...)
    pot_dict = Dict{Symbol, FewBodyPotential}(nameof(typeof(potential)) => potential)
    system = FewBodySystem(bodies, pot_dict)
    return system
end

function FewBodySystem(bodies, params, potential::Vector)#::T) where T <: FewBodyPotential
    # pot = potential(params...)
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
                        :callbacks => ["collision", "escape"], :showprogress => false,
                        :verbose => false, :ode_params => Dict(), :max_cpu_time => Inf
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
- `alg = DPRKN6()`: ODE solver to use for simulation. See [DiffEq website](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) for list of all solvers.
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
- `callbacks::Vector = ["collision", "escape"]`: callbacks to use in the integration. Can be used to define stopping conditions or other checks.
                                        Value should be an array containing a mix of strings, referencing a pre-defined callback,
                                         or a custom callback from the `DifferentialEquations.jl` ecosystem. See [`Event Handling and Callback Functions`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) for more.
- `showprogress = true`: whether to display the progress of the simulation.
...

The function also accepts all keyword arguments supported by the `CommonSolve.solve` interface from
`DifferentialEquations.jl`. See ['Common Solver Options'](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options) for more information.


# Example
```jldoctest
julia> triple = multibodysystem([1.0, 1.0, 1.0]u"Msun")
julia> sim = simulation(t_sim=10) # simulate for 10 outer orbits
julia> sim = simulation(t_sim=500u"kyr") # simulate for 500 kyr
julia> sim = simulation(t_sim=1, npoints=1000) # save 1000 datapoints
```
"""
function simulation(system::MultiBodySystem; kwargs...)


    kwargs = Dict{Symbol, Any}(kwargs)

    args = parse_arguments!(kwargs)

    particles = system.particles

    # Setup time step (only used if using symplectic integrator)
    periods = [bin.elements.P |> upreferred for bin in values(system.binaries)]
    Pᵢ, Pₒ = extrema(periods)

    args[:dt] *= Pᵢ.val # time step is multiple of inner period

    # Setup time span
    t0 = args[:t0]
    t0 = isnothing(t0) ? u"s"(system.time).val : u"s"(t0).val
    args[:t0] = t0
    t_sim = args[:t_sim]
    if t_sim isa Quantity
        t_sim *= 1.0
        t_sim = u"s"(t_sim).val + t0
        tspan = (t0, t_sim)
    else
        tspan = (t0, t0 + t_sim*Pₒ.val)
        t_sim = t_sim*Pₒ.val
    end
    
    args[:tspan] = tspan

    if !iszero(args[:npoints])
        args[:saveat] = range(tspan..., length=args[:npoints])
    end

    if any(x -> x == 14, [p.structure.type for p in values(system.particles)])
        if ! ("tidal_disruption" in args[:callbacks])
            if args[:verbose]
                @warn "Tidal disruption callback not included even though system contains black hole."
            end
        end
    end

    # Setup parameters for the system
    if args[:potential] isa Vector && any(x -> x isa EquilibriumTidalPotential, args[:potential])

        core_masses = typeof(upreferred(1.0u"Msun"))[]
        core_radii = typeof(upreferred(1.0u"Rsun"))[]
        ages = typeof(upreferred(1.0u"yr"))[]
        for particle_index in 1:system.n
            p = system.particles[particle_index]

            core_mass = p.structure.m_core |> upreferred
            core_radius = p.structure.R_core |> upreferred

            push!(core_masses, core_mass)
            push!(core_radii, core_radius)
            push!(ages, upreferred(system.time))
        end

        core_masses = SA[core_masses...]
        core_radii = SA[core_radii...]
        ages = SA[ages...]

        args[:ode_params][:core_masses] = core_masses
        args[:ode_params][:core_radii] = core_radii
        args[:ode_params][:ages] = ages
    end

    ode_params = setup_params(system.binaries, particles, args[:ode_params])
   
    simulation = FewBodySimulation(system, args[:tspan], 
                                          args[:potential_params], 
                                          args[:potential], ode_params,
                                          args, kwargs)

end


function setup_params(binaries, particles, extras)
    semi_major_axes = typeof(upreferred(1.0u"m"))[]
    masses = typeof(upreferred(1.0u"Msun"))[]
    luminosities = typeof(upreferred(1.0u"Lsun"))[]
    radii = typeof(upreferred(1.0u"Rsun"))[]
    spins = typeof(upreferred(1.0u"1/yr"))[]
    types = typeof(1.0u"stp")[]

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
        stellar_type = p.structure.type.index

        push!(masses, mass)
        push!(luminosities, luminosity)
        push!(radii, radius)
        push!(spins, spin)
        push!(types, stellar_type)
    end

    semi_major_axes = SA[semi_major_axes...]
    masses = SA[masses...]
    luminosities = SA[luminosities...]
    radii = SA[radii...]
    spins = SA[spins...]
    types = SA[types...]

    all_params = Dict(:a => semi_major_axes, :R => radii, 
                      :M => masses, :S => spins,
                      :L => luminosities,
                      :stellar_type => types)
    merge!(all_params, extras)

    tpes = [typeof(arr) for arr in values(all_params)] |> unique
    nmes = tuple(collect(keys(all_params))...)
    ode_params = @LVector Union{tpes...} nmes
    
    for (k, v) in all_params
        setproperty!(ode_params, k, v)
    end

    # ode_params = LVector(NamedTuple(all_params))

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

