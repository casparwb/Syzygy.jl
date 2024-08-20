using Printf, Unitful, UnitfulAstro

# abstract type AbstractSyzygyCallback end

# struct CollisionCB{T}         <: AbstractSyzygyCallback 
#     check_every::Int
#     grav_rad_multiple::T
#     CollisionCB(check_every=1, grav_rad_multiple=1000) = new{typeof(grav_rad_multiple)}(check_every)
# end

# struct EscapeCB{T}         <: AbstractSyzygyCallback 
#     max_a_factor::T
#     check_every::Int
#     EscapeCB(max_a_factor=100, check_every=100) = new{typeof(max_a_factor)}(max_a_factor, check_every)
# end

# struct RocheLobeOverflowCB <: AbstractSyzygyCallback 
#     check_every::Int
#     RocheLobeOverflowCB(check_every=1) = new(check_every)
# end

# struct CPUTimeCB           <: AbstractSyzygyCallback 
#     check_every::Int
#     CPUTimeCB(check_every=1) = new(check_every)
# end

# struct CentreOfMassCB      <: AbstractSyzygyCallback 
#     check_every::Int
#     CentreOfMassCB(check_every=1) = new(check_every)
# end

# struct HubbleTimeCB        <: AbstractSyzygyCallback 
#     check_every::Int
#     HubbleTimeCB(check_every=1) = new(check_every)
# end

# struct DemocraticCheckCB   <: AbstractSyzygyCallback 
#     check_every::Int
#     DemocraticCheckCB(check_every=1) = new(check_every)
# end

# struct IonizationCB{T}     <: AbstractSyzygyCallback 
#     check_every::Int
#     max_a_factor::T
#     IonizationCB(max_a_factor=100, check_every=100) = new{typeof(max_a_factor)}(max_a_factor, check_every)
# end

function bodies(system::T where T <: MultiBodyInitialConditions)

    positions  = [zeros(3) for i = 1:system.n]
    velocities = [zeros(3) for i = 1:system.n]
    spins      = [zeros(3) for i = 1:system.n]
    masses     =  zeros(system.n)

    for (key, particle) in system.particles
        positions[key]  .= upreferred.(particle.position) |> ustrip
        velocities[key] .= upreferred.(particle.velocity) |> ustrip
        spins[key]      .= upreferred.(particle.structure.S) |> ustrip
        masses[key]      = upreferred(particle.mass) |> ustrip
    end

    SA[[MassBody(SA[r...], SA[v...], SA[s...], m) for (r, v, s, m) in zip(positions, velocities, spins, masses)]...]
end

function bodies(positions, velocities, spins, masses)
    n = length(masses)

    positions  = [ustrip(upreferred(unit(p[1])), p) for p in positions]
    velocities = [ustrip(upreferred(unit(v[1])), v) for v in velocities]
    spins      = [ustrip(upreferred(unit(v[1])), v) for v in spins]
    masses     = [ustrip(upreferred(unit(m)), m) for m in masses]

    SA[[MassBody(SA[r...], SA[v...], SA[s...], m) for (r, v, s, m) in zip(positions, velocities, spins, masses)]...]
end

function get_potential_dict(potential::MultiBodyPotential)
    pot_dict = Dict{Symbol, MultiBodyPotential}(nameof(typeof(potential)) => potential)
    return pot_dict
end

function get_potential_dict(potential::Vector)
    pot_dict = Dict{Symbol, MultiBodyPotential}()
    for pot in potential
        pot_dict[nameof(typeof(pot))] = pot
    end
    return pot_dict
end


function multibodysimulation(system::T where T <: MultiBodyInitialConditions, tspan,
                             potential,
                             ode_params, args, diffeq_args)
    massbodies = bodies(system)
    pot_dict = get_potential_dict(potential)
    
    MultiBodySimulation(system, massbodies, pot_dict, tspan, 
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
                        :potential => PureGravitationalPotential(),
                        :callbacks => [CollisionCB()], :showprogress => false,
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
[`MultiBodySimulation`](@ref) object.

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
                                             of abstract type ``MultiBodyPotential``, or a ``Vector{MultiBodyPotential}``, in which case
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
    ode_params = setup_params(particles, system.time)
   
    simulation = multibodysimulation(system, args[:tspan], args[:potential], 
                                     ode_params, args, kwargs)

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
function simulation(masses, positions, velocities, spins=nothing; multibodysystem_args=(;), kwargs...)

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

    spins = isnothing(spins) ? [zeros(3)*upreferred(u"kg/m^2/s") for i = 1:n] : spins

    mass_bodies = bodies(positions, velocities, spins, masses)
    pot_dict = get_potential_dict(kwargs[:potential])
    ode_params = setup_params(multiple_system.particles, multiple_system.time)

    
    return MultiBodySimulation(multiple_system, mass_bodies, pot_dict, args[:tspan], ode_params, args, kwargs)
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


function setup_params(particles, time)

    # masses        = typeof(upreferred(1.0u"Msun"))[]
    # luminosities  = typeof(upreferred(1.0u"Lsun"))[]
    # radii         = typeof(upreferred(1.0u"Rsun"))[]
    # core_masses   = typeof(upreferred(1.0u"Msun"))[]
    # core_radii    = typeof(upreferred(1.0u"Rsun"))[]
    # ages          = typeof(upreferred(1.0u"yr"))[]

    masses        = Float64[]
    luminosities  = Float64[]
    radii         = Float64[]
    core_masses   = Float64[]
    core_radii    = Float64[]
    ages          = Float64[]
    stellar_types = Int[]

    particle_keys = keys(particles) |> collect |> sort

    for i in particle_keys
        p = particles[i]
        
        mass         = p.structure.m      |> upreferred |> ustrip 
        luminosity   = p.structure.L      |> upreferred |> ustrip
        radius       = p.structure.R      |> upreferred |> ustrip 
        core_mass    = p.structure.m_core |> upreferred |> ustrip
        core_radius  = p.structure.R_core |> upreferred |> ustrip
        stellar_type = p.structure.stellar_type.index 

        push!(core_masses,   core_mass)
        push!(core_radii,    core_radius)
        push!(ages,          upreferred(time).val)
        push!(masses,        mass)
        push!(luminosities,  luminosity)
        push!(radii,         radius)
        push!(stellar_types, stellar_type)
    end

    masses = MVector(masses...)
    luminosities = MVector(luminosities...)
    radii = MVector(radii...)
    stellar_types = MVector(stellar_types...)
    core_masses = MVector(core_masses...)
    core_radii = MVector(core_radii...)
    ages = MVector(ages...)

    all_params = Dict(:R            => radii, 
                      :M            => masses, 
                      :L            => luminosities,
                      :stellar_type => stellar_types,
                      :core_masses  => core_masses,
                      :core_radii   => core_radii,
                      :ages         => ages)

    ode_params = DefaultSimulationParams(all_params[:R], 
                                         all_params[:M], 
                                         all_params[:L], 
                                         all_params[:stellar_type],
                                         all_params[:core_masses], 
                                         all_params[:core_radii], 
                                         all_params[:ages])

    return ode_params
end
