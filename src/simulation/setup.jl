using Printf, Unitful, UnitfulAstro
using ArbNumerics, DoubleFloats

function bodies(system::T where T <: MultiBodyInitialConditions, dtype=Float64)

    positions  = [zeros(dtype, 3) for i = 1:system.n]
    velocities = [zeros(dtype, 3) for i = 1:system.n]
    spins      = [zeros(dtype, 3) for i = 1:system.n]
    masses     =  zeros(dtype, system.n)

    for (key, particle) in system.particles
        positions[key]  .= upreferred.(particle.position) |> ustrip
        velocities[key] .= upreferred.(particle.velocity) |> ustrip
        spins[key]      .= upreferred.(particle.structure.S) |> ustrip
        masses[key]      = upreferred(particle.mass) |> ustrip
    end

    SA[[MassBody(SA[r...], SA[v...], SA[s...], m) for (r, v, s, m) in zip(positions, velocities, spins, masses)]...]
end

function bodies(positions, velocities, spins, masses, dtype=Float64)

    positions  = [dtype.(ustrip(upreferred(unit(p[1])), p)) for p in positions]
    velocities = [dtype.(ustrip(upreferred(unit(v[1])), v)) for v in velocities]
    spins      = [dtype.(ustrip(upreferred(unit(v[1])), v)) for v in spins]
    masses     = [dtype(ustrip(upreferred(unit(m)), m)) for m in masses]

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
    massbodies = bodies(system, args[:dtype])
    pot_dict = get_potential_dict(potential)
    
    MultiBodySimulation(system, massbodies, pot_dict, tspan, 
                        ode_params, args, diffeq_args)
end


function get_kwarg!(kwargs::Dict, arg::Symbol, default_value)
    haskey(kwargs, arg) ? pop!(kwargs, arg) : default_value
end

function parse_arguments!(kwargs::Dict)

    default_args = Dict(
                        :t0        => nothing,  :dt     => 1/10, :t_sim => 1.0,
                        :alg       => DPRKN8(), :saveat => [], :npoints => 0,
                        :save_every => nothing,
                        :maxiters  => Inf,
                        :abstol    => 1.0e-10, :reltol  => 1.0e-10,
                        :potential => [PureGravitationalPotential()],
                        :callbacks => AbstractSyzygyCallback[CollisionCB()], 
                        :showprogress => false,
                        :params => DefaultSimulationParams,
                        :verbose   => false, :max_cpu_time => Inf,
                        :precision => :Float64, :stellar_evolution => false,
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
    simulation(system::MultiBodyInitialConditions; <simulation kwargs>)

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
- `stellar_evolution = false`: whether the stars will be evolved, meaning their structural properties will change. If true, the structural parameters of the stars will be types `MVector`, otherwise `SVector`.
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
function simulation(system::MultiBodyInitialConditions; kwargs...)


    kwargs = Dict{Symbol, Any}(kwargs)
    args = parse_arguments!(kwargs)

    dtype = get_datatype_from_precision(args[:precision])
    args[:dtype] = dtype
    particles = system.particles
    
    if any(x -> x isa SpinPotential, args[:potential])
        check_spins(system.particles.S)
    end

    # if any(x -> typoef(x) in SA[PN1Potential, PN2Potential, PN2p5Potential, PNPotential]) && 
    #    (kwargs[:alg] == DPRKN6() || kwargs[:alg] == DPRKN8() || kwargs[:alg] == DPRKN10() || kwargs[:alg] == DPRKN12())
    #    @warn "The chosen solver is not compatible with a post-Newtonian potential. Recommendations are: Feagin10, Tsit5, Vern6/7/8/9, FineRKN5"
    #    return nothing
    # end

    args[:callbacks] = AbstractSyzygyCallback[args[:callbacks]...]

    periods = if system isa NonHierarchicalSystem
        periods_ = typeof(1.0*unit_time)[]
        for pair in system.pairs
            i, j = pair
            r = particles[i].position - particles[j].position
            v = particles[i].velocity - particles[j].velocity

            d = norm(r)
            v² = norm(v)^2

            M = sum(particles.mass[[i, j]])
            a = semi_major_axis(d, v², M)
            if a < zero(a) # not a bound binary
                push!(periods_, NaN*unit_time)
            else
                push!(periods_, 2π*√(a^3/(GRAVCONST*M)))
            end
        end
        filter(!isnan, periods_)
    else
        [bin.elements.P |> upreferred for bin in values(system.binaries)]
    end

    P_in, P_out = isempty(periods) ? (Inf*unit_time, Inf*unit_time) : extrema(periods)
    
    # setup time step (only used if using symplectic integrator)
    if args[:dt] isa Real
        if isinf(P_in)
            throw(DomainError(args[:dt], "None of the bodies are bound, therefore giving dt as a multiple of the smallest period is not possible. Solution: give dt as a number with a time unit."))
        end
        args[:dt] *= P_in.val # time step is multiple of smallest period
    else
        args[:dt] = ustrip(unit_time, args[:dt])
    end

    # Setup time span
    t0 = args[:t0]
    t0 = isnothing(t0) ? ustrip(unit_time, system.time) : ustrip(unit_time, t0)
    args[:t0] = t0
    t_sim = args[:t_sim]
    # tspan = setup_timespan(t0, t_sim, P_out, dtype)
    t_final = get_final_time(t0, t_sim, P_out, dtype)
    args[:tspan] = (t0, t_final)

    # args[:callbacks] = AbstractSyzygyCallback[args[:callbacks]...]
    # push!(args[:callbacks], FinalTimeCB(t_final))

    ########################################### Set up saving ###########################################
    if haskey(args, :saveat)
        saveat = args[:saveat]
        if saveat isa Number
            saveat = saveat isa Quantity ? ustrip(unit_time, saveat) : saveat
            args[:saveat] = t0:saveat:t_final
        end
    else
        args[:saveat] = []
    end

    # Setup optional saving points
    if !iszero(args[:npoints])
        args[:saveat] = range(t0, t_final, length=args[:npoints])
    end

    if !isnothing(args[:save_every])
        save_every = args[:save_every]
        if isone(save_every)
            kwargs[:save_everystep] = true
        else
            save_every_cb = SavingCB(save_every)
            push!(args[:callbacks], save_every_cb)
            kwargs[:save_everystep] = false
        end
    end
    ####################################################################################################

    # Setup parameters
    if any(x -> any(pot -> x isa pot, [EquilibriumTidalPotential, TimeDependentEquilibriumTidalPotential, DynamicalTidalPotential]), args[:potential])
        if args[:params] == DefaultSimulationParams
            args[:params] = TidalSimulationParams
        end
    end
    ode_params = setup_params(args[:params], system, dtype)
   
    simulation = multibodysimulation(system, args[:tspan], args[:potential], 
                                     ode_params, args, kwargs)

    return simulation
end


function check_spins(spins)
    if isnothing(spins)
        @error "Please give spins if using a spin potential."
    elseif !(typeof(spins) <: AbstractVector{<:AbstractVector{<:Quantity}})
        @error "Please give spins as a Vector of Vectors"
    else
        𝐋, 𝐌, 𝐓 = Unitful.𝐋, Unitful.𝐌, Unitful.𝐓
        correct_dim = 𝐋^3*𝐌/𝐓^2
        if dimension(eltype(first(spins))) != correct_dim
            throw(DimensionMismatch("Given spins do not have correct dimensions. $(dimension(eltype(first(spins)))) was given, but must be $correct_dim"))
        end
    end
end

function get_final_time(t0, t_sim, P_out, datatype=Float64)
    if t_sim isa Quantity
        t_sim *= one(t_sim)
        return (ustrip(unit_time, t_sim) + t0)
    else
        P_out = ustrip(unit_time, P_out)
        return t0 + t_sim*P_out
    end
end

function setup_params(::Type{<:DefaultSimulationParams}, system, datatype=Float64)
    particles = system.particles

    masses        = datatype[]
    radii         = datatype[]
    stellar_types = AbstractStellarType[]
    stellar_type_nums = Int[]

    particle_keys = keys(particles) |> collect |> sort


    for i in particle_keys
        p = particles[i]
        
        mass         = p.structure.m      |> upreferred |> ustrip 
        radius       = p.structure.R      |> upreferred |> ustrip 
        stellar_type = p.structure.stellar_type 

        push!(masses,        mass)
        push!(radii,         radius)
        push!(stellar_types, stellar_type)
        push!(stellar_type_nums, stellar_type.number)
    end

    masses = SVector(masses...)
    radii = SVector(radii...)
    # stellar_types = SVector{length(particles), Union{unique(typeof.(stellar_types))...}}(stellar_types...)
    stellar_types = tuple(stellar_types...)
    stellar_type_nums = SVector(stellar_type_nums...)

    ode_params = DefaultSimulationParams(radii, masses, stellar_types, stellar_type_nums)

    return ode_params
end

function setup_params(::Type{<:TidalSimulationParams}, system, datatype=Float64)
    particles, time = system.particles, system.time

    masses        = datatype[]
    luminosities  = datatype[]
    radii         = datatype[]
    core_masses   = datatype[]
    core_radii    = datatype[]
    ages          = datatype[]
    stellar_types = AbstractStellarType[]
    stellar_type_nums = Int[]


    particle_keys = keys(particles) |> collect |> sort

    luminosity_unit = if unit_mass == u"Msun"
        u"Lsun"
    else
        upreferred(u"Lsun")
    end

    for i in particle_keys
        p = particles[i]
        
        mass         = p.structure.m      |> upreferred |> ustrip 
        luminosity   = p.structure.L      |> luminosity_unit |> ustrip
        radius       = p.structure.R      |> upreferred |> ustrip 
        core_mass    = p.structure.m_core |> upreferred |> ustrip
        core_radius  = p.structure.R_core |> upreferred |> ustrip
        stellar_type = p.structure.stellar_type 

        push!(core_masses,   core_mass)
        push!(core_radii,    core_radius)
        push!(ages,          ustrip(unit_time, time))
        push!(masses,        mass)
        push!(luminosities,  luminosity)
        push!(radii,         radius)
        push!(stellar_types, stellar_type)
        push!(stellar_type_nums, stellar_type.number)

    end

    masses = SVector(masses...)
    luminosities = SVector(luminosities...)
    radii = SVector(radii...)
    stellar_types = SVector(stellar_types...)
    core_masses = SVector(core_masses...)
    core_radii = SVector(core_radii...)
    ages = SVector(ages...)
    stellar_type_nums = SVector(stellar_type_nums...)


    ode_params = TidalSimulationParams(radii, 
                                       masses, 
                                       luminosities, 
                                       stellar_types,
                                       stellar_type_nums,
                                       core_masses, 
                                       core_radii, 
                                       ages)
    # ode_params = DefaultSimulationParams(all_params[:R], 
    #                                      all_params[:M], 
    #                                      all_params[:L], 
    #                                      all_params[:stellar_type],
    #                                      all_params[:core_masses], 
    #                                      all_params[:core_radii], 
    #                                      all_params[:ages])

    return ode_params
end

function get_datatype_from_precision(precision)

    if precision isa Int
        setworkingprecision(ArbFloat, precision)
        return ArbFloat
    elseif precision isa Symbol
        if precision == :Float64
            return Float64
        elseif !(precision ∈ propertynames(DoubleFloats))
            throw(ArgumentError("Given precision is not supported."))
        end

        try 
            return eval(precision)
        catch e
            throw(e)
        end
    else
        @info "Precision should be either an integer for the number of bits, or a symbol."
    end
end