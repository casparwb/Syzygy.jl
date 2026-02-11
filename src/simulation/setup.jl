using Printf
using ArbNumerics, DoubleFloats

function bodies(system::T where {T <: MultiBodyInitialConditions}, dtype = Float64)

    positions = [zeros(dtype, 3) for _ in 1:system.n]
    velocities = [zeros(dtype, 3) for _ in 1:system.n]
    masses = zeros(dtype, system.n)

    u_l, u_m, u_t = system.units.u_length, system.units.u_mass, system.units.u_time

    for (key, particle) in system.particles
        positions[key] .= ustrip.(u_l, particle.position)
        velocities[key] .= ustrip.(u_l / u_t, particle.velocity)
        masses[key] = ustrip.(u_m, particle.mass)
    end

    return SA[[MassBody(SA[r...], SA[v...], m) for (r, v, m) in zip(positions, velocities, masses)]...]
end

function bodies(positions, velocities, masses, units, dtype = Float64)

    u_l, u_m, u_t = units.u_length, units.u_mass, units.u_time

    positions = [dtype.(ustrip.(u_l, p)) for p in positions]
    velocities = [dtype.(ustrip.(u_l / u_t, v)) for v in velocities]
    masses = [dtype(ustrip.(u_m, m)) for m in masses]

    return SA[[MassBody(SA[r...], SA[v...], m) for (r, v, m) in zip(positions, velocities, masses)]...]
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


function multibodysimulation(
        system::T where {T <: MultiBodyInitialConditions},
        tspan, potential, ode_params, args, diffeq_args
    )

    massbodies = bodies(system, args[:dtype])
    pot_dict = get_potential_dict(potential)

    return MultiBodySimulation(
        system, massbodies, pot_dict, tspan,
        ode_params, args, diffeq_args
    )
end


function get_kwarg!(kwargs::Dict, arg::Symbol, default_value)
    return haskey(kwargs, arg) ? pop!(kwargs, arg) : default_value
end

function parse_arguments!(kwargs::Dict)

    default_args = Dict(
        :t0 => nothing, :dt => 0.0u"s", :t_sim => 1.0,
        :alg => ODESolvers.DPRKN8, :saveat => [], :npoints => 0,
        :save_every => nothing, :maxiters => Inf,
        :abstol => 1.0e-10, :reltol => 1.0e-10,
        :callbacks => AbstractSyzygyCallback[CollisionCB()],
        :showprogress => false,
        :params => DefaultSimulationParams,
        :verbose => false, :max_cpu_time => Inf,
        :precision => :Float64, :stellar_evolution => false,
        :param_options => Dict(),
        :softening => 0.0,
    )

    args = copy(default_args)
    for (k, _) in kwargs
        if haskey(args, k)
            args[k] = pop!(kwargs, k)
        end
    end

    # args = merge(args, kwargs)
    return args
end

function check_potentials_and_units(system, potentials)

    for pot in potentials
        !hasproperty(pot, :G) && continue
        G = getproperty(pot, :G)
        @assert G ≈ get_G_in_system_units(system) "Given potential is not in the correct unit system."
    end

    return
end


"""
    simulation(system::MultiBodyInitialConditions; <simulation kwargs>)

Setup a simulation with a given system and simulation arguments. Returns a 
`Syzygy.MultiBodySimulation` object.

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
                                             the total acceleration will be the sum of all acceleration functions for each potential.
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
```
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
    args[:potential] = pop!(kwargs, :potential, [PureGravitationalPotential(system, dtype = dtype, softening = args[:softening])])
    args[:dtype] = dtype
    particles = system.particles

    unit_length, unit_mass, unit_time = system.units.u_length, system.units.u_mass, system.units.u_time

    check_potentials_and_units(system, args[:potential])

    args[:callbacks] = AbstractSyzygyCallback[args[:callbacks]...]

    periods = if system isa NonHierarchicalSystem
        periods_ = Quantity[]
        for pair in system.pairs
            i, j = pair
            r = particles[i].position - particles[j].position
            v = particles[i].velocity - particles[j].velocity

            d = norm(r)
            v² = norm(v)^2

            M = sum(particles.mass[[i, j]])
            a = semi_major_axis(d, v², M)
            if a < zero(a) # not a bound binary
                push!(periods_, NaN * default_unit_time)
            else
                push!(periods_, 2π * √(a^3 / (GRAVCONST * M)))
            end
        end
        filter(!isnan, periods_)
    else
        [bin.elements.P for bin in values(system.binaries)]
    end

    P_in, P_out = isempty(periods) ? (Inf * default_unit_time, Inf * default_unit_time) : extrema(periods)

    # setup time step (only used if using symplectic integrator)
    if args[:dt] isa Real
        if isinf(P_in)
            throw(DomainError(args[:dt], "None of the bodies are bound, therefore giving dt as a multiple of the smallest period is not possible. Solution: give dt as a number with a time unit."))
        end
        args[:dt] *= ustrip(unit_time, P_in) # time step is multiple of smallest period
    else
        args[:dt] = ustrip(unit_time, args[:dt])
    end

    # Setup time span
    t0 = args[:t0]
    t0 = isnothing(t0) ? ustrip(unit_time, system.time) : ustrip(unit_time, t0)
    args[:t0] = t0
    t_sim = args[:t_sim]
    t_final = get_final_time(t0, t_sim, P_out, unit_time, dtype)
    args[:tspan] = dtype.((t0, t_final))

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
        args[:saveat] = range(t0, t_final, length = args[:npoints])
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
    ode_params = setup_params(args[:params], system, dtype, options = args[:param_options])

    simulation = multibodysimulation(
        system, args[:tspan], args[:potential],
        ode_params, args, kwargs
    )

    return simulation
end

function get_final_time(t0, t_sim, P_out, unit_time, datatype = Float64)
    if t_sim isa Quantity
        t_sim *= one(t_sim)
        return (ustrip(unit_time, t_sim) + t0)
    else
        P_out = ustrip(unit_time, P_out)
        if isinf(P_out)
            throw(ArgumentError("Bodies are not bound, so t_sim can not be given as a Real, but must be given as a Quantity of time."))
        end
        return t0 + t_sim * P_out
    end
end

function setup_params(::Type{<:DefaultSimulationParams}, system, datatype = Float64; options)
    unit_length, unit_mass = system.units.u_length, system.units.u_mass

    particles = system.particles

    masses = datatype[]
    radii = datatype[]
    stellar_types = StellarType[]
    stellar_type_nums = Int[]

    particle_keys = keys(particles) |> collect |> sort


    for i in particle_keys
        p = particles[i]

        mass = ustrip(unit_mass, p.structure.m)
        radius = ustrip(unit_length, p.structure.R)
        stellar_type = p.structure.stellar_type

        push!(masses, mass)
        push!(radii, radius)
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

function setup_params(::Type{<:TidalSimulationParams}, system, datatype = Float64; options)
    unit_length, unit_mass, unit_time = system.units.u_length, system.units.u_mass, system.units.u_time
    unit_power = unit_mass * unit_length^2 / unit_time^3
    # @warn "Tidal params expects units of R⊙, M⊙, and yr."
    # if unit_system != "Solar"
    #     @warn """Default unit system is not set to solar units, which is what the tidal prescription expects. Conversion is currently not supported. Set the units by calling `Syzygy.set_units("Solar")` """
    # end

    age = system.time
    n_bodies = system.n

    R_envs = datatype[]
    m_envs = datatype[]
    masses = datatype[]
    luminosities = datatype[]
    radii = datatype[]
    stellar_types = StellarType[]
    stellar_type_nums = Int[]

    metallicity = get(options, :Z, 0.0134)

    set_stellar_structure = get(options, :set_stellar_structure, false)


    for i in 1:n_bodies
        particle = system.particles[i]
        envelope_radius, envelope_mass = if set_stellar_structure
            if particle.structure.stellar_type isa Star && particle.mass < 1.25Msun
                envelope_structure(system.particles[i], age, Z)
            else
                0.0Rsun, 0.0Rsun
            end
        else
            particle.structure.R_env, particle.structure.m_env
        end

        push!(R_envs, ustrip(unit_length, envelope_radius))
        push!(m_envs, ustrip(unit_mass, envelope_mass))
    end

    envelope_radii = SA[R_envs...]
    envelope_masses = SA[m_envs...]

    lb_multiplier = get(options, :lb_multiplier, 1.1)
    ub_multiplier = get(options, :lb_multiplier, 1.1)
    supplied_apsidal_motion_constants = get(options, :supplied_apsidal_motion_constants, nothing)
    supplied_rotational_angular_velocities = get(options, :supplied_rotational_angular_velocities, nothing)
    set_spin = get(options, :set_spin, false)

    logk_interpolator = if isnothing(supplied_apsidal_motion_constants)
        order = get(options, :order, (3, 3))
        get_k_interpolator(Z = metallicity, lb_multiplier = lb_multiplier, ub_multiplier = ub_multiplier, order = order)
    else
        nothing
    end


    apsidal_motion_constants = Float64[]
    rotational_angular_velocities = Float64[]
    for i in 1:n_bodies
        particle = system.particles[i]

        mass = ustrip(unit_mass, particle.structure.m)
        luminosity = ustrip(unit_power, particle.structure.L)
        radius = ustrip(unit_length, particle.structure.R)
        stellar_type = particle.structure.stellar_type

        push!(masses, mass)
        push!(luminosities, luminosity)
        push!(radii, radius)
        push!(stellar_types, stellar_type)
        push!(stellar_type_nums, stellar_type.number)

        if !(particle.structure.stellar_type isa Star)
            push!(apsidal_motion_constants, 0.0)
            push!(rotational_angular_velocities, 0.0)
            continue
        else
            m, R = particle.mass, particle.radius

            if isnothing(supplied_apsidal_motion_constants)
                logg = log10(ustrip(u"cm/s^2", (GRAVCONST * m / R^2)))
                logm = log10(ustrip(Msun, m))

                logg = clamp(logg, -0.4617, 4.61961)
                logm = clamp(logm, 0.74565, 34.99814)
                try
                    logk = logk_interpolator(logm, logg)
                    k = 10^logk
                    k = ifelse(isinf(k) || isnan(k), 0.0, k)
                    push!(apsidal_motion_constants, k)
                catch e
                    throw(e)
                end
            end

            if isnothing(supplied_rotational_angular_velocities) && set_spin
                Ω = stellar_rotational_frequency(m, R)
                push!(rotational_angular_velocities, ustrip(1 / unit_time, 2π * Ω))
            end
        end
    end

    if !set_spin
        rotational_angular_velocities = 2π * ustrip.(unit_time^(-1), system.particles.spin)
    end

    apsidal_motion_constants = isnothing(supplied_apsidal_motion_constants) ? SA[apsidal_motion_constants...] : SA[supplied_apsidal_motion_constants...]
    rotational_angular_velocities = if isnothing(supplied_rotational_angular_velocities)
        SA[rotational_angular_velocities...]
    else
        SA[ustrip.(unit_time^-1, supplied_rotational_angular_velocities)...]
    end

    rotational_angular_velocities = if get(options, :evolve_spins, true)
        MVector(rotational_angular_velocities)
    else
        rotational_angular_velocities
    end

    kT_convs = let
        ms = system.particles.mass
        Rs = system.particles.radius
        M_envs = system.particles.envelope_mass
        R_envs = system.particles.envelope_radius
        Ls = system.particles.luminosity

        [Syzygy.k_over_T_convective(m, R, m_env, R_env, L) for (m, R, m_env, R_env, L) in zip(ms, Rs, M_envs, R_envs, Ls)]
    end

    fake_perturber_mass = 1.0 # M⊙
    kT_rad_factors = let
        ms = ustrip.(Msun, system.particles.mass)
        R²s = ustrip.(Rsun^2, system.particles.radius .^ 2)

        q₂s = fake_perturber_mass ./ ms

        tmp = [Syzygy.k_over_T_radiative(m, R², 1.0, fake_perturber_mass) for (m, R²) in zip(ms, R²s)]
        @. tmp / (1 + q₂s)^(5 / 6)
    end

    kT_convs = ustrip.(u"yr^-1", kT_convs)

    kT_convs = map(x -> ifelse(isnan(x) || isinf(x), 0.0, x), kT_convs)
    kT_rad_factors = map(x -> ifelse(isnan(x) || isinf(x), 0.0, x), kT_rad_factors)

    pertuber_mass_ratio_factors = let
        m_perturbers = system.particles.mass
        m_objects = system.particles.mass
        [(1 + m_perturbers[j] / m_objects[i])^(5 / 6) for j in 1:n_bodies, i in 1:n_bodies]
    end

    R³_over_Gms = let
        ms = ustrip.(unit_mass, system.particles.mass)
        R³s = ustrip.(Rsun^3, system.particles.radius .^ 3)

        [R³ / (UNITLESS_G * m) for (R³, m) in zip(R³s, ms)]
    end

    R³_over_Gms = SVector(R³_over_Gms...)
    pertuber_mass_ratio_factors = SMatrix{n_bodies, n_bodies}(pertuber_mass_ratio_factors)
    kT_rad_factors = SVector(kT_rad_factors...)
    kT_convs = SVector(kT_convs...)

    masses = SVector(masses...)
    luminosities = SVector(luminosities...)
    radii = SVector(radii...)
    stellar_types = SVector(stellar_types...)
    stellar_type_nums = SVector(stellar_type_nums...)

    k_over_T_conversion_factor = ustrip(u"yr^-1", cbrt(1.0 * unit_length^2 * unit_mass / unit_time^3 * unit_mass^-1 * unit_length^-2))

    ode_params = TidalSimulationParams(
        radii,
        masses,
        luminosities,
        stellar_types,
        stellar_type_nums,
        metallicity,
        envelope_masses,
        envelope_radii,
        apsidal_motion_constants,
        rotational_angular_velocities,
        R³_over_Gms,
        pertuber_mass_ratio_factors,
        kT_rad_factors,
        kT_convs,
        k_over_T_conversion_factor
    )

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
    elseif precision isa DataType
        return precision
    else
        @info "Precision should be either an integer for the number of bits, or a symbol."
    end
end
