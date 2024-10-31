using DiffEqBase, StaticArrays
using FunctionWranglers
using ArbNumerics

abstract type MultiBodyInitialConditions end
abstract type AbstractBinary end
abstract type AbstractParticle end

struct ParticleIndex
    i::Int
end

struct BinaryIndex
    i::Int
end

#################################### Multibody system setup #########################################
struct OrbitalElements{aT, PT, eT, ωT, iT, ΩT, νT}
    a::aT # semi-major axis
    P::PT # orbital period
    e::eT # eccentricity
    ω::ωT # argument of periapsis
    i::iT # inclination (with respect to xy-plane)
    Ω::ΩT # longitude of ascending node
    ν::νT # true anomaly
end

OrbitalElements(;a=0.0u"AU", P=0.0u"d", e=0.0, ω=0.0u"°", i=0.0u"°", Ω=0.0u"°", ν=0.0u"°") = OrbitalElements(a, P, e, ω, i, Ω, ν)

struct StellarStructure{tT, mT, RT, ST, LT}
    stellar_type::tT   
    m::mT      # total mass
    R::RT      # total radius
    S::ST      # total spin
    L::LT      # total luminosity
    R_core::RT # core radius
    m_core::mT # core mass
    R_env::RT  # envelope radius
    m_env::mT  # envelope mass
end

struct Particle{siblingType, massType, posType, velType, structType} <: AbstractParticle
    key::ParticleIndex
    parent::BinaryIndex
    sibling::siblingType
    mass::massType
    position::posType
    velocity::velType
    structure::structType
end

struct Binary{siblingType, childType, nchildType, posType, velType, massType, elType} <: AbstractBinary
    key::BinaryIndex
    level::Int
    parent::BinaryIndex
    sibling::siblingType
    children::childType
    nested_children::nchildType
    position::posType
    velocity::velType
    masses::massType
    elements::elType
end

struct HierarchicalMultiple{timeType, bodType, pairType, binType, hierType} <: MultiBodyInitialConditions
    n::Int
    time::timeType
    particles::bodType
    pairs::pairType
    binaries::binType
    levels::SVector{N, Int} where N
    root::Binary
    hierarchy::hierType
end

struct NonHierarchicalSystem{tT, T, pT} <: MultiBodyInitialConditions
    n::Int
    time::tT
    particles::T
    pairs::pT
end
####################################################################################################


####################################### Simulation  setup ##########################################
abstract type AbstractMassBody end

struct MassBody{T} <: AbstractMassBody
    position::SVector{3, T}
    velocity::SVector{3, T}
    spin::SVector{3, T}
    mass::T
end

struct MultiBodySimulation{tType, pType, aType, bType <: AbstractMassBody, potType <: MultiBodyPotential}
    ic::MultiBodyInitialConditions
    bodies::SVector{N, bType} where N
    potential::Dict{Symbol, potType}
    tspan::Tuple{tType, tType}
    params::pType 
    args::aType
    diffeq_args::aType
end

struct SimulationResult{cType, rType <: Quantity{T} where T <: Real, opType, aType}
    solution::DiffEqBase.AbstractTimeseriesSolution
    simulation::MultiBodySimulation
    retcode::cType
    runtime::rType
    ode_params::opType
    args::aType
end

struct MultiBodySolution{tT, rT, vT, ST, SvT, sT, oT, pT}
    ic::MultiBodyInitialConditions # initial conditions
    t::tT # time
    r::rT # positions
    v::vT # velocities
    S::ST # spin
    Sv::SvT # spin velocities
    structure::sT
    ode_system::oT 
    ode_params::pT
end

struct AccelerationFunctions{T}
    fs::T
    N::Int
end
########################################################################################################


################################ Framework for the different potentials ################################
function get_accelerating_function(potential::PureGravitationalPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> pure_gravitational_acceleration!(dvi, dvj, rs, pair, params)
end

function get_accelerating_function(potential::DynamicalTidalPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> dynamical_tidal_drag_force!(dvi, dvj, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::EquilibriumTidalPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> equilibrium_tidal_drag_force!(dvi, dvj, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::StaticEquilibriumTidalPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> equilibrium_tidal_drag_force!(dvi, dvj, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::PN1Potential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN1_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2Potential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN2_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2p5Potential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN2p5_acceleration!(dvi, dvj, rs, vs, pair, params)
end

# function get_accelerating_function(potential::PN3Potential)
#     (dv, u, v, p, t, i) -> PN3_acceleration!(dv, u, v, p, i, n, potential)
# end

# function get_accelerating_function(potential::PN3_5Potential)
#     (dv, u, v, p, t, i) -> PN3_5acceleration!(dv, u, v, p, i, n, potential)
# end

function get_accelerating_function(potential::PNPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN1p5SpinPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN1p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2SpinPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN2_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2p5SpinPotential)
    (dvi, dvj, rs, vs, pair, time, params) -> PN2p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::SpinPrecessionPotential)
    (dvi, dvj, dvs, rs, vs, pair, time, params) -> spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN1SpinPrecessionPotential)
    (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN1p5SpinPrecessionPotential)
    (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1p5_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2SpinPrecessionPotential)
    (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN2_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
end
######################################################################################################

function gather_accelerations_for_potentials(simulation::MultiBodySimulation)
    acceleration_functions = []

    for (potential_name, potential) in simulation.potential
        potential isa SpinPotential && continue
        push!(acceleration_functions, get_accelerating_function(potential))
    end

    AccelerationFunctions(SA[acceleration_functions...], length(acceleration_functions))
end

function gather_spin_accelerations_for_potentials(simulation::MultiBodySimulation)
    acceleration_functions = []

    for (potential_name, potential) in simulation.potential
        !(potential isa SpinPotential) && continue
        push!(acceleration_functions, get_accelerating_function(potential))
    end

    AccelerationFunctions(SA[acceleration_functions...], length(acceleration_functions))
end

###################################### The in-place ODE solver ######################################
function make_initial_conditions(us, vs, ss, spinvel, dtype::Type{<:AbstractFloat})
    n = length(us)
    u0 = MMatrix{6, n, dtype}([reduce(hcat, us); reduce(hcat, ss)])
    v0 = MMatrix{6, n, dtype}([reduce(hcat, vs); reduce(hcat, spinvel)])

    return u0, v0
end

function make_initial_conditions(us, vs, ss, spinvel, dtype::Type{ArbFloat})
    n = length(us)
    u0 = SizedMatrix{6, n, dtype}([reduce(hcat, us); reduce(hcat, ss)])
    v0 = SizedMatrix{6, n, dtype}([reduce(hcat, vs); reduce(hcat, spinvel)])

    return u0, v0
end

function make_initial_conditions(us, vs, dtype::Type{<:AbstractFloat})
    n = length(us)
    u0 = MMatrix{3, n, dtype}(reduce(hcat, us))
    v0 = MMatrix{3, n, dtype}(reduce(hcat, vs))

    return u0, v0
end

function make_initial_conditions(us, vs, dtype::Type{ArbFloat})
    n = length(us)
    u0 = SizedMatrix{3, n, dtype}(reduce(hcat, us))
    v0 = SizedMatrix{3, n, dtype}(reduce(hcat, vs))

    return u0, v0
end


function get_initial_conditions(simulation::MultiBodySimulation, dtype, ::Type{SpinPotential})
    bodies = simulation.bodies

    spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    if any(x -> x isa SpinPotential, values(simulation.potential))
        for pot in values(simulation.potential)
            !(pot isa SpinPotential) && continue

            system = simulation.ic
            for pair in system.pairs
                i, j = pair
                b1 = bodies[i]
                b2 = bodies[j]
                spinvelocity[i] += get_spin_precession_velocity(b1, b2, pot)
                spinvelocity[j] += get_spin_precession_velocity(b2, b1, pot)
            end
        end
    end

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]
    ss = [b.spin for b in bodies]

    u0, v0 = make_initial_conditions(us, vs, ss, spinvelocity, dtype)

    return u0, v0
end

function get_initial_conditions(simulation::MultiBodySimulation, dtype)
    bodies = simulation.bodies

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]

    u0, v0 = make_initial_conditions(us, vs, dtype)

    return u0, v0
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{ArbFloat})
                                          
    u0, v0 = get_initial_conditions(simulation, dtype)
    ai     = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj     = SizedVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0, ai, aj)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{<:AbstractFloat})

    u0, v0 = get_initial_conditions(simulation, dtype)
    ai     = MVector{3, dtype}(zeros(dtype, 3)...)
    aj     = MVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0, ai, aj)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          spin_acc_funcs::AccelerationFunctions, 
                                          dtype::Type{ArbFloat})
                                          
    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    ai     = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj     = SizedVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0, ai, aj)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          spin_acc_funcs::AccelerationFunctions, 
                                          dtype::Type{<:AbstractFloat})

    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    ai     = MVector{3, dtype}(zeros(dtype, 3)...)
    aj     = MVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0, ai, aj)
end


function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          u0, v0, ai, aj)
    pairs = simulation.ic.pairs

    N = acc_funcs.N
    accelerations = FunctionWrangler(acc_funcs.fs)
    output = Vector{Nothing}(undef, N)

    dtype = eltype(u0)
    dtype_0 = zero(dtype)
    function soode_system!(dv, v, u, p, t)
        fill!(dv, dtype_0)
        @inbounds for pair in pairs
            i, j = pair
            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(output, accelerations, ai, aj, u, v, pair, t, p)

            # @inbounds for k = 1:3
            #     dv[k, i] += ai[k]
            #     dv[k, j] += aj[k]
            # end

            dv[1, i] += ai[1]
            dv[1, j] += aj[1]

            dv[2, i] += ai[2]
            dv[2, j] += aj[2]

            dv[3, i] += ai[3]
            dv[3, j] += aj[3]

        end

    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          spin_acc_funcs::AccelerationFunctions,
                                          u0, v0, ai, aj)
    pairs = simulation.ic.pairs

    accelerations = FunctionWrangler(acc_funcs.fs)
    output = Vector{Nothing}(undef, acc_funcs.N)

    spin_accelerations = FunctionWrangler(spin_acc_funcs.fs)
    spin_out = Vector{Nothing}(undef, spin_acc_funcs.N)
 

    dtype = eltype(u0)
    dtype_0 = zero(dtype)
    function soode_system!(dv, v, u, p, t)
        fill!(dv, dtype_0)
        @inbounds for pair in pairs
            i, j = pair
            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(output, accelerations, ai, aj, u, v, pair, t, p)

            # @inbounds for k = 1:3
            #     dv[k, i] += ai[k]
            #     dv[k, j] += aj[k]
            # end

            dv[1, i] += ai[1]
            dv[1, j] += aj[1]

            dv[2, i] += ai[2]
            dv[2, j] += aj[2]

            dv[3, i] += ai[3]
            dv[3, j] += aj[3]

            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(spin_out, spin_accelerations, ai, aj, dv, u, v, pair, t, p)

            dv[4, i] += ai[4]
            dv[4, j] += aj[4]

            dv[5, i] += ai[5]
            dv[5, j] += aj[5]

            dv[6, i] += ai[6]
            dv[6, j] += aj[6]
        end

    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end
######################################################################################################


######################################## The static ODE solver #######################################
function get_initial_conditions_static(simulation::MultiBodySimulation)
    bodies = simulation.bodies;
    n = length(bodies)

    spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    if any(x -> x isa SpinPotential, values(simulation.potential))
        system = simulation.ic
        for pair in system.pairs
            i, j = pair
            b1 = bodies[i]
            b2 = bodies[j]
            dS = spin_precession_velocity(b1, b2)
            spinvelocity[i] += dS
        end
    end

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]
    ss = [b.spin for b in bodies]

    u0 = SMatrix{6, n}([reduce(hcat, us); reduce(hcat, ss)])
    v0 = SMatrix{6, n}([reduce(hcat, vs); reduce(hcat, spinvelocity)])

    (u0, v0, n)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Type{ArbFloat})
    (u0, v0, n) = get_initial_conditions_static(simulation)

    # dtype = eltype(u0)
    ai = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj = SizedVector{3, dtype}(zeros(dtype, 3)...)
    dv = SizedMatrix{6, n, dtype}(undef)
    a_spin = SizedVector{3, dtype}(zeros(dtype, 3)...)

    sodeprob_static(simulation, u0, v0, ai, aj, a_spin, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Type{<:AbstractFloat})
    (u0, v0, n) = get_initial_conditions_static(simulation)
    
    ai = MVector{3, dtype}(zeros(dtype, 3)...)
    aj = MVector{3, dtype}(zeros(dtype, 3)...)
    dv = MMatrix{6, n, dtype}(undef)
    a_spin = MVector{3, dtype}(zeros(dtype, 3)...)

    sodeprob_static(simulation, u0, v0, ai, aj, a_spin, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, u0, v0, ai, aj, a_spin, dv)
    acc_funcs = gather_accelerations_for_potentials(simulation)
    pairs = simulation.ic.pairs
    
    # spin_precession = false
    # d²Sdt²! = if any(x -> x isa SpinPotential, values(simulation.potential))
    #     spin_precession = true
    #     spinpot = first([k for (k, v) in simulation.potential if v isa SpinPotential])
    #     get_accelerating_function(simulation.potential[spinpot])
    # end


    N = acc_funcs.N
    fs = acc_funcs.fs

    dtype = eltype(u0)
    dtype_0 = zero(dtype)
    soode_system = let fs = fs
        function soode_system(v, u, p, t)
            fill!(dv, dtype_0)
            @inbounds for pair in pairs
                i, j = pair
                fill!(ai, dtype_0)
                fill!(aj, dtype_0)
    
                ntuple(i -> fs[i]((ai, aj, u, v, pair, t, p)...), N)

                dv[1:3, i] .= ai
                dv[1:3, j] .= aj
            end

            # if spin_precession
            #     @inbounds for pair in pairs
            #         i, j = pair
            #         fill!(ai, dtype_0)
            #         fill!(aj, dtype_0)

            #         d²Sdt²!(ai, aj, dv, u, v, pair, t, p)
            #         dv[4:6, i] = ai
            #         dv[4:6, j] = aj
            #     end
            # end

            SMatrix(dv)
        end
    end

    SecondOrderODEProblem(soode_system, v0, u0, simulation.tspan, simulation.params)
end

#################################################################################################


######################################## Helper functions #######################################
function Base.getproperty(particles::Dict{Int, Syzygy.Particle}, sym::Symbol)
    if sym in fieldnames(Particle)
        return [getfield(particles[i], sym) for i in keys(sort(particles))]
    elseif sym in fieldnames(StellarStructure)
        return [getfield(particles[i].structure, sym) for i in keys(sort(particles))]
    else
        getfield(particles, sym)
    end

end

function Base.getproperty(binaries::Dict{Int, Syzygy.Binary}, sym::Symbol)
    if sym in fieldnames(Binary)
        return [getfield(binaries[i], sym) for i in keys(sort(binaries))]
    elseif sym in fieldnames(OrbitalElements)
        return [getfield(binaries[i].elements, sym) for i in keys(sort(binaries))]
    else
        getfield(binaries, sym)
    end

end
#################################################################################################