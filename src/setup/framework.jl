using DiffEqBase, StaticArrays
using Unrolled

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

struct PhysicalQuantities{hT, ET}
    h::hT # (Specific) Angular momentum
    E::ET # total energy
    K::ET # kinetic energy
    U::ET # potential energy
end

struct Particle{siblingType, massType, posType, velType, structType, attType} <: AbstractParticle
    key::ParticleIndex
    parent::BinaryIndex
    sibling::siblingType
    mass::massType
    position::posType
    velocity::velType
    structure::structType
    attributes::attType
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

struct MultiBodySystem{timeType, bodType, pairType, binType, hierType, quanType} <: MultiBodyInitialConditions
    n::Int
    time::timeType
    particles::bodType
    pairs::pairType
    binaries::binType
    levels::SVector{N, Int} where N
    root::Binary
    hierarchy::hierType
    quantities::quanType
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
    ic::MultiBodySystem
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
    ic::MultiBodySystem # initial conditions
    t::tT # time
    r::rT # positions
    v::vT # velocities
    S::ST # spin
    Sv::SvT # spin velocities
    structure::sT
    ode_system::oT 
    ode_params::pT
end
########################################################################################################


################################ Framework for the different potentials ################################
function get_accelerating_function(potential::PureGravitationalPotential, n)
    (dvi, dvj, rs, vs, pair, time, params) -> pure_gravitational_acceleration!(dvi, dvj, rs, pair, params)
end

function get_accelerating_function(potential::DynamicalTidalPotential, n)
    (dv, u, v, p, t, i) -> dynamical_tidal_drag_force!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::EquilibriumTidalPotential, n)
    (dv, u, v, p, t, i) -> equilibrium_tidal_drag_force!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::StaticEquilibriumTidalPotential, n)
    (dv, u, v, p, t, i) -> equilibrium_tidal_drag_force!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PN1Potential, n)
    (dvi, dvj, rs, vs, pair, time, params) -> PN1_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2Potential, n)
    (dv, u, v, p, t, i) -> PN2_acceleration!(dv, u, v, p, i, n, potential)
end


function get_accelerating_function(potential::PN2_5Potential, n)
    (dv, u, v, p, t, i) -> PN2_5_acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PN3Potential, n)
    (dv, u, v, p, t, i) -> PN3_acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PN3_5Potential, n)
    (dv, u, v, p, t, i) -> PN3_5acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PNPotential, n)
    (dv, u, v, p, t, i) -> PN1_to_3_5_acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::SpinPrecessionPotential, n)
    (dv, dvs, u, v, p, t, i) -> spin_precession!(dv, dvs, u, v, p, i, n, potential)
end


######################################################################################################




function gather_accelerations_for_potentials(simulation::MultiBodySimulation)
    acceleration_functions = []

    n = simulation.ic.n
    for (potential, parameters) in simulation.potential
        push!(acceleration_functions, get_accelerating_function(parameters, n))
    end

    tuple(acceleration_functions...)
end

###################################### The in-place ODE solver ######################################

function get_initial_conditions(simulation::MultiBodySimulation)
    bodies = simulation.bodies
    n = length(bodies)

    L = 6*n
    u0 = MMatrix{6, n, eltype(bodies[1].position), L}(undef)
    v0 = MMatrix{6, n, eltype(bodies[1].velocity), L}(undef)

    # S0 = MMatrix{3, n, eltype(bodies[1].spin), L}(undef)
    # dS0 = MMatrix{3, n, eltype(bodies[1].spin), L}(undef)

    spinvelocity = []
    if :deSitterPotential in keys(simulation.potential)
        system = simulation.ic
        for i = 1:n
            particle = system.particles[i]
            parent_binary = system.binaries[particle.parent.i]
            dS = deSitter_spin_velocity(particle, parent_binary)
            push!(spinvelocity, ustrip.(upreferred(unit(dS[1])), dS))
        end
    else
        spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    end


    for i = 1:n
        u0[1:3, i] = bodies[i].position
        v0[1:3, i] = bodies[i].velocity

        u0[4:6, i] = bodies[i].spin
        v0[4:6, i] = spinvelocity[i]
    end

    u0, v0, n
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, acc_funcs::Tuple)
    u0, v0, n = get_initial_conditions(simulation)
    pairs = simulation.ic.pairs

    spin_precession = false
    if :SpinPrecessionPotential in keys(simulation.potential)
        spin_precession = true
        idx = findall(x -> x.parameters == SpinPrecessionPotential(), acc_funcs) |> only
        d²Sdt²! = acc_funcs[idx]
        acc_funcs = tuple([f for f in acc_funcs if !(f.parameters == SpinPrecessionPotential())]...)
    end

    fg = acc_funcs[1]

    ai = MVector{3, Float64}(zeros(3)...)
    aj = MVector{3, Float64}(zeros(3)...)
    a_spin = MVector{3, Float64}(zeros(3)...)
    function soode_system!(dv, v, u, p, t)

        @inbounds for pair in pairs
            i, j = pair
            fill!(ai, 0.0)
            fill!(aj, 0.0)
            
            fg(ai, aj, u, v, pair, t, p)

            dv[1:3, i] .= ai
            dv[1:3, j] .= aj
        end

        # if spin_precession
        #     @inbounds for i = 1:n
        #         fill!(a_spin, 0.0)
        #         d²Sdt²!(a_spin, dv, u, v, p, t, i)

        #         dv[4:6, i] = a_spin
        #     end
        # end
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end


######################################################################################################

######################################## The static ODE solver #######################################


function get_initial_conditions_static(simulation::MultiBodySimulation)
    bodies = simulation.bodies;
    n = length(bodies)

    u0 = MMatrix{6, n, eltype(bodies[1].position)}(undef)
    v0 = MMatrix{6, n, eltype(bodies[1].velocity)}(undef)

    spinvelocity = []
    if :deSitterPotential in keys(simulation.potential)
        system = simulation.ic
        for i = 1:n
            particle = system.particles[i]
            parent_binary = system.binaries[particle.parent.i]
            dS = deSitter_spin_velocity(particle, parent_binary)
            push!(spinvelocity, ustrip.(upreferred(unit(dS[1])), dS))
        end
    else
        spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    end

    for i = 1:n
        u0[1:3, i] = bodies[i].position
        v0[1:3, i] = bodies[i].velocity

        u0[4:6, i] = bodies[i].spin
        v0[4:6, i] = spinvelocity[i]
    end

    u0 = SMatrix{6, n}(u0)
    v0 = SMatrix{6, n}(v0)

    (u0, v0, n)
end

function sodeprob_static(simulation::MultiBodySimulation)
    (u0, v0, n) = get_initial_conditions_static(simulation)
    acc_funcs = gather_accelerations_for_potentials(simulation)
    pairs = simulation.ic.pairs
    
    spin_precession = false
    if :SpinPrecessionPotential in keys(simulation.potential)
        spin_precession = true
        idx = findall(x -> x.parameters == SpinPrecessionPotential(), acc_funcs) |> only
        d²Sdt²! = acc_funcs[idx]
        acc_funcs = tuple([f for f in acc_funcs if !(f.parameters == SpinPrecessionPotential())]...)
    end

    fg = acc_funcs[1]

    ai = MVector{3, Float64}(zeros(3)...)
    aj = MVector{3, Float64}(zeros(3)...)
    a_spin = MVector{3, Float64}(zeros(3)...)

    dv = MMatrix{6, n, Float64}(undef)
    soode_system = let acc_funcs = tuple(acc_funcs...) 
        function soode_system(v, u, p, t)
            fill!(dv, zero(dv[1]))
            @inbounds for pair in pairs
                i, j = pair
                fill!(ai, 0.0)
                fill!(aj, 0.0)
    
                fg(ai, aj, u, v, pair, t, p)
    
                dv[1:3, i] .= ai
                dv[1:3, j] .= aj
            end

            if spin_precession
                @inbounds for i = 1:n
                    fill!(a_spin, 0.0)
                    d²Sdt²!(a_spin, dv, u, v, p, t, i)
                    dv[4:6, i] = a_spin
                end
            end

            SMatrix(dv)
        end
    end

    SecondOrderODEProblem(soode_system, v0, u0, simulation.tspan, simulation.params)
end


#################################################################################################

"""
    apply_acc_funcs(state::Tuple, acc_funcs::Tuple)

Apply the acceleration functions to the given state using an unrolled loop.
Unrolling this loop makes it non-allocating when you have more than one acceleration
function.
"""
@unroll function apply_acc_funcs(state::Tuple, acc_funcs::Tuple)
    @unroll for i in 1:length(acc_funcs)
        # acc_funcs[i](state...)
        acc_funcs[i](state[1], state[2], state[3], state[4], state[5], state[6])
    end
end


######################################## The new ODE solver #######################################


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