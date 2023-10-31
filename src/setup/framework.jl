using DiffEqBase, StaticArrays
using Unrolled

abstract type FewBodyInitialConditions end
abstract type AbstractBinary end
abstract type AbstractParticle end

struct ParticleIndex
    i::Int
end

struct BinaryIndex
    i::Int
end

pI(i) = ParticleIndex(i)
bI(i) = BinaryIndex(i)

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
    type::tT   # stellar type
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

struct MultiBodySystem{timeType, bodType, binType, hierType, quanType} <: FewBodyInitialConditions
    n::Int
    time::timeType
    particles::bodType
    binaries::binType
    levels::SVector{N, Int} where N
    root::Binary
    hierarchy::hierType
    quantities::quanType
end

####################################################################################################


####################################### Simulation  setup ##########################################

abstract type CelestialBody end

struct MassBody{posType <: Real, velType <: Real, mType <: Real} <: CelestialBody
    position::SVector{3, posType}
    velocity::SVector{3, velType}
    mass::mType
end

struct FewBodySystem{bType <: CelestialBody, pType <: FewBodyPotential}
    bodies::SVector{N, bType} where N
    potential::Dict{Symbol, pType}
end

struct FewBodySimulation{tType, pType, aType}
    ic::MultiBodySystem
    system::FewBodySystem
    tspan::Tuple{tType, tType}
    params::pType 
    args::aType
    diffeq_args::aType
end

struct SimulationResult{cType, rType <: Quantity{T} where T <: Real, opType, aType}
    solution::DiffEqBase.AbstractTimeseriesSolution
    simulation::FewBodySimulation
    retcode::cType
    runtime::rType
    ode_params::opType
    args::aType
end

struct FewBodySolution{tT, rT, vT, eT, sT, qT, oT, pT}
    initial_conditions::MultiBodySystem
    t::tT
    r::rT
    v::vT
    elements::eT
    structure::sT
    quantities::qT
    ode_system::oT
    ode_params::pT
end
##################################################################################################


################################ Framework for the different potentials ################################
function get_accelerating_function(parameters::PureGravitationalPotential, n)
    (dv, u, v, p, t, i) -> pure_gravitational_acceleration!(dv, u, p, i, n, parameters)
end

function get_accelerating_function(parameters::DynamicalTidalPotential, n)
    (dv, u, v, p, t, i) -> dynamical_tidal_drag_force!(dv, u, v, p, i, n, parameters)
end

function get_accelerating_function(parameters::EquilibriumTidalPotential, n)
    (dv, u, v, p, t, i) -> equilibrium_tidal_drag_force!(dv, u, v, p, i, n, parameters)
end

function get_accelerating_function(parameters::StaticEquilibriumTidalPotential, n)
    (dv, u, v, p, t, i) -> equilibrium_tidal_drag_force!(dv, u, v, p, i, n, parameters)
end


function gather_accelerations_for_potentials(simulation::FewBodySimulation)
    acceleration_functions = []

    n = simulation.ic.n
    for (potential, parameters) in simulation.system.potential
        push!(acceleration_functions, get_accelerating_function(parameters, n))
    end

    tuple(acceleration_functions...)
    # SA[acceleration_functions...]
    # acceleration_functions
end

function get_initial_conditions(simulation::FewBodySimulation)
    system = simulation.system
    bodies = system.bodies;
    n = length(bodies)

    L = 3*n
    u0 = MMatrix{3, n, eltype(bodies[1].position), L}(undef)
    v0 = MMatrix{3, n, eltype(bodies[1].velocity), L}(undef)

    for i = 1:n
        u0[:, i] = bodies[i].position
        v0[:, i] = bodies[i].velocity
    end

    u0, v0, n
end

function DiffEqBase.SecondOrderODEProblem(simulation::FewBodySimulation, acc_funcs::Tuple)

    u0, v0, n = get_initial_conditions(simulation)

    a = MVector{3, Float64}(0.0, 0.0, 0.0)
    ids = 1:length(acc_funcs)
    function soode_system!(dv, v, u, p, t)
        @inbounds for i = 1:n
            fill!(a, 0.0)

            apply_acc_funcs((a, u, v, p, t, i), acc_funcs)
            
            dv[:, i] = a
        end
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end

@unroll function apply_acc_funcs(state::Tuple, acc_funcs::Tuple)
    @unroll for i in 1:length(acc_funcs)
        acc_funcs[i](state...)
    end
end


function get_initial_conditions_static(simulation::FewBodySimulation)
    system = simulation.system
    bodies = system.bodies;
    n = length(bodies)

    # u0 = zeros(eltype(bodies[1].position), 3, n)
    # v0 = zeros(eltype(bodies[1].velocity), 3, n)

    u0 = MMatrix{3, n, eltype(bodies[1].position)}(undef)
    v0 = MMatrix{3, n, eltype(bodies[1].velocity)}(undef)

    for i = 1:n
        u0[:, i] = bodies[i].position
        v0[:, i] = bodies[i].velocity
    end

    u0 = SMatrix{3, n}(u0)
    v0 = SMatrix{3, n}(v0)

    (u0, v0, n)
end

function sodeprob_static(simulation::FewBodySimulation)
    (u0, v0, n) = get_initial_conditions_static(simulation)

    acceleration_functions = gather_accelerations_for_potentials(simulation)

    a = MVector{3, Float64}(0.0, 0.0, 0.0)
    dv = MMatrix{3, n, Float64}(undef)
    soode_system = let acceleration_functions = tuple(acceleration_functions...) 
        function soode_system(v, u, p, t)
            fill!(dv, 0.0)
            @inbounds for i = 1:n
                fill!(a, 0.0)

                for acceleration! in acceleration_functions
                    acceleration!(a, u, v, p, t, i);
                end

                dv[:, i] = a
            end

            SMatrix(dv)
        end
    end

    SecondOrderODEProblem(soode_system, v0, u0, simulation.tspan, simulation.params)
end
#################################################################################################


################################### Potential setup functions ###################################



#################################################################################################



######################################## Helper functions #######################################

function Base.getproperty(particles::Dict{Int, Syzygy.Particle}, sym::Symbol)
    if sym in fieldnames(Particle)
        return [getfield(particles[i], sym) for i = 1:length(particles)]
    elseif sym in fieldnames(StellarStructure)
        return [getfield(particles[i].structure, sym) for i = 1:length(particles)]
    else
        getfield(particles, sym)
    end

end

function Base.getproperty(binaries::Dict{Int, Syzygy.Binary}, sym::Symbol)
    if sym in fieldnames(Binary)
        return [getfield(binaries[i], sym) for i = 1:length(binaries)]
    elseif sym in fieldnames(OrbitalElements)
        return [getfield(binaries[i].elements, sym) for i = 1:length(binaries)]
    else
        getfield(binaries, sym)
    end

end

#################################################################################################