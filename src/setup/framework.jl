using DiffEqBase, StaticArrays

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

struct TripleInitialConditions 
end

struct OrbitalElements{aT, PT, eT, ωT, iT, ΩT, νT}
    a::aT
    P::PT
    e::eT
    ω::ωT
    i::iT
    Ω::ΩT
    ν::νT
end

OrbitalElements(;a=0.0u"AU", P=0.0u"d", e=0.0, ω=0.0u"°", i=0.0u"°", Ω=0.0u"°", ν=0.0u"°") = OrbitalElements(a, P, e, ω, i, Ω, ν)

struct StellarStructure{tT, mT, RT, ST, LT}
    type::tT
    m::mT # mass
    R::RT # radius
    S::ST # spin
    L::LT # luminosity
end

struct PhysicalQuantities{hT, ET}
    h::hT # (Specific) Angular momentum
    E::ET # Total energy
    K::ET # Kinetic energy
    U::ET # Potential energy
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

# function Base.getproperty(obj::Union{Particle, Binary}, sym::Symbol)
#     if sym === :key
#         return getfield(obj, sym).i
#     else
#         return getfield(obj, sym)
#     end
# end

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

struct SimulationResult{cType, rType <: Quantity{T} where T <: Real, aType}
    solution::DiffEqBase.AbstractTimeseriesSolution
    simulation::FewBodySimulation
    retcode::cType
    runtime::rType
    args::aType
end

struct FewBodySolution{tT, rT, vT, eT, sT, qT, oT}
    initial_conditions::MultiBodySystem
    t::tT
    r::rT
    v::vT
    elements::eT
    structure::sT
    quantities::qT
    ode_system::oT
end


################################ Framework for the different potentials ################################
function get_accelerating_function(parameters::PureGravitationalPotential, simulation::FewBodySimulation)
    n = simulation.ic.n
    (dv, u, v, p, t, i) -> pure_gravitational_acceleration!(dv, u, p, i, n, parameters)
end

function get_accelerating_function(parameters::DynamicalTidalPotential, simulation::FewBodySimulation)
    n = simulation.ic.n
    (dv, u, v, p, t, i) -> dynamical_tidal_drag_force!(dv, u, v, p, i, n, parameters)
end

function get_accelerating_function(parameters::EquilibriumTidalPotential, simulation::FewBodySimulation)
    n = simulation.ic.n
    (dv, u, v, p, t, i) -> equilibrium_tidal_drag_force!(dv, u, v, p, i, n, parameters)
end


function gather_accelerations_for_potentials(simulation::FewBodySimulation)
    acceleration_functions = []

    for (potential, parameters) in simulation.system.potential
        push!(acceleration_functions, get_accelerating_function(parameters, simulation))
    end

    tuple(acceleration_functions...)
end

function get_initial_conditions(simulation::FewBodySimulation)
    system = simulation.system
    bodies = system.bodies;
    n = length(bodies)

    u0 = zeros(eltype(bodies[1].position), 3, n)
    v0 = zeros(eltype(bodies[1].velocity), 3, n)

    for i = 1:n
        u0[:, i] = bodies[i].position
        v0[:, i] = bodies[i].velocity
    end

    (u0, v0, n)
end

function DiffEqBase.SecondOrderODEProblem(simulation::FewBodySimulation)

    (u0, v0, n) = get_initial_conditions(simulation)

    acceleration_functions = gather_accelerations_for_potentials(simulation)

    # a = @SVector zeros(3)
    a = MVector(0.0, 0.0, 0.0)
    function soode_system!(dv, v, u, p, t)
        @inbounds for i = 1:n
            # a = MVector(0.0, 0.0, 0.0)
            fill!(a, 0.0)
            for acceleration! in acceleration_functions
                acceleration!(a, u, v, p, t, i);
            end
            dv[:, i] .= a
        end
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end
#################################################################################################


################################ Potential setup functions ################################


function EquilibriumTidalPotential(system::MultiBodySystem, τ, G=𝒢.val)

    m = [u"Msun"(p.mass) for p in values(system.particles)] |> ustrip

    if τ[1] isa Quantity
        τ = upreferred.(τ) |> ustrip
    else
        τ = τ .* minimum(upreferred(bin.elements.P) |> ustrip for bin in values(system.binaries))
    end

    return EquilibriumTidalPotential(m, G=G, τ=τ)
end

###########################################################################################