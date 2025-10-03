using DiffEqBase, StaticArrays
using FunctionWranglers
using ArbNumerics
using RecursiveArrayTools: ArrayPartition
import PreallocationTools

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

struct StellarStructure{tT, m1T, m2T, m3T, RT, ST, LT}
    stellar_type::tT   
    mass::m1T             # total mass
    radius::RT           # total radius
    spin::ST             # spin
    luminosity::LT       # total luminosity
    core_radius::RT      # core radius
    core_mass::m2T        # core mass
    envelope_radius::RT  # envelope radius
    envelope_mass::m3T    # envelope mass
end

struct Particle{siblingType, posType, velType, structType} <: AbstractParticle
    key::ParticleIndex
    parent::BinaryIndex
    sibling::siblingType
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
########################################################################################################


####################################### Simulation postprocess #########################################
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
########################################################################################################


struct AccelerationFunctions{T}
    fs::T
    N::Int
end


################################ Framework for the different potentials ################################
function get_accelerating_function(potential::PureGravitationalPotential)
    (dv, rs, vs, pair, time, params) -> pure_gravitational_acceleration(dv, rs, pair, params)
end

function get_accelerating_function(potential::DynamicalTidalPotential)
    (dv, rs, vs, pair, time, params) -> dynamical_tidal_acceleration(dv, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::TimeDependentEquilibriumTidalPotential)
    (dv, rs, vs, pair, time, params) -> equilibrium_tidal_acceleration(dv, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::EquilibriumTidalPotential)
    (dv, rs, vs, pair, time, params) -> equilibrium_tidal_acceleration(dv, rs, vs, pair, params, potential)
end

function get_accelerating_function(potential::PN1Potential)
    (dv, rs, vs, pair, time, params) -> PN1_acceleration(dv, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2Potential)
    (dv, rs, vs, pair, time, params) -> PN2_acceleration(dv, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN2p5Potential)
    (dv, rs, vs, pair, time, params) -> PN2p5_acceleration(dv, rs, vs, pair, params)
end

function get_accelerating_function(potential::PNPotential)
    (dv, rs, vs, pair, time, params) -> PN1_to_2p5_acceleration(dv, rs, vs, pair, params)
end

# function get_accelerating_function(potential::PN1p5SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN1p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN2_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2p5SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN2p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN1SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN1p5SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1p5_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN2_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end
######################################################################################################

function gather_accelerations_for_potentials(simulation::MultiBodySimulation)
    acceleration_functions = []

    for (potential_name, potential) in simulation.potential
        push!(acceleration_functions, get_accelerating_function(potential))
    end

    AccelerationFunctions(SA[acceleration_functions...], length(acceleration_functions))
end

###################################### The in-place ODE solver ######################################
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

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{<:AbstractFloat})

    u0, v0 = get_initial_conditions(simulation, dtype)

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          u0, v0)
    pairs = simulation.ic.pairs

    N = acc_funcs.N
    accelerations = FunctionWrangler(acc_funcs.fs)
    output = Vector{Nothing}(undef, N)

    function soode_system!(dv, v, u, p, t)
        fill!(dv, 0)
        @inbounds for pair in pairs
            smap!(output, accelerations, dv, u, v, pair, t, p)
        end
        nothing
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end

######################################################################################################


######################################## The static ODE solver #######################################
function get_initial_conditions_static(simulation::MultiBodySimulation)
    bodies = simulation.bodies;
    n = length(bodies)

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]

    u0 = SMatrix{3, n}(reduce(hcat, us))
    v0 = SMatrix{3, n}(reduce(hcat, vs))

    (u0, v0, n)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Type{ArbFloat})
    (u0, v0, n) = get_initial_conditions_static(simulation)

    dv = SizedMatrix{3, n, dtype}(undef)

    sodeprob_static(simulation, u0, v0, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Type{<:AbstractFloat})
    (u0, v0, n) = get_initial_conditions_static(simulation)
    
    dv = MMatrix{3, n, dtype}(undef)

    sodeprob_static(simulation, u0, v0, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, u0, v0, dv)
    acc_funcs = gather_accelerations_for_potentials(simulation)
    pairs = simulation.ic.pairs

    N = acc_funcs.N
    fs = acc_funcs.fs

    soode_system = let fs = fs
        function soode_system(v, u, p, t)            
            dv = similar(u)
            fill!(dv, 0)

            @inbounds for pair in pairs
                ntuple(i -> fs[i]((dv, u, v, pair, t, p)...), N)
            end

            return SMatrix(dv)
        end
    end

    SecondOrderODEProblem(soode_system, v0, u0, simulation.tspan, simulation.params)
end

#################################################################################################

function DiffEqBase.ODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{ArbFloat})
                                          
    u0, v0 = get_initial_conditions(simulation, dtype)

    ODEProblem(simulation, acc_funcs, u0, v0)
end

function DiffEqBase.ODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{<:AbstractFloat})

    u0, v0 = get_initial_conditions(simulation, dtype)

    ODEProblem(simulation, acc_funcs, u0, v0)
end

function DiffEqBase.ODEProblem(simulation::MultiBodySimulation, 
                               acc_funcs::AccelerationFunctions, 
                               r0, v0)

    pairs = simulation.ic.pairs

    n = size(r0, 2)
    N = acc_funcs.N
    accelerations = FunctionWrangler(acc_funcs.fs)
    output = Vector{Nothing}(undef, N)

    dtype = eltype(r0)
    dtype_0 = zero(dtype)
    function ode_system!(du, u, p, t)
        fill!(du, dtype_0)
        @inbounds for pair in pairs
            
            dv = du.x[1]
            r = u.x[2]
            v = u.x[1]

            smap!(output, accelerations, dv, r, v, pair, t, p)
        end

        du.x[2] .= u.x[1]
        return nothing
    end

    u0 = ArrayPartition(v0, r0)

    ODEProblem(ode_system!, u0, simulation.tspan, simulation.params)
end

######################################## Helper functions #######################################

const alternative_stellar_property_names = Dict{Symbol, Symbol}(:m => :mass, :R => :radius, :L => :luminosity, :m_core => :core_mass,
                                                                :R_core => :core_radius, :m_env => :envelope_mass, :R_env => :envelope_radius,
                                                                :S => :spin)

const alternative_orbital_element_names = Dict{Symbol, Symbol}(:semimajor_axis => :a, :semi_major_axis => :a, :sma => :a,
                                                               :eccentricity => :e, :ecc => :a,
                                                               :period => :P, 
                                                               :aop => :ω, :argument_of_periapsis => :ω, :argument_of_pericenter => :ω,
                                                               :inclination => :i, :mutual_inclination => :i, :i_mut => :i,
                                                               :loan => :Ω, :longitude_of_ascending_node => :Ω,
                                                               :true_anomaly => :ν)

function Base.getproperty(particles::Dict{Int, Syzygy.Particle}, sym::Symbol)
    sym = get(alternative_stellar_property_names, sym, sym)
    if sym in fieldnames(Particle)
        return [getfield(particles[i], sym) for i in keys(sort(particles))]
    elseif sym in fieldnames(StellarStructure)
        return [getfield(particles[i].structure, sym) for i in keys(sort(particles))]
    else
        getfield(particles, sym)
    end
end 


function Base.getproperty(particle::Particle, sym::Symbol)
    sym = get(alternative_stellar_property_names, sym, sym)
    if sym in fieldnames(Particle)
        return getfield(particle, sym)
    elseif sym in fieldnames(StellarStructure)
        return getfield(particle.structure, sym)
    else
        getfield(particle, sym)
    end
end

function Base.getproperty(structure::StellarStructure, sym::Symbol)
    sym = get(alternative_stellar_property_names, sym, sym)
    getfield(structure, sym)
end

function Base.getproperty(binaries::Dict{Int, Syzygy.Binary}, sym::Symbol)
    sym = get(alternative_orbital_element_names, sym, sym)
    if sym in fieldnames(Particle)
        return [getfield(binaries[i], sym) for i in keys(sort(binaries))]
    elseif sym in fieldnames(OrbitalElements)
        return [getfield(binaries[i].elements, sym) for i in keys(sort(binaries))]
    else
        getfield(binaries, sym)
    end
end 


function Base.getproperty(binary::Binary, sym::Symbol)
    sym = get(alternative_orbital_element_names, sym, sym)
    if sym in fieldnames(Particle)
        return getfield(binary, sym)
    elseif sym in fieldnames(OrbitalElements)
        return getfield(binary.elements, sym)
    else
        getfield(binary, sym)
    end
end
#################################################################################################