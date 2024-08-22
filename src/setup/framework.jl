using DiffEqBase, StaticArrays
using Unrolled
using FunctionWranglers

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

struct AccelerationFunctions{T}
    fs::T
    N::Int
end
########################################################################################################


################################ Framework for the different potentials ################################
# function get_accelerating_function(potential::PureGravitationalPotential, n)
#     (dvi, dvj, rs, vs, pair, time, params) -> pure_gravitational_acceleration!(dvi, dvj, rs, pair, params)
# end

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
    (dvi, dvj, rs, vs, pair, time, params) -> PN2_acceleration!(dvi, dvj, rs, vs, pair, params)
end


function get_accelerating_function(potential::PN2_5Potential, n)
    (dvi, dvj, rs, vs, pair, time, params) -> PN1_5_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::PN3Potential, n)
    (dv, u, v, p, t, i) -> PN3_acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PN3_5Potential, n)
    (dv, u, v, p, t, i) -> PN3_5acceleration!(dv, u, v, p, i, n, potential)
end

function get_accelerating_function(potential::PNPotential, n)
    (dvi, dvj, rs, vs, pair, time, params) -> PN_acceleration!(dvi, dvj, rs, vs, pair, params)
end

function get_accelerating_function(potential::SpinPrecessionPotential, n)
    (dvi, dvj, rs, vs, pair, time, params) -> spin_precession!(dvi, dvj, rs, vs, pair, params)
end


######################################################################################################




function gather_accelerations_for_potentials(simulation::MultiBodySimulation)
    acceleration_functions = []

    n = simulation.ic.n
    for (potential, parameters) in simulation.potential
        push!(acceleration_functions, get_accelerating_function(parameters, n))
    end

    # tuple(acceleration_functions...)
    AccelerationFunctions(SA[acceleration_functions...], length(acceleration_functions))
end

###################################### The in-place ODE solver ######################################

function make_initial_conditions(us, vs, ss, spinvel, dtype::Union{Type{Float64}, Type{<:MultipartFloat}})
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


function get_initial_conditions(simulation::MultiBodySimulation, dtype)
    bodies = simulation.bodies
    n = length(bodies)

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

    us = [b.position for b in bodies]
    ss = [b.spin for b in bodies]
    vs = [b.velocity for b in bodies]

    u0, v0 = make_initial_conditions(us, vs, ss, spinvelocity, dtype)

    return u0, v0
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Type{ArbFloat})
    u0, v0 = get_initial_conditions(simulation, dtype)
 
    ai = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj = SizedVector{3, dtype}(zeros(dtype, 3)...)
    a_spin = SizedVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0, ai, aj, a_spin)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          dtype::Union{Type{Float64}, Type{<:MultipartFloat}})

    u0, v0 = get_initial_conditions(simulation, dtype)
    ai = MVector{3, dtype}(zeros(dtype, 3)...)
    aj = MVector{3, dtype}(zeros(dtype, 3)...)
    a_spin = MVector{3, dtype}(zeros(dtype, 3)...)

    SecondOrderODEProblem(simulation, acc_funcs, u0, v0, ai, aj, a_spin)
end

function DiffEqBase.SecondOrderODEProblem(simulation::MultiBodySimulation, 
                                          acc_funcs::AccelerationFunctions, 
                                          u0, v0, ai, aj, a_spin)
    pairs = simulation.ic.pairs

    spin_precession = false
    if :SpinPrecessionPotential in keys(simulation.potential)
        spin_precession = true
        idx = findall(x -> x.parameters == SpinPrecessionPotential(), acc_funcs) |> only
        d²Sdt²! = acc_funcs[idx]
        acc_funcs = tuple([f for f in acc_funcs.fs if !(f.parameters == SpinPrecessionPotential())]...)
        acc_funcs = AccelerationFunctions(SA[acc_funcs...], length(acc_funcs))
    end

    N = acc_funcs.N
    w = FunctionWrangler(acc_funcs.fs)
    o = Vector{Nothing}(undef, N)

    dtype = eltype(u0)
    dtype_0 = zero(dtype)
    function soode_system!(dv, v, u, p, t)

        @inbounds for pair in pairs
            i, j = pair
            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(o, w, ai, aj, u, v, pair, t, p)

            dv[1:3, i] .= ai
            dv[1:3, j] .= aj
        end

        # if spin_precession
        #     @inbounds for i = 1:n
        #         fill!(a_spin, dtype_0)
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

    us = [b.position for b in bodies]
    ss = [b.spin for b in bodies]
    vs = [b.velocity for b in bodies]

    u0 = SMatrix{6, n}([reduce(hcat, us); reduce(hcat, ss)])
    v0 = SMatrix{6, n}([reduce(hcat, vs); reduce(hcat, spinvelocity)])


    (u0, v0, n)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Type{ArbFloat})
    (u0, v0, n) = get_initial_conditions_static(simulation)

    # dtype = eltype(u0)
    ai = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj = SizedVector{3, dtype}(zeros(dtype, 3)...)
    a_spin = SizedVector{3, dtype}(zeros(dtype, 3)...)
    dv = SizedMatrix{6, n, dtype}(undef)


    sodeprob_static(simulation::MultiBodySimulation, u0, v0, ai, aj, a_spin, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, dtype::Union{Type{Float64}, Type{<:MultipartFloat}})
    (u0, v0, n) = get_initial_conditions_static(simulation)
    
    ai = MVector{3, dtype}(zeros(dtype, 3)...)
    aj = MVector{3, dtype}(zeros(dtype, 3)...)
    a_spin = MVector{3, dtype}(zeros(dtype, 3)...)
    dv = MMatrix{6, n, dtype}(undef)

    sodeprob_static(simulation::MultiBodySimulation, u0, v0, ai, aj, a_spin, dv)
end

function sodeprob_static(simulation::MultiBodySimulation, u0, v0, ai, aj, a_spin, dv)
    acc_funcs = gather_accelerations_for_potentials(simulation)
    pairs = simulation.ic.pairs
    
    spin_precession = false
    if :SpinPrecessionPotential in keys(simulation.potential)
        spin_precession = true
        idx = findall(x -> x.parameters == SpinPrecessionPotential(), acc_funcs) |> only
        d²Sdt²! = acc_funcs[idx]
        acc_funcs = tuple([f for f in acc_funcs if !(f.parameters == SpinPrecessionPotential())]...)
    end

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

            if spin_precession
                @inbounds for i = 1:n
                    fill!(a_spin, dtype_0)
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
# @unroll function apply_acc_funcs(state::sT, acc_funcs::aT) where {sT, aT}
#     @unroll for i in 1:length(acc_funcs)
#         # acc_funcs[i](state...)
#         acc_funcs[i](state[1], state[2], state[3], state[4], state[5], state[6], state[7])
#     end
# end

# @unroll function apply_acc_funcs(state::sT, acc_funcs::aT, ::Val{N}()) where {sT, aT, N}
#     @unroll for i in 1:N
#         # acc_funcs[i](state...)
#         acc_funcs[i](state[1], state[2], state[3], state[4], state[5], state[6], state[7])
#     end
# end

@inline function apply_acc_funcs(state::Tuple, acc_funcs::Tuple, ::Val{1})
    acc_funcs[1](state[1], state[2], state[3], state[4], state[5], state[6], state[7])
end

@inline function apply_acc_funcs(state::Tuple, acc_funcs::T, ::Val{1}) where T
    acc_funcs[1](state[1], state[2], state[3], state[4], state[5], state[6], state[7])
end

# function apply_acc_funcs(ai, aj, u, v, pair, t, p, acc_funcs::Tuple, ::Val{1})
#     acc_funcs[1](ai, aj, u, v, pair, t, p)
# end



# @unroll function apply_acc_funcs(ai, aj, u, v, pair, t, p, acc_funcs::Tuple)
#     @unroll for i in 1:length(acc_funcs)
#         # acc_funcs[i](state...)
#         acc_funcs[i](ai, aj, u, v, pair, t, p)
#     end
# end

# @unroll function apply_acc_funcs(ai, aj, u, v, pair, t, p, acc_funcs, ::Val{N}) where N
#     @unroll for i in 1:N
#         # acc_funcs[i](state...)
#         acc_funcs[i](ai, aj, u, v, pair, t, p)
#     end
# end



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