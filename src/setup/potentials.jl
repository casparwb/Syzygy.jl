# abstract type Body end
using StaticArrays, JLD2, Printf
using LinearAlgebra: dot, norm, ×

include("../physics/tides.jl")

abstract type MultiBodyPotential end
abstract type SimulationParams end


struct DefaultSimulationParams{FloatVecType, IntVecType, stpVecType} <: SimulationParams
    radii::FloatVecType # radii
    masses::FloatVecType # masses
    stellar_types::stpVecType 
    stellar_type_numbers::IntVecType
end

struct TidalSimulationParams{FloatVecType, IntVecType, stpVecType, MutableFloatVecType, IntMatType, T} <: SimulationParams
    radii::FloatVecType # radii
    masses::FloatVecType # masses
    luminosities::FloatVecType # luminosities
    stellar_types::stpVecType 
    stellar_type_numbers::IntVecType
    Z::T 
    envelope_masses::FloatVecType 
    envelope_radii::FloatVecType 
    apsidal_motion_constants::FloatVecType 
    rotational_angular_velocities::MutableFloatVecType
    R³_over_Gm::FloatVecType
    pertuber_mass_ratio_factors::IntMatType
    k_over_T_rad_factors::FloatVecType
    k_over_T_conv::FloatVecType
    k_over_T_conversion_factor::T
end

# struct PNSimulationParams{RType, MType, M2_type, LType, SType, stpType, cMType, cRType, ageType} <: SimulationParams
#     R::RType # radii
#     M::MType # masses
#     M_squared::M2_type # masses squared
#     L::LType # luminosities
#     S::SType # spins
#     stellar_types::stpType 
# end

################################################ Potential structs ################################################

"""
    PureGravitationalPotential()


Newtonian gravitational potential. Corresponds to the acceleration function `Syzygy.pure_gravitational_acceleration!!`.
"""
struct PureGravitationalPotential{T <: Real} <: MultiBodyPotential
    G::T
end

function PureGravitationalPotential(system)
    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time

    G = get_G_in_system_units(system)
    return PureGravitationalPotential(G)
end


struct PN1Potential{T} <: MultiBodyPotential
    G::T
    G²::T
    c⁻²::T
end

function PN1Potential(system)
    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time

    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time
    G = get_G_in_system_units(system)
    c = get_c_in_system_units(system)

    G² = G*G
    c⁻² = 1/c^2

    return PN1Potential(G, G², c⁻²)
end

PN1Potential(G, c) = PN1Potential(G, G*G, 1/c^2)

struct PN2Potential{T} <: MultiBodyPotential
    G::T
    G²::T
    G³::T
    c⁻⁴::T
end

function PN2Potential(system)
    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time

    G = get_G_in_system_units(system)
    c = get_c_in_system_units(system)

    G² = G*G
    G³ = G²*G
    c⁻⁴ = 1/c^4

    return PN2Potential(G, G², G³, c⁻⁴)
end

PN2Potential(G, c) = PN2Potential(G, G^2, G^3, 1/c^4)

struct PN2p5Potential{T} <: MultiBodyPotential
    G::T
    G²::T
    G³::T
    c⁻⁵::T
end

function PN2p5Potential(system)
    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time

    G = get_G_in_system_units(system)
    c = get_c_in_system_units(system)

    G² = G*G
    G³ = G²*G
    c⁻⁵ = 1/c^5

    return PN2p5Potential(G, G², G³, c⁻⁵)
end

PN2p5Potential(G, c) = PN2p5Potential(G, G^2, G^3, 1/c^5)

struct PNPotential{T} <: MultiBodyPotential 
    G::T
    G²::T
    G³::T
    c⁻²::T
    c⁻⁴::T
    c⁻⁵::T
end

function PNPotential(system)
    u_length, u_mass, u_time = system.units.u_length, system.units.u_mass, system.units.u_time

    G = get_G_in_system_units(system)
    c = get_c_in_system_units(system)

    G² = G*G
    G³ = G²*G

    c⁻² = 1/c^2
    c⁻⁴ = 1/c^4
    c⁻⁵ = 1/c^5

    return PNPotential(G, G², G³, c⁻², c⁻⁴, c⁻⁵)
end

PNPotential(G, c) = PNPotential(G, G^2, G^3, 1/c^2, 1/c^4, 1/c^5)

"""
    DynamicalTidalPotential(n_t, polytropic_index)


Set up the dynamical tidal potential for a system as defined by Samsing, Leigh & Trani 2018. 
Corresponds to the acceleration function `Syzygy.dynamical_tidal_acceleration!!`.

# Arguments

- `n_t`: tidal force power index. Can be either 4 or 10.
- `polytropic_index`: vector of polytropic indices of each body in the system
"""
struct DynamicalTidalPotential{nType, fType <: Function} <: MultiBodyPotential
    n_t::Int  # Tidal force power constant
    polytropic_index::nType # Polytropic index of each star
    tidal_factor::fType

    function DynamicalTidalPotential(n_t, polytropic_index)

        if n_t == 4
            f = tidal_factor_n4
        elseif n_t == 10
            f = tidal_factor_n10
        else
            @warn "Tidal factor for n ≠ {4, 10} is not defined."
            f = x -> x
        end
    
        new{typeof(polytropic_index), typeof(f)}(n_t, polytropic_index, f)
    end
end

"""
    TimeDependentEquilibriumTidalPotential()


Equilibrium tidal potential, with the assumption that the stars are evolving over time, which means
the internal structure changes. The envelope structure is calculated throughout the simulation, rather than being fixed
from the start. 

Corresponds to the acceleration function `Syzygy.pure_gravitational_acceleration!!`.
"""
struct TimeDependentEquilibriumTidalPotential <: MultiBodyPotential 
    function TimeDependentEquilibriumTidalPotential()
        if unit_system != "Solar"
            @warn """Default unit system is not set to solar units, which is what the tidal prescription expects. Conversion is currently not supported. Set the units by calling `Syzygy.set_units("Solar")` """
        end

        return new()
    end
end


"""
    EquilibriumTidalPotential(system; Z=0.02, lb_multiplier=1.1, ub_multiplier=1.1,
                                      supplied_apsidal_motion_constants=nothing, 
                                      supplied_rotational_angular_velocities=nothing)


Set up the equilibrium tidal potential corresponding to the equilibrium tidal force prescription 
from Hut 1981. The apsidal motion constant `k` for each stellar body is acquired by interpolating
a grid of MESA models, using data from Claret 2023. If the structural properties of one or more of the 
stellar objects fall outside the grid, the interpolator can extrapolate to the given boundaries, which can
be adapted using the keyword arguments `lb_multiplier` for the lower bounds, and `ub_multiplier` for the upper bounds.

Alternatively, the apsidal motion constants and rotational angular velocities can be supplied using the keyword arguments
`supplied_apsidal_motion_constants` and `supplied_rotational_angular_velocities`.

Corresponds to the acceleration function `Syzygy.equilibrium_tidal_acceleration!!`.

# Arguments
- `system`: an instance of a `HierarchicalMultiple` or `NonHierarchichalSystem` type.

# Keyword arguments

- `Z`: metallicity
- `lb_multiplier`: multiplier for the lower bounds of the apsidal motion constant. 
- `ub_multiplier`: multiplier for the upper bounds of the apsidal motion constant. 
- `supplied_apsidal_motion_constants`
- `supplied_rotational_angular_velocities`
"""
struct EquilibriumTidalPotential <: MultiBodyPotential end

###################################################################################################################


############################################## Acceleration functions #############################################


"""
    pure_gravitational_acceleration!(dv, rs, pair::Tuple{Int, Int}, params::SimulationParams)

Gravitational acceleration on bodies i and j, with `(i, j) = pair`.
"""
function pure_gravitational_acceleration!(dv, rs,
                                         pair::Tuple{Int, Int},
                                         params::SimulationParams,
                                         pot::PureGravitationalPotential)
    
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]

    r̄ = r̄₁ - r̄₂
    r = norm(r̄)
    n̂ = r̄/r

    m₁ = params.masses[i]
    m₂ = params.masses[j]
    Gr⁻² = -pot.G/r^2

    a = Gr⁻²*n̂

    a₁ =  a*m₂
    a₂ = -a*m₁

    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]
    return nothing
end



"""
    dynamical_tidal_acceleration!!(dv, rs, vs, params::SimulationParams, i::Int, n::Int, potential::DynamicalTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Implementing Tidal and Gravitational Wave Energy Losses in Few-body Codes: A Fast and Easy Drag Force Model](https://arxiv.org/abs/1803.08215)
"""
function dynamical_tidal_acceleration!(dv, rs, vs,
                                     pair::Tuple{Int, Int},
                                     params::SimulationParams,
                                     pot::DynamicalTidalPotential)
    
    # by j on i -> j is (p)erturber, i is (t)idal object
    
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄)
    v = norm(v̄)

    ms = params.masses
    Rs = params.radii

    m₁, m₂ = ms[i], ms[j]
    M = m₁ + m₂
    a = semi_major_axis(r, v^2, M)
    e = eccentricity(r̄, v̄, a, M)
    rₚ = a*(1 - e)

    J = potential.tidal_factor(e)

    a₁ = let
        Rₜ = Rs[i]
        ΔE = tidal_ΔE(m₁, Rₜ, m₂, rₚ, potential.polytropic_index[i], UNITLESS_G)

        ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
        ε = drag_force_coefficient(ΔE, J, a, e, M, potential.n_t, UNITLESS_G)

        F₁₂ = @. (-ε*(v/r^potential.n_t)*v̄/v)
        F₁₂ / m₁
    end

    a₂ = let
        Rₜ = Rs[j]
        ΔE = tidal_ΔE(m₂, Rₜ, m₁, rₚ, potential.polytropic_index[j], UNITLESS_G)
        
        ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
        ε = drag_force_coefficient(ΔE, J, a, e, M, potential.n_t, UNITLESS_G)

        F₂₁ = @. (-ε*(v/r^potential.n_t)*(-v̄)/v)
        F₂₁ / m₂
    end

    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]

    return nothing
end



"""
Acceleration function from equilibrium tides using the Hut 1981 prescription.
"""
function equilibrium_tidal_acceleration!(dv, rs, vs,
                                        pair::Tuple{Int, Int},
                                        params::SimulationParams,
                                        pot::TimeDependentEquilibriumTidalPotential) 

    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    ms = params.masses
    Rs = params.radii

    m₁, m₂ = ms[i], ms[j]
    
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂
    
    r = norm(r̄)
    v = norm(v̄)
    
    r² = r^2
    r_hat = r̄/r

    θ_dot = (r̄ × v̄)/r²
    θ_dot_norm = norm(θ_dot)
    θ_hat = v̄/v

    a = semi_major_axis(r, v^2, m₂+m₁)

     # tidal force on 1 by 2
    a₁ = let k = i
        stellar_type = params.stellar_types[k]
        if !(stellar_type isa Star)
            SA[0.0, 0.0, 0.0]
        else
            S̄₁  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

            core_mass   = params.core_masses[k]
            core_radius = params.core_radii[k]
            luminosity  = params.luminosities[k]
            age_Myr     = params.ages[k]#*to_Myr_conversion_factor # change

            R = Rs[k]

            Ω = norm(S̄₁)
            
            logg = pot.logg[k]
            logm = pot.logm[k]
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = UNITLESS_G*m₂/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, age_Myr, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₂, a)

            kτ = R^3/(UNITLESS_G*m₁)*k_T

            @. -3μ*m₂/m₁*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end

    # tidal force on 2 by 1
    a₂ = let k = j
        stellar_type = params.stellar_types[k]

        if !(stellar_type isa Star)
            SA[0.0, 0.0, 0.0]
        else
            S̄₂  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

            core_mass   = params.core_masses[k]
            core_radius = params.core_radii[k]
            luminosity  = params.luminosities[k]
            age_Myr     = params.ages[k]#*to_Myr_conversion_factor # change

            R = Rs[k]

            Ω = norm(S̄₂)
            
            logg = pot.logg[k]
            logm = pot.logm[k]
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = UNITLESS_G*m₁/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, age_Myr, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₁, a)

            kτ = R^3/(UNITLESS_G*m₂)*k_T

            @. -3μ*m₁/m₂*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end
    
    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]
    nothing
end


function equilibrium_tidal_acceleration!(dv, rs, vs,
                                        pair::Tuple{Int, Int},
                                        params::SimulationParams,
                                        pot::EquilibriumTidalPotential) 

    i, j = pair
    stellar_type_1 = params.stellar_type_numbers[i]
    stellar_type_2 = params.stellar_type_numbers[j]

    if (stellar_type_1 > 9) && (stellar_type_2 > 9) # tides are (currently) only for stars
        return nothing
    end

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄)
    v = norm(v̄)
    
    r⁻¹ = 1/r
    r_hat = r̄*r⁻¹

    m₁ = params.masses[i]
    m₂ = params.masses[j]

    θ_dot = (r̄ × v̄)*r⁻¹*r⁻¹
    θ_dot_norm = norm(θ_dot)
    θ_hat = v̄/v

    separation⁻⁵ = 1/semi_major_axis(r, v^2, m₂+m₁m, pot.G)^5

    # tidal force on 1 by 2
    a₁ = let
        if stellar_type_1 > 9
            SA[0.0, 0.0, 0.0]
        else    
            μ = pot.G*m₂*r⁻¹*r⁻¹

            R = params.radii[i]
            Ω = params.rotational_angular_velocities[i]

            q_fac = params.pertuber_mass_ratio_factors[i,j]
            k_over_T = get_k_over_T(m₁, stellar_type_1, 
                                    params.k_over_T_conv[i], 
                                    params.k_over_T_rad_factors[i], 
                                    separation⁻⁵, q_fac)*params.k_over_T_conversion_factor
            k = params.apsidal_motion_constants[i]
            kτ = params.R³_over_Gm[i]*k_over_T#R^3/(UNITLESS_G*m₁)*k_T

            @. -3μ*m₂/m₁*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end

    # # tidal force on 2 by 1
    a₂ = let r_hat = -r_hat, θ_hat = -θ_hat
        if stellar_type_2 > 9 
            SA[0.0, 0.0, 0.0]
        else    
            μ = pot.G*m₁*r⁻¹*r⁻¹

            R = params.radii[j]
            Ω = params.rotational_angular_velocities[j]

            q_fac = params.pertuber_mass_ratio_factors[j,i]
            k_over_T = get_k_over_T(m₁, stellar_type_2, 
                                    params.k_over_T_conv[j],
                                    params.k_over_T_rad_factors[j], 
                                    separation⁻⁵, q_fac)*params.k_over_T_conversion_factor
            k = params.apsidal_motion_constants[j]
            kτ = params.R³_over_Gm[j]*k_over_T

            @. -3μ*m₁/m₂*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end

    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]
    nothing
end

@fastmath function PN1_acceleration!(dv, rs, vs,
                                    pair::Tuple{Int, Int},
                                    params::SimulationParams,
                                    pot::PN1Potential)
                           
    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    v₂ = norm(v̄₂)
    
    m₁ = params.masses[i]
    m₂ = params.masses[j]
        
    v₁² = v₁^2
    v₂² = v₂^2

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r
    n = r̄*r⁻¹

    v₁v₂ = v₂v₁ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    G = pot.G
    # G_r = G*r⁻¹
    # G_r² = G_r*r⁻¹

    m₁m₂ = m₂m₁ = m₁*m₂

    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3

    ai = (5G²_r³*m₁m₂ + 4G²_r³*m₂^2 + G_r²*m₂*(1.5nv₂^2 - v₁² + 4v₁v₂ - 2v₂²) )*n + G_r²*m₂*(4nv₁ - 3nv₂)*v̄

    aj = let n = -n, v̄ = -v̄, nv₁ = -nv₁, nv₂ = -nv₂ 
         (5G²_r³*m₂m₁ + 4G²_r³*m₁^2 + G_r²*m₁*(1.5nv₁^2 - v₂² + 4v₂v₁ - 2v₁²) )*n + G_r²*m₁*(4nv₂ - 3nv₁)*v̄
    end

    ai *= pot.c⁻²
    aj *= pot.c⁻²

    dv[1, i] += ai[1]
    dv[1, j] += aj[1]

    dv[2, i] += ai[2]
    dv[2, j] += aj[2]

    dv[3, i] += ai[3]
    dv[3, j] += aj[3]

    
    nothing
end



@fastmath function PN2_acceleration!(dv, rs, vs,
                                    pair::Tuple{Int, Int},
                                    params::SimulationParams,
                                    pot::PN2Potential)
                           
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)
    v₁² = v₁^2

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    v₂ = norm(v̄₂)
    v₂² = v₂^2

    m₁ = params.masses[i]
    m₂ = params.masses[j]
    
    # i = 1, j = 2
    # add @fastmath?

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂
    r⁻¹ = 1/r
    # r² = r^2
    n = r̄*r⁻¹

    v₁v₂ = v₂v₁ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    v₁v₂² = v₂v₁² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₂m₁ = m₁*m₂
    m₁²m₂ = m₂m₁² = m₁^2*m₂
    m₁m₂² = m₂²m₁ = m₁*m₂^2

    G_r = pot.G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = pot.G²*r⁻¹^3
    G³_r⁴ = pot.G³*r⁻¹^4

    # PN-2 acceleration:
    # expression is split up to avoid allocations that can appear in long expressions

    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57m₁²m₂*0.25 - 69m₁m₂²*0.5 - 9m₂^3) 
    a_num += G_r²*m₂*(-1.875nv₂⁴ + 1.5nv₂²*v₁² - 6nv₂²*v₁v₂ - 2v₁v₂² + 4.5nv₂²*v₂² + 
                        4v₁v₂*v₂² - 2v₂^4)
    a_num += G²_r³*m₁m₂*(19.5nv₁² - 39nv₁*nv₂ + 8.5nv₂² - 3.75v₁² - 2.5v₁v₂ + 1.25v₂²) 
    a_num += G²_r³*m₂^2*(2nv₁² - 4nv₁*nv₂ - 6nv₂² - 8v₁v₂ + 4v₂²) 
    a₂1 = n*a_num

    a_num = G²_r³*m₂^2*(-2nv₁ - 2nv₂) + G²_r³*m₁m₂*(-15.75nv₁ + 13.75nv₂) 
    a_num += G_r²*m₂*(-6nv₁*nv₂² + 4.5nv₂^3 + nv₂*v₁² - 4nv₁*v₁v₂ + 
                        4nv₂*v₁v₂ + 4nv₁*v₂² - 5nv₂*v₂²)
    a₂2 = v̄*a_num
    a₁ = a₂1 + a₂2


    # acceleration for body 2 (j)
    a₂ = let nv₁ = -nv₁, nv₂ = -nv₂
        a_num = G³_r⁴*(-57m₂²m₁*0.25 - 69m₂m₁²*0.5 - 9*m₁^3) 
        a_num += G_r²*m₁*(-1.875nv₁⁴ + 1.5nv₁²*v₂² - 6nv₁²*v₂v₁ - 2v₂v₁² + 4.5nv₁²*v₁² + 
                            4v₂v₁*v₁² - 2v₁^4)
        a_num += G²_r³*m₂m₁*(19.5nv₂² - 39nv₂*nv₁ + 8.5nv₁² - 3.75v₂² - 2.5v₂v₁ + 1.25v₁²) 
        a_num += G²_r³*m₁^2*(2nv₂² - 4nv₂*nv₁ - 6nv₁² - 8v₂v₁ + 4v₁²) 
        a₁1 = (-n)*a_num

        a_num = G²_r³*m₁^2*(-2nv₂ - 2nv₁) + G²_r³*m₂m₁*(-15.75nv₂ + 13.75nv₁) 
        a_num += G_r²*m₁*(-6nv₂*nv₁² + 4.5nv₁^3 + nv₁*v₂² - 4nv₂*v₂v₁ + 
                            4nv₁*v₂v₁ + 4nv₂*v₁² - 5nv₁*v₁²)
        a₂2 = (-v̄)*a_num
        a₁1 + a₂2
    end

    a₁ *= pot.c⁻⁴
    a₂ *= pot.c⁻⁴

    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]
    nothing
end



@fastmath function PN2p5_acceleration!(dv, rs, vs,
                                      pair::Tuple{Int, Int},
                                      params::SimulationParams,
                                      pot::PN2p5Potential)                            
    # i = 1, j = 2
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    m₁ = params.masses[i]
    m₂ = params.masses[j]

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂

    r⁻¹ = 1/r
    v² = v^2

    n = r̄*r⁻¹

    nv = dot(n, v̄)

    m₁m₂ = m₂m₁ = m₁*m₂
    m₁²m₂ = m₂m₁² = m₁^2*m₂
    m₁m₂² = m₂²m₁ = m₁*m₂^2

    G²_r³ = pot.G²*r⁻¹^3
    G³_r⁴ = pot.G³*r⁻¹^4

    #################### PN-2.5 acceleration ####################
    # acceleration for body 1 (i)
    a_num = div_208_15*G³_r⁴*m₁m₂²*nv - 24G³_r⁴*m₁²m₂*0.2nv + 12G²_r³*m₁m₂*0.2v²
    a1 = a_num*n
    a_num = 8G³_r⁴*m₁²m₂*0.2 - 32G³_r⁴*m₁m₂²*0.2 - 4G²_r³*m₁m₂*0.2v²
    a2 = a_num*v̄
    a₁ = a1 + a2


    # acceleration for body 2 (j)
    a_num = div_208_15*G³_r⁴*m₂m₁²*nv - 24G³_r⁴*m₂²m₁*0.2nv + 12G²_r³*m₂m₁*0.2v²
    a1 = a_num*(-n)
    a_num = 8G³_r⁴*m₂²m₁*0.2 - 32G³_r⁴*m₂m₁²*0.2 - 4G²_r³*m₂m₁*0.2v²
    a2 = a_num*(-v̄)
    a₂ = a1 + a2
    ###############################################################


    a₁ *= pot.c⁻⁵
    a₂ *= pot.c⁻⁵

    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]
    nothing
end

@fastmath function PN1_to_2p5_acceleration!(dv, rs, vs,
                                  pair::Tuple{Int, Int},
                                  params::SimulationParams,
                                  pot::PNPotential)                           
    i, j = pair

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

    m₁ = params.masses[i]
    m₂ = params.masses[j]

    v₁ = norm(v̄₁)
    v₂ = norm(v̄₂)

    v₁² = v₁^2
    v₂² = v₂^2

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂

    r⁻¹ = 1/r
    # r² = r^2
    v² = v^2
    n = r̄*r⁻¹

    v₁v₂ = v₂v₁ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)
    nv = dot(n, v̄)

    v₁v₂² = v₂v₁² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₂m₁ = m₁*m₂
    m₁²m₂ = m₂m₁² = m₁^2*m₂
    m₁m₂² = m₂²m₁ = m₁*m₂^2

    G = pot.G
    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = pot.G²*r⁻¹^3
    G³_r⁴ = pot.G³*r⁻¹^4

    ai_PN1, ai_PN2, ai_PN2p5 = let
        
        #################### PN-1 acceleration ##################
        # PN1 = (5G²*m₁*m₂*r⁻¹^3 + 4G²*m₂^2*r⁻¹^3 + G*m₂*r⁻¹^2*(1.5*nv₂^2 - v₁² + 4v₁v₂ - 2v₂²) )*n + G*m₂*r⁻¹^2*(4nv₁ - 3nv₂)*v̄
        PN1 = (5G²_r³*m₁m₂ + 4G²_r³*m₂^2 + G_r²*m₂*(1.5nv₂^2 - v₁² + 4v₁v₂ - 2v₂²) )*n + G_r²*m₂*(4nv₁ - 3nv₂)*v̄
        #########################################################

        #################### PN-2 acceleration ##################
        a_num = G³_r⁴*(-57m₁²m₂*0.25 - 69m₁m₂²*0.5 - 9m₂^3) 
        a_num += G_r²*m₂*(-1.875nv₂⁴ + 1.5nv₂²*v₁² - 6nv₂²*v₁v₂ - 2v₁v₂² + 4.5nv₂²*v₂² + 
                            4v₁v₂*v₂² - 2v₂^4)
        a_num += G²_r³*m₁m₂*(19.5nv₁² - 39nv₁*nv₂ + 8.5nv₂² - 3.75v₁² - 2.5v₁v₂ + 1.25v₂²) 
        a_num += G²_r³*m₂^2*(2nv₁² - 4nv₁*nv₂ - 6nv₂² - 8v₁v₂ + 4v₂²) 
        a₂1 = n*a_num
    
        a_num = G²_r³*m₂^2*(-2nv₁ - 2nv₂) + G²_r³*m₁m₂*(-15.75nv₁ + 13.75nv₂) 
        a_num += G_r²*m₂*(-6nv₁*nv₂² + 4.5nv₂^3 + nv₂*v₁² - 4nv₁*v₁v₂ + 
                            4nv₂*v₁v₂ + 4nv₁*v₂² - 5nv₂*v₂²)
        a₂2 = v̄*a_num
        PN2 = a₂1 + a₂2
        #########################################################


        ################### PN-2.5 acceleration #################
        a_num = div_208_15*G³_r⁴*m₁m₂²*nv - 24G³_r⁴*m₁²m₂*0.2nv + 12G²_r³*m₁m₂*0.2v²
        a1 = a_num*n
        a_num = 8G³_r⁴*m₁²m₂*0.2 - 32G³_r⁴*m₁m₂²*0.2 - 4G²_r³*m₁m₂*0.2v²
        a2 = a_num*v̄
        PN2p5 = a1 + a2
        #########################################################


        PN1, PN2, PN2p5
    end

    
    
    aj_PN1, aj_PN2, aj_PN2p5 = let n = -n, v̄ = -v̄, nv₁ = -nv₁, nv₂ = -nv₂

        #################### PN-1 acceleration ##################
        # PN1 = (5G²*m₂*m₁*r⁻¹^3 + 4G²*m₁^2*r⁻¹^3 + G*m₁*r⁻¹^2*(1.5*nv₁^2 - v₂² + 4v₂v₁ - 2v₁²) )*n + G*m₁*r⁻¹^2*(4nv₂ - 3nv₁)*v̄
        PN1 = (5G²_r³*m₂m₁ + 4G²_r³*m₁^2 + G_r²*m₁*(1.5nv₁^2 - v₂² + 4v₂v₁ - 2v₁²) )*n + G_r²*m₁*(4nv₂ - 3nv₁)*v̄
        #########################################################


        #################### PN-2 acceleration ##################
        a_num = G³_r⁴*(-57m₂²m₁*0.25 - 69m₂m₁²*0.5 - 9m₁^3) 
        a_num += G_r²*m₁*(-1.875nv₁⁴ + 1.5nv₁²*v₂² - 6nv₁²*v₂v₁ - 2v₂v₁² + 4.5nv₁²*v₁² + 
                            4v₂v₁*v₁² - 2v₁^4)
        a_num += G²_r³*m₂m₁*(19.5nv₂² - 39nv₂*nv₁ + 8.5nv₁² - 3.75v₂² - 2.5v₂v₁ + 1.25v₁²) 
        a_num += G²_r³*m₁^2*(2nv₂² - 4nv₂*nv₁ - 6nv₁² - 8v₂v₁ + 4v₁²) 
        a₁1 = n*a_num

        a_num = G²_r³*m₁^2*(-2nv₂ - 2nv₁) + G²_r³*m₂m₁*(-15.75nv₂ + 13.75nv₁) 
        a_num += G_r²*m₁*(-6nv₂*nv₁² + 4.5nv₁^3 + nv₁*v₂² - 4nv₂*v₂v₁ + 
                            4nv₁*v₂v₁ + 4nv₂*v₁² - 5*nv₁*v₁²)
        a₂2 = v̄*a_num
        PN2 = a₁1 + a₂2
        #########################################################


        ################## PN-2p5 acceleration ##################
        a_num = div_208_15*G³_r⁴*m₂m₁²*nv - 24G³_r⁴*m₂²m₁*0.2nv + 12G²_r³*m₂m₁*0.2v²
        a1 = a_num*n
        a_num = 8G³_r⁴*m₂²m₁*0.2 - 32G³_r⁴*m₂m₁²*0.2 - 4G²_r³*m₂m₁*0.2v²
        a2 = a_num*v̄
        PN2p5 = a1 + a2
        #########################################################

        PN1, PN2, PN2p5
    end

    a₁ = ai_PN1*pot.c⁻² + ai_PN2*pot.c⁻⁴ + ai_PN2p5*pot.c⁻⁵
    a₂ = aj_PN1*pot.c⁻² + aj_PN2*pot.c⁻⁴ + aj_PN2p5*pot.c⁻⁵


    dv[1, i] += a₁[1]
    dv[1, j] += a₂[1]

    dv[2, i] += a₁[2]
    dv[2, j] += a₂[2]

    dv[3, i] += a₁[3]
    dv[3, j] += a₂[3]

    nothing
end


###################################################################################################################