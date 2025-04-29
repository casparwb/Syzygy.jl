# abstract type Body end
using StaticArrays, JLD2, Printf
using LinearAlgebra: dot, norm, ×

include("../physics/tides.jl")

abstract type MultiBodyPotential end
abstract type SpinPotential <: MultiBodyPotential end
abstract type SimulationParams end

# struct DefaultSimulationParams{FloatVecType, stpVecType} <: SimulationParams
#     R::FloatVecType # radii
#     M::FloatVecType # masses
#     L::FloatVecType # luminosities
#     stellar_types::stpVecType 
#     M_cores::FloatVecType # core masses
#     R_cores::FloatVecType # core radii
#     ages::FloatVecType 
# end

struct DefaultSimulationParams{FloatVecType, IntVecType, stpVecType} <: SimulationParams
    radii::FloatVecType # radii
    masses::FloatVecType # masses
    stellar_types::stpVecType 
    stellar_type_numbers::IntVecType
end

struct TidalSimulationParams{FloatVecType, stpVecType} <: SimulationParams
    radii::FloatVecType # radii
    masses::FloatVecType # masses
    luminosities::FloatVecType # luminosities
    stellar_types::stpVecType 
    core_masses::FloatVecType # core masses
    core_radii::FloatVecType # core radii
    ages::FloatVecType 
end

# struct PNSimulationParams{RType, MType, M2_type, LType, SType, stpType, cMType, cRType, ageType} <: SimulationParams
#     R::RType # radii
#     M::MType # masses
#     M_squared::M2_type # masses squared
#     L::LType # luminosities
#     S::SType # spins
#     stellar_types::stpType 
#     M_cores::cMType # core masses
#     R_cores::cRType # core radii
#     ages::ageType 
# end

################################################ Potential structs ################################################

"""
    PureGravitationalPotential()


Newtonian gravitational potential. Corresponds to the acceleration function [`pure_gravitational_acceleration!`](@ref)
"""
struct PureGravitationalPotential <: MultiBodyPotential end


"""
    DynamicalTidalPotential(tidal_force_power_index, polytropic_index)


Set up the dynamical tidal potential for a system as defined by Samsing, Leigh & Trani 2018. 
Corresponds to the acceleration function [`dynamical_tidal_acceleration!`](@ref)

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

Corresponds to the acceleration function [`pure_gravitational_acceleration!`](@ref)
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

Corresponds to the acceleration function [`equilibrium_tidal_acceleration!`](@ref)

# Arguments
- `system`: an instance of a `HierarchicalMultiple` or `NonHierarchichalSystem` type.

# Keyword arguments

- `Z`: metallicity
- `lb_multiplier`: multiplier for the lower bounds of the apsidal motion constant. 
- `ub_multiplier`: multiplier for the upper bounds of the apsidal motion constant. 
- `supplied_apsidal_motion_constants`
- `supplied_rotational_angular_velocities`
"""
struct EquilibriumTidalPotential{T} <: MultiBodyPotential
    Z::Float64 
    M_env::T 
    R_env::T 
    apsidal_motion_constant::T 
    rotational_angular_velocity::T

    function EquilibriumTidalPotential(system; Z=0.02, lb_multiplier=1.1, ub_multiplier=1.1,
                                               supplied_apsidal_motion_constants=nothing, 
                                               supplied_rotational_angular_velocities=nothing,
                                               set_spin=false)
    
        if unit_system != "Solar"
            @warn """Default unit system is not set to solar units, which is what the tidal prescription expects. Conversion is currently not supported. Set the units by calling `Syzygy.set_units("Solar")` """
        end

        age = system.time
        n_bodies = system.n
        R_envs = Float64[]
        m_envs = Float64[]
    
        for i = 1:n_bodies
            particle = system.particles[i]
            envelope_radius, envelope_mass = if particle.structure.stellar_type isa Star && particle.mass < 1.25u"Msun"
                                                envelope_structure(system.particles[i], age, Z)
                                             else
                                                0.0u"Rsun", 0.0u"Rsun"
                                             end
    
            push!(R_envs, ustrip(envelope_radius))
            push!(m_envs, ustrip(envelope_mass))
        end
        
        R_envs = SA[R_envs...]
        m_envs = SA[m_envs...]
    
        logk_interpolator = get_k_interpolator(Z=Z, lb_multiplier=lb_multiplier, ub_multiplier=ub_multiplier)
        apsidal_motion_constants = Float64[]
        rotational_angular_velocities = Float64[]
        for i = 1:n_bodies
            particle = system.particles[i]
    
            if !(particle.structure.stellar_type isa Star)
                push!(apsidal_motion_constants, 0.0)
                push!(rotational_angular_velocities, 0.0)
                continue
            else
                m, R = particle.mass, particle.radius
    
                if isnothing(supplied_apsidal_motion_constants)
                    logg = log10(ustrip(u"cm/s^2", (GRAVCONST*m/R^2))) 
                    logm = log10(ustrip(u"Msun", m))
                    try
                        logk = logk_interpolator(logm, logg)
                        push!(apsidal_motion_constants, 10^logk)
                    catch e
                        throw(e)
                    end
                end
    
                if isnothing(supplied_rotational_angular_velocities) && set_spin
                    Ω = stellar_rotational_frequency(m, R)
                    push!(rotational_angular_velocities, ustrip(unit(1/unit_time), 2π*Ω))
                end
    
            end
        end

        if !set_spin
            rotational_angular_velocities = 2π*ustrip.(unit_time^(-1), system.particles.spin)
        end
    
        apsidal_motion_constants      = isnothing(supplied_apsidal_motion_constants)      ? SA[apsidal_motion_constants...]      : SA[supplied_apsidal_motion_constants...]
        rotational_angular_velocities = if isnothing(supplied_rotational_angular_velocities) 
            SA[rotational_angular_velocities...] 
            else
                SA[ustrip.(unit_time^-1, supplied_rotational_angular_velocities)...]
            end
    
        new{typeof(m_envs)}(Z, m_envs, R_envs, apsidal_motion_constants, rotational_angular_velocities)
    end
end



struct PN1Potential            <: MultiBodyPotential end

struct PN2Potential            <: MultiBodyPotential end

struct PN2p5Potential          <: MultiBodyPotential end

struct PN3Potential            <: MultiBodyPotential end

struct PN3_5Potential          <: MultiBodyPotential end

struct PNPotential             <: MultiBodyPotential end

###################################################################################################################


############################################## Acceleration functions #############################################


"""
    pure_gravitational_acceleration!(dvi, dvj, rs, pair::Tuple{Int, Int}, params::SimulationParams)

Gravitational acceleration on bodies i and j, with `(i, j) = pair`.
"""
function pure_gravitational_acceleration!(dvi, dvj, rs,
                                          pair::Tuple{Int, Int},
                                          params::SimulationParams)
    
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]

    r̄ = r̄₁ - r̄₂
    r = norm(r̄)
    n̂ = r̄/r

    m₁ = params.masses[i]
    m₂ = params.masses[j]
    G_r² = -UNITLESS_G/r^2

    a = G_r²*n̂

    a₁ =  a*m₂
    a₂ = -a*m₁

    dvi .+= a₁
    dvj .+= a₂
    nothing
end



"""
    dynamical_tidal_acceleration!(dv, rs, vs, params::SimulationParams, i::Int, n::Int, potential::DynamicalTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Implementing Tidal and Gravitational Wave Energy Losses in Few-body Codes: A Fast and Easy Drag Force Model](https://arxiv.org/abs/1803.08215)
"""
function dynamical_tidal_acceleration!(dvi, dvj, rs, vs,
                                     pair::Tuple{Int, Int},
                                     params::SimulationParams,
                                     potential::DynamicalTidalPotential)
    
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

    dvi .+= a₁
    dvj .+= a₂
    nothing
end



"""
Acceleration function from equilibrium tides using the Hut 1981 prescription.
"""
function equilibrium_tidal_acceleration!(dvi, dvj, rs, vs,
                                       pair::Tuple{Int, Int},
                                       params::SimulationParams,
                                       potential::TimeDependentEquilibriumTidalPotential) 

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
    θ_hat = θ_dot/θ_dot_norm

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
            age_Myr     = params.ages[k]*to_Myr_conversion_factor

            R = Rs[k]

            Ω = norm(S̄₁)
            
            logg = potential.logg[k]
            logm = potential.logm[k]
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
            age_Myr     = params.ages[k]*to_Myr_conversion_factor

            R = Rs[k]

            Ω = norm(S̄₂)
            
            logg = potential.logg[k]
            logm = potential.logm[k]
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = UNITLESS_G*m₁/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, age_Myr, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₁, a)

            kτ = R^3/(UNITLESS_G*m₂)*k_T

            @. -3μ*m₁/m₂*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end
    
    dvi .+= a₁
    dvj .+= a₂
    nothing
end


function equilibrium_tidal_acceleration!(dvi, dvj, rs, vs,
                                       pair::Tuple{Int, Int},
                                       params::SimulationParams,
                                       potential::EquilibriumTidalPotential) 

    i, j = pair
    stellar_type_1 = params.stellar_types[i]
    stellar_type_2 = params.stellar_types[j]

    if !(stellar_type_1 isa Star) && !(stellar_type_2 isa Star) # tides are (currently) only for stars
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
    θ_hat = θ_dot/θ_dot_norm

    a = semi_major_axis(r, v^2, m₂+m₁)

    # tidal force on 1 by 2
    a₁ = let
        if !(stellar_type_1 isa Star)
                SA[0.0, 0.0, 0.0]
            else    
                μ = UNITLESS_G*m₂*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[i]
                envelope_radius = potential.R_env[i]
                luminosity = params.luminosities[i]
                R = params.radii[i]
                Ω = potential.rotational_angular_velocity[i]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, envelope_mass, envelope_radius,
                                                                   stellar_type_1, luminosity, 
                                                                   m₂, a, potential.Z)
                k = potential.apsidal_motion_constant[i]
                kτ = R^3/(UNITLESS_G*m₁)*k_T
    
                @. -3μ*m₂/m₁*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    # # tidal force on 2 by 1
    a₂ = let 
        if !(stellar_type_2 isa Star)
                SA[0.0, 0.0, 0.0]
            else    
                μ = UNITLESS_G*m₁*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[j]
                envelope_radius = potential.R_env[j]
                luminosity  = params.luminosities[j]
                R = params.radii[j]
                Ω = potential.rotational_angular_velocity[j]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, envelope_mass, envelope_radius,
                                                                   stellar_type_2, luminosity, 
                                                                   m₁, a, potential.Z)
                k = potential.apsidal_motion_constant[j]
                kτ = R^3/(UNITLESS_G*m₂)*k_T
    
                @. -3μ*m₁/m₂*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    dvi .+= a₁
    dvj .+= a₂
    nothing
end

function PN1_acceleration!(dvi, dvj, rs, vs,
                           pair::Tuple{Int, Int},
                           params::SimulationParams)
                           
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

    G = UNITLESS_G
    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹

    m₁m₂ = m₂m₁ = m₁*m₂

    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3

    ai = (5G²_r³*m₁m₂ + 4G²_r³*m₂^2 + G_r²*m₂*(3/2*nv₂^2 - v₁² + 4v₁v₂ - 2v₂²) )*n + G_r²*m₂*(4nv₁ - 3nv₂)*v̄

    aj = let n = -n, v̄ = -v̄, nv₁ = -nv₁, nv₂ = -nv₂ 
         (5G²_r³*m₂m₁ + 4G²_r³*m₁^2 + G_r²*m₁*(3/2*nv₁^2 - v₂² + 4v₂v₁ - 2v₁²) )*n + G_r²*m₁*(4nv₂ - 3nv₁)*v̄
    end
    
    dvi .+= ai*c⁻²
    dvj .+= aj*c⁻²
    nothing
end



function PN2_acceleration!(dvi, dvj, rs, vs,
                           pair::Tuple{Int, Int},
                           params::SimulationParams)
                           
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

    G_r = UNITLESS_G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    # PN-2 acceleration:
    # expression is split up to avoid allocations that can appear in long expressions

    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
    a_num += G_r²*m₂*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
                        4*v₁v₂*v₂² - 2v₂^4)
    a_num += G²_r³*m₁m₂*(39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂² - 15/4*v₁² - 5/2*v₁v₂ + 5/4*v₂²) 
    a_num += G²_r³*m₂^2*(2*nv₁² - 4*nv₁*nv₂ - 6*nv₂² - 8*v₁v₂ + 4v₂²) 
    a₂1 = n*a_num

    a_num = G²_r³*m₂^2*(-2*nv₁ - 2*nv₂) + G²_r³*m₁m₂*(-63/4*nv₁ + 55/4*nv₂) 
    a_num += G_r²*m₂*(-6*nv₁*nv₂² + 9/2*nv₂^3 + nv₂*v₁² - 4*nv₁*v₁v₂ + 
                        4*nv₂*v₁v₂ + 4*nv₁*v₂² - 5*nv₂*v₂²)
    a₂2 = v̄*a_num
    ai = a₂1 + a₂2


    # acceleration for body 2 (j)
    aj = let nv₁ = -nv₁, nv₂ = -nv₂
        a_num = G³_r⁴*(-57*m₂²m₁/4 - 69*m₂m₁²/2 - 9*m₁^3) 
        a_num += G_r²*m₁*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₂v₁ - 2*v₂v₁² + 9/2*nv₁²*v₁² + 
                            4*v₂v₁*v₁² - 2v₁^4)
        a_num += G²_r³*m₂m₁*(39/2*nv₂² - 39*nv₂*nv₁ + 17/2*nv₁² - 15/4*v₂² - 5/2*v₂v₁ + 5/4*v₁²) 
        a_num += G²_r³*m₁^2*(2*nv₂² - 4*nv₂*nv₁ - 6*nv₁² - 8*v₂v₁ + 4v₁²) 
        a₁1 = (-n)*a_num

        a_num = G²_r³*m₁^2*(-2*nv₂ - 2*nv₁) + G²_r³*m₂m₁*(-63/4*nv₂ + 55/4*nv₁) 
        a_num += G_r²*m₁*(-6*nv₂*nv₁² + 9/2*nv₁^3 + nv₁*v₂² - 4*nv₂*v₂v₁ + 
                            4*nv₁*v₂v₁ + 4*nv₂*v₁² - 5*nv₁*v₁²)
        a₂2 = (-v̄)*a_num
        a₁1 + a₂2
    end


    dvi .+= ai*c⁻⁴
    dvj .+= aj*c⁻⁴
    nothing
end



function PN2p5_acceleration!(dvi, dvj, rs, vs,
                             pair::Tuple{Int, Int},
                             params::SimulationParams)                            
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

    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    #################### PN-2.5 acceleration ####################
    # acceleration for body 1 (i)
    a_num = 208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²
    a1 = a_num*n
    a_num = 8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²
    a2 = a_num*v̄
    ai = a1 + a2


    # acceleration for body 2 (j)
    a_num = 208G³_r⁴*m₂m₁²/15*nv - 24G³_r⁴*m₂²m₁/5*nv + 12G²_r³*m₂m₁/5*v²
    a1 = a_num*(-n)
    a_num = 8G³_r⁴*m₂²m₁/5 - 32G³_r⁴*m₂m₁²/5 - 4G²_r³*m₂m₁/5*v²
    a2 = a_num*(-v̄)
    aj = a1 + a2
    ###############################################################


    dvi .+= ai*c⁻⁵
    dvj .+= aj*c⁻⁵
    nothing
end

function PN1_to_2p5_acceleration!(dvi, dvj, rs, vs,
                                  pair::Tuple{Int, Int},
                                  params::SimulationParams)                           
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
    r² = r^2
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

    G = UNITLESS_G
    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    ai_PN1, ai_PN2, ai_PN2p5 = let
        
        #################### PN-1 acceleration ##################
        PN1 = (5G²*m₁*m₂*r⁻¹^3 + 4G²*m₂^2*r⁻¹^3 + G*m₂*r⁻¹^2*(3/2*nv₂^2 - v₁² + 4v₁v₂ - 2v₂²) )*n + G*m₂*r⁻¹^2*(4nv₁ - 3nv₂)*v̄
        #########################################################

        #################### PN-2 acceleration ##################
        a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
        a_num += G_r²*m₂*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
                            4*v₁v₂*v₂² - 2v₂^4)
        a_num += G²_r³*m₁m₂*(39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂² - 15/4*v₁² - 5/2*v₁v₂ + 5/4*v₂²) 
        a_num += G²_r³*m₂^2*(2*nv₁² - 4*nv₁*nv₂ - 6*nv₂² - 8*v₁v₂ + 4v₂²) 
        a₂1 = n*a_num
    
        a_num = G²_r³*m₂^2*(-2*nv₁ - 2*nv₂) + G²_r³*m₁m₂*(-63/4*nv₁ + 55/4*nv₂) 
        a_num += G_r²*m₂*(-6*nv₁*nv₂² + 9/2*nv₂^3 + nv₂*v₁² - 4*nv₁*v₁v₂ + 
                            4*nv₂*v₁v₂ + 4*nv₁*v₂² - 5*nv₂*v₂²)
        a₂2 = v̄*a_num
        PN2 = a₂1 + a₂2
        #########################################################


        ################### PN-2.5 acceleration #################
        a_num = 208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²
        a1 = a_num*n
        a_num = 8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²
        a2 = a_num*v̄
        PN2p5 = a1 + a2
        #########################################################


        PN1, PN2, PN2p5
    end

    
    
    aj_PN1, aj_PN2, aj_PN2p5 = let n = -n, v̄ = -v̄, nv₁ = -nv₁, nv₂ = -nv₂

        #################### PN-1 acceleration ##################
        PN1 = (5G²*m₂*m₁*r⁻¹^3 + 4G²*m₁^2*r⁻¹^3 + G*m₁*r⁻¹^2*(3/2*nv₁^2 - v₂² + 4v₂v₁ - 2v₁²) )*n + G*m₁*r⁻¹^2*(4nv₂ - 3nv₁)*v̄
        #########################################################


        #################### PN-2 acceleration ##################
        a_num = G³_r⁴*(-57*m₂²m₁/4 - 69*m₂m₁²/2 - 9*m₁^3) 
        a_num += G_r²*m₁*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₂v₁ - 2*v₂v₁² + 9/2*nv₁²*v₁² + 
                            4*v₂v₁*v₁² - 2v₁^4)
        a_num += G²_r³*m₂m₁*(39/2*nv₂² - 39*nv₂*nv₁ + 17/2*nv₁² - 15/4*v₂² - 5/2*v₂v₁ + 5/4*v₁²) 
        a_num += G²_r³*m₁^2*(2*nv₂² - 4*nv₂*nv₁ - 6*nv₁² - 8*v₂v₁ + 4v₁²) 
        a₁1 = n*a_num

        a_num = G²_r³*m₁^2*(-2*nv₂ - 2*nv₁) + G²_r³*m₂m₁*(-63/4*nv₂ + 55/4*nv₁) 
        a_num += G_r²*m₁*(-6*nv₂*nv₁² + 9/2*nv₁^3 + nv₁*v₂² - 4*nv₂*v₂v₁ + 
                            4*nv₁*v₂v₁ + 4*nv₂*v₁² - 5*nv₁*v₁²)
        a₂2 = v̄*a_num
        PN2 = a₁1 + a₂2
        #########################################################


        ################## PN-2p5 acceleration ##################
        a_num = 208G³_r⁴*m₂m₁²/15*nv - 24G³_r⁴*m₂²m₁/5*nv + 12G²_r³*m₂m₁/5*v²
        a1 = a_num*n
        a_num = 8G³_r⁴*m₂²m₁/5 - 32G³_r⁴*m₂m₁²/5 - 4G²_r³*m₂m₁/5*v²
        a2 = a_num*v̄
        PN2p5 = a1 + a2
        #########################################################

        PN1, PN2, PN2p5
    end

    dvi .+= ai_PN1*c⁻² + ai_PN2*c⁻⁴ + ai_PN2p5*c⁻⁵
    dvj .+= aj_PN1*c⁻² + aj_PN2*c⁻⁴ + aj_PN2p5*c⁻⁵

    nothing
end


###################################################################################################################