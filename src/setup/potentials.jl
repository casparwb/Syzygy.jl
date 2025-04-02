# abstract type Body end
using StaticArrays, JLD2, Printf
using LinearAlgebra: dot, norm, ×

include("../physics/tides.jl")

abstract type MultiBodyPotential end
abstract type SpinPotential <: MultiBodyPotential end
abstract type SimulationParams end

struct DefaultSimulationParams{FloatVecType, IntVecType} <: SimulationParams
    R::FloatVecType # radii
    M::FloatVecType # masses
    L::FloatVecType # luminosities
    stellar_types::IntVecType 
    M_cores::FloatVecType # core masses
    R_cores::FloatVecType # core radii
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

struct PureGravitationalPotential <: MultiBodyPotential end

"""
    DynamicalTidalPotential(tidal_force_power_index, polytropic_index)


Set up the dynamical tidal potential for a system as defined by Samsing, Leigh & Trani 2018. 

# Keyword arguments

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

struct EquilibriumTidalPotential <: MultiBodyPotential end

struct StaticEquilibriumTidalPotential{M_env_Type, R_env_Type, kT, ΩT} <: MultiBodyPotential
    M_env::M_env_Type
    R_env::R_env_Type
    apsidal_motion_constant::kT
    rotational_angular_velocity::ΩT
    # logg# # convert from R⊙/yr² to cm/s²
    # logm = log10(m₁)
end

function StaticEquilibriumTidalPotential(system; Z=0.02, lb_multiplier=1.1, ub_multiplier=1.1)

    age = system.time
    n_bodies = system.n
    R_envs = Float64[]
    m_envs = Float64[]

    for i = 1:n_bodies
        particle = system.particles[i]
        envelope_radius, envelope_mass = if particle.structure.stellar_type isa Star && particle.structure.m < 1.25u"Msun"
                                            envelope_structure(system.particles[i], age, Z)
                                         else
                                            0.0u"Rsun", 0.0u"Rsun"
                                         end

        push!(R_envs, ustrip(envelope_radius))
        push!(m_envs, ustrip(envelope_mass))
    end
    
    R_envs = SA[R_envs...]
    m_envs = SA[m_envs...]

    k_interpolator = get_k_interpolator(Z=Z, lb_multiplier=lb_multiplier, ub_multiplier=ub_multiplier)
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
            g = GRAVCONST*m/R^2

            logg = log10(ustrip(u"cm/s^2", g))
            logm = log10(ustrip(u"Msun", m))
            try
                k = k_interpolator(logm, logg)
                push!(apsidal_motion_constants, k)
            catch e
                # if e isa ArgumentError
                #     k = get_k_interpolator(Z=Z, lb_multiplier=1.1, ub_multiplier=1.1)
                # end
                throw(e)
            end

            Ω = stellar_spin(m, R)
            push!(rotational_angular_velocities, ustrip(unit(1/unit_time), Ω))

        end
    end

    apsidal_motion_constants = SA[apsidal_motion_constants...]
    rotational_angular_velocities = SA[rotational_angular_velocities...]

    StaticEquilibriumTidalPotential(m_envs, R_envs, apsidal_motion_constants, rotational_angular_velocities)
end


struct PN1Potential            <: MultiBodyPotential end

struct PN2Potential            <: MultiBodyPotential end

struct PN2p5Potential          <: MultiBodyPotential end

struct PN3Potential            <: MultiBodyPotential end

struct PN3_5Potential          <: MultiBodyPotential end

struct PNPotential             <: MultiBodyPotential end

##########################################################

struct PN1p5SpinPotential         <: SpinPotential end

struct PN2SpinPotential           <: SpinPotential end

struct PN2p5SpinPotential         <: SpinPotential end

struct PN1SpinPrecessionPotential <: SpinPotential end

struct PN1p5SpinPrecessionPotential <: SpinPotential end

struct PN2SpinPrecessionPotential <: SpinPotential end

struct SpinPrecessionPotential    <: SpinPotential end

###################################################################################################################


############################################## Acceleration functions #############################################


"""
    pure_gravitational_acceleration!(dvi, dvj, rs, pair::Tuple{Int, Int}, params::SimulationParams)

Gravitational acceleration on bodies i and j, with `(i, j) = pair`.
"""
function pure_gravitational_acceleration!(dvi,
                                          dvj,
                                          rs,
                                          pair::Tuple{Int, Int},
                                          params::SimulationParams)
    
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]

    r̄ = r̄₁ - r̄₂
    r = norm(r̄)
    n̂ = r̄/r

    m₁ = params.M[i]
    m₂ = params.M[j]
    G_r² = -UNITLESS_G/r^2

    a = G_r²*n̂

    a₁ =  a*m₂
    a₂ = -a*m₁

    dvi .+= a₁
    dvj .+= a₂
    nothing
end



"""
    dynamical_tidal_drag_force!(dv, rs, vs, params::SimulationParams, i::Int, n::Int, potential::DynamicalTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Implementing Tidal and Gravitational Wave Energy Losses in Few-body Codes: A Fast and Easy Drag Force Model](https://arxiv.org/abs/1803.08215)
"""
function dynamical_tidal_drag_force!(dvi,
                                     dvj,
                                     rs,
                                     vs,
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

    ms = params.M
    Rs = params.R

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
function equilibrium_tidal_drag_force!(dvi,
                                       dvj,
                                       rs,
                                       vs,
                                       pair::Tuple{Int, Int},
                                       params::SimulationParams,
                                       potential::EquilibriumTidalPotential) 

    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    ms = params.M
    Rs = params.R

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
        if !(stellar_types[stellar_type] isa Star)
            SA[0.0, 0.0, 0.0]
        else
            S̄₁  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

            core_mass   = params.M_cores[k]
            core_radius = params.R_cores[k]
            luminosity  = params.L[k]
            age         = params.ages[k]

            R = Rs[k]

            Ω = norm(S̄₁)
            
            logg = log10(surface_gravity_unit_conversion_factor*(UNITLESS_G*m₁/R^2)) # convert from default units to cm/s²
            logm = log10(m₁)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = UNITLESS_G*m₂/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, age, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₂, a)

            kτ = R^3/(UNITLESS_G*m₁)*k_T

            @. -3μ*m₂/m₁*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end

    # tidal force on 2 by 1
    a₂ = let k = j
        stellar_type = params.stellar_types[k]

        if !(stellar_types[stellar_type] isa Star)
            SA[0.0, 0.0, 0.0]
        else
            S̄₂  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

            core_mass   = params.M_cores[k]
            core_radius = params.R_cores[k]
            luminosity  = params.L[k]
            age         = params.ages[k]

            R = Rs[k]

            Ω = norm(S̄₂)
            
            logg = log10(surface_gravity_unit_conversion_factor*(UNITLESS_G*m₂/R^2)) # convert from R⊙/yr² to cm/s²
            logm = log10(m₂)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = UNITLESS_G*m₁/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, age, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₁, a)

            kτ = R^3/(UNITLESS_G*m₂)*k_T

            @. -3μ*m₁/m₂*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end
    
    dvi .+= a₁
    dvj .+= a₂
    nothing
end


function equilibrium_tidal_drag_force!(dvi,
                                       dvj,
                                       rs,
                                       vs,
                                       pair::Tuple{Int, Int},
                                       params::SimulationParams,
                                       potential::StaticEquilibriumTidalPotential) 
                               

    i, j = pair
    stellar_type_1 = params.stellar_types[i]
    stellar_type_2 = params.stellar_types[j]

    # println(i, " ", stellar_type_1)
    # println(j, " ", stellar_type_2)

    if (stellar_type_1 > 10) && (stellar_type_2 > 10) # tides are (currently) only for stars
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
    # r² = r^2
    r_hat = r̄*r⁻¹

    # ms = params.M
    m₁ = params.M[i]
    m₂ = params.M[j]
    # Rs = params.R
    # R₁ = params.R[i]
    # R₂ = params.R[j]

    # Ω₁ = potential.rotational_angular_velocity[i]
    # Ω₂ = potential.rotational_angular_velocity[j]

    # m₁, m₂ = ms[i], ms[j]

    θ_dot = (r̄ × v̄)*r⁻¹*r⁻¹
    θ_dot_norm = norm(θ_dot)
    θ_hat = θ_dot/θ_dot_norm

    a = semi_major_axis(r, v^2, m₂+m₁)

    # k₁ = potential.apsidal_motion_constant[i]
    # k₂ = potential.apsidal_motion_constant[j]

    # tidal force on 1 by 2
    a₁ = let
        if !(stellar_type_1 < 10)
                SA[0.0, 0.0, 0.0]
            else    
                μ = UNITLESS_G*m₂*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[i]
                envelope_radius = potential.R_env[i]
                luminosity = params.L[i]
                R = params.R[i]
                Ω = potential.rotational_angular_velocity[i]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, envelope_mass, envelope_radius,
                                                                   stellar_type_1, luminosity, 
                                                                   m₂, a)
                k = potential.apsidal_motion_constant[i]
                kτ = R^3/(UNITLESS_G*m₁)*k_T
    
                @. -3μ*m₂/m₁*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    # # tidal force on 2 by 1
    a₂ = let 
        if !(stellar_type_2 < 10)
                SA[0.0, 0.0, 0.0]
            else    
                μ = UNITLESS_G*m₁*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[j]
                envelope_radius = potential.R_env[j]
                luminosity  = params.L[j]
                R = params.R[j]
                Ω = potential.rotational_angular_velocity[j]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, envelope_mass, envelope_radius,
                                                                   stellar_type_1, luminosity, 
                                                                   m₁, a)
                k = potential.apsidal_motion_constant[j]
                kτ = R^3/(UNITLESS_G*m₂)*k_T
    
                @. -3μ*m₁/m₂*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    dvi .+= a₁
    dvj .+= a₂
    nothing
end

function PN1_acceleration!(dvi, 
                           dvj,
                           rs,
                           vs,
                           pair::Tuple{Int, Int},
                           params::SimulationParams)
                           
    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    v₂ = norm(v̄₂)
    
    m₁ = params.M[i]
    m₂ = params.M[j]
        
    v₁² = v₁^2
    v₂² = v₂^2

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r

    n = r̄*r⁻¹

    v₁v₂ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    G_r = UNITLESS_G*r⁻¹
    G_r² = G_r*r⁻¹

    # PN-1 acceleration
    ai = @. n*(G_r²*m₂)*(5*G_r*m₁ + 4*G_r*m₂ + 3/2*nv₂^2 - v₁² + 4*v₁v₂ - 2*v₂²) +
            (4*nv₁ - 3*nv₂)*v̄

    aj = @. (-n)*(G_r²*m₁)*(5*G_r*m₂ + 4*G_r*m₁ + 3/2*nv₁^2 - v₂² + 4*v₁v₂ - 2*v₁²) +
            (4*(-nv₂) - 3*(-nv₁))*(-v̄)
    
    dvi .+= ai*c⁻²
    dvj .+= aj*c⁻²
    nothing
end



function PN2_acceleration!(dvi, 
                           dvj,
                           rs,
                           vs,
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

    m₁ = params.M[i]
    m₂ = params.M[j]
    
    # i = 1, j = 2
    # add @fastmath?

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r

    r² = r^2

    n = r̄*r⁻¹

    v₁v₂ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    v₁v₂² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₁*m₂
    m₁²m₂ = m₁^2*m₂
    m₁m₂² = m₁*m₂^2

    G_r = UNITLESS_G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    # PN-2 acceleration:
    # expression is split up to avoid allocations that can appear in long expressions

    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
    a_num += UNITLESS_G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
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
    a_num = -57G³_r⁴*m₁m₂²/4 - 69G³_r⁴*m₁²m₂/2 - 9G³_r⁴*m₁^3 
    a_num += UNITLESS_G*m₁/r²*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₁v₂ - 2*v₁v₂² + 9/2*nv₁²*v₁² + 
                        4*v₁v₂*v₁² - 2v₂^4)
    a_num += G²_r³*m₁m₂*(39/2*nv₂² - 39*(-nv₂)*(-nv₁) + 17/2*nv₁² - 15/4*v₂² - 5/2*v₁v₂ + 5/4*v₁²) 
    a_num += G²_r³*m₁^2*(2*nv₂² - 4*(-nv₂)*(-nv₁) - 6*nv₁² - 8*v₁v₂ + 4v₁²) 
    a₂1 = (-n)*a_num

    a_num = G²_r³*m₁^2*(-2*(-nv₂) - 2*(-nv₁)) + G²_r³*m₁m₂*(-63/4*(-nv₂) + 55/4*(-nv₁)) 
    a_num += G_r²*m₁*(-6*(-nv₂)*nv₁² + 9/2*(-nv₁)^3 + (-nv₁)*v₂² - 4*(-nv₂)*v₁v₂ + 
                        4*(-nv₁)*v₁v₂ + 4*(-nv₂)*v₁² - 5*(-nv₁)*v₁²)
    a₂2 = (-v̄)*a_num
    aj = a₂1 + a₂2

    dvi .+= ai*c⁻⁴
    dvj .+= aj*c⁻⁴
    nothing
end



function PN2p5_acceleration!(dvi,
                             dvj,
                             rs,
                             vs,
                             pair::Tuple{Int, Int},
                             params::SimulationParams)                            
    # i = 1, j = 2
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    m₁ = params.M[i]
    m₂ = params.M[j]

    # add @fastmath?

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂

    r⁻¹ = 1/r

    v² = v^2

    n = r̄*r⁻¹

    nv = dot(n, v̄)

    m₁m₂ = m₁*m₂
    m₁²m₂ = m₁^2*m₂
    m₁m₂² = m₁*m₂^2

    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    # # PN-2.5 acceleration
    # acceleration for body 1 (i)
    a_num = 208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²
    a1 = a_num*n
    a_num = 8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²
    a2 = a_num*v̄
    ai = a1 + a2

    # acceleration for body 2 (j)
    a_num = 208G³_r⁴*m₁²m₂/15*nv - 24G³_r⁴*m₁m₂²/5*nv + 12G²_r³*m₁m₂/5*v²
    a1 = a_num*(-n)
    a_num = 8G³_r⁴*m₁m₂²/5 - 32G³_r⁴*m₁²m₂/5 - 4G²_r³*m₁m₂/5*v²
    a2 = a_num*(-v̄)
    aj = a1 + a2

    dvi .+= ai*c⁻⁵
    dvj .+= aj*c⁻⁵
    nothing
end


function PN1_to_2p5_acceleration!(dvi,
                                 dvj,
                                 rs,
                                 vs,
                                 pair::Tuple{Int, Int},
                                 params::SimulationParams)                           
    i, j = pair

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
        
    m₁ = params.M[i]
    m₂ = params.M[j]
    
    # m = m₁ + m₂

    # com = (m₁*r̄₁ + m₂*r̄₂)/m
    # v_com = (m₁*v̄₁ + m₂*v̄₂)/m

    # r̄₁ = r̄₁ - com
    # v̄₁ = v̄₁ - v_com

    # r̄₂ = r̄₂ - com
    # v̄₂ = v̄₂ - v_com

    # println(r̄₁)

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

    v₁v₂ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)
    nv = dot(n, v̄)

    v₁v₂² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₁*m₂
    m₁²m₂ = m₁^2*m₂
    m₁m₂² = m₁*m₂^2

    G_r = UNITLESS_G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    ################## PN-1 acceleration ##################
    ai = n*(G_r²*m₂)*(5*G_r*m₁ + 4*G_r*m₂ + 3/2*nv₂^2 - v₁² + 4*v₁v₂ - 2*v₂²) +
         (4*nv₁ - 3*nv₂)*v̄

    aj = (-n)*(G_r²*m₁)*(5*G_r*m₂ + 4*G_r*m₁ + 3/2*nv₁^2 - v₂² + 4*v₁v₂ - 2*v₁²) +
         (4*(-nv₂) - 3*(-nv₁))*(-v̄)
    
    dvi .+= ai*c⁻²
    dvj .+= aj*c⁻²
    #########################################################


    ################## PN-2 acceleration ##################
    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
    a_num += UNITLESS_G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
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
    a_num = -57G³_r⁴*m₁m₂²/4 - 69G³_r⁴*m₁²m₂/2 - 9G³_r⁴*m₁^3 
    a_num += UNITLESS_G*m₁/r²*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₁v₂ - 2*v₁v₂² + 9/2*nv₁²*v₁² + 
                        4*v₁v₂*v₁² - 2v₂^4)
    a_num += G²_r³*m₁m₂*(39/2*nv₂² - 39*(-nv₂)*(-nv₁) + 17/2*nv₁² - 15/4*v₂² - 5/2*v₁v₂ + 5/4*v₁²) 
    a_num += G²_r³*m₁^2*(2*nv₂² - 4*(-nv₂)*(-nv₁) - 6*nv₁² - 8*v₁v₂ + 4v₁²) 
    a₂1 = (-n)*a_num

    a_num = G²_r³*m₁^2*(-2*(-nv₂) - 2*(-nv₁)) + G²_r³*m₁m₂*(-63/4*(-nv₂) + 55/4*(-nv₁)) 
    a_num += G_r²*m₁*(-6*(-nv₂)*nv₁² + 9/2*(-nv₁)^3 + (-nv₁)*v₂² - 4*(-nv₂)*v₁v₂ + 
                        4*(-nv₁)*v₁v₂ + 4*(-nv₂)*v₁² - 5*(-nv₁)*v₁²)
    a₂2 = (-v̄)*a_num
    aj = a₂1 + a₂2

    dvi .+= ai*c⁻⁴
    dvj .+= aj*c⁻⁴
    #########################################################


    ################## PN-2.5 acceleration ##################
    # acceleration for body 1 (i)
    a_num = 208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²
    a1 = a_num*n
    a_num = 8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²
    a2 = a_num*v̄
    ai = a1 + a2

    # acceleration for body 2 (j)
    a_num = 208G³_r⁴*m₁²m₂/15*nv - 24G³_r⁴*m₁m₂²/5*nv + 12G²_r³*m₁m₂/5*v²
    a1 = a_num*(-n)
    a_num = 8G³_r⁴*m₁m₂²/5 - 32G³_r⁴*m₁²m₂/5 - 4G²_r³*m₁m₂/5*v²
    a2 = a_num*(-v̄)
    aj = a1 + a2

    dvi .+= ai*c⁻⁵
    dvj .+= aj*c⁻⁵
    #########################################################

    nothing
end


###################################################################################################################