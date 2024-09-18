# abstract type Body end
using StaticArrays, JLD2, Printf
using LinearAlgebra: dot, norm, ×

include("../physics/tides.jl")

abstract type MultiBodyPotential end
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


struct PureGravitationalPotential <: MultiBodyPotential end

struct DynamicalTidalPotential{nType, fType <: Function} <: MultiBodyPotential
    nₜ::Int  # Tidal force power constant
    γ::nType # Polytropic index of each star
    tidal_factor::fType
end

"""
    DynamicalTidalPotential(;G, n, γ)


Set up the dynamical tidal potential for a system. 

# Keyword arguments
- `G`: gravitational constant.
- `n`: tidal force power index
- `γ`: vector of polytropic indices of each body in the system
"""
function DynamicalTidalPotential(;n, γ)

    if n == 4
        f = tidal_factor_n4
    elseif n == 10
        f = tidal_factor_n10
    else
        @warn "Tidal factor for n ≠ (4, 10) is not defined."
        f = x -> x
    end

    DynamicalTidalPotential(n, γ, f)
end

struct EquilibriumTidalPotential <: MultiBodyPotential end

struct StaticEquilibriumTidalPotential{M_env_Type, R_env_Type} <: MultiBodyPotential
    M_env::M_env_Type
    R_env::R_env_Type
end

function StaticEquilibriumTidalPotential(system; Z=0.02)

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

    StaticEquilibriumTidalPotential(m_envs, R_envs)
end

"""
    PN1Potential{gType <: Real}

Potential for Post-Newtonian 1 (PN1) acceleration. 
"""
struct PN1Potential{gType <: Real} <: MultiBodyPotential
    G::gType
    PN1Potential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end

"""
    PN2Potential{gType <: Real}

Potential for Post-Newtonian 2 (PN2) acceleration.
"""
struct PN2Potential{gType <: Real} <: MultiBodyPotential
    G::gType
    PN2Potential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end

"""
    PN1Potential{gType <: Real}

Potential for Post-Newtonian 2.5 (PN2.5) acceleration.
"""
struct PN2_5Potential{gType <: Real} <: MultiBodyPotential
    G::gType
    PN2_5Potential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end

struct PN3Potential{gType <: Real} <: MultiBodyPotential
    G::gType
    PN3Potential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end

struct PN3_5Potential{gType <: Real} <: MultiBodyPotential
    G::gType
    PN3_5Potential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end


"""
    PN1Potential{gType <: Real}

Potential for Post-Newtonian 1 to 2.5 (PN1, PN2, PN2.5) acceleration
"""
struct PNPotential{gType <: Real} <: MultiBodyPotential
    G::gType
    PNPotential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)
end

"""
    SpinPrecessionPotential{gType <: Real}

Potential for spin precession.
"""
struct SpinPrecessionPotential <: MultiBodyPotential end



"""
    pure_gravitational_acceleration!(dv,, rs,, params::SimulationParams,, i::Integer,, n::Integer,, potential::PureGravitationalPotential)

Acceleration function from gravitational acceleration.
"""
# function pure_gravitational_acceleration!(dv,
#                                           rs,
#                                           params::SimulationParams,
#                                           i::Int,
#                                           n::Int,
#                                           potential::PureGravitationalPotential)
#     accel = SA[0.0, 0.0, 0.0];
#     ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
#     @inbounds for j = 1:n
#         if j != i
#             m_num::Float64 = ustrip(params.M[j])
#             rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
#             rij = ri - rj
#             accel -= G * m_num * rij / (norm(rij)^3)
#         end
#     end
#     @. dv += accel
# end

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
    G_r² = -G/r^2

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
    a = semi_major_axis(r, v^2, M, G)
    e = eccentricity(r̄, v̄, a, M, G)
    rₚ = a*(1 - e)

    J = potential.tidal_factor(e)

    a₁ = let
        Rₜ = Rs[i]
        ΔE = tidal_ΔE(m₁, Rₜ, m₂, rₚ, potential.γ[i], G)

        ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
        ε = drag_force_coefficient(ΔE, J, a, e, M, potential.nₜ, G)

        F₁₂ = @. (-ε*(v/r^potential.nₜ)*v̄/v)
        F₁₂ / m₁
    end

    a₂ = let
        Rₜ = Rs[j]
        ΔE = tidal_ΔE(m₂, Rₜ, m₁, rₚ, potential.γ[j], G)
        
        ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
        ε = drag_force_coefficient(ΔE, J, a, e, M, potential.nₜ, G)

        F₂₁ = @. (-ε*(v/r^potential.nₜ)*(-v̄)/v)
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

    a = semi_major_axis(r, v^2, m₂+m₁, G)

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
            
            logg = log10(6.985766564066957e-5*(G*m₁/R^2)) # convert from R⊙/yr² to cm/s²
            logm = log10(m₁)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = G*m₂/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, age, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₂, a)

            kτ = R^3/(G*m₁)*k_T

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
            
            logg = log10(6.985766564066957e-5*(G*m₂/R^2)) # convert from R⊙/yr² to cm/s²
            logm = log10(m₂)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = G*m₁/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, age, core_mass, core_radius, 
                                                               stellar_type, luminosity, m₁, a)

            kτ = R^3/(G*m₂)*k_T

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
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄)
    v = norm(v̄)
    
    r² = r^2
    r_hat = r̄/r

    ms = params.M
    Rs = params.R

    m₁, m₂ = ms[i], ms[j]

    θ_dot = (r̄ × v̄)/r²
    θ_dot_norm = norm(θ_dot)
    θ_hat = θ_dot/θ_dot_norm

    a = semi_major_axis(r, v^2, m₂+m₁, G)

    # tidal force on 1 by 2
    a₁ = let k = i
        stellar_type = params.stellar_types[k]

        if !(stellar_types[stellar_type] isa Star)
            SA[0.0, 0.0, 0.0]
        else
            S̄₁  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

            envelope_mass = potential.M_env[k]
            envelope_radius = potential.R_env[k]
            luminosity  = params.L[k]

            R = Rs[k]

            Ω = norm(S̄₁)
            
            logg = log10(6.985766564066957e-5*(G*m₁/R^2)) # convert from R⊙/yr² to cm/s²
            logm = log10(m₁)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = G*m₂/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, envelope_mass, envelope_radius,
                                                               stellar_type, luminosity, 
                                                               m₂, a)

            kτ = R^3/(G*m₁)*k_T

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

            envelope_mass = potential.M_env[k]
            envelope_radius = potential.R_env[k]
            luminosity  = params.L[k]

            R = Rs[k]

            Ω = norm(S̄₂)
            
            logg = log10(6.985766564066957e-5*(G*m₂/R^2)) # convert from R⊙/yr² to cm/s²
            logm = log10(m₂)
            k = asidal_motion_constant_interpolated(logm, logg)

            μ = G*m₁/r²
            
            k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, envelope_mass, envelope_radius,
                                                               stellar_type, luminosity, 
                                                               m₁, a)

            kτ = R^3/(G*m₂)*k_T

            @. -3μ*m₁/m₂*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
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

    G_r = G*r⁻¹
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
    r³ = r²*r

    n = r̄*r⁻¹

    v₁v₂ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    v₁v₂² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    # nv₁³ = nv₁^3
    nv₂³ = nv₂^3

    nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₁*m₂
    m₁²m₂ = m₁^2*m₂
    m₁m₂² = m₁*m₂^2

    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    # PN-2 acceleration:
    # expression is split up to avoid allocations that can appear in long expressions

    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
    a_num += G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
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
    a_num += G*m₁/r²*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₁v₂ - 2*v₁v₂² + 9/2*nv₁²*v₁² + 
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

function PN2_5_acceleration!(dvi,
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

function PN3_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN3Potential)

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

     
    
    m₁ = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    a₂ =  @MVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @inbounds for j = 1:n
        if j != i   
            r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂ = norm(v̄₂)

            v₂² = v₂^2

            r̄ = r̄₁ - r̄₂
            v̄ = v̄₁ - v̄₂

            r = norm(r̄) # r₁₂
            v = norm(v̄) # v₁₂

            r⁻¹ = 1/r

            r² = r^2
            r³ = r²*r
            r⁴ = r³*r
            r⁵ = r⁴*r
            
            v² = v^2
            n = r̄*r⁻¹

            v₁v₂ = dot(v̄₁, v̄₂) 
            nv₁ = dot(n, v̄₁)
            nv₂ = dot(n, v̄₂)
            nv = dot(n, v̄)
            # nv² = nv^2

            v₁v₂² = v₁v₂^2 

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            # nv₁³ = nv₁^3
            nv₂³ = nv₂^3

            # nv₁⁴ = nv₁^4
            nv₂⁴ = nv₂^4

            m₂ = params.M[j]
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2
            m₁²m₂² = m₁^2*m₂^2

            G_r = G*r⁻¹
            G_r² = G_r*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            r₁′ = r₂′ = 1.0
      
            # TO-DO: split up to avoid allocations; calculate the gauge constants r′
            a = @. n*(G_r²*m₂*(35/16*nv₂^6 - 15/8*nv₂⁴*v₁² + 15/2*nv₂⁴*v₁v₂ + 3*nv₂²*v₁v₂² -
                                     15/2*nv₂⁴*v₂² + 3/2*nv₂²*v₁²*v₂² - 12*nv₂²*v₁v₂*v₂² - 2*v₁v₂²*v₂² + 
                                     15/2*nv₂²*v₂^4 + 4*v₁v₂*v₂^4 - 2v₂^6
                                    ) +
                      G²_r³*m₁m₂*(-171/8*nv₁^4 + 171/2*nv₁³*nv₂ - 723/4*nv₁²*nv₂² +
                                          383/2*nv₁*nv₂³ - 455/8*nv₂⁴ + 229/4*nv₁²*v₁² - 
                                          205/2*nv₁*nv₁*v₁² + 191/4*nv₂²*v₁² - 91/8*v₁^4 - 229/2*nv₁²*v₁v₂ +
                                          244*nv₁*nv₂*v₁v₂ - 225/2*nv₂²*v₁v₂ + 91/2*v₁²*v₁v₂ -
                                          177/4*v₁v₂² + 229/4*nv₁²*v₂² - 283/2*nv₁*nv₂*v₂² +
                                          259/4*nv₂²*v₂² - 91/4*v₁²*v₂² + 43*v₁v₂*v₂² - 81/8*v₂^4
                                         ) +
                     G²_r³*m₂^2*(-6*nv₁²*nv₂² + 12*nv₁*nv₂³ + 6*nv₂⁴ + 
                                         4*nv₁*nv₂*v₁v₂ + 12*nv₂³*v₁v₂ + 4*v₁v₂ -
                                         4*nv₁*nv₂*v₂² - 12*nv₂²*v₂² - 8*v₁v₂*v₂² + 4v₂^4
                                       ) +
                     G³_r⁴*m₂^3*(-nv₁² + 2*nv₁*nv₂ + 43/2*nv₂² + 18*v₁v₂ - 9v₂²) +
                     G³_r⁴*m₁m₂²*(415/8*nv₁² - 375/4*nv₁*nv₂ + 1113/8*nv₂² - 615/64*nv²*π² +
                                           18v₁² + 123/64*π²*v² + 33*v₁v₂ - 33/2*v₂²) + 
                     G³_r⁴*m₁²m₂*(-45887/168*nv₁² + 24025/42*nv₁*nv₂ - 10469/42*nv₂² + 48197/840*v₁² -
                                           36227/420*v₁v₂ + 36227*v₂² + 110*nv²*log(r̄/r₁′) - 22*v²*log(r̄/r₁′)) + 
                     16G⁴*m₂^4/r⁵ + G⁴*m₁²m₂²/r⁵*(175 - 41/16*π² - 44/3*log(r̄/r₂′))) +
                     (G_r²*m₂*(15/2*nv₁*nv₂⁴ - 45/8*nv₂^5 - 3/2*nv₂³*v₁² + 6*nv₁*nv₂²*v₁v₂ -
                                     6*nv₂³*v₁v₂ - 2*nv₂*v₁v₂² - 12*nv₁*nv₂²*v₂² + 12*nv₂³*v₂² +
                                     nv₂*v₁²*v₂² - 4*nv₁*v₁v₂*v₂² + 8*nv₂*v₁v₂*v₂² + 4*nv₁*v₂^4 -
                                     7*nv₂*v₂^4) +
                      G²_r³*m₂^2*(-2*nv₁²*nv₂ + 8*nv₁*nv₂² + 2*nv₂³ + 2*nv₁*v₁v₂ + 
                                         4*nv₂*v₁v₂ - 2*nv₁*v₂² - 4*nv₂*v₂²) +
                      G²_r³*m₁m₂*(-243/4*nv₁³ + 565/4*nv₁²*nv₂ - 269/4*nv₁*nv₂² -
                                          95/12*nv₂³ + 207/8*nv₁*v₁² - 137/8*nv₂*v₁² - 36*nv₁*v₁v₂ + 
                                          27/4*nv₂*v₁v₂ + 81/8*nv₁*v₂² + 83/8*nv₂*v₂²) + 
                      G³_r⁴*m₂^3*(4*nv₁ + 5*nv₂) + 
                      G³_r⁴*m₁m₂²*(-307/8*nv₁ - 479/8*nv₂ + 123/32*nv*π²) + 
                      G³_r⁴*m₁²m₂*(31397/420*nv₁ - 36227/427*nv₂ - 44*nv*log(r/r₁′)))*v̄

            accel += a
        end
    end

    @. dv += accel * c⁻⁶
end

function PN3_5_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN3_5Potential)

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

     
    
    m₁ = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @inbounds for j = 1:n
        if j != i   
            r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂ = norm(v̄₂)

            v₂² = v₂^2

            r̄ = r̄₁ - r̄₂
            v̄ = v̄₁ - v̄₂

            r = norm(r̄) # r₁₂
            v = norm(v̄) # v₁₂

            r⁻¹ = 1/r

            r² = r^2
            r³ = r²*r
            r⁴ = r³*r
            r⁵ = r⁴*r
         
            v² = v^2


            n = r̄*r⁻¹

            v₁v₂ = dot(v̄₁, v̄₂) 
            nv₁ = dot(n, v̄₁)
            nv₂ = dot(n, v̄₂)
            nv = dot(n, v̄)
            # nv² = nv^2

            v₁v₂² = v₁v₂^2 

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            # nv₁³ = nv₁^3
            nv₂³ = nv₂^3

            # nv₁⁴ = nv₁^4

            m₂ = params.M[j]
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2
            m₁²m₂² = m₁^2*m₂^2

            G_r = G*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            # r₁′ = r₂′ = 1.0
      

            a = @. n*(G⁴*m₁^3*m₂/r⁵*(3992/105*nv₁ - 4328/105*nv₂) + 
                      G⁴*m₁²m₂²/r⁶*(-13576/105*nv₁ + 2872/21*nv₂) - 3172/21*G⁴*m₁*m₂^3/r⁶*nv +
                      G³_r⁴*m₁²m₂*(48*nv₁³ - 696/5*nv₁²*nv₂ + 744/5*nv₁*nv₂² - 288/5*nv₂³ -
                                            4888/105*nv₁*v₁² + 5056*nv₂*v₁² + 2056/21*nv₁*v₁v₂ -
                                            2224/21*nv₂*v₁v₂ - 1028/21*nv₁*v₂² + 5812/105*nv₂*v₂²) + 
                      G³_r⁴*m₁m₂²*(-582/5*nv₂³ + 1746/5*nv₁²*nv₂ - 1954/5*nv₁*nv₂² +
                                            158*nv₂³ + 3568/105*nv*v₁² - 2864/35*nv₁*v₁v₂ +
                                            10048/105*nv₂*v₁v₂ + 1432/35*nv₁*v₂² - 5752/105*nv₂*v₂²) +
                      G²_r³*m₁m₂*(-56*nv^5 + 60*nv₁³*v² - 180*nv₁²*nv₂*v² + 
                                          174*nv₁*nv₂²*v² - 54*nv₂³*v² - 246/35*nv*v₁^4 +
                                          1068/35*nv₁*v₁²*v₁v₂ - 984/35*nv₂*v₂²*v₁v₂ - 1068/35*nv₁*v₁v₂² +
                                          180/7*nv₂*v₁v₂² - 534/35*nv₁*v₁²*v₂² + 90/7*nv₂*v₁²*v₂² +
                                          984/35*nv₁*v₁v₂*v₂² - 732/35*nv₂*v₁v₂*v₂² - 204/35*nv₁*v₂^4 + 
                                          24/7*nv₂*v₂^4)) + 
                   v*(-184/21*G⁴*m₁^3*m₂/r⁵ + 6224/105*G⁴*m₁²m₂²/r⁶ + 6388/105*G⁴*m₁*m₂^3/r⁶ +
                      G³*m₁²m₂*(52/15*nv₁² - 56/15*nv₁*nv₂ - 44/15*nv₂² - 132/35*v₁² + 152/35*v₁v₂ - 48/35*v₂²) +
                      G³_r⁴*m₁m₂²*(454/15*nv₁² - 372/5*nv₁*nv₂ + 854/15*nv₂² - 152/21*v₁² + 
                                            2864/105*v₁v₂ - 1768/105*v₂²) +
                      G²_r³*m₁m₂*(60*nv^4 - 348/5*nv₁²*v² + 684/5*nv₁*nv₁*v² -
                                          66*nv₂²*v² + 334/35*v₁^4 - 1336/35*v₁²*v₁v₂ + 1308/35*v₁v₂² + 654/35*v₁²*v₂² -
                                          1252/35*v₁v₂*v₂² + 292/35*v₂^4))

            accel += a
        end
    end

    @. dv += accel * c⁻⁷
end

function PN1_to_3_5_acceleration!(dvi,
                                 dvj,
                                 rs,
                                 vs,
                                 pair::Tuple{Int, Int},
                                 params::SimulationParams)                           
    i, j = pair
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    v₂ = norm(v̄₂)
    
    v₁² = v₁^2
    v₂² = v₂^2

    m₁ = params.M[i]
    m₂ = params.M[j]

    # a₂ =  @MVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂

    r⁻¹ = 1/r

    r² = r^2
    r³ = r²*r
    # r⁴ = r³*r
    # r⁵ = r⁴*r
    # r⁶ = r⁵*r
    
    v² = v^2

    n = r̄*r⁻¹

    v₁v₂ = dot(v̄₁, v̄₂) 
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)
    nv = dot(n, v̄)
    # nv² = nv^2

    v₁v₂² = v₁v₂^2 

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    # nv₁³ = nv₁^3
    nv₂³ = nv₂^3

    # nv₁⁴ = nv₁^4
    nv₂⁴ = nv₂^4

    m₁m₂ = m₁*m₂
    m₁²m₂ = m₁^2*m₂
    m₁m₂² = m₁*m₂^2
    # m₁²m₂² = m₁^2*m₂^2

    G_r = G*r⁻¹
    G_r² = G_r*r⁻¹
    G²_r³ = G²*r⁻¹^3
    G³_r⁴ = G³*r⁻¹^4

    # r₁′ = r₂′ = 1.0

    ################## PN-1 acceleration ##################
    ai = @. n*(G_r²*m₂)*(5*G_r*m₁ + 4*G_r*m₂ + 3/2*nv₂^2 - v₁² + 4*v₁v₂ - 2*v₂²) +
            (4*nv₁ - 3*nv₂)*v̄

    aj = @. (-n)*(G_r²*m₁)*(5*G_r*m₂ + 4*G_r*m₁ + 3/2*nv₁^2 - v₂² + 4*v₁v₂ - 2*v₁²) +
            (4*(-nv₂) - 3*(-nv₁))*(-v̄)
    
    dvi .+= ai*c⁻²
    dvj .+= aj*c⁻²

    ################## PN-2 acceleration ##################
    # acceleration for body 1 (i)
    a_num = G³_r⁴*(-57*m₁²m₂/4 - 69*m₁m₂²/2 - 9*m₂^3) 
    a_num += G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
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
    a_num += G*m₁/r²*(-15/8*nv₁⁴ + 3/2*nv₁²*v₂² - 6*nv₁²*v₁v₂ - 2*v₁v₂² + 9/2*nv₁²*v₁² + 
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

    # a₄ = @. n*(G_r²*m₂*(35/16*nv₂^6 - 15/8*nv₂⁴*v₁² + 15/2*nv₂⁴*v₁v₂ + 3*nv₂²*v₁v₂² -
    #                          15/2*nv₂⁴*v₂² + 3/2*nv₂²*v₁²*v₂² - 12*nv₂²*v₁v₂*v₂² - 2*v₁v₂²*v₂² + 
    #                          15/2*nv₂²*v₂^4 + 4*v₁v₂*v₂^4 - 2v₂^6
    #                         ) +
    #           G²_r³*m₁m₂*(-171/8*nv₁^4 + 171/2*nv₁³*nv₂ - 723/4*nv₁²*nv₂² +
    #                               383/2*nv₁*nv₂³ - 455/8*nv₂⁴ + 229/4*nv₁²*v₁² - 
    #                               205/2*nv₁*nv₁*v₁² + 191/4*nv₂²*v₁² - 91/8*v₁^4 - 229/2*nv₁²*v₁v₂ +
    #                               244*nv₁*nv₂*v₁v₂ - 225/2*nv₂²*v₁v₂ + 91/2*v₁²*v₁v₂ -
    #                               177/4*v₁v₂² + 229/4*nv₁²*v₂² - 283/2*nv₁*nv₂*v₂² +
    #                               259/4*nv₂²*v₂² - 91/4*v₁²*v₂² + 43*v₁v₂*v₂² - 81/8*v₂^4
    #                              ) +
    #          G²_r³*m₂^2*(-6*nv₁²*nv₂² + 12*nv₁*nv₂³ + 6*nv₂⁴ + 
    #                              4*nv₁*nv₂*v₁v₂ + 12*nv₂³*v₁v₂ + 4*v₁v₂ -
    #                              4*nv₁*nv₂*v₂² - 12*nv₂²*v₂² - 8*v₁v₂*v₂² + 4v₂^4
    #                            ) +
    #          G³_r⁴*m₂^3*(-nv₁² + 2*nv₁*nv₂ + 43/2*nv₂² + 18*v₁v₂ - 9v₂²) +
    #          G³_r⁴*m₁m₂²*(415/8*nv₁² - 375/4*nv₁*nv₂ + 1113/8*nv₂² - 615/64*nv²*π² +
    #                                18v₁² + 123/64*π²*v² + 33*v₁v₂ - 33/2*v₂²) + 
    #          G³_r⁴*m₁²m₂*(-45887/168*nv₁² + 24025/42*nv₁*nv₂ - 10469/42*nv₂² + 48197/840*v₁² -
    #                                36227/420*v₁v₂ + 36227*v₂² + 110*nv²*log(r̄/r₁′) - 22*v²*log(r̄/r₁′)) + 
    #          16G⁴*m₂^4/r⁵ + G⁴*m₁²m₂²/r⁵*(175 - 41/16*π² - 44/3*log(r̄/r₂′))) +
    #          (G_r²*m₂*(15/2*nv₁*nv₂⁴ - 45/8*nv₂^5 - 3/2*nv₂³*v₁² + 6*nv₁*nv₂²*v₁v₂ -
    #                          6*nv₂³*v₁v₂ - 2*nv₂*v₁v₂² - 12*nv₁*nv₂²*v₂² + 12*nv₂³*v₂² +
    #                          nv₂*v₁²*v₂² - 4*nv₁*v₁v₂*v₂² + 8*nv₂*v₁v₂*v₂² + 4*nv₁*v₂^4 -
    #                          7*nv₂*v₂^4) +
    #           G²_r³*m₂^2*(-2*nv₁²*nv₂ + 8*nv₁*nv₂² + 2*nv₂³ + 2*nv₁*v₁v₂ + 
    #                              4*nv₂*v₁v₂ - 2*nv₁*v₂² - 4*nv₂*v₂²) +
    #           G²_r³*m₁m₂*(-243/4*nv₁³ + 565/4*nv₁²*nv₂ - 269/4*nv₁*nv₂² -
    #                               95/12*nv₂³ + 207/8*nv₁*v₁² - 137/8*nv₂*v₁² - 36*nv₁*v₁v₂ + 
    #                               27/4*nv₂*v₁v₂ + 81/8*nv₁*v₂² + 83/8*nv₂*v₂²) + 
    #           G³_r⁴*m₂^3*(4*nv₁ + 5*nv₂) + 
    #           G³_r⁴*m₁m₂²*(-307/8*nv₁ - 479/8*nv₂ + 123/32*nv*π²) + 
    #           G³_r⁴*m₁²m₂*(31397/420*nv₁ - 36227/427*nv₂ - 44*nv*log(r/r₁′)))*v

    # a₅ = @. n*(G⁴*m₁^3*m₂/r⁵*(3992/105*nv₁ - 4328/105*nv₂) + 
    #           G⁴*m₁²m₂²/r⁶*(-13576/105*nv₁ + 2872/21*nv₂) - 3172/21*G⁴*m₁*m₂^3/r⁶*nv +
    #           G³_r⁴*m₁²m₂*(48*nv₁³ - 696/5*nv₁²*nv₂ + 744/5*nv₁*nv₂² - 288/5*nv₂³ -
    #                                 4888/105*nv₁*v₁² + 5056*nv₂*v₁² + 2056/21*nv₁*v₁v₂ -
    #                                 2224/21*nv₂*v₁v₂ - 1028/21*nv₁*v₂² + 5812/105*nv₂*v₂²) + 
    #           G³_r⁴*m₁m₂²*(-582/5*nv₂³ + 1746/5*nv₁²*nv₂ - 1954/5*nv₁*nv₂² +
    #                                 158*nv₂³ + 3568/105*nv*v₁² - 2864/35*nv₁*v₁v₂ +
    #                                 10048/105*nv₂*v₁v₂ + 1432/35*nv₁*v₂² - 5752/105*nv₂*v₂²) +
    #           G²_r³*m₁m₂*(-56*nv^5 + 60*nv₁³*v² - 180*nv₁²*nv₂*v² + 
    #                               174*nv₁*nv₂²*v² - 54*nv₂³*v² - 246/35*nv*v₁^4 +
    #                               1068/35*nv₁*v₁²*v₁v₂ - 984/35*nv₂*v₂²*v₁v₂ - 1068/35*nv₁*v₁v₂² +
    #                               180/7*nv₂*v₁v₂² - 534/35*nv₁*v₁²*v₂² + 90/7*nv₂*v₁²*v₂² +
    #                               984/35*nv₁*v₁v₂*v₂² - 732/35*nv₂*v₁v₂*v₂² - 204/35*nv₁*v₂^4 + 
    #                               24/7*nv₂*v₂^4)) + 
    #        v*(-184/21*G⁴*m₁^3*m₂/r⁵ + 6224/105*G⁴*m₁²m₂²/r⁶ + 6388/105*G⁴*m₁*m₂^3/r⁶ +
    #           G³*m₁²m₂*(52/15*nv₁² - 56/15*nv₁*nv₂ - 44/15*nv₂² - 132/35*v₁² + 152/35*v₁v₂ - 48/35*v₂²) +
    #           G³_r⁴*m₁m₂²*(454/15*nv₁² - 372/5*nv₁*nv₂ + 854/15*nv₂² - 152/21*v₁² + 
    #                                 2864/105*v₁v₂ - 1768/105*v₂²) +
    #           G²_r³*m₁m₂*(60*nv^4 - 348/5*nv₁²*v² + 684/5*nv₁*nv₁*v² -
    #                               66*nv₂²*v² + 334/35*v₁^4 - 1336/35*v₁²*v₁v₂ + 1308/35*v₁v₂² + 654/35*v₁²*v₂² -
    #                               1252/35*v₁v₂*v₂² + 292/35*v₂^4))

    nothing
end

# function deSitter_precession!(dv,
#                                dvs,
#                                rs,
#                                vs,
#                                params::SimulationParams,
#                                i::Int,
#                                n::Int,
#                                potential::deSitterPotential)

#     ā₁  = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
#     r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
#     v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]
#     S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
#     dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]
#     v₁  = norm(v̄₁)

#     v₁² = v₁^2
    
#     m₁ = params.M[i]
#     accel = @SVector [0.0, 0.0, 0.0]
    
#     # i = 1, j = 2
#     # add @fastmath?
#     @inbounds for j = 1:n
#         if j != i   
#             m₂ = params.M[j]
#             μ = reduced_mass(m₁, m₂)

#             ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
#             r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
#             v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
#             v₂ = norm(v̄₂)

#             ā = ā₁ - ā₂
#             r̄ = r̄₁ - r̄₂
#             v̄ = v̄₁ - v̄₂

#             r = norm(r̄) # r₁₂
#             # v = norm(v̄) # v₁₂

#             n = r̄/r

#             L̄ = angular_momentum(r̄, μ*v̄)
#             L = norm(L̄)
#             L̂ = L̄/L

#             τ = μ/L*(r̄ × ā)

#             # Ωds = 3G*n*(m₂ + μ/3)/(2*c²*a*(1 - e^2))

#             accel += (S̄₁ × τ) + (L̂ × dS̄₁)
#         end 

#     end

#     @. dv = accel
# end

function spin_precession!(dvi,
                          dvj,
                          dvs,
                          rs,
                          vs,
                          pair::Tuple{Int, Int},
                          params::SimulationParams)

    # i = 1, j = 2
    i, j = pair
    ā₁  = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂ = @SVector [rs[1,  j], rs[2,  j], rs[3,  j]]
    v̄₂ = @SVector [vs[1,  j], vs[2,  j], vs[3,  j]]

    S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂  = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]
    
    # v₁  = norm(v̄₁)
    # v₁² = v₁^2
    
    m₁ = params.M[i]
    m₂ = params.M[j]
    δm = m₁ - m₂
    
    # add @fastmath?

    # v₂ = norm(v̄₂)
    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂
    # rxv = r̄ × v̄

    rv = dot(r̄, v̄)

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂

    r² = r*r
    r⁻¹ = 1/r
    r⁻² = 1/r²

    n    = r̄/r
    # nv   = dot(n, v̄)
    # nv₁  = dot(n, v̄₁)
    # nv₂  = dot(n, v̄₂)
    # vv₂  = dot(v̄, v̄₂)
    # nS₁  = dot(n, S̄₁)
    # vS₁  = dot(v̄, S̄₁)
    # v₁S₁ = dot(v̄₁, S̄₁)
    # v₂S₁ = dot(v̄₂, S̄₁)

    # nv₂² = nv₂^2
    
    # dr_dt    = rv*r⁻¹
    # dn_dt    = r⁻¹*(v̄ - rv .* r⁻²*r̄)
    # dnv_dt   = n  .* ā   .+ v  .* dn_dt
    # dnv₁_dt  = n  .* ā₁  .+ v̄₁ .* dn_dt
    # dnv₂_dt  = n  .* ā₂  .+ v̄₂ .* dn_dt
    # dnS₁_dt  = n  .* dS̄₁ .+ S̄₁ .* dn_dt
    # dvS₁_dt  = v̄  .* dS̄₁ .+ S̄₁ .* ā
    # dv₁S₁_dt = v̄₁ .* dS̄₁ .+ S̄₁ .* ā₁
    # dv₂S₁_dt = v̄₂ .* dS̄₁ .+ S̄₁ .* ā₂
    # dvv₂_dt  = v̄  .* ā₂  .+ v̄₂ .* ā
    
    Gm₁ = G*m₁
    Gm₂ = G*m₂

    # dT2PN_dt = @. -2*Gm₂*dr_dt*r⁻²*r⁻¹*((Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 
    #                                  2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁)*n + 
    #                                 (-5*G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv)*v̄₂ - 
    #                                 (-G*(6*δm)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + nv₂*vS₁)*v̄₁ + 
    #                                 (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ - 3*nv*nv₂²/2 + nv₂*vv₂)*S̄₁
    #                                    ) + 
    #             Gm₂*r⁻²*((Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁)*dn_dt + 
    #                      (-5*G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv)*ā₂ - 
    #                      (-G*(6*δm)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + nv₂*vS₁)*ā₁ + 
    #                      (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ -3*nv*nv₂²/2 + nv₂*vv₂)*dS̄₁ + 
    #                      (5*G*δm*nS₁*dr_dt*r⁻² - 5*G*δm*dnS₁_dt*r⁻¹ + (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*nS₁ + 
    #                       (3*nv₂² + 2*vv₂)*dnS₁_dt + 2*(v₁S₁ + v₂S₁)*dnv_dt +2*(dv₁S₁_dt + dv₂S₁_dt)*nv)*v̄₂ - 
    #                      (G*(6*δm)*nS₁*dr_dt*r⁻² - G*(6*δm)*dnS₁_dt*r⁻¹ + 
    #                       3*nS₁*nv₂*dnv₂_dt + 3*nv₂²*dnS₁_dt/2 + nv₂*dvS₁_dt + vS₁*dnv₂_dt)*v̄₁ + 
    #                      (-Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*dr_dt*r⁻² + 
    #                       Gm₁*(-16*nS₁*dnv_dt - 16*nv*dnS₁_dt + 
    #                       3*dv₁S₁_dt - 7*dv₂S₁_dt)*r⁻¹ - 
    #                       2*Gm₂*nS₁*nv*dr_dt*r⁻² + 2*Gm₂*nS₁*dnv_dt*r⁻¹ + 
    #                       2*Gm₂*nv*dnS₁_dt*r⁻¹ + (6*nv₂*dnv₂_dt + 
    #                       2*dvv₂_dt)*vS₁ + (3*nv₂² + 2*vv₂)*dvS₁_dt)*n + 
    #                      (-Gm₁*nv₁*dr_dt*r⁻² + Gm₁*dnv₁_dt*r⁻¹ + 
    #                       Gm₂*nv*dr_dt*r⁻² - Gm₂*dnv_dt*r⁻¹ - 
    #                       3*nv*nv₂*dnv₂_dt - 3*nv₂²*dnv_dt/2 + 
    #                       nv₂*dvv₂_dt + vv₂*dnv₂_dt)*S̄₁
    #                     )

    dT1PN_dt₁, dT2PN_dt₁ = let 

        nv   = dot(n, v̄)
        nv₁  = dot(n, v̄₁)
        nv₂  = dot(n, v̄₂)
        vv₂  = dot(v̄, v̄₂)
        nS₁  = dot(n, S̄₁)
        vS₁  = dot(v̄, S̄₁)
        v₁S₁ = dot(v̄₁, S̄₁)
        v₂S₁ = dot(v̄₂, S̄₁)

        dr_dt    = rv*r⁻¹
        dn_dt    = r⁻¹*(v̄ .- (rv*r⁻²*r̄)) 
        dnv_dt   = dot(n,  ā)   + dot(v̄,  dn_dt)
        dnv₁_dt  = dot(n,  ā₁)  + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)  + dot(v̄₂, dn_dt)
        dnS₁_dt  = dot(n,  dS̄₁) + dot(S̄₁, dn_dt)
        dvS₁_dt  = dot(v̄,  dS̄₁) + dot(S̄₁, ā)
        dv₁S₁_dt = dot(v̄₁, dS̄₁) + dot(S̄₁, ā₁)
        dv₂S₁_dt = dot(v̄₂, dS̄₁) + dot(S̄₁, ā₂)
        dvv₂_dt  = dot(v̄,  ā₂)  + dot(v̄₂, ā)

        nv₂² = nv₂^2

        dT1PN = -2*Gm₂*((v̄₁ - 2*v̄₂)*nS₁ + S̄₁*nv - 2*n*vS₁)*v*r⁻²*r⁻¹ + 
                   Gm₂*((v̄₁ - 2*v̄₂)*dnS₁_dt + (ā₁ - 2*ā₂)*nS₁ + 
                   S̄₁*dnv_dt - 2*n*dvS₁_dt + nv*dS̄₁ - 2*vS₁*dn_dt)*r⁻²

        num = (Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁)*n
        num += (-5*G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv)*v̄₂
        num += -(-G*(6*δm)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + nv₂*vS₁)*v̄₁
        num += (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ - 3*nv*nv₂²/2 + nv₂*vv₂)*S̄₁
        num *= -2*Gm₂*dr_dt*r⁻²*r⁻¹
        num += (Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁)*dn_dt
        num += (-5*G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv)*ā₂ 
        num += -(-G*(6*δm)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + nv₂*vS₁)*ā₁ 
        num += (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ -3*nv*nv₂²/2 + nv₂*vv₂)*dS̄₁ 
        num += (5*G*δm*nS₁*dr_dt*r⁻² - 5*G*δm*dnS₁_dt*r⁻¹ + (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*nS₁ + 
                (3*nv₂² + 2*vv₂)*dnS₁_dt + 2*(v₁S₁ + v₂S₁)*dnv_dt +2*(dv₁S₁_dt + dv₂S₁_dt)*nv)*v̄₂ 
        num += -(G*(6*δm)*nS₁*dr_dt*r⁻² - G*(6*δm)*dnS₁_dt*r⁻¹ +
                3*nS₁*nv₂*dnv₂_dt + 3*nv₂²*dnS₁_dt/2 + nv₂*dvS₁_dt + vS₁*dnv₂_dt)*v̄₁
        num += (-Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*dr_dt*r⁻² +
                Gm₁*(-16*nS₁*dnv_dt - 16*nv*dnS₁_dt + 
                    3*dv₁S₁_dt - 7*dv₂S₁_dt)*r⁻¹ - 
                2*Gm₂*nS₁*nv*dr_dt*r⁻² + 2*Gm₂*nS₁*dnv_dt*r⁻¹ + 
                2*Gm₂*nv*dnS₁_dt*r⁻¹ + 
                (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*vS₁ + 
                (3*nv₂² + 2*vv₂)*dvS₁_dt)*n
        num += (-Gm₁*nv₁*dr_dt*r⁻² + Gm₁*dnv₁_dt*r⁻¹ + 
                Gm₂*nv*dr_dt*r⁻² - Gm₂*dnv_dt*r⁻¹ - 
                3*nv*nv₂*dnv₂_dt - 3*nv₂²*dnv_dt/2 + 
                nv₂*dvv₂_dt + vv₂*dnv₂_dt)*S̄₁
        num *= Gm₂*r⁻²
        dT2PN = num *  Gm₂*r⁻²
        
        dT1PN, dT2PN
    end

    dT1PN_dt₂, dT2PN_dt₂ = let n = -n, v̄ = -v̄, r̄ = -r̄, ā = -ā, δm = -δm

        nv   = dot(n, v̄)
        nv₁  = dot(n, v̄₁)
        nv₂  = dot(n, v̄₂)
        vv₁  = dot(v̄, v̄₁)
        nS₂  = dot(n, S̄₂)
        vS₂  = dot(v̄, S̄₂)
        v₂S₂ = dot(v̄₂, S̄₂)
        v₁S₂ = dot(v̄₁, S̄₂)        

        dr_dt    = rv*r⁻¹
        dn_dt    = r⁻¹*(v̄ - (rv*r⁻²*r̄))
        dnv_dt   = dot(n,  ā)     + dot(v̄, dn_dt)
        dnv₁_dt  = dot(n,  ā₁)    + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)    + dot(v̄₂, dn_dt)
        dnS₂_dt  = dot(n,  dS̄₂)   + dot(S̄₂, dn_dt)
        dvS₂_dt  = dot(v̄,  dS̄₂)   + dot(S̄₂, ā)
        dv₂S₂_dt = dot(v̄₂, dS̄₂)   + dot(S̄₂, ā₂)
        dv₁S₂_dt = dot(v̄₁, dS̄₂)   + dot(S̄₂, ā₁)
        dvv₁_dt  = dot(v̄,   ā₁)   + dot(v̄₁, ā)
    
        nv₁² = nv₁^2

        dT1PN = -2*Gm₁*((v̄₂ - 2*v̄₁)*nS₂ + S̄₂*nv - 2*n*vS₂)*v*r⁻²*r⁻¹ + 
                   Gm₁*((v̄₂ - 2*v̄₁)*dnS₂_dt + (ā₂ - 2*ā₁)*nS₂ + 
                   S̄₂*dnv_dt - 2*n*dvS₂_dt + nv*dS̄₂ - 2*vS₂*dn_dt)*r⁻²

        num =  (Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*r⁻¹ + 2*Gm₁*nS₂*nv*r⁻¹ + (3*nv₁² + 2*vv₁)*vS₂)*n
        num += (-5*G*δm*nS₂*r⁻¹ + (3*nv₁² + 2*vv₁)*nS₂ + 2*(v₂S₂ + v₁S₂)*nv)*v̄₁
        num += -(-G*(6*δm)*nS₂*r⁻¹ + 3*nS₂*nv₁²/2 + nv₁*vS₂)*v̄₂
        num += (Gm₂*nv₂*r⁻¹ - Gm₁*nv*r⁻¹ - 3*nv*nv₁²/2 + nv₁*vv₁)*S̄₂
        num *= -2*Gm₁*dr_dt*r⁻²*r⁻¹
        num += (Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*r⁻¹ + 2*Gm₁*nS₂*nv*r⁻¹ + (3*nv₁² + 2*vv₁)*vS₂)*dn_dt
        num += (-5*G*δm*nS₂*r⁻¹ + (3*nv₁² + 2*vv₁)*nS₂ + 2*(v₂S₂ + v₁S₂)*nv)*ā₁ 
        num += -(-G*(6*δm)*nS₂*r⁻¹ + 3*nS₂*nv₁²/2 + nv₁*vS₂)*ā₂ 
        num += (Gm₂*nv₂*r⁻¹ - Gm₁*nv*r⁻¹ -3*nv*nv₁²/2 + nv₁*vv₁)*dS̄₂ 
        num += (5*G*δm*nS₂*dr_dt*r⁻² - 5*G*δm*dnS₂_dt*r⁻¹ + (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*nS₂ + 
                (3*nv₁² + 2*vv₁)*dnS₂_dt + 2*(v₂S₂ + v₁S₂)*dnv_dt +2*(dv₂S₂_dt + dv₁S₂_dt)*nv)*v̄₁ 
        num += -(G*(6*δm)*nS₂*dr_dt*r⁻² - G*(6*δm)*dnS₂_dt*r⁻¹ +
                3*nS₂*nv₁*dnv₁_dt + 3*nv₁²*dnS₂_dt/2 + nv₁*dvS₂_dt + vS₂*dnv₁_dt)*v̄₂
        num += (-Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*dr_dt*r⁻² +
                Gm₂*(-16*nS₂*dnv_dt - 16*nv*dnS₂_dt + 
                    3*dv₂S₂_dt - 7*dv₁S₂_dt)*r⁻¹ - 
                2*Gm₁*nS₂*nv*dr_dt*r⁻² + 2*Gm₁*nS₂*dnv_dt*r⁻¹ + 
                2*Gm₁*nv*dnS₂_dt*r⁻¹ + 
                (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*vS₂ + 
                (3*nv₁² + 2*vv₁)*dvS₂_dt)*n
        num += (-Gm₂*nv₂*dr_dt*r⁻² + Gm₂*dnv₂_dt*r⁻¹ + 
                Gm₁*nv*dr_dt*r⁻² - Gm₁*dnv_dt*r⁻¹ - 
                3*nv*nv₁*dnv₁_dt - 3*nv₁²*dnv_dt/2 + 
                nv₁*dvv₁_dt + vv₁*dnv₁_dt)*S̄₂
        num *= Gm₁*r⁻²
        dT2PN = num *  Gm₁*r⁻²

        dT1PN, dT2PN
    end

    
    dvi .+= dT1PN_dt₁*c⁻² + dT2PN_dt₁*c⁻⁴ 
    dvj .+= dT1PN_dt₂*c⁻² + dT2PN_dt₂*c⁻⁴ 
    nothing
end

function spin_precession_COM!(dvi,
                              dvj,
                              dvs,
                              rs,
                              vs,
                              pair::Tuple{Int, Int},
                              params::SimulationParams)

    i, j = pair

    ā₁  = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂ = @SVector [rs[1,  j], rs[2,  j], rs[3,  j]]
    v̄₂ = @SVector [vs[1,  j], vs[2,  j], vs[3,  j]]
    
    # S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    # dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]
    
    v₁  = norm(v̄₁)
    # v₂ = norm(v̄₂)
    # v₁² = v₁^2
    
    m₁ = params.M[i]
    m₂ = params.M[j]
    
    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂
    # rxv = r̄ × v̄

    rv = dot(r̄, v̄)
    # av = dot(ā, v̄)

    r = norm(r̄) # r₁₂
    v = norm(v̄) # v₁₂
    v² = v*v
    v³ = v²*v

    n = r̄/r

    M = m₁ + m₂
    X1 = m₁/M
    X2 = m₂/M

    GM = G*M

    Δ = X1 - X2
    δm = m₁ - m₂
    ν = X1*X2
    ν² = ν*ν
    ν³ = ν²*ν

    n̄ = r̄/r

    nxv = n̄ × v̄
    nv = dot(n̄, v̄)
    nxv_norm = norm(nxv)

    # G² G³
    dnxv_dt = 1/r*(v̄ - rv/r^2*r̄)

    # dΩ_dt =  -4*G³*M^3*(ν^3/2 - 9*ν^2/8 - 9*ν/4 + (7/16) + dm*(-ν^2/8 - ν/8 + (-7/16))/M)*dr_dt/r^5 - 
    #          3*G²*M^2*((-9*ν^3/8 + 75*ν^2/32 + 27*ν/4 + (3/16) + dm*(35*ν^2/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
    #                    (13*ν^3/4 - 159*ν^2/16 - 525*ν/32 + (1/4) + dm*(-87*ν^2/16 - 75*ν/32 + (-1/4))/M)*nv^2)*dr_dt/r^4 + 
    #          G²*M^2*(2*(-9*ν^3/8 + 75*ν^2/32 + 27*ν/4 + (3/16) + dm*(35*ν^2/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 2*(13*ν^3/4 - 
    #                     159*ν^2/16 - 525*ν/32 + (1/4) + dm*(-87*ν^2/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3 - 
    #          2*GM*((-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv^2*v^2 + 
    #                (15*ν^3/16 - 195*ν^2/32 + 15*ν/8 + dm*(-75*ν^2/32 + 15*ν/8)/M)*nv^4 +
    #                (17*ν^3/16 - 31*ν^2/8 + 19*ν/16 + (1/32) + dm*(-11*ν^2/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3 + 
    #          GM*(2*(-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv^2*v*dv_dt + 
    #              2*(-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv*v^2*dnv_dt + 
    #              4*(15*ν^3/16 - 195*ν^2/32 + 15*ν/8 + dm*(-75*ν^2/32 + 15*ν/8)/M)*nv^3*dnv_dt + 
    #              4*(17*ν^3/16 - 31*ν^2/8 + 19*ν/16 + (1/32) + dm*(-11*ν^2/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt)/r^2 - 
    #           2*GM*(ν/2 + (3/4) - 3*dm/(4*M))*dr_dt/(c^2*r^3) + (-3*G²*M^2*(ν^2/2 - 3*ν/8 + (-1/4) + dm*((1/4) - ν/8)/M)*dr_dt/r^4 + 
    #           2*GM*(-3*ν/4 - 3*dm*ν/(2*M))*nv*dnv_dt/r^2 - 2*GM*(ν^2*(-3*ν^2/8 + 11*ν/8 + (1/16) + dm*(ν/2 + (-1/16))/M) + 
    #           (-3*ν/4 - 3*dm*ν/(2*M))*nv^2)*dr_dt/r^3)/c^4

    dΩ₁_dt = let δm = δm
        num   = -2*GM*(ν/2 + (3/4) - 3*δm/(4*M))*dr_dt/(c^2*r^3)
        numm  = -3*G^2*M^2*(ν²/2 - 3*ν/8 + (-1/4) + δm*((1/4) - ν/8)/M)*dr_dt/r^4 
        numm += 2*GM*(-3*ν/4 - 3*δm*ν/(2*M))*nv*dnv_dt/r^2 
        numm -= 2*GM*(ν²*(-3*ν²/8 + 11*ν/8 + (1/16) + δm*(ν/2 + (-1/16))/M) + (-3*ν/4 - 3*δm*ν/(2*M))*nv^2)*dr_dt/r^3
        num  += c⁻⁴*numm

        numm  = -4*G^3*M^3*(ν³/2 - 9*ν²/8 - 9*ν/4 + (7/16) + δm*(-ν²/8 - ν/8 + (-7/16))/M)*dr_dt/r^5
        numm -= 3*G^2*M^2*((-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
                (13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + 
                (-1/4))/M)*nv^2)*dr_dt/r^4
        numm += G^2*M^2*(2*(-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 
                2*(13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3
        numm -= 2*GM*((-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v^2 + 
                    (15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^4 + 
                    (17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3

        nummm  = 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v*dv_dt
        nummm += 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv*v^2*dnv_dt
        nummm += 4*(15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^3*dnv_dt
        nummm += 4*(17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt

        numm += nummm/r^2 
        num  += numm*c⁻⁶
        num
    end

    dΩ₂_dt = let δm = -δm
        num   = -2*GM*(ν/2 + (3/4) - 3*δm/(4*M))*dr_dt/(c^2*r^3)
        numm  = -3*G^2*M^2*(ν²/2 - 3*ν/8 + (-1/4) + δm*((1/4) - ν/8)/M)*dr_dt/r^4 
        numm += 2*GM*(-3*ν/4 - 3*δm*ν/(2*M))*nv*dnv_dt/r^2 
        numm -= 2*GM*(ν²*(-3*ν²/8 + 11*ν/8 + (1/16) + δm*(ν/2 + (-1/16))/M) + (-3*ν/4 - 3*δm*ν/(2*M))*nv^2)*dr_dt/r^3
        num  += c⁻⁴*numm

        numm  = -4*G^3*M^3*(ν³/2 - 9*ν²/8 - 9*ν/4 + (7/16) + δm*(-ν²/8 - ν/8 + (-7/16))/M)*dr_dt/r^5
        numm -= 3*G^2*M^2*((-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
                (13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + 
                (-1/4))/M)*nv^2)*dr_dt/r^4
        numm += G^2*M^2*(2*(-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 
                2*(13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3
        numm -= 2*GM*((-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v^2 + 
                    (15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^4 + 
                    (17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3

        nummm  = 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v*dv_dt
        nummm += 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv*v^2*dnv_dt
        nummm += 4*(15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^3*dnv_dt
        nummm += 4*(17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt

        numm += nummm/r^2 
        num  += numm*c⁻⁶
        num
    end

    ai = nxv*dΩ₁_dt + Ω*dnxv_dt
    aj = nxv*dΩ₂_dt + Ω*dnxv_dt

    dvi .+= ai
    dvj .+= aj
    nothing
end