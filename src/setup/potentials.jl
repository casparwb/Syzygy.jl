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

# struct UnitfulSimulationParams{LengthVecType, MassVecType, WattVecType, IntVecType, TimeVecType} <: SimulationParams
#     R::LengthVecType # radii
#     M::MassVecType # masses
#     L::WattVecType # luminosities
#     stellar_types::IntVecType 
#     M_cores::MassVecType # core masses
#     R_cores::LengthVecType # core radii
#     ages::TimeVecType 
# end

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
        @warn "Tidal factor for n ≠ {4, 10} is not defined."
        f = x -> x
    end

    DynamicalTidalPotential(n, γ, f)
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
            
            logg = log10(surface_gravity_unit_conversion_factor*(G*m₁/R^2)) # convert from R⊙/yr² to cm/s²
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
            
            logg = log10(surface_gravity_unit_conversion_factor*(G*m₂/R^2)) # convert from R⊙/yr² to cm/s²
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

    a = semi_major_axis(r, v^2, m₂+m₁, G)

    # k₁ = potential.apsidal_motion_constant[i]
    # k₂ = potential.apsidal_motion_constant[j]

    # tidal force on 1 by 2
    a₁ = let
        if !(stellar_type_1 < 10)
                SA[0.0, 0.0, 0.0]
            else    
                μ = G*m₂*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[i]
                envelope_radius = potential.R_env[i]
                luminosity = params.L[i]
                R = params.R[i]
                Ω = potential.rotational_angular_velocity[i]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, envelope_mass, envelope_radius,
                                                                   stellar_type_1, luminosity, 
                                                                   m₂, a)
                k = potential.apsidal_motion_constant[i]
                kτ = R^3/(G*m₁)*k_T
    
                @. -3μ*m₂/m₁*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    # # tidal force on 2 by 1
    a₂ = let 
        if !(stellar_type_2 < 10)
                SA[0.0, 0.0, 0.0]
            else    
                μ = G*m₁*r⁻¹*r⁻¹
    
                envelope_mass = potential.M_env[j]
                envelope_radius = potential.R_env[j]
                luminosity  = params.L[j]
                R = params.R[j]
                Ω = potential.rotational_angular_velocity[j]
                k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, envelope_mass, envelope_radius,
                                                                   stellar_type_1, luminosity, 
                                                                   m₁, a)
                k = potential.apsidal_motion_constant[j]
                kτ = R^3/(G*m₂)*k_T
    
                @. -3μ*m₁/m₂*(R*r⁻¹)^5*((k + 3v̄*r⁻¹*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
            end

        end

    # a₁ = let k = i
    #     stellar_type = params.stellar_types[k]

    #     if !(stellar_types[stellar_type] isa Star)
    #         SA[0.0, 0.0, 0.0]
    #     else
    #         S̄₁  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

    #         envelope_mass = potential.M_env[k]
    #         envelope_radius = potential.R_env[k]
    #         luminosity  = params.L[k]

    #         R = Rs[k]

    #         Ω = norm(S̄₁)
            
    #         logg = log10(6.985766564066957e-5*(G*m₁/R^2)) # convert from R⊙/yr² to cm/s²
    #         logm = log10(m₁)
    #         # k = asidal_motion_constant_interpolated(logm, logg)
    #         # k = potential.k_interpolator(logm, logg)
    #         k = potential.apsidal_motion_constant[]

    #         μ = G*m₂/r²
            
    #         k_T = apsidal_motion_constant_over_tidal_timescale(m₁, R, envelope_mass, envelope_radius,
    #                                                            stellar_type, luminosity, 
    #                                                            m₂, a)

    #         kτ = R^3/(G*m₁)*k_T

    #         @. -3μ*m₂/m₁*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
    #     end
    # end

    # a₂ = let k = j
    #     stellar_type = params.stellar_types[k]

    #     if !(stellar_types[stellar_type] isa Star)
    #         SA[0.0, 0.0, 0.0]
    #     else
    #         S̄₂  = @SVector [rs[4, k], rs[5, k], rs[6, k]]

    #         envelope_mass = potential.M_env[k]
    #         envelope_radius = potential.R_env[k]
    #         luminosity  = params.L[k]

    #         R = Rs[k]

    #         Ω = norm(S̄₂)
            
    #         logg = log10(6.985766564066957e-5*(G*m₂/R^2)) # convert from R⊙/yr² to cm/s²
    #         logm = log10(m₂)
    #         # k = asidal_motion_constant_interpolated(logm, logg)
    #         # k = potential.k_interpolator(logm, logg)
    #         k = potential.apsidal_motion_constant[]

    #         μ = G*m₁/r²
            
    #         k_T = apsidal_motion_constant_over_tidal_timescale(m₂, R, envelope_mass, envelope_radius,
    #                                                            stellar_type, luminosity, 
    #                                                            m₁, a)

    #         kτ = R^3/(G*m₂)*k_T

    #         @. -3μ*m₁/m₂*(R/r)^5*((k + 3v̄/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
    #     end
    # end
    # println(i, " ", a₁)
    # println(j, " ", a₂)
    
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

function PN1p5_spin_acceleration!(dvi, 
                                 dvj,
                                 rs,
                                 vs,
                                 pair::Tuple{Int, Int},
                                 params::SimulationParams)
                           
    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    
    m₁ = params.M[i]
    m₂ = params.M[j]
        
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r
    r⁻³ = r⁻¹/r^2   

    n = r̄*r⁻¹
    nv = dot(n, v̄)
    nxv = n × v̄

    # @show (n × S̄₂) (v̄ × S̄₁)

    G_r = G*r⁻¹
    G_r³ = G_r*r⁻³

    # PN-1.5 acceleration from spin contribution
    ai = G_r³*m₂*( (6*dot(S̄₁, nxv)/m₁+ 6*dot(S̄₂, nxv)/m₂)*n + 3nv*(n × S̄₁)/m₁ + 
                          6nv*(n × S̄₂)/m₂ - 3*(v̄ × S̄₁)/m₁ - 4*(v̄ × S̄₂)/m₂
                         )

    aj = G_r³*m₁*( (6*dot(S̄₂, nxv)/m₂ + 6*dot(S̄₁, nxv)/m₁)*n + 3nv*(n × S̄₂)/m₂ + 
                          6nv*(n × S̄₁)/m₁ - 3*(v̄ × S̄₂)/m₂ - 4*(v̄ × S̄₁)/m₁
                         )
    
    dvi .+= ai*c⁻³
    dvj .+= aj*c⁻³
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

function PN2_spin_acceleration!(dvi, 
                                dvj,
                                rs,
                                vs,
                                pair::Tuple{Int, Int},
                                params::SimulationParams)
                           
    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    # v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    # v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    
    m₁ = params.M[i]
    m₂ = params.M[j]
        
    r̄ = r̄₁ - r̄₂
    # v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r

    n = r̄*r⁻¹
    S₁S₂ = dot(S̄₁, S̄₂)
    S₂S₁ = S₁S₂
    nS₁  = dot(n, S̄₁)
    nS₂  = dot(n, S̄₂)

    G_r = G*r⁻¹
    G_r⁴ = G_r*r⁻¹*r⁻¹*r⁻¹

    # PN-2 spin-spin interaction 
    ai = -G_r⁴*3/m₁*( n*S₁S₂ + S̄₁*nS₁ + S̄₂*nS₁ - 5n*nS₁*nS₂ )

    aj = -G_r⁴*3/m₂*( n*S₂S₁ + S̄₂*nS₂ + S̄₁*nS₂ - 5n*nS₂*nS₁ )
    
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

function PN2p5_spin_acceleration!(dvi, 
                                  dvj,
                                  rs,
                                  vs,
                                  pair::Tuple{Int, Int},
                                  params::SimulationParams)
                               
    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    
    m₁ = params.M[i]
    m₂ = params.M[j]
        
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    v₁xv₂ = v̄₁ × v̄₂
    v₂xv₁ = v̄₂ × v̄₁

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r

    n = r̄*r⁻¹
    nv = dot(n, v̄)
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    vv₁ = dot(v̄, v̄₁)
    vv₂ = dot(v̄, v̄₂)

    nS₁ = dot(n, S̄₁)
    nS₂ = dot(n, S̄₂)

    nxv  = n × v̄
    nxv₁ = n × v̄₁
    nxv₂ = n × v̄₂

    nxS₁ = n × S̄₁
    nxS₂ = n × S̄₂

    vxS₁ = v̄ × S̄₁
    vxS₂ = v̄ × S̄₂

    v₁S₁ = dot(v̄₁, S̄₁)
    v₂S₂ = dot(v̄₂, S̄₂)

    nv₁v₂ = dot(n, v₁xv₂)
    nv₂v₁ = dot(n, v₂xv₁)
    S₁nv  = dot(S̄₁, nxv)
    S₂nv  = dot(S̄₂, nxv)
    S₁nv₁ = dot(S̄₁, nxv₁)
    S₁nv₂ = dot(S̄₁, nxv₂)
    S₂nv₁ = dot(S̄₂, nxv₁)
    S₂nv₂ = dot(S̄₂, nxv₂)

    G_r = G*r⁻¹
    G_r³ = G_r*r⁻¹*r⁻¹

    Gm₁ = G*m₁
    Gm₂ = G*m₂ 

    ai = n*(-6*nv₁v₂*(v₁S₁/m₁ + v₂S₂/m₂) - 
            S₁nv/m₁*(15nv₂^2 + 6vv₂ + 26Gm₁/r + 18Gm₂/r)  - 
            S₂nv/m₂*(15nv₂^2 + 6vv₂ + 49/2*Gm₁/r + 20Gm₂/r))
    ai += v̄₁*(-3*S₁nv₁/m₁*(nv₁ + nv₂) + 6nv₁*S₁nv₂/m₁ - 3*dot(S̄₁, v₁xv₂)/m₁  - 
              6nv₁*S₂nv₁/m₂ + S₂nv₂/m₂*(12nv₁ - 6nv₂)  - 4*dot(S̄₂, v₁xv₂)/m₂)
    ai += v̄₂*(6nv₁*S₁nv/m₁ + 6nv₁*S₂nv/m₂)
    ai -= nxv₁*(3nv*v₁S₁/m₁ + 4Gm₁/r*nS₂/m₂)  
    ai -= nxv₂*(6nv*v₂S₂/m₂ - 4Gm₁/r*nS₂/m₂) 
    ai += v₁xv₂*(3v₁S₁/m₁ + 4v₂S₂/m₂) 
    ai += nxS₁/m₁*(-15/2*nv*nv₂^2 + 3nv₂*vv₂  - 
                       14Gm₁/r*nv - 9Gm₂/r*nv) 
    ai += nxS₂/m₂*(-15nv*nv₂^2 - 6nv₁*vv₂ + 12nv₂*vv₂ + Gm₁/r*(-35/2*nv₁ + 39/2*nv₂) - 16Gm₂/r*nv)
    ai += vxS₁/m₁*(-3nv₁*nv₂ + 15/2*nv₂^2 + G/r*(14m₁ + 9m₂) + 3vv₂)
    ai += vxS₂/m₂*(6nv₂^2 + 4vv₂ + 23/2*Gm₁/r + 12Gm₂/r)
    ai *= G_r³*m₂

    aj = let n = -n

        nv₁ = -nv₁
        nv₂ = -nv₂
        vv₁ = -vv₁
        nS₁ = -nS₁
    
        nxv₁ = -nxv₁
        nxv₂ = -nxv₂
    
        nxS₁ = -nxS₁
        nxS₂ = -nxS₂
    
        vxS₁ = -vxS₁
        vxS₂ = -vxS₂
    
        nv₂v₁ = -nv₂v₁
        S₁nv₁ = -S₁nv₁
        S₁nv₂ = -S₁nv₂
        S₂nv₁ = -S₂nv₁
        S₂nv₂ = -S₂nv₂

        a = n*(-6*nv₂v₁*(v₂S₂/m₂ + v₁S₁/m₁) - 
                S₂nv/m₂*(15nv₁^2 + 6vv₁ + 26Gm₂/r + 18Gm₁/r)  - 
                S₁nv/m₁*(15nv₁^2 + 6vv₁ + 49/2*Gm₂/r + 20Gm₁/r))
        a += v̄₂*(-3*S₂nv₂/m₂*(nv₂ + nv₁) + 6nv₂*S₂nv₁/m₂ - 3*dot(S̄₂, v₂xv₁)/m₂  - 
                6nv₂*S₁nv₂/m₁ + S₁nv₁/m₁*(12nv₂ - 6nv₁)  - 4*dot(S̄₁, v₂xv₁)/m₁)
        a += v̄₁*(6nv₂*S₂nv/m₂ + 6nv₂*S₁nv/m₁) 
        a -= nxv₂*(3nv*v₂S₂/m₂ + 4Gm₂/r*nS₁/m₁) 
        a -= nxv₁*(6nv*v₁S₁/m₁  - 4Gm₂/r*nS₁/m₁) 
        a += v₂xv₁*(3v₂S₂/m₂ + 4v₁S₁/m₁) 
        a += nxS₂/m₂*(-15/2*nv*nv₁^2 + 3nv₁*vv₁  - 
                        14Gm₂/r*nv - 9Gm₁/r*nv) 
        a += nxS₁/m₁*(-15nv*nv₁^2 - 6nv₂*vv₁ + 12nv₁*vv₁ + Gm₂/r*(-35/2*nv₂ + 39/2*nv₁) - 16Gm₁/r*nv) 
        a += vxS₂/m₂*(-3nv₂*nv₁ + 15/2*nv₁^2 + G/r*(14m₂ + 9m₁) + 3vv₁) 
        a += vxS₁/m₁*(6nv₁^2 + 4vv₁ + 23/2*Gm₂/r + 12Gm₁/r)
        a *= G_r³*m₁

        a
    end
    
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

    G_r = G*r⁻¹
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

function PN1_spin_precession!(dvi,
                              dvj,
                              dvs,
                              rs,
                              vs,
                              pair::Tuple{Int, Int},
                              params::SimulationParams)
    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

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
    

    m₁ = params.M[i]
    m₂ = params.M[j]
    
    # add @fastmath?

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r² = r*r
    r⁻¹ = 1/r
    r⁻² = 1/r²

    n    = r̄/r
    nv   = dot(n, v̄)

    dr_dt = nv
    dn_dt = (r*v̄ - r̄*dr_dt)*r⁻²
    
    Gm₁ = G*m₁
    Gm₂ = G*m₂

    dT1PN_dt₁ = let 

        nS₁  = dot(n, S̄₁)
        vS₁  = dot(v̄, S̄₁)

        dnS₁_dt  = dot(n,  dS̄₁) + dot(S̄₁, dn_dt)
        dvS₁_dt  = dot(v̄,  dS̄₁) + dot(S̄₁, ā)
        dnv_dt   = dot(n,  ā)   + dot(v̄,  dn_dt)

        Gm₂r⁻² = Gm₂*r⁻²
        dGm₂r⁻²_dt = -2*Gm₂*dr_dt*r⁻¹*r⁻²
        fac = (v̄₁ - 2*v̄₂)*nS₁ + S̄₁*nv - 2*n*vS₁
        dfac_dt = (v̄₁ - 2*v̄₂)*dnS₁_dt + (ā₁ - 2*ā₂)*nS₁ + 
                  S̄₁*dnv_dt + nv*dS̄₁ - 
                  2*n*dvS₁_dt - 2*vS₁*dn_dt
        fac*dGm₂r⁻²_dt + Gm₂r⁻²*dfac_dt
    end

    dT1PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā,

        dn_dt = -dn_dt
        nS₂  = dot(n, S̄₂)
        vS₂  = dot(v̄, S̄₂)

        dnv_dt   = dot(n,  ā)     + dot(v̄, dn_dt)
        dnS₂_dt  = dot(n,  dS̄₂)   + dot(S̄₂, dn_dt)
        dvS₂_dt  = dot(v̄,  dS̄₂)   + dot(S̄₂, ā)


        Gm₁r⁻² = Gm₁*r⁻²
        dGm₁r⁻²_dt = -2*Gm₁*dr_dt*r⁻¹*r⁻²
        fac = (v̄₂ - 2*v̄₁)*nS₂ + S̄₂*nv - 2*n*vS₂
        dfac_dt = (v̄₂ - 2*v̄₁)*dnS₂_dt + (ā₂ - 2*ā₁)*nS₂ + 
                   S̄₂*dnv_dt + nv*dS̄₂ - 
                   2*n*dvS₂_dt - 2*vS₂*dn_dt
        fac*dGm₁r⁻²_dt + Gm₁r⁻²*dfac_dt
    end
    
    dvi .+= dT1PN_dt₁*c⁻² 
    dvj .+= dT1PN_dt₂*c⁻² 
    nothing
end

function PN1p5_spin_precession!(dvi,
                              dvj,
                              dvs,
                              rs,
                              vs,
                              pair::Tuple{Int, Int},
                              params::SimulationParams)
    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    r̄₂ = @SVector [rs[1,  j], rs[2,  j], rs[3,  j]]
    v̄₂ = @SVector [vs[1,  j], vs[2,  j], vs[3,  j]]

    S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂  = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]
    
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1/r
    r⁻² = r⁻¹*r⁻¹
    r⁻³ = r⁻²*r⁻¹

    n    = r̄/r
    nv   = dot(n, v̄)

    dr_dt = nv
    dn_dt = (r*v̄ - r̄*dr_dt)*r⁻²
    
    Gr⁻³ = -G*r⁻³
    dGr⁻³_dt = 3*G*dr_dt*r⁻²*r⁻²

    dT1p5PN_dt₁ = let 
        nS₂ = dot(n, S̄₂)
        dnS₂_dt = dot(n,  dS̄₂)   + dot(S̄₂, dn_dt)

        # F1p5PN = S̄₂ - 3nS₂*n
        # dF1p5PN_dt = dS̄₂ - 3*n*dnS₂_dt - 3nS₂*dn_dt
        # dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₁) + (F1p5PN × dS̄₁))

        F = (3nS₂*n × S̄₁) - (S̄₂ × S̄₁)
        dF_dt = 3*((nS₂*dn_dt + n*dnS₂_dt) × S̄₁ + nS₂*n × dS̄₁) - (S̄₁ × dS̄₂) + (S̄₂ × dS̄₁)

        Gr⁻³*dF_dt + F*dGr⁻³_dt
    end

    dT1p5PN_dt₂ = let n = -n, dn_dt = -dn_dt
        nS₁ = dot(n, S̄₁)
        dnS₁_dt  = dot(n,  dS̄₁) + dot(S̄₁, dn_dt)

        # F1p5PN = S̄₁ - 3nS₁*n
        # dF1p5PN_dt = dS̄₁ - 3*n*dnS₁_dt - 3nS₁*dn_dt
        # dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₂) + (F1p5PN × dS̄₂))
        # ₁ ₂ ₃

        F = (3nS₁*n × S̄₂) - (S̄₁ × S̄₂)
        dF_dt = 3*((nS₁*dn_dt + n*dnS₁_dt) × S̄₂ + nS₁*n × dS̄₂) - (S̄₂ × dS̄₁) + (S̄₁ × dS̄₂)

        Gr⁻³*dF_dt + F*dGr⁻³_dt
    end
    
    dvi .+= dT1p5PN_dt₁*c⁻³
    dvj .+= dT1p5PN_dt₂*c⁻³
    nothing
end

function PN2_spin_precession!(dvi,
                              dvj,
                              dvs,
                              rs,
                              vs,
                              pair::Tuple{Int, Int},
                              params::SimulationParams)

    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

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
    
    m₁ = params.M[i]
    m₂ = params.M[j]
    δm = m₁ - m₂    

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂
    r = norm(r̄) # r₁₂

    r² = r*r
    r⁻¹ = 1/r
    r⁻² = 1/r²

    n    = r̄/r
    nv   = dot(n, v̄)
    nv₁  = dot(n, v̄₁) 
    nv₂  = dot(n, v̄₂)

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    dr_dt = nv
    dn_dt = (r*v̄ - r̄*dr_dt)*r⁻²
    
    Gm₁ = G*m₁
    Gm₂ = G*m₂
    
    Gm₁r⁻² = Gm₁*r⁻²
    dGm₁r⁻²_dt = -2*Gm₁*dr_dt*r⁻¹*r⁻²

    Gm₂r⁻² = Gm₂*r⁻²
    dGm₂r⁻²_dt = -2*Gm₂*dr_dt*r⁻¹*r⁻²

    dT2PN_dt₁ = let 

        vv₂  = dot(v̄, v̄₂)
        nS₁  = dot(n, S̄₁)
        vS₁  = dot(v̄, S̄₁)
        v₁S₁ = dot(v̄₁, S̄₁)
        v₂S₁ = dot(v̄₂, S̄₁)


        dnS₁_dt  = dot(n,  dS̄₁) + dot(S̄₁, dn_dt)
        dvS₁_dt  = dot(v̄,  dS̄₁) + dot(S̄₁, ā)
        dv₁S₁_dt = dot(v̄₁, dS̄₁) + dot(S̄₁, ā₁)
        dv₂S₁_dt = dot(v̄₂, dS̄₁) + dot(S̄₁, ā₂)
        dvv₂_dt  = dot(v̄,  ā₂)  + dot(v̄₂, ā)
        dnv_dt   = dot(n,  ā)   + dot(v̄,  dn_dt)
        dnv₁_dt  = dot(n,  ā₁)  + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)  + dot(v̄₂, dn_dt)


        # T2PN = m₂*r⁻²*(S̄₁*(nv₂*vv₂ - 3/2*nv₂²*nv + Gm₁*r⁻¹*nv₁ - Gm₂*r⁻¹*nv) + 
        #                 n*(vS₁*(3*nv₂² + 2*vv₂) + Gm₁*r⁻¹*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁) +
        #                     2nS₁*Gm₂*r⁻¹*nv) - v̄₁*(3/2*nS₁*nv₂² + vS₁*nv₂ -
        #                     nS₁*G*r⁻¹*(6m₁ - m₂)) + v̄₂*(nS₁*(2vv₂ + 3nv₂²) +
        #                     2nv*(v₁S₁ + v₂S₁) - 5nS₁*G*r⁻¹*δm)
        #                 )

        # dF2PN_dt = dn_dt*(Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁) + 
        #            ā₂*(-5*G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv) - 
        #            (-G*(6*m₁ - m₂)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + v̄₂*vS₁)*ā₁ + 
        #            (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ - 3*nv*nv₂²/2 + nv₂*vv₂)*dS₁_dt + 
        #            (5*G*δm*nS₁*dr_dt*r⁻² - 5*G*δm*dnS₁_dt*r⁻¹ + 
        #             (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*nS₁ + (3*nv₂² + 2*vv₂)*dnS₁_dt + 
        #             2*(v₁S₁ + v₂S₁)*dnv_dt + 2*(dv₁S₁_dt + dv₂S₁_dt)*nv)*v̄₂ - 
        #            (G*(6*m₁ - m₂)*nS₁*dr_dt*r⁻² - G*(6*m₁ - m₂)*dnS₁_dt*r⁻¹ + 
        #             3*nS₁*nv₂*dnv₂_dt + 3*nv₂²*dnS₁_dt/2 + nv₂*dvS₁_dt + vS₁*dnv₂_dt)*v̄₁ + 
        #             (-Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*dr_dt*r⁻² + Gm₁*(-16*nS₁*dnv_dt - 
        #              16*nv*dnS₁_dt + 3*dv₁S₁_dt - 7*dv₂S₁_dt)*r⁻¹ - 2*Gm₂*nS₁*nv*dr_dt*r⁻² + 
        #              2*Gm₂*nS₁*dnv_dt*r⁻¹ + 2*Gm₂*nv*dnS₁_dt*r⁻¹ + (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*vS₁ + 
        #             (3*nv₂² + 2*vv₂)*dvS₁_dt)*n + 
        #             (-Gm₁*nv₁*dr_dt*r⁻² + Gm₁*dnv₁_dt*r⁻¹ + Gm₂*nv*dr_dt*r⁻² - 
        #              Gm₂*dnv_dt*r⁻¹ - 3*nv*nv₂*dnv₂_dt - 3*nv₂²*dnv_dt/2 + nv₂*dvv₂_dt + vv₂*dnv₂_dt)*S₁
        
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
        
        # Gm₁r⁻²*dF2PN_dt + dGm₁r⁻²_dt*T2PN
        num*Gm₂*r⁻²
    end

    dT2PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā, δm = -δm

        dn_dt = -dn_dt
        nv₁  = -nv₁
        nv₂  = -nv₂
        vv₁  = dot(v̄, v̄₁)
        nS₂  = dot(n, S̄₂)
        vS₂  = dot(v̄, S̄₂)
        v₂S₂ = dot(v̄₂, S̄₂)
        v₁S₂ = dot(v̄₁, S̄₂)        

        dnv_dt   = dot(n,  ā)     + dot(v̄, dn_dt)
        dnv₁_dt  = dot(n,  ā₁)    + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)    + dot(v̄₂, dn_dt)
        dnS₂_dt  = dot(n,  dS̄₂)   + dot(S̄₂, dn_dt)
        dvS₂_dt  = dot(v̄,  dS̄₂)   + dot(S̄₂, ā)
        dv₂S₂_dt = dot(v̄₂, dS̄₂)   + dot(S̄₂, ā₂)
        dv₁S₂_dt = dot(v̄₁, dS̄₂)   + dot(S̄₂, ā₁)
        dvv₁_dt  = dot(v̄,   ā₁)   + dot(v̄₁, ā)

        # T2PN = m₁*r⁻²*(S̄₂*(nv₁*vv₁ - 3/2*nv₁²*nv + Gm₂*r⁻¹*nv₂ - Gm₁*r⁻¹*nv) + 
        #                 n*(vS₂*(3*nv₁² + 2*vv₁) + Gm₂*r⁻¹*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂) +
        #                     2nS₂*Gm₁*r⁻¹*nv) - v̄₂*(3/2*nS₂*nv₁² + vS₂*nv₁ -
        #                     nS₂*G*r⁻¹*(6m₂ - m₁)) + v̄₁*(nS₂*(2vv₁ + 3nv₁²) +
        #                     2nv*(v₂S₂ + v₁S₂) - 5nS₂*G*r⁻¹*δm)
        #                 )

        # dF2PN_dt = dn_dt*(Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*r⁻¹ + 2*Gm₁*nS₂*nv*r⁻¹ + (3*nv₁² + 2*vv₁)*vS₂) + 
        #            ā₁*(-5*G*δm*nS₂*r⁻¹ + (3*nv₁² + 2*vv₁)*nS₂ + 2*(v₂S₂ + v₁S₂)*nv) - 
        #            (-G*(6*m₂ - m₁)*nS₂*r⁻¹ + 3*nS₂*nv₁²/2 + v̄₁*vS₂)*ā₂ + 
        #            (Gm₂*nv₂*r⁻¹ - Gm₁*nv*r⁻¹ - 3*nv*nv₁²/2 + nv₁*vv₁)*dS₂_dt + 
        #            (5*G*δm*nS₂*dr_dt*r⁻² - 5*G*δm*dnS₂_dt*r⁻¹ + 
        #             (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*nS₂ + (3*nv₁² + 2*vv₁)*dnS₂_dt + 
        #             2*(v₂S₂ + v₁S₂)*dnv_dt + 2*(dv₂S₂_dt + dv₁S₂_dt)*nv)*v̄₁ - 
        #            (G*(6*m₂ - m₁)*nS₂*dr_dt*r⁻² - G*(6*m₂ - m₁)*dnS₂_dt*r⁻¹ + 
        #             3*nS₂*nv₁*dnv₁_dt + 3*nv₁²*dnS₂_dt/2 + nv₁*dvS₂_dt + vS₂*dnv₁_dt)*v̄₂ + 
        #             (-Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*dr_dt*r⁻² + Gm₂*(-16*nS₂*dnv_dt - 
        #              16*nv*dnS₂_dt + 3*dv₂S₂_dt - 7*dv₁S₂_dt)*r⁻¹ - 2*Gm₁*nS₂*nv*dr_dt*r⁻² + 
        #              2*Gm₁*nS₂*dnv_dt*r⁻¹ + 2*Gm₁*nv*dnS₂_dt*r⁻¹ + (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*vS₂ + 
        #             (3*nv₁² + 2*vv₁)*dvS₂_dt)*n + 
        #             (-Gm₂*nv₂*dr_dt*r⁻² + Gm₂*dnv₂_dt*r⁻¹ + Gm₁*nv*dr_dt*r⁻² - 
        #              Gm₁*dnv_dt*r⁻¹ - 3*nv*nv₁*dnv₁_dt - 3*nv₁²*dnv_dt/2 + nv₁*dvv₁_dt + vv₁*dnv₁_dt)*S₂
        # Gm₂r⁻²*dF2PN_dt + dGm₂r⁻²_dt*T2PN

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
        
        num*Gm₁*r⁻²
    end
    
    dvi .+= dT2PN_dt₁*c⁻⁴ 
    dvj .+= dT2PN_dt₂*c⁻⁴ 
    nothing
end

function spin_precession!(dvi,
                          dvj,
                          dvs,
                          rs,
                          vs,
                          pair::Tuple{Int, Int},
                          params::SimulationParams)

    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    ā₁  = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂  = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂  = @SVector [rs[1,  j], rs[2,  j], rs[3,  j]]
    v̄₂  = @SVector [vs[1,  j], vs[2,  j], vs[3,  j]]

    S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂  = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]
    

    m₁ = params.M[i]
    m₂ = params.M[j]
    δm = m₁ - m₂

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r² = r*r
    r⁻¹ = 1/r
    r⁻² = 1/r²

    n    = r̄/r
    nv   = dot(n, v̄)
    nv₁  = dot(n, v̄₁)
    nv₂  = dot(n, v̄₂)

    nS₁  = dot(n, S̄₁)
    nS₂  = dot(n, S̄₂)

    dr_dt = nv
    dn_dt = (r*v̄ - r̄*dr_dt)*r⁻²
    
    Gm₁ = G*m₁
    Gm₂ = G*m₂


    Gr⁻³ = G*r⁻²*r⁻¹
    dGr⁻³_dt = 3*dr_dt*Gr⁻³*r⁻¹

    dT1PN_dt₁, dT1p5PN₁, dT2PN_dt₁ = let 

        vv₂  = dot(v̄, v̄₂)
        vS₁  = dot(v̄, S̄₁)
        v₁S₁ = dot(v̄₁, S̄₁)
        v₂S₁ = dot(v̄₂, S̄₁)

        dnv_dt   = dot(n,  ā)   + dot(v̄,  dn_dt)
        dnv₁_dt  = dot(n,  ā₁)  + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)  + dot(v̄₂, dn_dt)
        dnS₁_dt  = dot(n,  dS̄₁) + dot(S̄₁, dn_dt)
        dnS₂_dt  = dot(n,  dS̄₂) + dot(S̄₂, dn_dt)        
        dvS₁_dt  = dot(v̄,  dS̄₁) + dot(S̄₁, ā)
        dv₁S₁_dt = dot(v̄₁, dS̄₁) + dot(S̄₁, ā₁)
        dv₂S₁_dt = dot(v̄₂, dS̄₁) + dot(S̄₁, ā₂)
        dvv₂_dt  = dot(v̄,  ā₂)  + dot(v̄₂, ā)

        nv₂² = nv₂^2

        ######################## PN-1 ########################
        Gm₂r⁻² = Gm₂*r⁻²
        dGm₂r⁻²_dt = -2*Gm₂*dr_dt*r⁻¹*r⁻²
        fac = (v̄₁ - 2*v̄₂)*nS₁ + S̄₁*nv - 2*n*vS₁
        dfac_dt = (v̄₁ - 2*v̄₂)*dnS₁_dt + (ā₁ - 2*ā₂)*nS₁ + 
                  S̄₁*dnv_dt + nv*dS̄₁ - 
                  2*n*dvS₁_dt - 2*vS₁*dn_dt
        dT1PN = fac*dGm₂r⁻²_dt + Gm₂r⁻²*dfac_dt
        #######################################################

        ####################### PN-1.5 ########################
        # F1p5PN = S̄₂ .- 3nS₂*n
        # dF1p5PN_dt = dS̄₂ .- 3*n*dnS₂_dt .- 3nS₂*dn_dt
        # dT1p5PN = dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₁) + (F1p5PN × dS̄₁))
        F = (3nS₂*n × S̄₁) - (S̄₂ × S̄₁)
        dF_dt = 3*((nS₂*dn_dt + n*dnS₂_dt) × S̄₁ + nS₂*n × dS̄₁) - (S̄₁ × dS̄₂) + (S̄₂ × dS̄₁)

        dT1p5PN = Gr⁻³*dF_dt + F*dGr⁻³_dt
 
        #######################################################

        ######################## PN-2 #########################
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
        #######################################################
        
        dT1PN, dT1p5PN, dT2PN
    end

    dT1PN_dt₂, dT1p5PN₂, dT2PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā, δm = -δm

        dn_dt = -dn_dt
        nv₁  = -nv₁
        nv₂  = -nv₂
        nS₁  = -nS₁
        nS₂  = -nS₂

        vv₁  = dot(v̄, v̄₁)
        vS₂  = dot(v̄, S̄₂)
        v₂S₂ = dot(v̄₂, S̄₂)
        v₁S₂ = dot(v̄₁, S̄₂)        

        dnv_dt   = dot(n,  ā)     + dot(v̄, dn_dt)
        dnv₁_dt  = dot(n,  ā₁)    + dot(v̄₁, dn_dt)
        dnv₂_dt  = dot(n,  ā₂)    + dot(v̄₂, dn_dt)
        dnS₁_dt  = dot(n,  dS̄₁)   + dot(S̄₁, dn_dt)
        dnS₂_dt  = dot(n,  dS̄₂)   + dot(S̄₂, dn_dt)
        dvS₂_dt  = dot(v̄,  dS̄₂)   + dot(S̄₂, ā)
        dv₂S₂_dt = dot(v̄₂, dS̄₂)   + dot(S̄₂, ā₂)
        dv₁S₂_dt = dot(v̄₁, dS̄₂)   + dot(S̄₂, ā₁)
        dvv₁_dt  = dot(v̄,   ā₁)   + dot(v̄₁, ā)


        nv₁² = nv₁^2

        ######################## PN-1 ########################
        Gm₁r⁻² = Gm₁*r⁻²
        dGm₁r⁻²_dt = -2*Gm₁*dr_dt*r⁻¹*r⁻²
        fac = (v̄₂ - 2*v̄₁)*nS₂ + S̄₂*nv - 2*n*vS₂
        dfac_dt = (v̄₂ - 2*v̄₁)*dnS₂_dt + (ā₂ - 2*ā₁)*nS₂ + 
                   S̄₂*dnv_dt + nv*dS̄₂ - 
                   2*n*dvS₂_dt - 2*vS₂*dn_dt
        dT1PN = fac*dGm₁r⁻²_dt + Gm₁r⁻²*dfac_dt
        #######################################################

        ###################### PN-1.5 #########################
        # F1p5PN = S̄₁ .- 3nS₁*n
        # dF1p5PN_dt = dS̄₁ - 3*n*dnS₁_dt - 3nS₁*dn_dt
        # dT1p5PN = dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₂) + (F1p5PN × dS̄₂))
        
        F = (3nS₁*n × S̄₂) - (S̄₁ × S̄₂)
        dF_dt = 3*((nS₁*dn_dt + n*dnS₁_dt) × S̄₂ + nS₁*n × dS̄₂) - (S̄₂ × dS̄₁) + (S̄₁ × dS̄₂)

        dT1p5PN = Gr⁻³*dF_dt + F*dGr⁻³_dt
        #######################################################


        ###################### PN-2 #########################
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
        #######################################################

        dT1PN, dT1p5PN, dT2PN
    end
    
    dvi .+= dT1PN_dt₁*c⁻² + dT1p5PN₁*c⁻³ + dT2PN_dt₁*c⁻⁴ 
    dvj .+= dT1PN_dt₂*c⁻² + dT1p5PN₂*c⁻³ + dT2PN_dt₂*c⁻⁴ 
    nothing
end

# function spin_precession_COM!(dvi,
#                               dvj,
#                               dvs,
#                               rs,
#                               vs,
#                               pair::Tuple{Int, Int},
#                               params::SimulationParams)

#     i, j = pair

#     ā₁  = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
#     r̄₁  = @SVector [rs[1, i], rs[2, i], rs[3, i]]
#     v̄₁  = @SVector [vs[1, i], vs[2, i], vs[3, i]]

#     ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
#     r̄₂ = @SVector [rs[1,  j], rs[2,  j], rs[3,  j]]
#     v̄₂ = @SVector [vs[1,  j], vs[2,  j], vs[3,  j]]
    
#     # S̄₁  = @SVector [rs[4, i], rs[5, i], rs[6, i]]
#     # dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]
    
#     v₁  = norm(v̄₁)
#     # v₂ = norm(v̄₂)
#     # v₁² = v₁^2
    
#     m₁ = params.M[i]
#     m₂ = params.M[j]
    
#     ā = ā₁ - ā₂
#     r̄ = r̄₁ - r̄₂
#     v̄ = v̄₁ - v̄₂
#     # rxv = r̄ × v̄

#     rv = dot(r̄, v̄)
#     # av = dot(ā, v̄)

#     r = norm(r̄) # r₁₂
#     v = norm(v̄) # v₁₂
#     v² = v*v
#     v³ = v²*v

#     n = r̄/r

#     M = m₁ + m₂
#     X1 = m₁/M
#     X2 = m₂/M

#     GM = G*M

#     Δ = X1 - X2
#     δm = m₁ - m₂
#     ν = X1*X2
#     ν² = ν*ν
#     ν³ = ν²*ν
    
#     n = r̄/r

#     nxv = n × v̄
#     nv = dot(n, v̄)
#     nxv_norm = norm(nxv)

#     # G² G³
#     dnxv_dt = 1/r*(v̄ - rv/r^2*r̄)

#     # dΩ_dt =  -4*G³*M^3*(ν^3/2 - 9*ν^2/8 - 9*ν/4 + (7/16) + dm*(-ν^2/8 - ν/8 + (-7/16))/M)*dr_dt/r^5 - 
#     #          3*G²*M^2*((-9*ν^3/8 + 75*ν^2/32 + 27*ν/4 + (3/16) + dm*(35*ν^2/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
#     #                    (13*ν^3/4 - 159*ν^2/16 - 525*ν/32 + (1/4) + dm*(-87*ν^2/16 - 75*ν/32 + (-1/4))/M)*nv^2)*dr_dt/r^4 + 
#     #          G²*M^2*(2*(-9*ν^3/8 + 75*ν^2/32 + 27*ν/4 + (3/16) + dm*(35*ν^2/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 2*(13*ν^3/4 - 
#     #                     159*ν^2/16 - 525*ν/32 + (1/4) + dm*(-87*ν^2/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3 - 
#     #          2*GM*((-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv^2*v^2 + 
#     #                (15*ν^3/16 - 195*ν^2/32 + 15*ν/8 + dm*(-75*ν^2/32 + 15*ν/8)/M)*nv^4 +
#     #                (17*ν^3/16 - 31*ν^2/8 + 19*ν/16 + (1/32) + dm*(-11*ν^2/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3 + 
#     #          GM*(2*(-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv^2*v*dv_dt + 
#     #              2*(-45*ν^3/16 + 291*ν^2/32 - 3*ν + dm*(177*ν^2/32 - 3*ν)/M)*nv*v^2*dnv_dt + 
#     #              4*(15*ν^3/16 - 195*ν^2/32 + 15*ν/8 + dm*(-75*ν^2/32 + 15*ν/8)/M)*nv^3*dnv_dt + 
#     #              4*(17*ν^3/16 - 31*ν^2/8 + 19*ν/16 + (1/32) + dm*(-11*ν^2/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt)/r^2 - 
#     #           2*GM*(ν/2 + (3/4) - 3*dm/(4*M))*dr_dt/(c^2*r^3) + (-3*G²*M^2*(ν^2/2 - 3*ν/8 + (-1/4) + dm*((1/4) - ν/8)/M)*dr_dt/r^4 + 
#     #           2*GM*(-3*ν/4 - 3*dm*ν/(2*M))*nv*dnv_dt/r^2 - 2*GM*(ν^2*(-3*ν^2/8 + 11*ν/8 + (1/16) + dm*(ν/2 + (-1/16))/M) + 
#     #           (-3*ν/4 - 3*dm*ν/(2*M))*nv^2)*dr_dt/r^3)/c^4

#     dΩ₁_dt = let δm = δm
#         num   = -2*GM*(ν/2 + (3/4) - 3*δm/(4*M))*dr_dt/(c^2*r^3)
#         numm  = -3*G^2*M^2*(ν²/2 - 3*ν/8 + (-1/4) + δm*((1/4) - ν/8)/M)*dr_dt/r^4 
#         numm += 2*GM*(-3*ν/4 - 3*δm*ν/(2*M))*nv*dnv_dt/r^2 
#         numm -= 2*GM*(ν²*(-3*ν²/8 + 11*ν/8 + (1/16) + δm*(ν/2 + (-1/16))/M) + (-3*ν/4 - 3*δm*ν/(2*M))*nv^2)*dr_dt/r^3
#         num  += c⁻⁴*numm

#         numm  = -4*G^3*M^3*(ν³/2 - 9*ν²/8 - 9*ν/4 + (7/16) + δm*(-ν²/8 - ν/8 + (-7/16))/M)*dr_dt/r^5
#         numm -= 3*G^2*M^2*((-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
#                 (13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + 
#                 (-1/4))/M)*nv^2)*dr_dt/r^4
#         numm += G^2*M^2*(2*(-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 
#                 2*(13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3
#         numm -= 2*GM*((-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v^2 + 
#                     (15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^4 + 
#                     (17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3

#         nummm  = 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v*dv_dt
#         nummm += 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv*v^2*dnv_dt
#         nummm += 4*(15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^3*dnv_dt
#         nummm += 4*(17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt

#         numm += nummm/r^2 
#         num  += numm*c⁻⁶
#         num
#     end

#     dΩ₂_dt = let δm = -δm
#         num   = -2*GM*(ν/2 + (3/4) - 3*δm/(4*M))*dr_dt/(c^2*r^3)
#         numm  = -3*G^2*M^2*(ν²/2 - 3*ν/8 + (-1/4) + δm*((1/4) - ν/8)/M)*dr_dt/r^4 
#         numm += 2*GM*(-3*ν/4 - 3*δm*ν/(2*M))*nv*dnv_dt/r^2 
#         numm -= 2*GM*(ν²*(-3*ν²/8 + 11*ν/8 + (1/16) + δm*(ν/2 + (-1/16))/M) + (-3*ν/4 - 3*δm*ν/(2*M))*nv^2)*dr_dt/r^3
#         num  += c⁻⁴*numm

#         numm  = -4*G^3*M^3*(ν³/2 - 9*ν²/8 - 9*ν/4 + (7/16) + δm*(-ν²/8 - ν/8 + (-7/16))/M)*dr_dt/r^5
#         numm -= 3*G^2*M^2*((-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v^2 + 
#                 (13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + 
#                 (-1/4))/M)*nv^2)*dr_dt/r^4
#         numm += G^2*M^2*(2*(-9*ν³/8 + 75*ν²/32 + 27*ν/4 + (3/16) + δm*(35*ν²/32 + 9*ν/8 + (-3/16))/M)*v*dv_dt + 
#                 2*(13*ν³/4 - 159*ν²/16 - 525*ν/32 + (1/4) + δm*(-87*ν²/16 - 75*ν/32 + (-1/4))/M)*nv*dnv_dt)/r^3
#         numm -= 2*GM*((-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v^2 + 
#                     (15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^4 + 
#                     (17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^4)*dr_dt/r^3

#         nummm  = 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv^2*v*dv_dt
#         nummm += 2*(-45*ν³/16 + 291*ν²/32 - 3*ν + δm*(177*ν²/32 - 3*ν)/M)*nv*v^2*dnv_dt
#         nummm += 4*(15*ν³/16 - 195*ν²/32 + 15*ν/8 + δm*(-75*ν²/32 + 15*ν/8)/M)*nv^3*dnv_dt
#         nummm += 4*(17*ν³/16 - 31*ν²/8 + 19*ν/16 + (1/32) + δm*(-11*ν²/8 + 3*ν/4 + (-1/32))/M)*v^3*dv_dt

#         numm += nummm/r^2 
#         num  += numm*c⁻⁶
#         num
#     end

#     ai = nxv*dΩ₁_dt + Ω*dnxv_dt
#     aj = nxv*dΩ₂_dt + Ω*dnxv_dt

#     dvi .+= ai
#     dvj .+= aj
#     nothing
# end

###################################################################################################################