# abstract type Body end
using StaticArrays, JLD2, Printf

include("../physics/tides.jl")

abstract type MultiBodyPotential end
abstract type SimulationParams end

struct DefaultSimulationParams{aType, RType, MType, LType, SType, stpType, cMType, cRType, ageType} <: SimulationParams
    a::aType # initial semi-major axes
    R::RType # radii
    M::MType # masses
    L::LType # luminosities
    S::SType # spins
    stellar_types::stpType 
    M_cores::cMType # core masses
    R_cores::cRType # core radii
    ages::ageType 
end

# struct GWSimulationParams{aType, RType, MType, M2_type, LType, SType, stpType, cMType, cRType, ageType} <: SimulationParams
#     a::aType # initial semi-major axes
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


struct PureGravitationalPotential{gType} <: MultiBodyPotential
    G::gType
    PureGravitationalPotential(G=upreferred(GRAVCONST).val) = new{typeof(G)}(G)

end

# PureGravitationalPotential() = PureGravitationalPotential(upreferred(GRAVCONST).val)

struct DynamicalTidalPotential{gType <: Real, nType, fType <: Function} <: MultiBodyPotential
    G::gType # Gravitational constant
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
function DynamicalTidalPotential(;G, n, γ)

    if n == 4
        f = tidal_factor_n4
    elseif n == 10
        f = tidal_factor_n10
    else
        f = x -> x
    end

    DynamicalTidalPotential(G, n, γ, f)
end

struct EquilibriumTidalPotential{gType <: Real} <: MultiBodyPotential
    G::gType
end

struct StaticEquilibriumTidalPotential{gType <: Real, M_env_Type, R_env_Type} <: MultiBodyPotential
    G::gType
    M_env::M_env_Type
    R_env::R_env_Type
end


function StaticEquilibriumTidalPotential(system, G=ustrip(upreferred(GRAVCONST)); Z=0.02)

    age = system.time
    n_bodies = system.n
    R_envs = typeof(1.0u"Rsun")[]
    m_envs = typeof(1.0u"Msun")[]

    for i = 1:n_bodies
        
        particle = system.particles[i]
        envelope_radius, envelope_mass = if particle.structure.stellar_type isa Star && particle.structure.m < 1.25u"Msun"
                                             envelope_structure(system.particles[i], age, Z)
                                         else
                                            0.0u"Rsun", 0.0u"Msun"
                                         end
        push!(R_envs, envelope_radius)
        push!(m_envs, envelope_mass)
    end
    
    R_envs = SA[R_envs...]
    m_envs = SA[m_envs...]

    StaticEquilibriumTidalPotential(G, m_envs, R_envs)
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
    pure_gravitational_acceleration!(dv,, rs,, params::SimulationParams,, i::Integer,, n::Integer,, potential::PureGravitationalPotential)

Acceleration function from gravitational acceleration.
"""
function pure_gravitational_acceleration!(dv,
                                          rs,
                                          params::SimulationParams,
                                          i::Int,
                                          n::Int,
                                          potential::PureGravitationalPotential)
    accel = SA[0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
   @inbounds for j = 1:n
        if j != i
            m_num::Float64 = ustrip(params.M[j])
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            accel -= potential.G * m_num * rij / (norm(rij)^3)
        end
    end
    @. dv += accel

end



"""
    dynamical_tidal_drag_force!(dv, rs, vs, params::SimulationParams, i::Int, n::Int, potential::DynamicalTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Implementing Tidal and Gravitational Wave Energy Losses in Few-body Codes: A Fast and Easy Drag Force Model](https://arxiv.org/abs/1803.08215)
"""
function dynamical_tidal_drag_force!(dv,
                           rs,
                           vs,
                           params::SimulationParams,
                           i::Int,
                           n::Int,
                           potential::DynamicalTidalPotential)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ms = params.M
    Rs = params.R

    Rₜ = ustrip(Rs[i])
    # by j on i -> j is (p)erturber, i is (t)idal object
    @inbounds for j = 1:n
        if j != i
            M = ustrip(ms[i]) + ustrip(ms[j])

            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]

            vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            rij, vij = ri - rj, vi - vj

            d = norm(rij)
            v = norm(vij)

            a = semi_major_axis(d, v^2, M, potential.G)
            e = eccentricity(rij, vij, a, M, potential.G)
            rₚ = a*(1 - e)


            J = potential.tidal_factor(e)
            ΔE::Float64 = tidal_ΔE(ustrip(ms[i]), Rₜ, ustrip(ms[j]), rₚ, 
                                   potential.γ[i], potential.G)

            ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
            ε = drag_force_coefficient(ΔE, J, a, e, M, potential.nₜ, potential.G)


            Fij = @. (-ε*(v/d^potential.nₜ)*vij/v)
            tidal_acc = Fij / ustrip(ms[i])
            accel += tidal_acc
        end
    end

    @. dv += accel
end



"""
equilibrium_tidal_drag_force!(dv, rs, vs, params::SimulationParams, i::Integer, n::Integer, potential::EquilibriumTidalPotential)

Acceleration function from equilibrium tides using the Hut 1981 prescription.
"""
function equilibrium_tidal_drag_force!(dv,
                               rs,
                               vs,
                               params::SimulationParams,
                               i::Int,
                               n::Int,
                               potential::EquilibriumTidalPotential) 

    stellar_type = ustrip(params.stellar_types[i]) |> Int
    accel = @SVector [0.0, 0.0, 0.0]

    if !(stellar_types[stellar_type] isa Star)
        return
    end

    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    Rs = params.R
    ms = params.M
    S = params.S
    
    M = ms[i]
    M_num = ustrip(M)
    R = Rs[i]
    R_num = ustrip(R)
    Ω = ustrip(S[i])
    G = unit(upreferred(GRAVCONST))*potential.G
    logg = log10(ustrip(u"cm/s^2", G*M/R^2))
    logm = log10(ustrip(u"Msun", M))
    k = asidal_motion_constant_interpolated(logm, logg)

    core_mass = params.M_cores[i]
    core_radius = params.R_cores[i]
    luminosity = params.L[i]
    age = params.ages[i]
    

    # tidal force on i by j
    @inbounds for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            rij = ri - rj
            vij = vi - vj

            r = norm(rij)
            r² = r^2
            r_hat = rij/r

            θ_dot = (rij × vij)/r²# × rij
            θ_dot_norm = norm(θ_dot)
            θ_hat = θ_dot/θ_dot_norm

            m = ms[j]
            m_num = ustrip(m)

            μ = potential.G*m_num/r²

            a = semi_major_axis(r, norm(vij)^2, m_num+M_num, potential.G)
            a_quant = a*upreferred(u"m")
            k_T = apsidal_motion_constant_over_tidal_timescale(M, R, age, core_mass, core_radius, 
                                                               stellar_type, luminosity, 
                                                               m, a_quant)# * upreferred(1.0u"yr^-1").val

            kτ = R_num^3/(potential.G*M_num)*k_T

            accel += @. -μ*3m_num/M_num*(R_num/r)^5*((k + 33vij/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end
    @. dv += accel
end


function equilibrium_tidal_drag_force!(dv,
                               rs,
                               vs,
                               params::SimulationParams,
                               i::Int,
                               n::Int,
                               potential::StaticEquilibriumTidalPotential) 

    stellar_type = ustrip(params.stellar_types[i]) |> Int
    accel = @SVector [0.0, 0.0, 0.0]

    if !(stellar_types[stellar_type] isa Star)
        return
    end

    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    Rs = params.R
    ms = params.M
    S = params.S
    
    M = ms[i]
    M_num = ustrip(M)
    R = Rs[i]
    R_num = ustrip(R)
    Ω = ustrip(S[i])
    G = unit(upreferred(GRAVCONST))*potential.G
    logg = log10(ustrip(u"cm/s^2", G*M/R^2))
    logm = log10(ustrip(u"Msun", M))
    k = asidal_motion_constant_interpolated(logm, logg)

    envelope_mass = potential.M_env[i]
    envelope_radius = potential.R_env[i]

    luminosity = params.L[i]    

    # tidal force on i by j
    @inbounds for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            rij = ri - rj
            vij = vi - vj

            r = norm(rij)
            r² = r^2
            r_hat = rij/r

            θ_dot = (rij × vij)/r²# × rij
            θ_dot_norm = norm(θ_dot)
            θ_hat = θ_dot/θ_dot_norm

            m = ms[j]
            m_num = ustrip(m)

            μ = potential.G*m_num/r²

            a = semi_major_axis(r, norm(vij)^2, m_num+M_num, potential.G)
            a_quant = a*upreferred(u"m")
            k_T::Float64 = apsidal_motion_constant_over_tidal_timescale(M, R,
                                                               envelope_mass, envelope_radius,
                                                               stellar_type, luminosity, 
                                                               m, a_quant) * upreferred(u"yr^-1").val

            kτ = R_num^3/(potential.G*M_num)*k_T

            accel += @. -μ*3m_num/M_num*(R_num/r)^5*((k + 33vij/r*kτ)*r_hat - (Ω - θ_dot_norm)*kτ*θ_hat)
        end
    end
    @. dv += accel
end

function PN1_acceleration!(dv,
                           rs,
                           vs,
                           params::SimulationParams,
                           i::Int,
                           n::Int,
                           potential::PN1Potential)
                           
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

    G = potential.G 
    
    m₁ = params.M[i].val
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

            r⁻¹ = 1/r

            n = r̄*r⁻¹

            v₁v₂ = dot(v̄₁, v̄₂) 
            nv₁ = dot(n, v̄₁)
            nv₂ = dot(n, v̄₂)
            nv₂² = nv₂^2

            m₂ = params.M[j].val
            G_r = G*r⁻¹
            G_r² = G_r*r⁻¹
      
            # PN-1 acceleration
            a = @. n*(G_r²*m₂)*(5*G_r*m₁ + 4*G_r*m₂ + 3/2*nv₂² - v₁² + 4*v₁v₂ - 2*v₂²) +
                    (4*nv₁ - 3*nv₂)*v̄

            accel += a
        end

    end
    # println(accel)
    @. dv += accel * c⁻²
end


function PN2_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN2Potential)
                           
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

    G = potential.G 
    
    m₁ = params.M[i].val
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

            # nv₁⁴ = nv₁^4
            nv₂⁴ = nv₂^4

            m₂ = params.M[j].val
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2

            G_r = G*r⁻¹
            G_r² = G_r*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

      
            # PN-2.5 acceleration
            # expression is split up to avoid allocations that can appear in long expressions
            a_num = -57G³_r⁴*m₁²m₂/4 - 69G³_r⁴*m₁m₂²/2 - 9G³_r⁴*m₂^3 
            a_num += G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
                                4*v₁v₂*v₂² - 2v₂^4)
            a_num += G²_r³*m₁m₂*(39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂² - 15/4*v₁² - 5/2*v₁v₂ + 5/4*v₂²) 
            a_num += G²_r³*m₂^2*(2*nv₁² - 4*nv₁*nv₂ - 6*nv₂² - 8*v₁v₂ + 4v₂²) 
            a₂1 = n*a_num

            a_num = G²_r³*m₂^2*(-2*nv₁ - 2*nv₂) + G²*m₁m₂/r³*(-63/4*nv₁ + 55/4*nv₂) 
            a_num += G_r²*m₂*(-6*nv₁*nv₂² + 9/2*nv₂³ + nv₂*v₁² - 4*nv₁*v₁v₂ + 
                                4*nv₂*v₁v₂ + 4*nv₁*v₂² - 5*nv₂*v₂²)
            a₂2 = v̄*a_num
            a = a₂1 + a₂2

            accel += @. a 
        end
    end

    @. dv += accel * c⁻⁴

end

function PN2_5_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN2_5Potential)
                           
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    
    m₁ = params.M[i].val
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @inbounds for j = 1:n
        if j != i   
            r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            r̄ = r̄₁ - r̄₂
            v̄ = v̄₁ - v̄₂

            r = norm(r̄) # r₁₂
            v = norm(v̄) # v₁₂

            r⁻¹ = 1/r

            v² = v^2

            n = r̄*r⁻¹

            nv = dot(n, v̄)

            m₂ = params.M[j].val
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2

            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            # # PN-2.5 acceleration
            a_num = 208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²
            a1 = a_num*n
            a_num = 8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²
            a2 = a_num*v̄
            a = a1 + a2

            # a = @. (208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²)*n + 
            #          (8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²)*v̄

            accel += a
        end
    end
    @. dv += accel * c⁻⁵
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

    G = potential.G 
    
    m₁ = params.M[i].val
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

            m₂ = params.M[j].val
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

    G = potential.G 
    
    m₁ = params.M[i].val
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

            m₂ = params.M[j].val
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

function PN1_to_3_5_acceleration!(dv,
                                rs,
                                vs,
                                params::SimulationParams,
                                i::Int,
                                n::Int,
                                potential::PNPotential)
                           

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

    G = potential.G 
    
    m₁ = params.M[i].val
    accel = @SVector [0.0, 0.0, 0.0]
    a₂ =  @MVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @fastmath @inbounds for j = 1:n
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

            m₂ = params.M[j].val
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2
            # m₁²m₂² = m₁^2*m₂^2

            G_r = G*r⁻¹
            G_r² = G_r*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            # r₁′ = r₂′ = 1.0
      
            # PN-1 acceleration
            a₁ = @. n*(G_r²*m₂)*(5*G_r*m₁ + 4*G_r*m₂ + 3/2*nv₂² - v₁² + 4*v₁v₂ - 2*v₂²) +
                    (4*nv₁ - 3*nv₂)*v̄

            # PN-2 acceleration
            # expression is split up to avoid allocations that can appear in long expressions
            a_num = -57G³_r⁴*m₁²m₂/4 - 69G³_r⁴*m₁m₂²/2 - 9G³_r⁴*m₂^3 
            a_num += G*m₂/r²*(-15/8*nv₂⁴ + 3/2*nv₂²*v₁² - 6*nv₂²*v₁v₂ - 2*v₁v₂² + 9/2*nv₂²*v₂² + 
                                4*v₁v₂*v₂² - 2v₂^4)
            a_num += G²_r³*m₁m₂*(39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂² - 15/4*v₁² - 5/2*v₁v₂ + 5/4*v₂²) 
            a_num += G²_r³*m₂^2*(2*nv₁² - 4*nv₁*nv₂ - 6*nv₂² - 8*v₁v₂ + 4v₂²) 
            a₂1 = n*a_num

            a_num = G²_r³*m₂^2*(-2*nv₁ - 2*nv₂) + G²*m₁m₂/r³*(-63/4*nv₁ + 55/4*nv₂) 
            a_num += G_r²*m₂*(-6*nv₁*nv₂² + 9/2*nv₂³ + nv₂*v₁² - 4*nv₁*v₁v₂ + 
                                4*nv₂*v₁v₂ + 4*nv₁*v₂² - 5*nv₂*v₂²)
            a₂2 = v̄*a_num
            a₂ = a₂1 + a₂2

            # PN-2.5 acceleration
            a₃ = @. (208G³_r⁴*m₁m₂²/15*nv - 24G³_r⁴*m₁²m₂/5*nv + 12G²_r³*m₁m₂/5*v²)*n + 
                     (8G³_r⁴*m₁²m₂/5 - 32G³_r⁴*m₁m₂²/5 - 4G²_r³*m₁m₂/5*v²)*v̄

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
            accel += @. a₁*c⁻² + a₂*c⁻⁴ + a₃*c⁻⁵# + a₄*c⁻⁶# + a₅*c⁻⁷


        end

    end
    # println(accel)
    @. dv += accel
end