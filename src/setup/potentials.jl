# abstract type Body end
using StaticArrays, JLD2, Printf

include("../physics/tides.jl")

abstract type FewBodyPotential end
abstract type SimulationParams end

struct DefaultSimulationParams{aType, RType, MType, LType, SType, stpType, cMType, cRType, ageType} <: SimulationParams
    a::aType
    R::RType
    M::MType
    L::LType
    S::SType
    stellar_types::stpType
    M_cores::cMType
    R_cores::cRType
    ages::ageType
end


struct PureGravitationalPotential{gType <: Real} <: FewBodyPotential
    G::gType
end

PureGravitationalPotential() = PureGravitationalPotential(upreferred(𝒢).val)

struct DynamicalTidalPotential{gType <: Real, nType, fType <: Function} <: FewBodyPotential
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

struct EquilibriumTidalPotential{gType <: Real} <: FewBodyPotential
    G::gType
end

struct StaticEquilibriumTidalPotential{gType <: Real, M_env_Type, R_env_Type} <: FewBodyPotential
    G::gType
    M_env::M_env_Type
    R_env::R_env_Type
end


function StaticEquilibriumTidalPotential(system, G=ustrip(upreferred(𝒢)); Z=0.02)

    age = system.time
    n_bodies = system.n
    R_envs = typeof(1.0u"Rsun")[]
    m_envs = typeof(1.0u"Msun")[]

    for i = 1:n_bodies
        
        particle = system.particles[i]
        envelope_radius, envelope_mass = if particle.structure.type isa Star && particle.structure.m < 1.25u"Msun"
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

Potential for Post-Newtonian 1 (PNA1) acceleration
"""
struct PN1Potential{gType <: Real} <: FewBodyPotential
    G::gType
end

"""
    PN2Potential{gType <: Real}

Potential for Post-Newtonian 2 (PNA2) acceleration
"""
struct PN2Potential{gType <: Real} <: FewBodyPotential
    G::gType
end

"""
    PN1Potential{gType <: Real}

Potential for Post-Newtonian 2.5 (PNA2.5) acceleration
"""
struct PN2_5Potential{gType <: Real} <: FewBodyPotential
    G::gType
end

"""
    PN1Potential{gType <: Real}

Potential for Post-Newtonian 1 to 2.5 (PNA1, PNA2, PNA2.5) acceleration
"""
struct PNPotential{gType <: Real} <: FewBodyPotential
    G::gType
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

    accel = @SVector [0.0, 0.0, 0.0];
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
    G = unit(upreferred(𝒢))*potential.G
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
    G = unit(upreferred(𝒢))*potential.G
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
                                                               m, a_quant) * upreferred(1.0u"yr^-1").val

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
                           
    r₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    G = potential.G

    m₁ = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    @inbounds for j = 1:n
        if j != i                 
            r₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            m₂ = params.M[j]

            r = r₂ - r₁
            n = r/norm(r)

            G_r = G/r

            a = n*(-v₁^2 - 2v₂^2 + 4*dot(v₁, v₂) + 3/2*(dot(n, v₂)^2) + 5*G_r*m₁ + 4*G_r*m₂) +
                (v₁ - v₂)*(4*dot(n, v₁) - 3*dot(n, v₂))

            accel += @. G_r*m₂/r*a*c⁻²
        end

    end
end


function PN2_acceleration!(dv,
                                       rs,
                                       vs,
                                       params::SimulationParams,
                                       i::Int,
                                       n::Int,
                                       potential::PN2Potential)
                           
    r₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁² = v₁^2

    G = potential.G
    
    mi = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    @inbounds for j = 1:n
        if j != i                 
            r₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂² = v₂^2

            r = r₁ - r₂
            v = v₁ - v₂
            n = r/norm(r)

            v₁v₂ = dot(v₁, v₂)
            nv₁ = dot(n, v₁)
            nv₂ = dot(n, v₂)

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            G_r = G/r

            a = n*(-2*v₂²^2 + 4v₂²*v₁v₂ - 2*v₁v₂^2 + 3/2*v₁²*nv₂² +
                   9/2*v₂²*nv₂² - 6*v₁v₂*nv₂² - 15/8*nv₂^4 + 
                   (G_r*mi)*(-15/4*v₁² + 5/4*v₂² - 5/2*v₁v₂ +
                             39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂²
                            ) + 
                   (G_r*mj)*(4*v₂² - 8*v₁v₂ + 2*nv₁² - 4*nv₁*nv₂ - 6*nv₂²)
                  ) + 
                v*(v₁²*nv₂ + 4*v₂²*nv₁ -5v₂²*nv₂^3 +
                   (G_r*mi)*(-63/4*nv₁ + 55/4*nv₂) + (G_r*mj)*(-2*nv₁ - 2*nv₂) 
                  ) + 
                G^3*mj/r^4*n*(-57/4*mi^2 - 9*mj^2 - 69/2*mi*mj)


            accel += @. G*mj/r^2*a * c⁻⁴
        end

    end
end

function PN2_5_acceleration!(dv,
                                    rs,
                                    vs,
                                    params::SimulationParams,
                                    i::Int,
                                    n::Int,
                                    potential::PN2_5Potential)
                           
    r₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    G = potential.G

    mi = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    @inbounds for j = 1:n
        if j != i                 
            r₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            r = r₁ - r₂
            v = v₁ - v₂
            n = r/norm(r)

            G_r = G/r

            a = 4/5*G^2*mi*mj/r^3*(v*(-v^2 + 2*(G_r*mi) - 8*(G_r*mj))) + 
                dot(n, dot(n, v))*(3*v₂ - 6*(G_r*mi) + 52/3*(G_r*mj))

            accel += a*c⁻⁵
        end

    end
end

function PN1_to_2_5_acceleration(dv,
                                       rs,
                                       vs,
                                       params::SimulationParams,
                                       i::Int,
                                       n::Int,
                                       potential::PNPotential)
                           
    r₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁² = v₁^2

    G = potential.G
    
    mi = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    @inbounds for j = 1:n
        if j != i                 
            r₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂² = v₂^2

            r = r₁ - r₂
            v = v₁ - v₂
            n = r/norm(r)

            v₁v₂ = dot(v₁, v₂) 
            nv₁ = dot(n, v₁)
            nv₂ = dot(n, v₂)

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            G_r = G/r

            a₂ = n*(-v₁^2 - 2v₂^2 + 4*v₁v₂ + 3/2*(nv₂^2) + 5*G_r*m₁ + 4*G_r*m₂) +
                 v*(4*nv₁ - 3*nv₂)

            a₄ = n*(-2*v₂²^2 + 4v₂²*v₁v₂ - 2*v₁v₂^2 + 3/2*v₁²*nv₂² +
                   9/2*v₂²*nv₂² - 6*v₁v₂*nv₂² - 15/8*nv₂^4 + 
                   (G_r*mi)*(-15/4*v₁² + 5/4*v₂² - 5/2*v₁v₂ +
                             39/2*nv₁² - 39*nv₁*nv₂ + 17/2*nv₂²
                            ) + 
                   (G_r*mj)*(4*v₂² - 8*v₁v₂ + 2*nv₁² - 4*nv₁*nv₂ - 6*nv₂²)
                  ) + 
                v*(v₁²*nv₂ + 4*v₂²*nv₁ -5v₂²*nv₂^3 +
                   (G_r*mi)*(-63/4*nv₁ + 55/4*nv₂) + (G_r*mj)*(-2*nv₁ - 2*nv₂) 
                  ) + 
                G^3*mj/r^4*n*(-57/4*mi^2 - 9*mj^2 - 69/2*mi*mj)

            a₅ = 4/5*G^2*mi*mj/r^3*(v*(-v^2 + 2*(G_r*mi) - 8*(G_r*mj))) + 
                 dot(n, dot(n, v))*(3*v₂ - 6*(G_r*mi) + 52/3*(G_r*mj))

            a₂ *= G_r*m₂/r
            a₄ *= G_r*m₂/r

            accel += a₂*c⁻² + a₄*c⁻⁴ + a₅*c⁻⁵
        end

    end
end

function Base.show(io::IO, params::DefaultSimulationParams)
    print(nameof(typeof(params)))
    print(":")
    println()
    for prop in propertynames(params)
        val = getproperty(params, prop) 
        un = unit(val[1])
        val = ustrip(val)

        @printf(io, "   %-16s %s %s", "$prop", "$val", "$un")
        println(io)
    end
end