# abstract type Body end
using StaticArrays, JLD2, LabelledArrays, FastChebInterp

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

function get_k_interpolator(order=(5,5))
    k_data_location = joinpath(@__DIR__, "..", "..", "deps", "tidal_evolution_constants", "grid.jld2")

    masses = JLD2.load(k_data_location, "Mass")#[1:10:end]
    logg = JLD2.load(k_data_location, "logg")#[1:10:end]
    logk2 = JLD2.load(k_data_location, "logk2")#[1:10:end]

    # unique_mass_ids = unique(i -> masses[i], eachindex(masses))

    logm = masses .|> u"Msun" |> ustrip .|> log10
    logg = logg .|> u"cm/s^2" |> ustrip
    logk2 = logk2

    coordinates = [SA[col...] for col in (eachcol([logm logg]'))]

    lb = [minimum(logm), minimum(logg)]
	ub = [maximum(logm), maximum(logg)]
    interpolator = chebregression(coordinates, logk2, lb, ub, order)
    k_itp(logm, logg) = interpolator(SA[logm, logg])

    return k_itp
end

const k_interpolator = get_k_interpolator()

function asidal_motion_constant_interpolated(logm::Float64, logg::Float64)
    return k_interpolator(logm, logg)
end

function StaticEquilibriumTidalPotential(system, G=ustrip(upreferred(𝒢)); Z=0.02)

    age = system.time
    n_bodies = system.n
    R_envs = Float64[]
    m_envs = Float64[]

    for i = 1:n_bodies
        
        particle = system.particles[i]
        envelope_radius, envelope_mass = if particle.structure.type isa Star && particle.structure.m < 1.25u"Msun"
                                             envelope_structure(system.particles[i], age, Z)
                                         else
                                            0.0, 0.0
                                         end
        push!(R_envs, envelope_radius)
        push!(m_envs, envelope_mass)
    end
    
    R_envs = SA[R_envs...]
    m_envs = SA[m_envs...]

    StaticEquilibriumTidalPotential(G, m_envs, R_envs)
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
    logg = log10(ustrip(u"cm/s^2", (potential.G*M_num/R_num^2 * upreferred(u"cm/s^2"))))
    logm = log10(M/DynamicQuantities.Constants.M_sun |> ustrip)
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
