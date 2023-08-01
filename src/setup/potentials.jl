# abstract type Body end
using StaticArrays, JLD2, LabelledArrays, DataInterpolations

include("../physics/tides.jl")

abstract type FewBodyPotential end

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


struct EquilibriumTidalPotential{gType <: Real, kType, k_TType} <: FewBodyPotential
    G::gType
    k::kType
    k_over_T::k_TType
end

function EquilibriumTidalPotential(system, G=upreferred(𝒢).val; Z=0.0122)
    n_particles = system.n

    k_itp = LinearInterpolation(JLD2.load("../../deps/tidal_evolution_constants/grid.jld2", "logk2"),
                                JLD2.load("../../deps/tidal_evolution_constants/grid.jld2", "logg"))

    k_over_Ts = Function[]
    for i ∈ 1:n_particles
        particle = system.particles[i]
        
        mass = particle.structure.m
        radius = particle.structure.R
        core_mass = particle.structure.m_core
        core_radius = particle.structure.R_core
        stellar_type = particle.structure.type
        luminosity = particle.structure.L
        stellar_type = particle.structure.type
        
        logg = log10(u"cm/s^2"(G*mass/radius^2).val)*u"cm/s^2"
        k = k_itp(logg)

        kT(mass_perturber) = k_over_T(mass, radius, core_mass, 
                                      core_radius, stellar_type, spin,
                                      luminosity, orbital_period,
                                      mass_perturber, semi_major_axis, Z)
        push!(k_over_Ts, kT)
    end

end



"""
    pure_gravitational_acceleration!(dv,, rs,, params::T where T <: LArray,, i::Integer,, n::Integer,, potential::PureGravitationalPotential)

Acceleration function from gravitational acceleration.
"""
function pure_gravitational_acceleration!(dv,
                                          rs,
                                          params::T where T <: LArray,
                                          i::Integer,
                                          n::Integer,
                                          potential::PureGravitationalPotential)

    ms = params.M
    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    @inbounds for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            accel -= potential.G * ms[j] * rij / (norm(rij)^3)
        end
    end
    @. dv += accel

    # nothing
    # accel
end


"""
    dynamical_tidal_drag_force!(dv, rs, vs, params::T where T <: LArray, i::Integer, n::Integer, potential::DynamicalTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Implementing Tidal and Gravitational Wave Energy Losses in Few-body Codes: A Fast and Easy Drag Force Model](https://arxiv.org/abs/1803.08215)
"""
function dynamical_tidal_drag_force!(dv,
                           rs,
                           vs,
                           params::T where T <: LArray,
                           i::Integer,
                           n::Integer,
                           potential::DynamicalTidalPotential)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ms = params.M
    Rs = params.R

    Rₜ = Rs[i]
    # by j on i -> j is (p)erturber, i is (t)idal object
    @inbounds for j = 1:n
        if j != i
            M = ms[i] + ms[j]

            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]

            vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]

            rij, vij = ri - rj, vi - vj

            d = norm(rij)
            v = norm(vij)

            a = semi_major_axis(d, v^2, M, potential.G)
            e = eccentricity(rij, vij, a, M, potential.G)
            rₚ = a*(1 - e)


            J = potential.tidal_factor(e)
            ΔE::Float64 = tidal_ΔE(ms[i], Rₜ, ms[j], rₚ, 
                                   potential.γ[i], potential.G)

            ΔE = ifelse(isinf(ΔE), 0.0, ΔE)
            ε = drag_force_coefficient(ΔE, J, a, e, M, potential.nₜ, potential.G)


            Fij = @. (-ε*(v/d^potential.nₜ)*vij/v)
            tidal_acc = Fij / ms[i]
            accel += tidal_acc
        end
    end

    @. dv += accel
end


"""
    dynamical_tidal_drag_force!(dv, rs, vs, params::T where T <: LArray, i::Integer, n::Integer, potential::EquilibriumTidalPotential)

Acceleration function from dynamical tides. This model is adapted from 
[Tidal evolution in close binary systems.](https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H)
"""
function equilibrium_tidal_drag_force!(dv,
                               rs,
                               vs,
                               params::T where T <: LArray,
                               i::Integer,
                               n::Integer,
                               potential::EquilibriumTidalPotential)

    accel = @SVector [0.0, 0.0, 0.0];

    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    Rs = params.R
    ms = params.M
    S = params.S

    τ = potential.τ[i]
    k = potential.k[i]

    M = ms[i]
    R = Rs[i]
    Ω = S[i]
    # by j on i => i is Primary
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

            μ = potential.G*m/r²
            accel += @. μ*3m/M*(R/r)^5*k*((1 + 3vij/r*τ)*r_hat - (Ω - θ_dot_norm)*τ*θ_hat)
        end
    end
    # @show accel
    @. dv += accel
end



