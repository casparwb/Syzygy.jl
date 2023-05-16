# abstract type Body end
using StaticArrays, JLD2, LabelledArrays#, LoopVectorization

include("../physics/tides.jl")

abstract type FewBodyPotential end

struct PureGravitationalPotential{gType <: Real} <: FewBodyPotential
    G::gType
end

struct DynamicalTidalPotential{gType <: Real, nType, fType <: Function} <: FewBodyPotential
    G::gType # Gravitational constant
    nₜ::Int  # Tidal force power constant
    γ::nType # Polytropic index of each star
    tidal_factor::fType
end

struct EquilibriumTidalPotential{gType <: Real, τType, kType, βType} <: FewBodyPotential
    G::gType
    τ::τType # time lag
    k::kType # apsidal motion constant
    β::βType # radius of gyration
end

struct EquilibriumTidalPotential2{gType <: Real, kType} <: FewBodyPotential
    G::gType
    k_over_T::kType
end

function PureGravitationalPotential(G=𝒢.val)
    PureGravitationalPotential(G)
end


"""
    DynamicalTidalPotential(;G, n, γ)


Set up the dynamical tidal potential for a system. 

# Keyword arguments
- `G`: gravitational constant.
- `n`: tidal force power index
- `γ`: vector of polytropic indixes of each body in the system
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

function EquilibriumTidalPotential(M; G, τ)

    tidal_evolution_constants = JLD2.load(joinpath(@__DIR__, "..", "..", "deps", 
                                          "tidal_evolution_constants", "z=0.0134.jld2"),
                                          "result")

    masses = keys(tidal_evolution_constants) |> collect

    βs = Float64[]
    ks = Float64[]
    for m in M
        idx = masses[argmin(abs.(m .- masses))]
        β = tidal_evolution_constants[idx].beta
        k = 10 ^ tidal_evolution_constants[idx].logk2

        push!(βs, β)
        push!(ks, k)

    end

    return EquilibriumTidalPotential(G, τ, ks, βs)
end

function EquilibriumTidalPotential2(mass, radius, core_mass, 
                                    core_radius, stellar_type, spin,
                                    luminosity, orbital_period,
                                    mass_perturber, semi_major_axis; 
                                    G=upreferred(𝒢))
    if (stellar_type == 1 && mass < 1.25u"Msun") && any(stellar_type .== (0, 2, 3, 5, 6, 8, 9))

        R_env = envelope_radius(mass, radius, core_radius, stellar_type)
        M_env = mass - core_mass
        Pspin = 2π/spin
        Ptid = 1.0/abs(1/orbital_period - 1/Pspin)

        τ_conv = 0.4311*cbrt((M_env*R_env*(R - 0.5*R_env))/(3*luminosity))
        f_conv = min(1, (Ptid/(2τ_conv))^2)

        k_over_T = 2/21*f_conv/τ_conv*M_env/mass

    elseif (stellar_type == 1 && mass > 1.25u"Msun") || stellar_type == 4 || stellar_type == 7
        E₂ = 1,592e-9*mass^2.84 # second-order tidal coefficient
        q₂ = mass_perturber/mass

        k_over_T = 1.9782e4*mass*radius/semi_major_axis^5*(1 - q₂)^(5/6)*E₂
    else
        k_over_T = 0.0
    end

    return k_over_T
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



