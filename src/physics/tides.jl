

######################################
# Dynamical Tides
######################################

"""
    tidal_ΔE(Mₜ, Rₜ, Mₚ, rₚ, γ, G)

Energy dissipated in a tidal interaction between a tidal object with mass
Mₜ and radius Rₜ and a perturber with mass Mₜ, where γ is the polytropic index
of the tidal object.
"""
function tidal_ΔE(Mₜ, Rₜ, Mₚ, rₚ, γ, G)
    rₚ_Rₜ = rₚ/Rₜ
    
    η = √(Mₜ/(Mₜ + Mₚ)*rₚ_Rₜ^3)
    # 0.08 > η > 6.3 && return 0.0
    if η < 1.0 || η > 10.0
        return 0.0
    end
    T₂η, T₃η = tidal_structure_function(η, γ)
    # if any(isinf, (T₂η, T₃η))
    #     @show rₚ_Rₜ, η, T₂η, T₃η
    # end

    ΔE = G*Mₚ^2/Rₜ*((Rₜ/rₚ)^6*T₂η) #+ (Rₜ/rₚ)^8*T₃η)
    return ΔE
end 

function tidal_factor_n4(e)::Float64
    return π/2*(2 + 7e^2 + e^4)
end

function tidal_factor_n10(e)::Float64
    return π/128*(128 + 2944e^2 + 10528e^4 + 
                  8960e^6 + 1715e^8 + 35e^10)
end

function drag_force_coefficient(ΔE, J, a, e, M, n, G)
    # @show a, e
    return ΔE*0.5*(((a*(1 - e^2))^(n - 0.5)))/(√(G*M)*J)
end


function tidal_structure_function(η, γ)
    # @show γ
    l2 = tidal_structure_coefficients[(γ, 2)]# .* 1e-4
    l3 = tidal_structure_coefficients[(γ, 3)]# .* 1e-4

    all(iszero(l2)) && return 0.0, 0.0

    x = log10(η)
    # @show x
    logT₂ = l2[1] + l2[2]*x + l2[3]*x^2 + l2[4]*x^3 + l2[5]*x^4 + l2[6]*x^5 
    logT₃ = l3[1] + l3[2]*x + l3[3]*x^2 + l3[4]*x^3 + l3[5]*x^4 + l3[6]*x^5 
    return 10^(logT₂), 10^(logT₃)
end


"""
    tertiary_tide_inner_inspiral_rate(a₁, a₂, masses, radii, τ)

Estimate the inspiral rate of an inner binary in a hierarchichal triple due to
tertiary tide effects. Formula is an empirical fit to viscoelastic simulations, taken
from [ref!](https://arxiv.org/abs/1910.12852)

τ is the viscoelastic relaxation time
"""
# ``a \frac{da\_1}{dt} = (2.22 \times 10^{-8} yrs^{-1}) \frac{4q}{1 + q}^2 \left(\frac{R_3}{100R_\odot}\right)^5.2`` 
function tertiary_tide_inner_inspiral_rate(a₁, a₂, masses, radii, τ=0.5u"yr")

    q = masses[1]/masses[2]
    R₃ = radii[3]

    return a₁*(2.22e-8u"1/yr")*4q/(1 + q)^2*(R₃/100u"Rsun")^(5.2)*
           (a₁/0.2u"AU")^4.8*(a₂/2.0u"AU")^(-10.2)*(τ/0.534u"yr")^(-1.0)
end

# function tertiary_tide_inner_inspiral_rate(triple::TripleInitialConditions, τ=0.5u"yr")
#     tertiary_tide_inner_inspiral_rate(triple.elements[1].a, triple.elements[2].a, 
#                                       triple.structure.m, triple.structure.R, τ)
# end


######################################
# Equilibrium Tides
######################################

function apsidal_motion_constant_over_tidal_timescale(mass, radius, core_mass, core_radius, 
                                                      stellar_type, spin, luminosity, 
                                                      mass_perturber,
                                                      orbital_period, semi_major_axis)

    return k_over_T(mass, radius, core_mass, 
                    core_radius, stellar_type, spin,
                    luminosity, orbital_period,
                    mass_perturber, semi_major_axis)
end

function k_over_T_convective(mass, radius, core_mass, core_radius, 
                             spin, luminosity, orbital_period)
        R_env = envelope_radius(mass, radius, core_radius, stellar_type)
        M_env = mass - core_mass
        Pspin = 2π/spin
        Ptid = 1.0/abs(1/orbital_period - 1/Pspin)

        τ_conv = 0.4311*cbrt((M_env*R_env*(R - 0.5*R_env))/(3*luminosity))
        f_conv = min(1, (Ptid/(2τ_conv))^2)

        return 2/21*f_conv/τ_conv*M_env/mass
end

function k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
    E₂ = 1,592e-9*mass^2.84 # second-order tidal coefficient
    q₂ = mass_perturber/mass

    return 1.9782e4*mass*radius/semi_major_axis^5*(1 - q₂)^(5/6)*E₂
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::MainSequence, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    if mass < 1.25u"Msun"
        return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
    else
        return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
    end
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::DeeplyOrFullyConvectiveLowMassMainSequence, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::HertzsprungGap, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end


function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::FirstGiantBranch, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end


function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::FirstAsymptoticGiantBranch, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end


function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::SecondAsymptoticGiantBranch, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end


function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::HertzsprungGapNakedHelium, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end


function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::GiantBranchNakedHelium, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::CoreHeliumBurning, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type::MainSequenceNakedHelium, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)

    return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
end

function k_over_T(mass, radius, core_mass, 
                  core_radius, stellar_type<:CompactObject, spin,
                  luminosity, orbital_period,
                  mass_perturber, semi_major_axis)
    return 0.0
end