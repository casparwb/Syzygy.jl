

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

######################################
# Equilibrium Tides
######################################
"""
apsidal_motion_constant_over_tidal_timescale(mass, radius, core_mass, core_radius, 
                                                    stellar_type, spin, luminosity, 
                                                    mass_perturber,
                                                    orbital_period, semi_major_axis)

Calculate and return the fraction of the apsidal motion constant `k` over the 
tidal timescale `T`. If the given stellar type is a type with a radiative envelope, such as 
a core helium burning star, or a massive main sequence star, the prescription from 
Hurley et al. 2002 will be used (Eq. 42). Otherwise, the the prescription from Preece et al. 2022 is used.                   
"""
function apsidal_motion_constant_over_tidal_timescale(mass, radius, age, core_mass, core_radius, 
                                                      stellar_type, luminosity, 
                                                      mass_perturber,
                                                      semi_major_axis)::typeof(1.0u"s^-1")
    return k_over_T(mass, radius, core_mass, age,
                    core_radius, stellar_type,
                    luminosity,
                    mass_perturber, semi_major_axis) |> upreferred
end


function k_over_T(mass, radius, core_mass, age,
                  core_radius, stellar_type,
                  luminosity,
                  mass_perturber, semi_major_axis, Z=0.0122)

    mass = u"Msun"(mass)
    radius = u"Rsun"(radius)
    core_mass = u"Msun"(core_mass)
    core_radius = u"Rsun"(core_radius)
    luminosity = u"Lsun"(luminosity)
    age = u"Myr"(age)

    mass_perturber = u"Msun"(mass_perturber)
    semi_major_axis = u"Rsun"(semi_major_axis)

    if !(stellar_types[stellar_type] isa Star)
        return 0.0u"1/yr"
    end

    if mass < 1.25u"Msun"
        envelope_radius, envelope_mass = envelope_structure(mass, radius, core_mass, 
                                                            core_radius, stellar_type, age)

        return k_over_T_convective(mass, radius, envelope_mass, envelope_radius, luminosity, Z)
    else
        return k_over_T_radiative(mass.val, radius.val, mass_perturber.val, semi_major_axis.val)
    end
end


function k_over_T_convective(mass, radius, envelope_mass, envelope_radius,
                             luminosity, Z=0.0122)
                             

        R_conv, M_conv = envelope_radius, envelope_mass
        t_conv = 0.4311*u"yr"(cbrt((3M_conv*R_conv^2)/luminosity)) # Preece et al. 2022

        log10z = log10(Z)
        a =  0.630*log10z + 2.72  # a(z) for envelope
        b = -0.219*log10z + 0.68  # b(z) for envelope
        c = -0.023*log10z + 0.220 # c(z) for envelope


        return (R_conv/radius)^a*(M_conv/mass)^b*c/t_conv
end

function k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
    E₂ = 1.592e-9*mass^2.84 # second-order tidal coefficient
    q₂ = mass_perturber/mass
    return 1.9782e4*mass*radius/semi_major_axis^5*(1 + q₂)^(5/6)*E₂ * u"1/yr"
end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::DeeplyOrFullyConvectiveLowMassMainSequence, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end

# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::HertzsprungGap, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::FirstGiantBranch, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::FirstAsymptoticGiantBranch, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::SecondAsymptoticGiantBranch, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::HertzsprungGapNakedHelium, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end


# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::GiantBranchNakedHelium, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_convective(mass, radius, core_mass, core_radius, spin, luminosity, orbital_period)
# end

# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::CoreHeliumBurning, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
# end

# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::MainSequenceNakedHelium, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#     return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
# end

# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::T where T <: CompactObject, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)
#     return 0.0
# end

# function k_over_T(mass, radius, core_mass, 
#                   core_radius, stellar_type::Int, spin,
#                   luminosity, orbital_period,
#                   mass_perturber, semi_major_axis)

#         if (stellar_type == 1 && mass < 1.25u"Msun") && (
#             stellar_type == 0 || stellar_type == 2 || stellar_type == 3 || 
#             stellar_type == 5 || stellar_type == 6 || stellar_type == 8 || stellar_type == 9)

#         R_env = envelope_radius(mass, radius, core_radius, stellar_type)
#         M_env = mass - core_mass
#         Pspin = 2π/spin
#         Ptid = 1.0/abs(1/orbital_period - 1/Pspin)

#         τ_conv = 0.4311*cbrt((M_env*R_env*(R - 0.5*R_env))/(3*luminosity))
#         f_conv = min(1, (Ptid/(2τ_conv))^2)

#         k_over_T = 2/21*f_conv/τ_conv*M_env/mass

#     elseif (stellar_type == 1 && mass > 1.25u"Msun") || stellar_type == 4 || stellar_type == 7
#         E₂ = 1.592e-9*mass^2.84 # second-order tidal coefficient
#         q₂ = mass_perturber/mass

#         k_over_T = 1.9782e4*mass*radius/semi_major_axis^5*(1 - q₂)^(5/6)*E₂
#     else
#         k_over_T = 0.0
#     end

#     return k_over_T
# end