using FastChebInterp

######################################
# Dynamical Tides
######################################

"""
    tidal_ΔE(Mₜ, Rₜ, Mₚ, rₚ, γ, G)

Energy dissipated in a tidal interaction between a tidal object with mass
Mₜ and radius Rₜ and a perturber with mass Mₚ where γ is the polytropic index
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

function tidal_factor_n4(e)
    return π/2*(2 + 7e^2 + e^4)
end

function tidal_factor_n10(e)
    return π/128*(128 + 2944e^2 + 10528e^4 + 
                  8960e^6 + 1715e^8 + 35e^10)
end

function drag_force_coefficient(ΔE, J, a, e, M, n, G)
    return ΔE*0.5*(((a*(1 - e^2))^(n - 0.5)))/(√(G*M)*J)
end


function tidal_structure_function(η, γ)
    l2 = tidal_structure_coefficients[(γ, 2)]# .* 1e-4
    l3 = tidal_structure_coefficients[(γ, 3)]# .* 1e-4

    all(iszero(l2)) && return 0.0, 0.0

    x = log10(η)
    logT₂ = l2[1] + l2[2]*x + l2[3]*x^2 + l2[4]*x^3 + l2[5]*x^4 + l2[6]*x^5 
    logT₃ = l3[1] + l3[2]*x + l3[3]*x^2 + l3[4]*x^3 + l3[5]*x^4 + l3[6]*x^5 
    return 10^(logT₂), 10^(logT₃)
end

######################################
# Equilibrium Tides
######################################
"""
apsidal_motion_constant_over_tidal_timescale(mass::Unitful.Mass, radius,
                                             envelope_mass, envelope_radius,
                                             stellar_type, luminosity, 
                                             mass_perturber,
                                             semi_major_axis, Z=0.02)

For a stellar object with the given properties `mass`, `radius`, `envelope_mass`, and `envelope_radius` being 
perturbed by an object with mass `mass_perturber` in a binary with the given `semi_major_axis`,
calculate and return the fraction of the apsidal motion constant `k` over the 
tidal timescale `T`. If the stellar type of the tidal object is a type with a radiative envelope, such as 
a core helium burning star, or a massive main sequence star, the prescription from 
Hurley et al. 2002 will be used (Eq. 42). Otherwise, the the prescription from Preece et al. 2022 is used.
"""
function apsidal_motion_constant_over_tidal_timescale(mass::Unitful.Mass, radius,
                                                      envelope_mass, envelope_radius,
                                                      stellar_type, luminosity, 
                                                      mass_perturber,
                                                      semi_major_axis, Z=0.02)::Float64

    if stellar_type isa Int
        if stellar_type > 9
            return 0.0
        end
    else
        if !(stellar_type isa Star)
            return 0.0
        end
    end

    mass = ustrip(u"Msun", mass)
    radius = ustrip(u"Rsun", radius)
    envelope_mass = ustrip(u"Msun", envelope_mass)
    envelope_radius = ustrip(u"Rsun", envelope_radius)
    luminosity = ustrip(u"Lsun", luminosity)

    mass_perturber = ustrip(u"Msun", mass_perturber)
    semi_major_axis = ustrip(u"Rsun", semi_major_axis)

    apsidal_motion_constant_over_tidal_timescale(mass, radius,
                                                 envelope_mass, envelope_radius,
                                                 stellar_type, luminosity, 
                                                 mass_perturber,
                                                 semi_major_axis, Z)
end

"""
apsidal_motion_constant_over_tidal_timescale(mass::Real, radius, core_mass, core_radius, 
                                             stellar_type, spin, luminosity, 
                                             mass_perturber,
                                             orbital_period, semi_major_axis)

If the parameters are given without units, solar units are assumed.
"""
function apsidal_motion_constant_over_tidal_timescale(mass::Real, radius,
                                                      envelope_mass, envelope_radius,
                                                      stellar_type::Int, luminosity, 
                                                      mass_perturber,
                                                      semi_major_axis, Z=0.02)::Float64

    if stellar_type > 9
        return 0.0
    end

    if mass < 1.25
        return k_over_T_convective(mass, radius, envelope_mass, envelope_radius, luminosity, Z)
    else
        return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
    end
end


# """
# apsidal_motion_constant_over_tidal_timescale(mass, radius, age, core_mass, core_radius, 
#                                              stellar_type, luminosity, 
#                                              mass_perturber,
#                                              semi_major_axis, Z=0.02)

# If the envelope structure of the star is not fixed, the age, core_mass, and core_radius 
# can be supplied instead, which will calculate the envelope mass and radius using `Syzygy.envelope_structure`.
# """
# function apsidal_motion_constant_over_tidal_timescale_core(mass, radius, age, core_mass, core_radius, 
#                                                            stellar_type, luminosity, 
#                                                            mass_perturber,
#                                                            semi_major_axis, Z=0.02)::Float64
#     return k_over_T(mass, radius, core_mass, age,
#                     core_radius, stellar_type,
#                     luminosity,
#                     mass_perturber, semi_major_axis, Z)
# end




# function k_over_T(mass::Unitful.Mass, radius, core_mass, age,
#                   core_radius, stellar_type, luminosity,
#                   mass_perturber, semi_major_axis, Z=0.02)

#     if !(stellar_type isa Star)
#         return 0.0
#     end

#     mass = ustrip(u"Msun", mass) 
#     radius = ustrip(u"Rsun", radius) 
#     core_mass = ustrip(u"Msun", core_mass) 
#     core_radius = ustrip(u"Rsun", core_radius) 
#     luminosity = ustrip(u"Lsun", luminosity) 
#     age = ustrip(u"Myr", age)

#     mass_perturber = ustrip(u"Msun", mass_perturber)
#     semi_major_axis = ustrip(u"Rsun", semi_major_axis)


#     k_over_T(mass, radius, core_mass, age,
#             core_radius, stellar_type,
#             luminosity,
#             mass_perturber, semi_major_axis, Z)
# end

# """

# 1/yr
# """
# function k_over_T(mass::Real, radius, core_mass, age,
#                   core_radius, stellar_type, luminosity,
#                   mass_perturber, semi_major_axis, Z=0.02)

#     if !(stellar_type isa Star)
#         return 0.0
#     end

#     if mass < 1.25
#         envelope_radius, envelope_mass = envelope_structure(mass, radius, core_mass, 
#                                                             core_radius, stellar_type, age)

#         return k_over_T_convective(mass, radius, envelope_mass, envelope_radius, luminosity, Z)
#     else
#         return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
#     end
# end

function k_over_T_convective(mass::Unitful.Mass, radius, envelope_mass, envelope_radius,
                             luminosity, Z=0.02)
                             
     return k_over_T_convective_(mass, radius, envelope_mass, envelope_radius, luminosity, Z)
end

function k_over_T_convective(mass::Float64, radius, envelope_mass, envelope_radius,
                             luminosity, Z=0.02)
    return k_over_T_convective_(mass, radius, envelope_mass, envelope_radius, luminosity, Z)*k_over_T_conversion_factor
end

function k_over_T_convective_(mass, radius, envelope_mass, envelope_radius,
                             luminosity, Z)
                             
    # @info " " mass radius envelope_mass envelope_radius luminosity
    R_conv, M_conv = envelope_radius, envelope_mass
    t_conv = 0.4311*(cbrt((3M_conv*R_conv^2)/luminosity)) # Preece et al. 2022

    log10z = log10(Z)
    a =  0.630*log10z + 2.72  # a(z) for envelope
    b = -0.219*log10z + 0.68  # b(z) for envelope
    c = -0.023*log10z + 0.220 # c(z) for envelope)

    return (R_conv/radius)^a*(M_conv/mass)^b*c/t_conv
end

function k_over_T_radiative(mass_tidal_object, radius_tidal_object, mass_perturber, semi_major_axis)
    E₂ = 1.592e-9*mass_tidal_object^2.84 # second-order tidal coefficient
    q₂ = mass_perturber/mass_tidal_object
    return 1.9782e4*mass_tidal_object*radius_tidal_object/semi_major_axis^5*(1 + q₂)^(5/6)*E₂
end

function get_k_interpolator(;order=(5,5), Z=0.0134, lb_multiplier=1, ub_multiplier=1)
    Z = Z == 0.02 ? 0.0134 : Z
    k_data_location = joinpath(@__DIR__, "..", "..", "deps", "tidal_evolution_constants", "grid_Z=$Z.jld2")
    
    if !isfile(k_data_location)
        @error "Tidal evolution constant data file $(split(k_data_location, "/")[end]) not found."
    end
    
    masses = JLD2.load(k_data_location, "Mass")
    logg = JLD2.load(k_data_location, "logg")
    logk2 = JLD2.load(k_data_location, "logk2")

    logm = ustrip.(u"Msun", masses) .|> log10
    logg = ustrip.(u"cm/s^2", logg)
    logk2 = logk2

    coordinates = [SA[col...] for col in (eachcol([logm logg]'))]

    lb = [minimum(logm), minimum(logg)] .* lb_multiplier
	ub = [maximum(logm), maximum(logg)] .* ub_multiplier
    interpolator = chebregression(coordinates, logk2, lb, ub, order)
    k_itp(logm, logg) = interpolator(SA[logm, logg])

    return k_itp
end

function asidal_motion_constant_interpolated(logm::Float64, logg::Float64)
    return k_interpolator(logm, logg)
end