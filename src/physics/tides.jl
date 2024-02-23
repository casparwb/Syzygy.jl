using FastChebInterp

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
                                                      semi_major_axis)::Float64
    return k_over_T(mass, radius, core_mass, age,
                    core_radius, stellar_type,
                    luminosity,
                    mass_perturber, semi_major_axis)
end

function apsidal_motion_constant_over_tidal_timescale(mass::Unitful.Mass, radius,
                                                      envelope_mass, envelope_radius,
                                                      stellar_type, luminosity, 
                                                      mass_perturber,
                                                      semi_major_axis)::Float64

    if !(stellar_types[stellar_type] isa Star)
        return 0.0
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
                                                 semi_major_axis)
end


"""
apsidal_motion_constant_over_tidal_timescale(mass, radius, core_mass, core_radius, 
                                                    stellar_type, spin, luminosity, 
                                                    mass_perturber,
                                                    orbital_period, semi_major_axis)

"""
function apsidal_motion_constant_over_tidal_timescale(mass::Real, radius,
                                                      envelope_mass, envelope_radius,
                                                      stellar_type, luminosity, 
                                                      mass_perturber,
                                                      semi_major_axis, Z=0.02)::Float64

    if !(stellar_types[stellar_type] isa Star)
        return 0.0
    end

    if mass < 1.25
        return k_over_T_convective(mass, radius, envelope_mass, envelope_radius, luminosity, Z) |> ustrip
    else
        return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis) |> ustrip
    end
end


function k_over_T(mass::Unitful.Mass, radius, core_mass, age,
                  core_radius, stellar_type,
                  luminosity,
                  mass_perturber, semi_major_axis, Z=0.02)

    if !(stellar_types[stellar_type] isa Star)
        return 0.0
    end

    mass = ustrip(u"Msun", mass) 
    radius = ustrip(u"Rsun", radius) 
    core_mass = ustrip(u"Msun", core_mass) 
    core_radius = ustrip(u"Rsun", core_radius) 
    luminosity = ustrip(u"Lsun", luminosity) 
    age = ustrip(u"Myr", age)

    mass_perturber = ustrip(u"Msun", mass_perturber)
    semi_major_axis = ustrip(u"Rsun", semi_major_axis)


    k_over_T(mass, radius, core_mass, age,
            core_radius, stellar_type,
            luminosity,
            mass_perturber, semi_major_axis, Z)
end

"""

1/yr
"""
function k_over_T(mass::Real, radius, core_mass, age,
                  core_radius, stellar_type,
                  luminosity,
                  mass_perturber, semi_major_axis, Z=0.02)

    if !(stellar_types[stellar_type] isa Star)
        return 0.0
    end

    if mass < 1.25
        envelope_radius, envelope_mass = envelope_structure(mass, radius, core_mass, 
                                                            core_radius, stellar_type, age)

        return k_over_T_convective(mass, radius, envelope_mass, envelope_radius, luminosity, Z) |> ustrip
    else
        return k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis) |> ustrip
    end
end


function k_over_T_convective(mass, radius, envelope_mass, envelope_radius,
                             luminosity, Z=0.02)
                             

        R_conv, M_conv = envelope_radius, envelope_mass
        t_conv = 0.4311*(cbrt((3M_conv*R_conv^2)/luminosity)) # Preece et al. 2022

        log10z = log10(Z)
        a =  0.630*log10z + 2.72  # a(z) for envelope
        b = -0.219*log10z + 0.68  # b(z) for envelope
        c = -0.023*log10z + 0.220 # c(z) for envelope

        return (R_conv/radius)^a*(M_conv/mass)^b*c/t_conv
end

function k_over_T_radiative(mass, radius, mass_perturber, semi_major_axis)
    E₂ = 1.592e-9*mass^2.84 # second-order tidal coefficient
    q₂ = mass_perturber/mass
    return 1.9782e4*mass*radius/semi_major_axis^5*(1 + q₂)^(5/6)*E₂
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