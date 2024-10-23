function envelope_structure(mass::Real, radius, core_mass, core_radius, stellar_type, age, Z=0.02)
    tMS, tBGB = main_sequence_lifetime(mass, Z)
    envelope_radius = convective_envelope_radius(mass, radius, core_radius, stellar_type, age, tMS, tBGB)
    envelope_mass = convective_envelope_mass(mass, core_mass, stellar_type, age, tMS, tBGB)

    return envelope_radius, envelope_mass
end

function envelope_structure(mass::Unitful.Mass, radius, core_mass, core_radius, stellar_type, age, Z=0.02)
    
    mass = ustrip(u"Msun", mass)
    radius = ustrip(u"Rsun", radius)
    core_mass = ustrip(u"Msun", core_mass)
    core_radius = ustrip(u"Rsun", core_radius)
    # stellar_type = ustrip(u"stp", stellar_type)
    age = ustrip(u"Myr", age)

    R_env, M_env = envelope_structure(mass, radius, core_mass, core_radius, stellar_type, age, Z)
    return R_env*u"Rsun", M_env*u"Msun"
end


function envelope_structure(star::Particle, age, Z=0.02)
    @assert star.structure.stellar_type isa Star "Envelope structure only relevant for stars."

    envelope_structure(star.structure.m, star.structure.R, 
                       star.structure.m_core, star.structure.R_core, 
                       star.structure.stellar_type.index, age)
end


"""
Radius of a zero-age main-sequence star. From Tout et al 1996.
"""
function zero_age_main_sequence_radius(M::Real)
    θ = 1.71535900
    ι = 6.59778800
    κ = 10.08855000
    λ = 1.01249500
    μ = 0.07490166
    ν = 0.01077422
    ξ = 3.08223400
    o = 17.84778000
    Π = 0.00022582


    (θ*M^2.5 + ι*M^6.5 + κ*M^11 + λ*M^19 + μ*M^19.5)/(ν + ξ*M^2 + o*M^8.5 + M^18.5 + Π*M^19.5)
end

function zero_age_main_sequence_radius(mass::Unitful.Mass)
    zero_age_main_sequence_radius(ustrip(u"Msun", mass))
end



"""
    main_sequence_radius_035_msun(τ, Z=0.02)

Radius of a 0.35 M⊙ main-sequence star at a time τ = t/tMS.  
"""
function main_sequence_radius_035_msun(τ::Real, Z=0.02)
    M = 0.35
    ζ = log10(Z/0.02)
    ζ² = ζ^2    
    aₙ(α, β=0.0, γ=0.0, η=0.0, μ=0.0) = α + β*ζ + γ*ζ² + η*ζ^3 + μ*ζ²^2

    a₁₈′ = aₙ(2.187715e-1, -2.154437e+0, -3.768678e+0, -1.975518e+0, -3.021475e-1)
    a₁₉′ = aₙ(1.466440e+0, 1.839725e+0, 6.442199e+0, 4.023635e+0, 6.957529e-1)
    a₂₀ = aₙ(2.652091e+1, 8.178458e+1, 1.156058e+2, 7.633811e+1, 1.950698e+1)

    a₁₈ = a₁₈′*a₂₀
    a₁₉ = a₁₉′*a₂₀
    a₂₁ = aₙ(1.472103e+0, -2.947609e+0, -3.312828e+0, -9.945065e-1)
    a₂₂ = aₙ(3.071048e+0, -5.679941e+0, -9.745523e+0, -3.594543e+0)

    R_zams = zero_age_main_sequence_radius(M)
    R_tms = (a₁₈ + a₁₉*M^a₂₁)/(a₂₀ + M^a₂₂) # Hurley et al 2000 eq. 9

    Mhook = 1.0185 + 0.16015ζ + 0.0892ζ² 
    @assert Mhook >= M "Only valid for M <= Mhook right now."


    a₆₂ = aₙ(8.4300e-2, -4.7500e-2, -3.5200e-2)
    a₇₆ = aₙ(1.192334e-2, 1.083057e-2, 1.230969e+0, 1.551656e+0)
    a₇₇ = aₙ(-1.668868e-1, 5.818123e-1, -1.105027e+1, -1.668070e+1)
    a₇₈ = aₙ(7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2)
    a₇₉ = aₙ(9.409838e+0, 1.522928e+0)

    a₇₆ = max(a₇₆, -0.1015564 - 0.2161264*ζ - 0.05182516*ζ²) 
    a₇₇ = max(-0.3868776 - 0.5457078*ζ - 0.1463472*ζ², min(0.0, a₇₇))
    a₇₈ = max(0.0, min(a₇₈, 7.454 + 9.046*ζ)) 
    a₇₉ = min(a₇₉, max(2.0, -13.3 - 18.6*ζ)) 

    αR = a₆₂
    βR = 1.06

    γ = a₇₆ + a₇₇*(M - a₇₈)^a₇₉

    logRMS_over_RZAMS = αR*τ + βR * τ^10 + γ*τ^40 + 
                            (log10(R_tms/R_zams) - αR - βR - γ)*τ^3


    10^logRMS_over_RZAMS*R_zams
end


"""
    envelope_radius(mass, radius, core_radius, stellar_type)

Calculate the radius of the envelope with given mass, radius, core radius, and stellar type.
Quantities must be in units of solar mass and solar radii.
Reference Hurley et al. 2002 - DOI: 10.1046/j.1365-8711.2002.05038.x
"""
function convective_envelope_radius(mass, radius, core_radius, stellar_type, age, tMS, tBGB)

    if any(stellar_type .== (3, 5, 6, 8, 9)) # giant-like stars
        return radius - core_radius
    elseif any(stellar_type .== (1, 7))   # main sequence stars
        τ = age/tMS 
        
        R_env₀ = if mass > 1.25
                    0.0
                elseif mass < 0.35
                    radius
                else

                    R′ = main_sequence_radius_035_msun(τ)
                    # R′ is the radius of a MS star with M = 0.35 M⊙ at τ
                    return R′*sqrt(1.25 - mass)/0.9
                end

        return R_env₀*(1 - τ)^0.25
    elseif any(stellar_type .== (2, 8)) # Hertzsprung gap stars
        τ = (age - tMS)/(tBGB - tMS)
        return sqrt(τ)*(radius - core_radius)
    end

end

"""
convective_envelope_mass(mass, radius, core_radius, stellar_type)

Calculate the mass of the envelope with given stellar mass, core mass, stellar age, 
stellar main sequence lifetime, stellar base giant branch (BHG) lifetime and stellar type.
Quantities must be in units of solar mass and solar radii.
Reference Hurley et al. 2000 - https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H
"""
function convective_envelope_mass(mass, core_mass, stellar_type, age, tMS, tBGB)
    @assert stellar_types[stellar_type] isa Star "Only stars have envelopes."

    if any(stellar_type .== (1, 7)) 
        M_env₀ = if mass < 0.35
                     mass
                 elseif mass > 1.25
                     0.0
                 else
                   ( 0.35*((1.25 - mass)/0.9)^2 )
                 end
        
        τ = age/tMS
        return M_env₀*(1 - τ)^0.25
    elseif any(stellar_type .== (2, 8))
        τ = (age - tMS)/(tBGB - tMS)
        return τ*(mass - core_mass)
    else 
        return mass - core_mass
    end
end 

"""

main_sequence_lifetime(M::Real, Z)

Return the main sequence lifetime of a star with mass M [M⊙] in Myr.
Reference Hurley et al. 2000 - https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H
"""
function main_sequence_lifetime(M::Real, Z=0.02)

    ζ = log10(Z/0.02) # Hurley et al 2000 page 5

    aₙ(α, β=0.0, γ=0.0, η=0.0, μ=0.0) = α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4

    a₁ = aₙ(1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2)
    a₂ = aₙ(2.706708e3, 1.483131e3, 5.772723e2, 7.411230)
    a₃ = aₙ(1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1)
    a₄ = aₙ(4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3)
    a₅ = aₙ(3.426349e-1)
    a₆ = aₙ(1.949814e1, 1.758178, -6.008212, -4.470533)
    a₇ = aₙ(4.903830)
    a₈ = aₙ(5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3)
    a₉ = aₙ(1.312179, -3.294936e-1, 9.231860e-2, 2.610989e-2)
    a₁₀ = aₙ(8.073972e-1)

    M⁷ = M^7
    μ = max(0.5, 1.0 - 0.01*max(a₆/M^a₇, a₈ + a₉/M^a₁₀))
    x = max(0.95, min(0.95 - 0.03*(ζ + 0.30103), 0.99))

    tBGB = (a₁ + a₂*M^4 + a₃*M^5.5 + M⁷)/(a₄*M^2 + a₅*M⁷)
    t_hook = μ*tBGB

    tMS = max(t_hook, x*tBGB)

    return tMS, tBGB
end

main_sequence_lifetime(M::Unitful.Mass, Z=0.02) = main_sequence_lifetime(ustrip(u"Msun", M), Z)


function mass_luminosity_relation(M)

    # M = u"Msun"(M).val

    local F, a
    if 0.2u"Msun" < M <= 0.85u"Msun"
        a = -141.7*M^4 + 232.4*M^3 - 129.1*M^2 + 33.29*M + 0.215
        F = 1.0
    elseif 0.85u"Msun" < M <= 2.0u"Msun"
        a = 4
        F = 1
    elseif 2.0u"Msun" < M <= 55.0u"Msun"
        a = 3.5
        F = 1.4
    elseif M > 55.0u"Msun"
        a = 1
        F = 32_000
    end


    return (F*M^a)u"Lsun"
end
