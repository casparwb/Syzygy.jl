using LinearAlgebra, StaticArrays


function centre_of_mass(positions::AbstractVector, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)

    com = @MVector zeros(eltype(positions[1]), 3)

    f = oneunit(one_over_m)
    @inbounds for (m, pos) in zip(masses, positions)
        com .+= m .* pos .* f
    end
    return SVector(com*ustrip(one_over_m))
end


function centre_of_mass(positions::AbstractMatrix, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)
    com = @MVector zeros(eltype(positions), 3)
    
    f = oneunit(one_over_m)
    @inbounds for (m, pos) in zip(masses, eachcol(positions))
        com .+= m .* pos .* f
    end
    return SVector(com*ustrip(one_over_m))
end

function centre_of_mass(positions::AbstractMatrix, masses::AbstractMatrix)
    one_over_m = 1.0/sum(masses)

    com = @MVector zeros(eltype(positions), 3) 

    f = oneunit(one_over_m)    
    @inbounds for i ∈ axes(masses, 2)
        com .+= masses[:,i] .* positions[:,i] .* f
    end
    return SVector(com*ustrip(one_over_m))
end

function centre_of_mass(binary::Binary)

    particles = get_particles_recursive(binary)

    positions = [p.position for p in particles]
    masses = [p.mass for p in particles]
    centre_of_mass(positions, masses)
end


function centre_of_mass(sol::FewBodySolution, bodies=eachindex(sol.initial_conditions.particles); 
                        tspan=nothing)
    time = sol.t
    tspan = isnothing(tspan) ? extrema(time) : tspan
    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:indices[2]
    com = Matrix{typeof(upreferred(1.0u"m"))}(undef, 3, length(indices))
    for (c, idx) in enumerate(indices)
        r = [sol.r[:,body,idx] for body in bodies]
        ms = [sol.structure.m[body,ifelse(idx == 1, 1, 2)] for body in bodies]
        com[:,c] = centre_of_mass(r, ms)
    end

    return com
end

function centre_of_mass_velocity(sol::FewBodySolution, bodies=eachindex(sol.initial_conditions.particles); 
                        tspan=nothing)
    time = sol.t
    tspan = isnothing(tspan) ? extrema(time) : tspan
    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:indices[2]
    com = Matrix{typeof(upreferred(1.0u"m/s"))}(undef, 3, length(indices))
    
    for (c, idx) in enumerate(indices)
        v = [sol.v[:,body,idx] for body in bodies]
        ms = [sol.structure.m[body,idx] for body in bodies]
        com[:,c] = centre_of_mass_velocity(v, ms)
    end

    return com
end


"""
Returns the center of mass velocity of the particles set.
The center of mass velocity is defined as the average
of the velocities of the particles, weighted by their masses.

"""
function centre_of_mass_velocity(velocities::AbstractVector, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)
    f = oneunit(one_over_m)

    vel_com = @MVector zeros(eltype(velocities[1]), 3)
    @inbounds for i in eachindex(velocities)
        vel_com .+= masses[i].*velocities[i] .* f
    end
    
    return SVector(vel_com*ustrip(one_over_m))
end

function centre_of_mass_velocity(velocities::AbstractMatrix, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)
    f = oneunit(one_over_m)
    vel_com = @MVector zeros(eltype(velocities), 3)
    @inbounds for (m, vel) in zip(masses, eachcol(velocities))
        vel_com .+= m .* vel .* f
    end
    return SVector(vel_com*ustrip(one_over_m))
end

function centre_of_mass_velocity(binary::Binary)

    particles = get_particles_recursive(binary)

    velocities = [p.velocity for p in particles]
    masses = [p.mass for p in particles]
    centre_of_mass(velocities, masses)
end

function distances(positions::AbstractMatrix)
    [norm(r) for r in eachcol(positions)]
end

function distances(pos1::AbstractMatrix, pos2::AbstractMatrix)
    [norm(r) for r in eachcol(pos2 .- pos1)]
end

function potential_energy(positions::AbstractVector, masses::AbstractVector, G)
    # U = 0.0u"J"
    U = zero((masses[1]*masses[2])/norm(positions[1] - positions[2]))
    n = length(positions)
    @inbounds for i in 1:n
        for j in i:n
            if i != j
                U -= masses[i]*masses[j]/(norm(positions[i] - positions[j]))
            end
        end
    end

    U*ifelse(U isa Quantity, upreferred(G), ustrip(upreferred(G)))
end

function potential_energy(positions::AbstractMatrix, masses, G)
    # U = 0.0
    U = zero((masses[1]*masses[2])/(norm(positions[:,1] - positions[:,2])))
    n = size(positions, 2)
    @inbounds for i in 1:n
        for j in i:n
            if i != j
                # @show i, j
                # @show U
                U -= masses[i]*masses[j]/(norm(positions[:,i] - positions[:,j]))
            end
        end
    end

    U*ifelse(U isa Quantity, upreferred(G), ustrip(upreferred(G)))
end

potential_energy(positions, masses) = potential_energy(positions, masses, 𝒢)

function potential_energy(sol::FewBodySolution)
    masses = sol.structure.m
    
    pot_energy = zeros(typeof(1.0u"J"), length(sol.t))
    @inbounds for i ∈ eachindex(sol.t)
        m = masses[:,i]
        pe = potential_energy(sol.r[:,:,i], m)
        pot_energy[i] = pe
    end

    pot_energy
end


"""
    kinetic_energy(velocity::AbstractVector, mass::Number)
    
Kinetic energy of a body with given velocity vector and mass.
"""
function kinetic_energy(velocity::AbstractVector, mass::Number)
    return 0.5*mass*norm(velocity)^2
end

"""
    kinetic_energy(velocity::AbstractVector{<:AbstractVector}, mass::Number)
    
Total kinetic energy of a bodies with velocity vectors and masses, where
velocities[1] is a Vector with the velocity components of particle 1. 
"""
function kinetic_energy(velocities::AbstractVector{<:AbstractVector}, masses::AbstractVector)
    T = 0.5*sum(masses .* norm.(velocities).^2)
end

"""
    kinetic_energy(velocity::AbstractMatrix, mass::AbstractVector)
    
Total kinetic energy of a bodies with velocity vectors and masses, where 
each column of velocities corresponds to the velocity components of each particle. 
"""
function kinetic_energy(velocities::AbstractMatrix, masses::AbstractVector)
    T = 0.5*sum(masses .* norm.(eachcol(velocities)).^2)
end

function kinetic_energy(sol::FewBodySolution)

    masses = sol.structure.m
    
    kin_energy = Vector{typeof(1.0u"J")}(undef, length(sol.t))
    @inbounds for i ∈ eachindex(sol.t)
        # v = eachcol(sol.v[:,:,i])
        m = masses[:,i]
        te = kinetic_energy(sol.v[:,:,i], m)
        kin_energy[i] = te
    end

    kin_energy
end


function total_energy(positions, velocities, masses)
    return potential_energy(positions, masses) + kinetic_energy(velocities, masses)
end


function total_energy(sol::FewBodySolution)
    return potential_energy(sol) .+ kinetic_energy(sol)
end


function specific_orbital_energy(r, v², μ, G)
    return v²/2 - μ/r
end

function total_angular_momentum(sol::FewBodySolution; step=1)
    indices = 1:step:length(sol.t)
    htot = Array{typeof(u"kg"*sol.quantities.h[1,1,1]), 2}(undef, 3, length(indices))
    for (i, idx) in enumerate(indices)
        htot[:,i] .= sum(sol.structure.m[:,idx]) .* sum(sol.quantities.h[:,:,idx], dims=2)
    end
    htot
end


function roche_radius(a, M₁, M₂)
    return a*roche_radius_fraction(M₁, M₂)
end

function roche_radius(a, q)
    return a*roche_radius_fraction(q)
end

function roche_radius_fraction(M₁, M₂)
    q = M₁/M₂

    q²³ = cbrt(q)^2

    return 0.49q²³/(0.6q²³ + log(1 + q^(1/3)))
end

function roche_radius_fraction(q::Real)
    cbrt_q² = cbrt(q)^2
    return 0.49cbrt_q²/(0.6cbrt_q² + log(1 + q^(1/3)))
end

""" 
    stellar_spin(m::T, R::T)

Return the stellar rotation of a star with mass 'm [M⊙]' and radius 'R [R⊙]', as
described by Hurley, Pols, & Tout 2000, eq 107-108.
"""
function stellar_spin(m::Quantity{<:Real, mS}, R::Quantity{<:Real, RS}) where {mS, RS}
    stellar_spin(u"Msun"(m).val, u"Rsun"(R).val)u"1/yr"
end

function stellar_spin(m::T, R::T) where T <: Real
    vᵣₒₜ = 330m^3.3/(15 + m^3.45)
    Ω = (45.35vᵣₒₜ/R)
end



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
    stellar_type = ustrip(u"stp", stellar_type)
    age = ustrip(u"Myr", age)

    return envelope_structure(mass, radius, core_mass, core_radius, stellar_type, age, Z)
end


function envelope_structure(star::Particle, age, Z=0.02)
    @assert star.structure.type isa Star "Envelope structure only relevant for stars."

    envelope_structure(star.structure.m, star.structure.R, 
                       star.structure.m_core, star.structure.R_core, 
                       star.structure.type.index, age)
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

zero_age_main_sequence_radius(M::DynamicQuantities.Quantity) = zero_age_main_sequence_radius(M / DynamicQuantities.Constants.M_sun |> dustrip)


# const main_sequence_radius_035_msun_factors = begin
#     M = 0.35
#     Z = 0.02
#     ζ = log10(Z/0.02)
#     aₙ(α, β, γ, η, μ) = α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4

#     a₁₈′ = aₙ(2.187715, -2.154437, -3.768678, -1.975518, -3.021475)
#     a₁₉′ = aₙ(1.466440, 1.839725, 6.442199, 4.023635, 6.957529)
#     a₂₀ = aₙ(2.652091, 8.178458, 1.156058, 7.633811, 1.950698)

#     a₁₈ = a₁₈′*a₂₀
#     a₁₉ = a₁₉′*a₂₀
#     a₂₁ = aₙ(1.472103, -2.947609, -3.312828, -9.945065, 0)
#     a₂₂ = aₙ(3.071048, -5.679941, -9.745523, -3.594543, 0)

#     R_zams = zero_age_main_sequence_radius(M)
#     R_tms = (a₁₈ + a₁₉*M^a₂₁)/(a₂₀ + M^a₂₂) # Hurley et al 2000 eq. 9

#     Mhook = 1.0185 + 0.16015ζ + 0.0892ζ^2
#     @assert Mhook >= M "Only valid for M <= Mhook right now."


#     a₆₂ = aₙ(8.4300, -4.7500, -3.5200, 0, 0)
#     a₇₆ = aₙ(1.192334, 1.083057, 1.230969, 1.551656, 0)
#     a₇₇ = aₙ(-1.668868, 5.818123, -1.105027, -1.668070, 0)
#     a₇₈ = aₙ(7.615495, 1.068243, -2.011333, -9.371415, 0)
#     a₇₉ = aₙ(9.409838, 1.522928, 0, 0, 0)

#     a₇₆ = max(a₇₆, -0.1015564 - 0.2161264ζ - 0.05182516*ζ^2) 
#     a₇₇ = max(-0.3868776 - 0.5457078*ζ - 0.1463472ζ^2, min(0.0, a₇₇))
#     a₇₈ = max(0.0, min(a₇₈, 7.454 + 9.046*ζ)) 
#     a₇₉ = min(a₇₉, max(2.0, -13.3 - 18.6*ζ)) 

#     αR = a₆₂
#     βR = 1.06

#     γ = a₇₆ + a₇₇*(M - a₇₈)^a₇₉
#     Dict{String, Float64}("αR" => αR, "βR" => βR, "γ" => γ, "R_tms" => R_tms,
#                           "R_zams" => R_zams)
# end

"""
    main_sequence_radius_035_msun(τ, Z=0.02)

Radius of a 0.35 M⊙ main-sequence star at a time τ = t/tMS.  
"""
function main_sequence_radius_035_msun(τ, Z=0.02)
    M = 0.35

    ζ = log10(Z/0.02)
    aₙ(α, β, γ, η, μ) = α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4

    a₁₈′ = aₙ(2.187715, -2.154437, -3.768678, -1.975518, -3.021475)
    a₁₉′ = aₙ(1.466440, 1.839725, 6.442199, 4.023635, 6.957529)
    a₂₀ = aₙ(2.652091, 8.178458, 1.156058, 7.633811, 1.950698)

    a₁₈ = a₁₈′*a₂₀
    a₁₉ = a₁₉′*a₂₀
    a₂₁ = aₙ(1.472103, -2.947609, -3.312828, -9.945065, 0)
    a₂₂ = aₙ(3.071048, -5.679941, -9.745523, -3.594543, 0)

    R_zams = zero_age_main_sequence_radius(M)
    R_tms = (a₁₈ + a₁₉*M^a₂₁)/(a₂₀ + M^a₂₂) # Hurley et al 2000 eq. 9

    Mhook = 1.0185 + 0.16015ζ + 0.0892ζ^2
    @assert Mhook >= M "Only valid for M <= Mhook right now."


    a₆₂ = aₙ(8.4300, -4.7500, -3.5200, 0, 0)
    a₇₆ = aₙ(1.192334, 1.083057, 1.230969, 1.551656, 0)
    a₇₇ = aₙ(-1.668868, 5.818123, -1.105027, -1.668070, 0)
    a₇₈ = aₙ(7.615495, 1.068243, -2.011333, -9.371415, 0)
    a₇₉ = aₙ(9.409838, 1.522928, 0, 0, 0)

    a₇₆ = max(a₇₆, -0.1015564 - 0.2161264ζ - 0.05182516*ζ^2) 
    a₇₇ = max(-0.3868776 - 0.5457078*ζ - 0.1463472ζ^2, min(0.0, a₇₇))
    a₇₈ = max(0.0, min(a₇₈, 7.454 + 9.046*ζ)) 
    a₇₉ = min(a₇₉, max(2.0, -13.3 - 18.6*ζ)) 

    αR = a₆₂
    βR = 1.06

    γ = a₇₆ + a₇₇*(M - a₇₈)^a₇₉

    logRMS_over_RZAMS = αR*τ + βR * τ^10 + γ*τ^40 + 
                            (log10(R_tms/R_zams) - αR - βR - γ)*τ^3


    10^logRMS_over_RZAMS*R_zams

    # let f = main_sequence_radius_035_msun_factors
    #     logRMS_over_RZAMS = f["αR"]*τ + f["βR"] * τ^10 + f["γ"]*τ^40 + 
    #                         (log10(f["R_tms"]/f["R_zams"]) - f["αR"] - f["βR"] - f["γ"])*τ^3


    #     10^logRMS_over_RZAMS*R_zams
    # end
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

    aₙ(α, β, γ, η, μ) = α + β*ζ + γ*ζ^2 + η*ζ^3 + μ*ζ^4

    a₁ = aₙ(1.593890, 2.053038, 1.231226, 2.327785, 0.0)
    a₂ = aₙ(2.706708, 1.483131, 5.772723, 7.411230, 0.0)
    a₃ = aₙ(1.466143, -1.048442, -6.795374, -1.391127, 0.0)
    a₄ = aₙ(4.141960, 4.564888, 2.958542, 5.571483, 0.0)
    a₅ = aₙ(3.426349, 0.0, 0.0, 0.0, 0.0)
    a₆ = aₙ(1.949814, 1.758178, -6.008212, -4.470533, 0.0)
    a₇ = aₙ(4.903830, 0.0, 0.0, 0.0, 0.0)
    a₈ = aₙ(5.212154, 3.166411, -2.750074, -2.271549, 0.0)
    a₉ = aₙ(1.312179, -3.294936, 9.231860, 2.610989, 0.0)
    a₁₀ = aₙ(8.073972, 0.0, 0.0, 0.0, 0.0)

    μ = max(0.5, 1.0 - 0.01*max(a₆/M^a₇, a₈ + a₉/M^a₁₀))
    x = max(0.95, min(0.95 - 0.03*(ζ + 0.30103), 0.99))

    tBGB = (a₁ + a₂*M^4 + a₃*M^5.5 + M^7)/(a₄*M^2 + a₅*M^7)
    t_hook = μ*tBGB

    tMS = max(t_hook, x*tBGB)

    return tMS, tBGB
end

main_sequence_lifetime(M::Unitful.Quantity, Z=0.02) = main_sequence_lifetime(ustrip(u"Msun", M), Z)
main_sequence_lifetime(M::DynamicQuantities.Quantity, Z=0.02) = main_sequence_lifetime(M / DynamicQuantities.Constants.M_sun |> dustrip, Z)


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


"""
    stability_criterion_ma01(m1, m2, m3, i, eout)

Return the critical semi-major axis ratio (a_out/a_in)_crit as defined in
Mardling & Aarset 1999.
"""
function stability_criterion_ma01(m1, m2, m3, i, eout)
    qout = m3/(m1 + m2)
    return 2.8/(1 - eout)*(1 - 0.3i/π)*((1 + qout)*(1 + eout)/√(1 - eout))^(2/5)
end

function stability_criterion_ma01(p::MultiBodySystem)
    @assert p.n == 3 "System must be a triple."
    eₒ = p.binaries[2].elements.e
    i = p.binaries[1].elements.i

    m1, m2, m3 = [p.particles[i].mass for i = 1:3]

    return stability_criterion_ma01(m1, m2, m3, i, eₒ)    
end

function get_a_out_on_stability_limit(aᵢ, m₁, m₂, m₃, i, eₒ; ϵ=1e-4)
    stability = stability_criterion_ma01(m₁, m₂, m₃, i, eₒ)
    aₒ = stability*aᵢ*(1-ϵ)
end


function is_unstable(p::MultiBodySystem; criterion="ma01")
    α_crit = stability_criterion_ma01(p)
    α = p.binaries[2].elements.a / p.binaries[1].elements.a
    return α < α_crit
end

function is_unstable(aₒ, aᵢ, m₁, m₂, m₃, i, eₒ; criterion="ma01")
    α_crit = stability_criterion_ma01(m₁, m₂, m₃, i, eₒ)
    α = aₒ/aᵢ
    return α < α_crit
end

""" 
Octupole term ϵₒ. Plays important role when 

`` |ϵₒ| ≥ {0.001, 0.01} `` 
"""
function octupole_parameter(triple::MultiBodySystem)
    @assert triple.n == 3 "Octupole parameter only valid for triple system."
    m₁, m₂ = [triple.particles[i].structure.m for i in 1:2]#[triple.mass[1], triple.mass[2]] |> sort |> reverse
    aᵢₙ  = triple.binaries[1].elements.a
    aₒᵤₜ = triple.binaries[2].elements.a
    eₒᵤₜ = triple.binaries[2].elements.e

    (m₁ - m₂)/(m₁ + m₂)*aᵢₙ/aₒᵤₜ*eₒᵤₜ/(1 - eₒᵤₜ^2)
end

function quadrupole_timescale(nbody::MultiBodySystem)
    m = nbody.m .|> u"kg"
    Pi = nbody.elements[1].P |> u"s"
    Po = nbody.elements[2].P |> u"s"
    eo = nbody.elements[2].e
    return 16/30π*sum(m)/m[3]*Po^2/Pi*cbrt(1 - eo^2)^2
end
