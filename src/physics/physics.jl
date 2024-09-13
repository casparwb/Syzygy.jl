using LinearAlgebra: norm, ×
using StaticArrays


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


function centre_of_mass(sol::MultiBodySolution, bodies=eachindex(sol.ic.particles); 
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

function centre_of_mass_velocity(sol::MultiBodySolution, bodies=eachindex(sol.ic.particles); 
                        tspan=nothing)
    time = sol.t
    tspan = isnothing(tspan) ? extrema(time) : tspan
    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:indices[2]
    com = Matrix{typeof(upreferred(1.0u"m/s"))}(undef, 3, length(indices))
    
    for (c, idx) in enumerate(indices)
        v = [sol.v[:,body,idx] for body in bodies]
        ms = [sol.structure.m[body,ifelse(idx == 1, 1, 2)] for body in bodies]
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

potential_energy(positions, masses) = potential_energy(positions, masses, GRAVCONST)

function potential_energy(sol::MultiBodySolution)
    masses = sol.structure.m
    
    pot_energy = zeros(typeof(1.0u"J"), length(sol.t))
    @inbounds for i ∈ eachindex(sol.t)
        m = ifelse(i > 1, masses[:,2], masses[:,1])
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

function kinetic_energy(sol::MultiBodySolution)

    masses = sol.structure.m
    
    kin_energy = Vector{typeof(1.0u"J")}(undef, length(sol.t))
    @inbounds for i ∈ eachindex(sol.t)
        m = ifelse(i > 1, masses[:,2], masses[:,1])
        te = kinetic_energy(sol.v[:,:,i], m)
        kin_energy[i] = te
    end

    kin_energy
end


function total_energy(positions, velocities, masses)
    return potential_energy(positions, masses) + kinetic_energy(velocities, masses)
end


function total_energy(sol::MultiBodySolution)
    return potential_energy(sol) .+ kinetic_energy(sol)
end


function specific_orbital_energy(r, v², μ)
    return v²/2 - μ/r
end

# function angular_momentum(sol::MultiBodySolution; step=1)
#     indices = 1:step:length(sol.t)
#     htot = Array{typeof(u"m^2/s"), 2}(undef, 3, length(indices))
#     for (i, idx) in enumerate(indices)
#         htot[:,i] .= angular_momentum(sol.r[])
#     end
#     htot
# end

function reduced_mass(m1, m2)
    return (m1*m2)/(m1 + m2)
end

function gravitational_radius(M::Unitful.Mass)
    2*GRAVCONST*M/(c²*unit(c^2))
end

function gravitational_radius(M::Real)
    2*G*M/c²
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
function stellar_spin(m::Unitful.Mass, R::Unitful.Length)
    stellar_spin(ustrip(u"Msun", m), ustrip(u"Rsun", R))*upreferred(1.0u"1/yr")
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

main_sequence_lifetime(M::Unitful.Quantity, Z=0.02) = main_sequence_lifetime(ustrip(u"Msun", M), Z)


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

function quadrupole_timescale(system::MultiBodySystem)
    m = system.particles.mass .|> u"kg"
    P_in = system.binaries[1].elements.P |> u"s"
    P_out = system.binaries[2].elements.P |> u"s"
    e_out = system.binaries[2].elements.e
    return 16/30π*sum(m)/m[3]*P_out^2/P_in*cbrt(1 - e_out^2)^2
end

"""
    PN1_energy(r1, r2, v1, v2, m1, m2; G=GRAVCONST)

    
Total energy of body 1 in a gravitational + PN1 potential. From Blanchet 2014.
"""
function PN1_energy(r1, r2, v1, v2, m1, m2; G=GRAVCONST)

    r = r1 - r2

    v1_norm = norm(v1)

    r_norm = norm(r)

    n = r/r_norm

    E = G^2*m1^2*m2/(2*r_norm^2) + 3*m1*v1_norm^4/8 + G*m1*m2/r_norm*(-0.25*dot(n, v1)*dot(n, v2) + 3/2*v1_norm^2 - 7/4*dot(v1, v2))
    Ekin = m1/2*v1_norm^2
    Epot = -G*m1*m2/(2*r_norm)
    return E*c⁻² + Ekin + Epot
end

function PN1_energy(sol::MultiBodySolution)

    n_bodies = sol.ic.n

    Etot = Vector{typeof(1.0u"J")}(undef, length(sol.t))

    @inbounds for idx in eachindex(sol.t)
        E = 0.0u"J"
        for i = 1:n_bodies
            # ri = sol.r[particle=i][:,idx]
            # vi = sol.v[particle=i][:,idx]

            ri = @SVector [sol.r[1, i, idx], sol.r[2, i, idx], sol.r[3, i, idx]]
            vi = @SVector [sol.v[1, i, idx], sol.v[2, i, idx], sol.v[3, i, idx]]

            mi = sol.structure.m[i,2]
            for j = 1:n_bodies
                if j != i
                    # rj = sol.r[particle=j][:,idx]
                    # vj = sol.v[particle=j][:,idx]

                    rj = @SVector [sol.r[1, j, idx], sol.r[2, j, idx], sol.r[3, j, idx]]
                    vj = @SVector [sol.v[1, j, idx], sol.v[2, j, idx], sol.v[3, j, idx]]

                    mj = sol.structure.m[j,2]
                    E += PN1_energy(ri, rj, vi, vj, mi, mj)
                end
            end
        end
        Etot[idx] = E
    end

    return Etot
end

function deSitter_factor(binary)
    a = binary.elements.a
    e = binary.elements.e
    m1, m2 = [p.mass for p in binary.children]
    μ = reduced_mass(m1, m2)
    n = √(G*(m1 + m2)/a^3)

    c2 = c²*unit(c)^2
    return 3G*n*(m2 + μ/3)/(2*c2*a*(1 - e^2))
end

function deSitter_spin_velocity(particle, parent_binary)
    if particle.sibling isa BinaryIndex
        return zeros(eltype(particle.structure.S)/oneunit(upreferred(1.0u"s")), 3)
    end
    
    Ωds = deSitter_factor(parent_binary)

    sibling = parent_binary.children[particle.sibling.i]
    r1 = particle.position
    r2 = sibling.position

    v1 = particle.velocity
    v2 = sibling.velocity

    r = r2 - r1
    v = v2 - v1

    m1, m2 = particle.mass, sibling.mass
    μ = reduced_mass(m1, m2)

    L = angular_momentum(r, μ*v)
    L̂ = L/norm(L)
    S = particle.structure.S
    # Ŝ = S/norm(S)

    return Ωds*L̂ × S#Ŝ

end

function precession_vector(r1, r2, v1, v2, m1, m2)

    M = m1 + m2
    X1 = m1/M
    X2 = m2/M

    Δ = X1 - X2
    ν = X1*X2

    r̄ = r1 - r2
    v̄ = v1 - v2

    r = norm(r̄)
    v² = norm(v̄)^2

    GM = GRAVCONST*M
    aDen = 2GM - v²*r
    a = GM*r/aDen
    # @show 
    GM_a³ = GM/a^3
    Ω = √GM_a³
    x = (GM*Ω/c^3)^(2/3)

    n̄ = r̄/r

    nxv = n̄ × v̄
    𝓁 = nxv/norm(nxv)


    # num = (0.75 + 0.5ν - 0.75*Δ).val 
    # num += x*(9/16 + 5/4*ν - 1/24*ν^2 + Δ*(-9/16 + 5/8*ν))
    # num += x^2*(27/32 + 3/16*ν - 105/32*ν^2 - 1/48*ν^3 + 
    #             Δ*(-27/32 + 39/8*ν - 5/32*ν^2))
    # Ω₁ = c³*x^(5/2)/(G*M)*𝓁*num

    Ω₁ = c^3*x^(5/2)/(G*M)*𝓁*(0.75 + 0.5ν - 0.75*Δ + 
                              x*(9/16 + 5/4*ν - 1/24*ν^2 + Δ*(-9/16 + 5/8*ν)) +
                              x^2*(27/32 + 3/16*ν - 105/32*ν^2 - 1/48*ν^3 + 
                                   Δ*(-27/32 + 39/8*ν - 5/32*ν^2)
                                  )
                             )

    return Ω₁
end

function spin_precession_velocity(particle1::Particle, particle2::Particle)
    S1 = particle1.structure.S
    r1 = particle1.position
    r2 = particle2.position
    v1 = particle1.velocity
    v2 = particle2.velocity

    m1 = particle1.mass
    m2 = particle2.mass
    spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
end

function spin_precession_velocity(body1::MassBody, body2::MassBody)
    S1 = body1.spin
    r1 = body1.position
    r2 = body2.position
    v1 = body1.velocity
    v2 = body2.velocity

    m1 = body1.mass
    m2 = body2.mass
    spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
end

function spin_precession_velocity(S1, r1, r2, v1, v2, m1::Quantity, m2::Quantity)
    T1PN, T2PN = spin_precession_velocity_factor(S1, r1, r2, v1, v2, m1, m2)
    return GRAVCONST*(T1PN/c^2 + T2PN/c^4)
end

function spin_precession_velocity(S1, r1, r2, v1, v2, m1::AbstractFloat, m2::AbstractFloat)
    T1PN, T2PN = spin_precession_velocity_factor(S1, r1, r2, v1, v2, m1, m2)
    return G*(T1PN/c² + T2PN/c⁴)
end

function spin_precession_velocity_factor(S1, r1, r2, v1, v2, m1, m2)
    
    r̄ = r1 - r2
    v̄ = v1 - v2

    r = norm(r̄)

    n = r̄/r
    nS1 = dot(n, S1)
    nv = dot(n, v̄)
    nv1 = dot(n, v1)
    nv2 = dot(n, v2)
    vv2 = dot(v̄, v2)
    vS1  = dot(v̄, S1)
    v1S1 = dot(v1, S1)
    v2S1 = dot(v2, S1)

    T1PN = m2/r^2*((v1 - 2*v2)*nS1 + S1*nv - 2*n*vS1)

    T2PN = m2/r^2*((G*m1*(-16*nS1*nv + 3*v1S1 - 7*v2S1)/r + 
                      2*G*m2*nS1*nv/r + (3*nv2^2 + 2*vv2)*vS1)*n + 
                     (-5*G*(m1 - m2)*nS1/r + (3*nv2^2 + 2*vv2)*nS1 + 
                      2*(v1S1 + v2S1)*nv)*v2 - 
                     (-G*(6*m1 - m2)*nS1/r + 3*nS1*nv2^2/2 + nv2*vS1)*v1 + 
                     (G*m1*nv1/r - G*m2*nv/r - 3*nv*nv2^2/2 + nv2*vv2)*S1
                    )

    T1PN, T2PN
end