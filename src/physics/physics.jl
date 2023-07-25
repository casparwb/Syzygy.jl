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
        ms = [sol.structure.m[body,idx] for body in bodies]
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

"""
    envelope_radius(mass, radius, core_radius, stellar_type)

Calculate the radius of the envelope with given mass, radius, core radius, and stellar type.
Reference Hurley et al. 2002 - DOI: 10.1046/j.1365-8711.2002.05038.x
"""
function envelope_radius(mass, radius, core_radius, stellar_type)

    if any(stellar_type .== (3, 5, 6, 8, 9))
        return radius - core_radius
    else
        if mass > 1.25u"Msun"
            return 0.0u"Rsun"
        elseif mass < 0.35u"Msun"
            return radius
        else
            @warn "Envelope radius for 0.35 < M < 1.25 not yet implemented." 
            return 0.0u"Rsun"
        end
    end
end

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
