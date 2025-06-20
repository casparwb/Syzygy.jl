using LinearAlgebra: norm, ×, ⋅
using StaticArrays

include("stellar_physics.jl")

"""
    centre_of_mass(positions, masses)

Return the centre of mass position given a collection of `positions` and `masses`. 

`positions` can either be a `Vector` of `Vector`s, or a `Matrix` in which each column
represents the position of each body.

# Examples
```julia-repl
julia> r = [rand(3) for i =1:3]
julia> m = rand(3)

julia> com = Syzygy.centre_of_mass(r, m)
```
"""
function centre_of_mass(positions::AbstractVector, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)

    com = @MVector zeros(eltype(positions[1]), 3)

    f = oneunit(one_over_m)
    @inbounds for i in eachindex(masses)
        m = masses[i]
        pos = positions[i]
        for k = 1:3
            com[k] += m*pos[k]*f
        end
    end
    return SVector(com*ustrip(one_over_m))
end


function centre_of_mass(positions::AbstractMatrix, masses::AbstractVector)
    one_over_m = 1.0/sum(masses)
    com = @MVector zeros(eltype(positions), 3)
    
    f = oneunit(one_over_m)
    @inbounds for i in eachindex(masses)# in zip(masses, eachcol(positions))
        for k = 1:3
            com[k] += masses[i]*positions[k,i]*f
        end
    end
    return SVector(com*ustrip(one_over_m))
end

function centre_of_mass(positions::AbstractMatrix, masses::AbstractMatrix)
    one_over_m = 1.0/sum(masses)

    com = @MVector zeros(eltype(positions), 3) 

    f = oneunit(one_over_m)    
    @inbounds for i ∈ axes(masses, 2)
        for k = 1:3
            com[k] += masses[k,i]*positions[k,i]*f
        end
    end
    return SVector(com*ustrip(one_over_m))
end

function center_of_mass(positions, masses)
    return centre_of_mass(positions, masses)
end

"""
    centre_of_mass(sol::MultiBodySolution, 
                   bodies=eachindex(sol.ic.particles); 
                   tspan=extrema(sol.t))

Return the centre of mass (COM) position given a `MultiBodySolution` object. 

By default, the COM position calculated for each body for each time step in `sol` is returned.
Specific bodies can be specified with their indices, and a specific time span can
be set with a tuple of `(t_start, t_end)`.

# Examples
```julia-repl
julia> com = Syzygy.centre_of_mass(sol)
julia> com12 = Syzygy.centre_of_mass(sol, [1, 2])
julia> # just get the COM position for the first two time steps:
julia> com_tspan = Syzygy.centre_of_mass(sol, tspan=(sol.t[1], sol.t[2])) 
```
"""
function centre_of_mass(sol::MultiBodySolution, 
                        bodies=eachindex(sol.ic.particles); 
                        tspan=extrema(sol.t))
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


"""
 centre_of_mass_velocity(velocities, masses)

Return the center of mass velocity.
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

function centre_of_mass_velocity(sol::MultiBodySolution, 
                                 bodies=eachindex(sol.ic.particles); 
                                 tspan=extrema(sol.t))
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
    potential_energy(positions, masses)

Return the total potential energy of bodies with `positions` and `masses`. 

`positions` can either be a `Vector` of `Vector`s, or a `Matrix` in which each column
represents the position of each body.

# Examples
```julia-repl
julia> r = [rand(3) for i =1:3]
julia> m = rand(3)

julia> com = Syzygy.potential_energy(r, m)
```
"""
function potential_energy(positions::AbstractVector, masses)
    U = zero((masses[1]*masses[2])/norm(positions[1] - positions[2]))
    n = length(positions)
    @inbounds for i in 1:n
        for j in i:n
            if i != j
                U -= masses[i]*masses[j]/(norm(positions[i] - positions[j]))
            end
        end
    end

    U*ifelse(U isa Quantity, GRAVCONST, UNITLESS_G)
end

function potential_energy(positions::AbstractMatrix, masses)
    U = zero((masses[1]*masses[2])/(norm(positions[:,1] - positions[:,2])))
    n = size(positions, 2)
    @inbounds for i in 1:n
        for j in i:n
            if i != j
                U -= masses[i]*masses[j]/(norm(positions[:,i] - positions[:,j]))
            end
        end
    end

    U*ifelse(U isa Quantity, GRAVCONST, UNITLESS_G)
end

"""
    potential_energy(sol::MultiBodySolution, 
                     bodies=keys(sol.ic.particles); 
                     tspan=extrema(sol.t))

Return the total potential energy for the given `bodies` and `tspan`.
"""
function potential_energy(sol::MultiBodySolution, 
                          bodies=keys(sol.ic.particles); 
                          tspan=extrema(sol.t))
    
    masses = sol.structure.m
    time = sol.t
    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:indices[2]
    
    pot_energy = zeros(typeof(upreferred(1.0u"J")), length(sol.t))
    for i ∈ indices
        m = ifelse(i > 1, masses[:,2], masses[:,1])
        r = [SVector{3}(sol.r[:,b,i]) for b in bodies]
        pe = potential_energy(r, m)
        pot_energy[i] = pe
    end

    pot_energy
end


"""
    kinetic_energy(velocity::AbstractVector, mass::Number)
    
Return the kinetic energy of a body with given velocity vector and mass.
"""
function kinetic_energy(velocity::AbstractVector, mass::Number)
    return 0.5*mass*norm(velocity)^2
end

"""
    kinetic_energy(velocity::AbstractVector{<:AbstractVector}, mass::Number)
    
Total kinetic energy of bodies with velocity vectors and `masses`, where
`velocities[1]` is a `Vector` with the velocity components of particle 1. 
"""
function kinetic_energy(velocities::AbstractVector{<:AbstractVector}, masses::AbstractVector)
    T = 0.5*sum(masses .* norm.(velocities).^2)
end

"""
    kinetic_energy(velocity::AbstractMatrix, mass::AbstractVector)
    
`velocity` can also be a `Matrix` where each column represent the velocity vector of each body. 
"""
function kinetic_energy(velocities::AbstractMatrix, masses::AbstractVector)
    T = 0.5*sum(masses .* norm.(eachcol(velocities)).^2)
end

"""
    kinetic_energy(sol::MultiBodySolution, 
                   bodies=keys(sol.ic.particles); 
                   tspan=extrema(sol.t))

Return the total kinetic energy for the given `bodies` and `tspan` from a simulation solution.
"""
function kinetic_energy(sol::MultiBodySolution, 
                        bodies=keys(sol.ic.particles); 
                        tspan=extrema(sol.t))

    masses = sol.structure.m
    time = sol.t

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:indices[2]
    
    kin_energy = Vector{typeof(1.0u"J")}(undef, length(sol.t))
    for i ∈ indices
        m = ifelse(i > 1, masses[:,2], masses[:,1])
        v = [SVector{3}(sol.v[:,b,i]) for b in bodies]
        te = kinetic_energy(v, m)
        kin_energy[i] = te
    end

    kin_energy
end

"""
    total_energy(positions, velocities, masses)
    
Total energy (kinetic + potential) of bodies with given `positions`, `velocities`, and `masses`.
"""
function total_energy(positions, velocities, masses)
    return potential_energy(positions, masses) + kinetic_energy(velocities, masses)
end

function total_energy(sol::MultiBodySolution, 
                      bodies=eachindex(sol.ic.particles); 
                      tspan=extrema(sol.t))
    return potential_energy(sol, bodies, tspan=tspan) .+ kinetic_energy(sol, bodies, tspan=tspan)
end

"""
    specific_orbital_energy(r, v², μ)

Return the specific orbital energy of two orbiting bodies with relative position `r`, relative velocity squared `v²`, and reduced mass `μ`

The specific orbital energy is sum of the mutual potential and kinetic energy divided by the reduced mass.
"""
function specific_orbital_energy(r, v², μ)
    return v²/2 - μ/r
end

@doc raw"""
    reduced_mass(m₁, m₂)

Return the reduced mass given the two masses `m₁` and `m₂`.

The reduced mass is defined as

```math
\frac{m_1 m_2}{m_1 + m_2}
```
"""
function reduced_mass(m₁, m₂)
    return (m₁*m₂)/(m₁ + m₂)
end

"""
    gravitational_radius(mass::Unitful.Mass)

Return the gravitational radius of an object with the given `mass`.

The gravitational radius is defined as GM/c².
"""
function gravitational_radius(mass::Unitful.Mass)
    GRAVCONST*mass/speed_of_light^2
end

function gravitational_radius(mass::Real)
    UNITLESS_G*mass/c²
end

"""
    schwarzschild_radius(mass::Unitful.Mass)

Return the Schwarzschild radius of an object with the given `mass`.

The Schwarzschild radius is defined as 2GM/c².
"""
function schwarzschild_radius(mass)
    return 2*gravitational_radius(mass)
end

@doc raw"""
    roche_radius(a, m₁, m₂)

Return the volume-equivalent Eggleton Roche lobe radius of two bodies
with masses `m₁` and `m₂`, and separation `a`, where `m₁ > m₂`.

The volume-equivalent Eggleton Roche radius is defined as

```math
R_L = a \frac{0.49q^{2/3}}{0.6q^{2/3} + \ln{(1 + q^{1/3})}}
```
"""
function roche_radius(a, m₁, m₂)
    return a*roche_radius_fraction(m₁, m₂)
end

"""
    roche_radius(a, q)

The mass ratio `q` can also be given directly, where `q = m₁/m₂`.
"""
function roche_radius(a, q)
    return a*roche_radius_fraction(q)
end

function roche_radius_fraction(m₁, m₂)
    q = m₁/m₂
    return roche_radius_fraction(q)
end

function roche_radius_fraction(q::Real)
    cbrt_q = cbrt(q)
    cbrt_q² = cbrt_q^2
    return 0.49cbrt_q²/(0.6cbrt_q² + log(1 + cbrt_q))
end

""" 
    stellar_rotational_frequency(m::Unitful.Mass, R::Unitful.Length)

Return the stellar rotation frequency of a star with mass 'm [M⊙]' and radius 'R [R⊙]', as
described by Hurley, Pols, & Tout 2000, eq 107-108.
"""
function stellar_rotational_frequency(m::Unitful.Mass, R::Unitful.Length)
    m = ustrip(u"Msun", m)
    v_rot = ((330m^3.3)/(15 + m^3.45))#u"km/s"
    Ω = (45.35v_rot/ustrip(u"Rsun", R))
    return u"yr^-1"*(Ω)
end

""" 
    stellar_rotational_frequency(m::Real, R::Real)

If `m` and `R` are given without units, they must be [`M⊙`] and [`R⊙`] respectively.
"""
function stellar_rotational_frequency(m::Real, R::Real)
    v_rot = (330m^3.3)/(15 + m^3.45) # km/s
    Ω = (45.35v_rot/R)
    return Ω
end

function stellar_angular_momentum(star::Particle, black_hole_spin_parameter=0.9)
    star_type = star.stellar_type
    mass = star.mass
    radius = star.radius

    k₂ = 0.1
    k₃ = 3/4
    if star_type isa Star
        Ω_zams = stellar_rotational_frequency(mass, radius)
        core_mass = star.core_mass
        return (k₂*(mass - core_mass)*radius^2 + k₃*core_mass*core_radius^2)*Ω_zams
    elseif (star_type isa WhiteDwarf) || (star_type isa NeutronStar)
        return k₃*mass*radius^2*Ω_zams
    elseif star_type isa BlackHole
        return (GRAVCONST*mass^2/speed_of_light)*black_hole_spin_parameter
    end
end


""" 
    compact_object_spin(m::Unitful.Mass, χ=0.1)

Return the spin angular momentum of a compact object with mass 'm [M⊙]' and spinning with a
fraction `χ` of its critical spin velocity.
"""
function compact_object_spin(m::Unitful.Mass, χ=1.0)
    return (GRAVCONST*m^2/speed_of_light)*χ
end

function compact_object_spin(m::Real, χ=1.0)
    return (UNITLESS_G*m^2/UNITLESS_c)*χ
end

"""
    stability_criterion_ma01(m₂, m₂, m₃, i, eout)

Return the critical semi-major axis ratio (a_out/a_in)_crit as defined in
Mardling & Aarset 1999, given a triple with masses `m₁`, `m₂`, `m₃`, mutual
inclination `i`, and outer eccentricity `eout`.
"""
function stability_criterion_ma01(m₁, m₂, m₃, i, eout)
    qout = m₃/(m₁ + m₂)
    return 2.8/(1 - eout)*(1 - 0.3i/π)*((1 + qout)*(1 + eout)/√(1 - eout))^(2/5)
end

"""
    stability_criterion_ma01(triple::HierarchicalMultiple)

The stability criterion can also be calculated for a triple system directly.
"""
function stability_criterion_ma01(triple::HierarchicalMultiple)
    @assert triple.n == 3 "System must be a triple."
    eₒ = triple.binaries[2].elements.e
    i = triple.binaries[1].elements.i

    m₁, m₂, m₃ = [triple.particles[i].mass for i = 1:3]

    return stability_criterion_ma01(m₁, m₂, m₃, i, eₒ)    
end

"""
    is_unstable(triple::HierarchicalMultiple)

Check if a given `triple` is unstable, i.e. whether aout/ain < (aout/ain)_crit.
"""
function is_unstable(triple::HierarchicalMultiple; criterion="ma01")
    α_crit = stability_criterion_ma01(triple)
    α = triple.binaries[2].elements.a / triple.binaries[1].elements.a
    return α < α_crit
end

function is_unstable(aₒ, aᵢ, m₁, m₂, m₃, i, eₒ; criterion="ma01")
    α_crit = stability_criterion_ma01(m₁, m₂, m₃, i, eₒ)
    α = aₒ/aᵢ
    return α < α_crit
end

@doc raw""" 
    octupole_parameter(triple::HierarchicalMultiple)

Return the octupole parameter ϵₒ for a given triple.

"""
function octupole_parameter(triple::HierarchicalMultiple)
    @assert triple.n == 3 "Octupole parameter only valid for triple system."
    m₁, m₂ = triple.particles.mass[[1, 2]]
    a_in  = triple.binaries[1].elements.a
    a_out = triple.binaries[2].elements.a
    e_out = triple.binaries[2].elements.e

    (m₁ - m₂)/(m₁ + m₂)*a_in/a_out*e_out/(1 - e_out^2)
end

function quadrupole_timescale(triple::HierarchicalMultiple)
    m = triple.particles.mass 
    P_in = triple.binaries[1].elements.P 
    P_out = triple.binaries[2].elements.P
    e_out = triple.binaries[2].elements.e
    return 16/30π*sum(m)/m[3]*P_out^2/P_in*cbrt(1 - e_out^2)^2 |> unit_time
end

function total_energy_2PN(r1, r2, v1, v2, m1, m2, G, c)

    r = r1 - r2
    r12 = norm(r1 - r2)
    n12 = r/r12

    vv1 = norm(v1)
    vv2 = norm(v2)

    n12v1 = dot(n12, v1)
    n12v2 = dot(n12, v2)
    v1v2 = dot(v1, v2)

    E_Newt = 0.5*m1*vv1^2 - G*m1*m2/2r12
    E_1PN = 1/c²*(G^2*m1^2*m2/2r12^2 + 3m1*vv1^4/8 + G*m1*m2/r12*(-0.25*n12v1*n12v2 + 1.5*vv1^2 - 7/4*v1v2))
    E_2PN = 1/c⁴*(-G^3*m1^3*m2/2r12^3 - 19*G^3*m1^2*m2^2/8r12^3 + 5m1*vv1^6/16 + 
                  G*m1*m2/r12*(3/8*n12v1^3*n12v2 + 3/16*n12v1^2*n12v2^2 - 9/8*n12v1*n12v2*vv1^2 -
                               13/8*n12v2^2*vv1^2 + 21/8*vv1^4 + 13/8*n12v1^2*v1v2 + 3/4*n12v1*n12v2*v1v2 -
                               -55/8*vv1^2*v1v2 + 17/8*v1v2^2 + 31/16*vv1^2*vv2^2) +
                  G^2*m1^2*m2/r12^2*(29/4*n12v1^2 - 13/4*n12v1*n12v2 + 0.5*(n12v2^2 - 3/2*vv1^2 + 7/4*vv2^2))
                  )
    E = E_Newt + E_1PN + E_2PN
    E_2p5PN = E + 4*G^2*m1^2*m2/(5*c^5*r12^2)*n12v1*(norm(v1 - v2)^2 - 2*G*(m1 - m2)/r12)
    return E + E_2p5PN
end

total_energy_2PN(r1, r2, v1, v2, m1, m2::Quantity) = total_energy_2PN(r1, r2, v1, v2, m1, m2, GRAVCONST, speed_of_light)
total_energy_2PN(r1, r2, v1, v2, m1, m2::Real) = total_energy_2PN(r1, r2, v1, v2, m1, m2, UNITLESS_G, UNITLESS_c)

function total_energy_2PN(res::SimulationResult, time, bodies=[1, 2])
    masses = res.ode_params.masses
    idx = findmin(x -> abs(x - time), res.solution.t)[2]
    E = 0.0
    i, j = bodies
    E += total_energy_2PN(SVector{3}(res.solution.u[idx].x[2][1:3, i]),
                            SVector{3}(res.solution.u[idx].x[2][1:3, j]), 
                            SVector{3}(res.solution.u[idx].x[1][1:3, i]),
                            SVector{3}(res.solution.u[idx].x[1][1:3, j]),
                            masses[i], masses[j])

    E += total_energy_2PN(SVector{3}(res.solution.u[idx].x[2][1:3, j]),
                            SVector{3}(res.solution.u[idx].x[2][1:3, i]), 
                            SVector{3}(res.solution.u[idx].x[1][1:3, j]),
                            SVector{3}(res.solution.u[idx].x[1][1:3, i]),
                            masses[j], masses[i])
        
    E
end


# """
#     PN1_energy(r1, r2, v1, v2, m1, m2; G=GRAVCONST)

    
# Total energy of body 1 in a gravitational + PN1 potential. From Blanchet 2014.
# """
# function PN1_energy(r1, r2, v1, v2, m1, m2; G=GRAVCONST)

#     r = r1 - r2

#     v1_norm = norm(v1)

#     r_norm = norm(r)

#     n = r/r_norm

#     E = G^2*m1^2*m2/(2*r_norm^2) + 3*m1*v1_norm^4/8 + G*m1*m2/r_norm*(-0.25*dot(n, v1)*dot(n, v2) + 3/2*v1_norm^2 - 7/4*dot(v1, v2))
#     Ekin = m1/2*v1_norm^2
#     Epot = -G*m1*m2/(2*r_norm)
#     return E*c⁻² + Ekin + Epot
# end

# function PN1_energy(sol::MultiBodySolution)

#     n_bodies = sol.ic.n

#     Etot = Vector{typeof(1.0u"J")}(undef, length(sol.t))

#     @inbounds for idx in eachindex(sol.t)
#         E = 0.0u"J"
#         for i = 1:n_bodies
#             # ri = sol.r[particle=i][:,idx]
#             # vi = sol.v[particle=i][:,idx]

#             ri = @SVector [sol.r[1, i, idx], sol.r[2, i, idx], sol.r[3, i, idx]]
#             vi = @SVector [sol.v[1, i, idx], sol.v[2, i, idx], sol.v[3, i, idx]]

#             mi = sol.structure.m[i,2]
#             for j = 1:n_bodies
#                 if j != i
#                     # rj = sol.r[particle=j][:,idx]
#                     # vj = sol.v[particle=j][:,idx]

#                     rj = @SVector [sol.r[1, j, idx], sol.r[2, j, idx], sol.r[3, j, idx]]
#                     vj = @SVector [sol.v[1, j, idx], sol.v[2, j, idx], sol.v[3, j, idx]]

#                     mj = sol.structure.m[j,2]
#                     E += PN1_energy(ri, rj, vi, vj, mi, mj)
#                 end
#             end
#         end
#         Etot[idx] = E
#     end

#     return Etot
# end

