####################################################
# This file contains the setup functions for various
# callbacks used in the simulations, such as collision, ejection, and
# roche-lobe overflows checks.
####################################################

using DiffEqCallbacks: ManifoldProjection

abstract type AbstractSyzygyCallback end

struct FinalTimeCB{T} <: AbstractSyzygyCallback 
    t_final::T
end

struct CollisionCB         <: AbstractSyzygyCallback 
    check_every::Int
    grav_rad_multiple::Int
    CollisionCB(check_every=1, grav_rad_multiple=1000) = new(check_every, grav_rad_multiple)
end

struct NewCollisionCB         <: AbstractSyzygyCallback 
    check_every::Int
    grav_rad_multiple::Int
    NewCollisionCB(check_every=1, grav_rad_multiple=1000) = new(check_every, grav_rad_multiple)
end

struct EscapeCB{T}         <: AbstractSyzygyCallback 
    max_a_factor::T
    check_every::Int
    EscapeCB(check_every=100, max_a_factor=100) = new{typeof(max_a_factor)}(max_a_factor, check_every)
end

struct RocheLobeOverflowCB <: AbstractSyzygyCallback 
    check_every::Int
    RocheLobeOverflowCB(check_every=1) = new(check_every)
end

struct CPUTimeCB           <: AbstractSyzygyCallback 
    check_every::Int
    CPUTimeCB(check_every=1) = new(check_every)
end

struct CentreOfMassCB      <: AbstractSyzygyCallback 
    check_every::Int
    CentreOfMassCB(check_every=1) = new(check_every)
end

struct HubbleTimeCB        <: AbstractSyzygyCallback 
    check_every::Int
    HubbleTimeCB(check_every=1) = new(check_every)
end

struct DemocraticCheckCB   <: AbstractSyzygyCallback 
    check_every::Int
    DemocraticCheckCB(check_every=1) = new(check_every)
end

struct IonizationCB{T}     <: AbstractSyzygyCallback 
    check_every::Int
    max_a_factor::T
    IonizationCB(check_every=100, max_a_factor=100) = new{typeof(max_a_factor)}(max_a_factor, check_every)
end

function setup_callbacks(conditions, system, p, retcodes, args)
    cbs = []

    for condition in conditions
        cb = get_callback(condition, system, retcodes, args)
        if cb isa Tuple
            append!(cbs, cb)
        else
            push!(cbs, cb)
        end
    end
    
    return cbs
end

function get_callback(cb::FinalTimeCB, system, retcodes, args)
    condition_final_time(u, t, integrator) = integrator.t >= cb.t_final
    
    function affect!(integrator) 
        retcodes[:FinalTime] = true
        terminate!(integrator)
    end

    return DiscreteCallback(condition_final_time, affect!, save_positions=(false, false))
end

function get_callback(cb::CollisionCB, system, retcodes, args)
    n = system.n
    condition_collision(u, t, integrator) = true
    affect_collision!(integrator) = collision_callback!(integrator, n, retcodes, cb.grav_rad_multiple)
    callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
    
    return callback_collision
end

function get_callback(cb::NewCollisionCB, system, retcodes, args)
    condition_collision(u, t, integrator) = true
    pairs = system.pairs
    affect_collision!(integrator) = collision_callback_idiomatic!(integrator, pairs, retcodes, cb.grav_rad_multiple)
    callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
    
    return callback_collision
end

function get_callback(cb::EscapeCB, system, retcodes, args)
    @assert system.n == 3 "Escape check callback is currently only available for triples."
    
    function condition_escape(u, t, integrator)
        iszero(integrator.iter % cb.check_every)  # check for escape every 'check_every' iteration
    end

    affect_escape!(integrator) = unbound_callback!(integrator, retcodes, max_a_factor=cb.max_a_factor)
    callback_escape = DiscreteCallback(condition_escape, affect_escape!, save_positions=(false, false))
    
    return callback_escape
end

function get_callback(cb::RocheLobeOverflowCB, system, retcodes, args)

    function condition_rlof(u, t, integrator)
        (integrator.iter % cb.check_every) == 0  # check for rlof every 'check_every' iteration
    end

    n = system.n

    rlof_rcodes = [Symbol(:RLOF_, i) for i = 1:n]
    rlof_rcodes = SA[rlof_rcodes...]
    affect_rlof_hier!(integrator) =  rlof_callback_hierarchical!(integrator, retcodes, 
                                                                system.particles, system.binaries, 
                                                                n, rlof_rcodes)
    affect_rlof_demo!(integrator) =  rlof_callback_democratic!(integrator, retcodes, n, rlof_rcodes)
    callback_rlof_demo = DiscreteCallback(condition_rlof, affect_rlof_demo!, save_positions=(false, false))
    callback_rlof_hier = DiscreteCallback(condition_rlof, affect_rlof_hier!, save_positions=(false, false))

    return callback_rlof_demo, callback_rlof_hier
end

function get_callback(cb::CPUTimeCB, system, retcodes, args)
    condition_cpu_time(u, t, integrator) = true
    
    t_start = time()
    affect_cpu_time!(integrator) = max_cpu_time_callback!(integrator, retcodes, t_start, args[:max_cpu_time])
    
    callback_cpu_time = DiscreteCallback(condition_cpu_time, affect_cpu_time!, save_positions=(false, false))
    
    return callback_cpu_time
end

function get_callback(cb::CentreOfMassCB, system, retcodes, args)

    affect_com!(integrator) = move_to_com_callback!(integrator)
    
    condition_com(u, t, integrator) = iszero(integrator.iter % cb.check_every)

    callback_com = DiscreteCallback(condition_com, affect_com!, save_positions=(false, false))
    return callback_com
end

function get_callback(cb::HubbleTimeCB, system, retcodes, args)

    affect_hubble!(integrator) = hubble_time_callback!(integrator, retcodes)
    
    condition_hubble(u, t, integrator) = true
    
    callback_hubble = DiscreteCallback(condition_hubble, affect_hubble!, save_positions=(false, false))
    
    return callback_hubble
end

function get_callback(cb::DemocraticCheckCB, system, retcodes, args)
    @assert system.n == 3 "Democratic interaction callback is currently only available for triples."
    
    pair = system.pairs
    affect_democratic1!(integrator) = democratic_check_callback_distance!(integrator, retcodes, system)
    affect_democratic2!(integrator) = democratic_check_callback_binary!(integrator, pair, retcodes)
    affect_democratic3!(integrator) = democratic_check_callback_hyperbolic(integrator, retcodes)
    
    condition_democratic(u, t, integrator) = true
    
    callback_democratic1 = DiscreteCallback(condition_democratic, affect_democratic1!, save_positions=(false, false))
    callback_democratic2 = DiscreteCallback(condition_democratic, affect_democratic2!, save_positions=(false, false))
    callback_democratic3 = DiscreteCallback(condition_democratic, affect_democratic3!, save_positions=(false, false))

    return callback_democratic1, callback_democratic2, callback_democratic3
end

function get_callback(cb::IonizationCB, system, retcodes, args)

    condition_ionization(u, t, integrator) = true
    max_distance =  cb.max_a_factor*upreferred(system.binaries[2].elements.a).val
    affect_ionization!(integrator) = ionization_callback!(integrator, retcodes, max_distance)
    callback_ionization = DiscreteCallback(condition_ionization, affect_ionization!, save_positions=(false, false))

    return callback_ionization
end
"""

Returns a callback for checking if collision has occured in system.
If the two objects are stars, the callback checks for overlapping radii,
if one of the objects is a compact object and the other is a star, the tidal
radius of the CO is used, and finally if both objects are COs, we use 100 × gravitational radius.
"""
function collision_callback_old!(integrator, n, retcode)
    # k = 1
    @inbounds for i ∈ 1:n
        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        for j ∈ i:n
            if i != j
                rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
                d = norm(ri - rj)

                if isless(d - ustrip(integrator.p.R[j]), ustrip(integrator.p.R[i]))
                    t = integrator.t * upreferred(u"s")
                    retcode[:Collision] = (SA[i, j], t)
                    terminate!(integrator)
                end
                # k += 1
            end
        end
    end
end

"""

Returns a callback for checking if collision has occured in system.
If the two objects are stars, the callback checks for overlapping radii,
if one of the objects is a compact object and the other is a star, the tidal
radius of the CO is used, and finally if both objects are COs, we use 100 × gravitational radius.
"""
function collision_callback!(integrator, n, retcode, grav_rad_multiple)
    # k = 1
    @inbounds for i ∈ 1:n
        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        Ri = integrator.p.R[i]
        Mi = integrator.p.M[i]

        stellar_type_i = integrator.p.stellar_types[i]
        for j ∈ i:n
            if i != j
                rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
                d = norm(ri - rj)
                
                Rj = integrator.p.R[j]
                Mj = integrator.p.M[j]

                stellar_type_j = integrator.p.stellar_types[j]

                collision = collision_check(d, Ri, Rj, Mi, Mj, stellar_type_i, stellar_type_j, 
                                                  grav_rad_multiple)::Bool
                if collision
                    t = integrator.t * upreferred(u"s")
                    retcode[:Collision] = (SA[i, j], t)
                    terminate!(integrator)
                end
                # k += 1
            end
        end
    end
end

function collision_check(d, R1, R2, m1, m2, stellar_type1::Int, stellar_type2::Int, grav_rad_multiple)
    
    if (0 <= stellar_type1 <= 9) && (0 <= stellar_type2 <= 9) # two stars
        return collision_check_stars(d, R1, R2)
    elseif  (0 <= stellar_type1 <= 9) && (10 <= stellar_type2 <= 14) # one star, one CO
        return collision_check_star_compact_object(d, R1, m2, m1)
    elseif (10 <= stellar_type1 <= 14) && (0 <= stellar_type2 <= 9) # one star, one CO
        return collision_check_star_compact_object(d, R2, m1, m2)
    elseif (10 <= stellar_type1 <= 14) && (10 <= stellar_type2 <= 14) # two COs
        return collision_check_compact_objects(d, m1, m2, grav_rad_multiple)
    end
end

function collision_check_stars(d, R1, R2)
    d <= (R1 + R2)
end

function collision_check_star_compact_object(d, R_star, m_CO, m_star)
    tidal_disruption_radius = R_star*cbrt(m_CO/m_star)
    if d <= (tidal_disruption_radius + R_star)
        println(d, " ", tidal_disruption_radius)
        return true
    end

    return false
end

function collision_check_compact_objects(d, m1, m2, grav_rad_multiple) 
    rg = G*(m1 + m2)*c⁻² # mutual gravitational radius

    return d <= rg*grav_rad_multiple
end


"""

An attempt to dispatch on the stellar types. Is faster than the conditional version, but allocates.
"""
function collision_callback_idiomatic!(integrator, pairs, retcode, grav_rad_multiple)
    # k = 1
    @inbounds for pair in pairs
        i, j = pair
        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        Ri = integrator.p.R[i]
        Mi = integrator.p.M[i]

        stellar_type_i = stellar_types[integrator.p.stellar_types[i]]

        rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
        d = norm(ri - rj)
        
        Rj = integrator.p.R[j]
        Mj = integrator.p.M[j]

        stellar_type_j = stellar_types[integrator.p.stellar_types[j]]

        collision::Bool = collision_check(d, Ri, Rj, Mi, Mj, 
                                          stellar_type_i, stellar_type_j, 
                                          grav_rad_multiple)

        if collision#collision_check(d, Ri, Rj, Mi, Mj, stellar_type_i, stellar_type_j, grav_rad_multiple)
            t = integrator.t * upreferred(u"s")
            retcode[:Collision] = (SA[i, j], t)
            terminate!(integrator)
        end
    end
    nothing
end


function collision_check(d, R1, R2, m1, m2, stellar_type1::Star, stellar_type2::Star, _)
    d <= (R1 + R2)
end

function collision_check(d, R1, R2, m1, m2, stellar_type1::CompactObject, stellar_type2::Star, _)
    tidal_disruption_radius = R2*cbrt(m1/m2) 

    d <= (tidal_disruption_radius + R2)
end

function collision_check(d, R1, R2, m1, m2, stellar_type1::Star, stellar_type2::CompactObject)
    return collision_check(d, R2, R1, m2, m1, stellar_type2, stellar_type1)  
end

function collision_check(d, R1, R2, m1, m2, stellar_type1::CompactObject, stellar_type2::CompactObject, grav_rad_multiple)
    rg = G*(m1 + m2)*c⁻² # mutual gravitational radius
    d <= grav_rad_multiple*rg
end


@inline function total_mass(masses, sibling_ids::SVector{N, Int}) where N
    M = zero(masses[1])
    @inbounds for k in sibling_ids
        M += masses[k]
    end

    M
end

@inline function get_positions(positions, masses, total_mass, sibling::T where T <: ParticleIndex, sibling_ids::SVector{N, Int}) where N
    return SA[positions[1,sibling.i], positions[2,sibling.i], positions[3,sibling.i]]
end

@inline function get_positions(positions, masses, total_mass, sibling::T where T <: BinaryIndex, sibling_ids::SVector{N, Int}) where N
    mapreduce(+, sibling_ids) do k
        r = SA[positions[1,k], positions[2,k], positions[3,k]]
        r * masses[k] / total_mass
    end
end


"""
Returns a callback for checking if system has become unbound.
"""
function unbound_callback!(integrator, retcode; max_a_factor=100)

    u = integrator.u
    combinations = SA[(1, SA[2, 3]), (2, SA[1, 3]), (3, SA[1, 2])]
    @inbounds for (particle, binary_ids) in combinations
        
        # escape candidate
        r_part = SA[u.x[2][1,particle], u.x[2][2,particle], u.x[2][3,particle]]
        v_part = SA[u.x[1][1,particle], u.x[1][2,particle], u.x[1][3,particle]]

        # remaining components
        r_comp1 = SA[u.x[2][1,binary_ids[1]], u.x[2][2,binary_ids[1]], u.x[2][3,binary_ids[1]]]
        r_comp2 = SA[u.x[2][1,binary_ids[2]], u.x[2][2,binary_ids[2]], u.x[2][3,binary_ids[2]]]

        v_comp1 = SA[u.x[1][1,binary_ids[1]], u.x[1][2,binary_ids[1]], u.x[1][3,binary_ids[1]]]
        v_comp2 = SA[u.x[1][1,binary_ids[2]], u.x[1][2,binary_ids[2]], u.x[1][3,binary_ids[2]]]

        r_rel_bin = r_comp1 - r_comp2
        v_rel_bin = v_comp1 - v_comp2

        # binary properties of remaining components
        M_bin = SA[integrator.p.M[binary_ids[1]], integrator.p.M[binary_ids[2]]]
        a_bin = semi_major_axis(norm(r_rel_bin), norm(v_rel_bin)^2, M_bin[1] + M_bin[2], G)
        a_bin < zero(a_bin) && continue

        r_bin = centre_of_mass(SA[r_comp1, r_comp2], M_bin)

        r_rel = r_part - r_bin

        d = norm(r_rel)

        new_position = r_part + v_part .* integrator.dt

        criteria_1 = d > (max_a_factor * a_bin)          # Body is a certain distance from sibling's centre of mass
        criteria_2 = norm(new_position - r_bin) > d      # Body is moving away from centre of mass
        
        escape = criteria_1 && criteria_2
        if escape
            escapee = particle
            M = ustrip(integrator.p.M[particle])
            T = kinetic_energy(v_part, M) 
            
            U = -(G*M)*(M_bin[1]/norm(r_part - r_comp1) + M_bin[2]/norm(r_part - r_comp2))
            
            Etot = T + U
            if Etot > zero(Etot)
                retcode[:Escape] = escapee
                terminate!(integrator)
            elseif d >= upreferred(1.0u"pc").val      # Body is 1 parsec away from binary
                retcode[:Drifter] = escapee
                terminate!(integrator)
            end
        end

    end

end


"""
The hierarchical RLOF check uses the siblings (the inner binary for the tertiary, and the companion)
"""
function rlof_callback_hierarchical!(integrator, retcode, particles, binaries, n, rlof_rcodes)
    u = integrator.u
    @inbounds for i ∈ 1:n
        rcode = rlof_rcodes[i]
        haskey(retcode, rcode) && continue
        
        if !(stellar_types[integrator.p.stellar_types[i]] isa Star)
            continue
        end
        
        particle = particles[i]

        position = SA[u.x[2][1,i], u.x[2][2,i], u.x[2][3,i]]

        sibling = particle.sibling
        sibling_ids = if sibling isa ParticleIndex
            SA[particles[sibling.i].key.i]
        else
            binaries[sibling.i].nested_children
        end 

        M₁ = integrator.p.M[i]
        M₂ = total_mass(integrator.p.M, sibling_ids)

        com = get_positions(u.x[2], integrator.p.M, M₂, sibling, sibling_ids)
       
        r_rel = position - com

        d = norm(r_rel)

        R_roche = roche_radius(d, M₁/M₂)
        rlof = isless(R_roche, integrator.p.R[i])
        if rlof
            retcode[rcode] = upreferred(u"s")*integrator.t
        end
    end
end

"""
The democratic RLOF check uses the nearest particle.
"""
function rlof_callback_democratic!(integrator, retcode, n, rlof_rcodes)
    u = integrator.u
    @inbounds for i ∈ 1:n
        rcode = rlof_rcodes[i]
        haskey(retcode, rcode) && continue

        if !(stellar_types[integrator.p.stellar_types[i]] isa Star)
            continue
        end

        position = SA[u.x[2][1,i], u.x[2][2,i], u.x[2][3,i]]

        dist = Inf
        sibling = i
        for k ∈ 1:n
            if k != i
                r = SA[u.x[2][1,k], u.x[2][2,k], u.x[2][3,k]]
                d = norm(position - r)
                if d < dist
                    dist = d
                    sibling = k
                end
            end
        end
        
        r = SA[u.x[2][1,sibling], u.x[2][2,sibling], u.x[2][3,sibling]]
        
        M₁ = integrator.p.M[i]
        M₂ = integrator.p.M[sibling]

        r_rel = position - r

        d = norm(r_rel)

        R_roche = roche_radius(d, M₁/M₂)#*(1 - e)

        rlof = R_roche <= integrator.p.R[i]
        if rlof
            retcode[rcode] = integrator.t * upreferred(u"s")
        end
    end
    nothing
end

function move_to_com_callback!(integrator)
    com = centre_of_mass(integrator.u.x[2], integrator.p.M)
    com_vel = centre_of_mass_velocity(integrator.u.x[1], integrator.p.M)

    integrator.u.x[2] .-= com
    integrator.u.x[1] .-= com_vel
end

function max_cpu_time_callback!(integrator, retcode, start_time, max_cpu_time)

    if (time() - start_time) >= max_cpu_time
        retcode[:MaxCPUTime] = true
        terminate!(integrator)
    end

end

function hubble_time_callback!(integrator, retcode)
    if upreferred(u"s")*integrator.t > 13.8u"Gyr"
        retcode[:HubbleTime] = true
        terminate!(integrator)
    end
end


"""
Check if the pair with the smallest distance is no longer the initial inner binary.
"""
function democratic_check_callback_distance!(integrator, retcode, system)

    if haskey(retcode, :Democratic_dist) # only need to raise flag once
        return
    end

    smallest_distance = Inf

    pair = (1, 2)
    @inbounds for i ∈ 1:system.n
        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        for j ∈ i:system.n
            if i != j
                rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
                d = norm(ri - rj)
                if d < smallest_distance
                    smallest_distance = d
                    pair = (i, j)
                end
            end

        end
    end

    if pair != (1, 2)
        retcode[:Democratic_dist] = (true, integrator.t)
    end

end

"""
Check if the original inner binary is no longer the one with the smallest semi-major axis.
"""
function democratic_check_callback_binary!(integrator, pairs, retcodes)
    if haskey(retcodes, :Democratic_sma) # only need to raise flags once
        return
    end

    sma = Inf
    bound_pair = (1, 2)
    @inbounds for pair in pairs
        i, j = pair

        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]

        vi = SA[integrator.u.x[1][1, i], integrator.u.x[1][2, i], integrator.u.x[1][3, i]]
        vj = SA[integrator.u.x[1][1, j], integrator.u.x[1][2, j], integrator.u.x[1][3, j]]

        Mi, Mj = integrator.p.M[i], integrator.p.M[j]
        M = Mi + Mj

        r_rel = rj - ri
        v_rel = vj - vi
        
        K = 0.5*(Mi*norm(vi)^2 + Mj*norm(vj)^2)
        U = -G*Mi*Mj/norm(r_rel)
    
        (K + U) > 0 && continue # if not bound, skip

        d = norm(r_rel)
        v² = norm(v_rel)^2
    
        a = semi_major_axis(d, v², M, G)
        a < zero(a) && continue

        if a < sma
            sma = a
            bound_pair = pair
        end

    end

    if bound_pair != (1, 2)
        retcodes[:Democratic_sma] = (true, integrator.t)
    end


end

"""
Check if any of the binaries have hyperbolic orbits (e > 1).
"""
function democratic_check_callback_hyperbolic(integrator, retcodes)
    if haskey(retcodes, :Democratic_ecc) # only need to raise flags once
        return
    end

    u = integrator.u

    r1 = SA[u.x[2][1,1], u.x[2][2,1], u.x[2][3,1]]
    r2 = SA[u.x[2][1,2], u.x[2][2,2], u.x[2][3,2]]

    v1 = SA[u.x[1][1,1], u.x[1][2,1], u.x[1][3,1]]
    v2 = SA[u.x[1][1,2], u.x[1][2,2], u.x[1][3,2]]

    m1, m2 = integrator.p.M[1], integrator.p.M[2]
    m12 = m1 + m2

    r_rel = r2 - r1
    v_rel = v2 - v1

    d = norm(r_rel)

    e = eccentricity(r_rel, v_rel, d, m12, G)

    if e >= 1
        retcodes[:Democratic_ecc] = (true, integrator.t)
    end

end

"""
Check if triple system has ionised (all three stars have become unbound).
"""
function ionization_callback!(integrator, retcodes, max_distance)

    u = integrator.u
			
    r1 = SA[u.x[2][1,1], u.x[2][2,1], u.x[2][3,1]]
    r2 = SA[u.x[2][1,2], u.x[2][2,2], u.x[2][3,2]]
    r3 = SA[u.x[2][1,3], u.x[2][2,3], u.x[2][3,3]]

    v1 = SA[u.x[1][1,1], u.x[1][2,1], u.x[1][3,1]]
    v2 = SA[u.x[1][1,2], u.x[1][2,2], u.x[1][3,2]]
    v3 = SA[u.x[1][1,3], u.x[1][2,3], u.x[1][3,3]]

    m = SA[integrator.p.M[1], integrator.p.M[2], integrator.p.M[3]]

    r12 = r2 .- r1
    r23 = r3 .- r2
    r13 = r3 .- r1
    
    com = Syzygy.centre_of_mass(SA[r1, r2, r3], m)

    d1 = norm(r1 .- com)
    d2 = norm(r2 .- com)
    d3 = norm(r3 .- com)

    d12 = norm(r12)
    d23 = norm(r23)
    d13 = norm(r13)

    distances_from_COM = SA[d1, d2, d3]
    distances_now = SA[d12, d23, d13]

    K = 0.5*m .* SA[norm(v1)^2, norm(v2)^2, norm(v3)^2]
    U = -G*m .* SA[m[2]/d12 + m[3]/d13,
                   m[1]/d12 + m[3]/d23,
                   m[1]/d13 + m[2]/d23]

    r1 = r1 + v1 .* integrator.dt
    r2 = r2 + v2 .* integrator.dt
    r3 = r3 + v3 .* integrator.dt

    r12 = r2 .- r1
    r23 = r3 .- r2
    r13 = r3 .- r1
    
    distances_next = SA[norm(r13), norm(r23), norm(r13)]

    criteria_1 = all(distances_next .> distances_now)
    criteria_2 = all(distances_from_COM .>= max_distance)
    criteria_3 = all((K .+ U) .> 0)

    if criteria_1 && criteria_2 && criteria_3
        retcodes[:Ionization] = true
        Syzygy.OrdinaryDiffEq.terminate!(integrator)        
    end

end

