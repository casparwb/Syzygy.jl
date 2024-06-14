####################################################
# This file contains the setup functions for various
# callbacks used in the simulations, such as collision, ejection, and
# roche-lobe overflows checks.
####################################################

using DiffEqCallbacks: ManifoldProjection

function setup_callbacks(conditions, system, p, retcodes, G, args; start_time=0)
    cbs = []

    n = system.n
    for condition in conditions
        cb = get_callback(condition, n, system, retcodes, G, args)
        if cb isa Tuple
            append!(cbs, cb)
        else
            push!(cbs, cb)
        end
    end
    return cbs
end

function get_callback(cb::CollisionCB, n, system, retcodes, G, args)
    condition_collision(u, t, integrator) = true
    affect_collision!(integrator) = collision_callback!(integrator, n, retcodes)
    callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
    
    return callback_collision
end

function get_callback(cb::EscapeCB, n, system, retcodes, G, args)
    @assert system.n == 3 "Escape check callback is currently only available for triples."
    
    function condition_escape(u, t, integrator)
        iszero(integrator.iter % cb.check_every)  # check for escape every 'check_every' iteration
    end

    affect_escape!(integrator) = unbound_callback!(integrator, retcodes, max_a_factor=cb.max_a_factor, G=G)
    callback_escape = DiscreteCallback(condition_escape, affect_escape!, save_positions=(false, false))
    
    return callback_escape
end

function get_callback(cb::RocheLobeOverflowCB, n, system, retcodes, G, args)

    function condition_rlof(u, t, integrator)
        (integrator.iter % cb.check_every) == 0  # check for rlof every 'check_every' iteration
    end

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

function get_callback(cb::CPUTimeCB, n, system, retcodes, G, args)
    condition_cpu_time(u, t, integrator) = true
    
    affect_cpu_time!(integrator) = max_cpu_time_callback!(integrator, retcode, start_time, args[:max_cpu_time])
    
    callback_cpu_time = DiscreteCallback(condition_cpu_time, affect_cpu_time!, save_positions=(false, false))
    
    return callback_cpu_time
end

function get_callback(cb::CentreOfMassCB, n, system, retcodes, G, args)

    affect_com!(integrator) = move_to_com_callback!(integrator)
    
    condition_com(u, t, integrator) = iszero(integrator.iter % cb.check_every)

    callback_com = DiscreteCallback(condition_com, affect_com!, save_positions=(false, false))
    return callback_com
end

function get_callback(cb::HubbleTimeCB, n, system, retcodes, G, args)

    affect_hubble!(integrator) = hubble_time_callback!(integrator, retcodes)
    
    condition_hubble(u, t, integrator) = true
    
    callback_hubble = DiscreteCallback(condition_hubble, affect_hubble!, save_positions=(false, false))
    
    return callback_hubble
end

function get_callback(cb::DemocraticCheckCB, n, system, retcodes, G, args)
    @assert system.n == 3 "Democratic interaction callback is currently only available for triples."
    
    affect_democratic1!(integrator) = democratic_check_callback_distance!(integrator, retcodes, system)
    affect_democratic2!(integrator) = democratic_check_callback_binary!(integrator, retcodes, system)
    affect_democratic3!(integrator) = democratic_check_callback_hyperbolic(integrator, retcodes, system)
    
    condition_democratic(u, t, integrator) = true
    
    callback_democratic1 = DiscreteCallback(condition_democratic, affect_democratic1!, save_positions=(false, false))
    callback_democratic2 = DiscreteCallback(condition_democratic, affect_democratic2!, save_positions=(false, false))
    callback_democratic3 = DiscreteCallback(condition_democratic, affect_democratic3!, save_positions=(false, false))

    return callback_democratic1, callback_democratic2, callback_democratic3
end

function get_callback(cb::IonizationCB, n, system, retcodes, G, args)

    condition_ionization(u, t, integrator) = true
    max_distance =  cb.max_a_factor*upreferred(system.binaries[2].elements.a).val
    affect_ionization!(integrator) = ionization_callback!(integrator, retcodes, max_distance)
    callback_ionization = Syzygy.OrdinaryDiffEq.DiscreteCallback(condition_ionization, affect_ionization!, save_positions=(false, false))

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
                    t = integrator.t * upreferred(1.0u"s")
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
function collision_callback!(integrator, n, retcode)
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

                collision::Bool = collision_check(d, Ri, Rj, Mi, Mj, stellar_type_i, stellar_type_j)::Bool
                if collision
                    t = integrator.t * upreferred(1.0u"s")
                    retcode[:Collision] = (SA[i, j], t)
                    terminate!(integrator)
                end
                # k += 1
            end
        end
    end
end

function collision_check(d, R1, R2, M1, M2, stellar_type1::Int, stellar_type2::Int)
    
    if (0 <= stellar_type1 <= 9) && (0 <= stellar_type2 <= 9) # two stars
        return collision_check_stars(d, R1, R2)
    elseif  (0 <= stellar_type1 <= 9) && (10 <= stellar_type2 <= 14) # one star, one CO
        return collision_check_star_compact_object(d, R1, M2, M1)
    elseif (10 <= stellar_type1 <= 14) && (0 <= stellar_type2 <= 9) # one star, one CO
        collision_check_star_compact_object(d, R2, M1, M2)
    elseif (10 <= stellar_type1 <= 14) && (10 <= stellar_type2 <= 14) # two COs
        return collision_check_compact_objects(d, M1, M2)
    end
end

function collision_check_stars(d, R1, R2)
    if isless(d, ustrip(R1 + R2))
        return true
    else
        return false
    end
end

function collision_check_star_compact_object(d, R2, M1, M2)
    tidal_disruption_radius = R2*cbrt(M1/M2) |> ustrip
    if d < (tidal_disruption_radius + ustrip(R2))
        return true
    else
        return false
    end
end

function collision_check_compact_objects(d, M1, M2) 
    rg = GRAVCONST*(M1 + M2)/c² # mutual gravitational radius
    if d <= 1000*upreferred(rg).val
        return true
    else
        return false
    end
end

# function collision_check(d, R1, R2, M1, M2, stellar_type1::T, stellar_type2::T) where {T}
#     if isless(d, ustrip(R1 + R2))
#         return true
#     else
#         return false
#     end
# end

# function collision_check(d, R1, R2, M1, M2, stellar_type1::TCO, stellar_type2::TS) where {TCO <: CompactObject, TS <: Star}
#     tidal_disruption_radius = R2*cbrt(M1/M2) |> ustrip
#     if d < (tidal_disruption_radius + ustrip(R2))
#         return true
#     else
#         return false
#     end
# end

# function collision_check(d, R1, R2, M1, M2, stellar_type1::T1, stellar_type2::T2) where {T1 <: CompactObject, T2 <: CompactObject}
#     rg = GRAVCONST*(M1 + M2)/c² # mutual gravitational radius
#     if d <= 100*upreferred(rg).val
#         return true
#     else
#         return false
#     end
# end

# function collision_check(d, R1, R2, M1, M2, stellar_type1::TS, stellar_type2::TCO)  where {TS <: Star, TCO <: CompactObject}
#     return collision_check(d, R2, R1, M2, M1, stellar_type2, stellar_type1)  
# end

@inline function get_masses(masses, sibling_ids::SVector{N, Int}) where N
    masses[sibling_ids]
end 

@inline function total_mass(masses, sibling_ids::SVector{N, Int}) where N
    M = zero(masses[1])
    @inbounds for k in sibling_ids
        M += masses[k]
    end

    M
end

function get_state_vectors(positions, velocities, masses, sibling::Particle, sibling_ids::SVector{N, Int}) where N
    return SA[SA[positions[1,particle.key.i], positions[2,particle.key.i], positions[3,particle.key.i]],
              SA[velocities[1,particle.key.i], velocities[2,particle.key.i], velocities[3,particle.key.i]]]
end

function get_state_vectors(positions, velocities, masses, binary::Binary, sibling_ids::SVector{N, Int}) where N
    m = get_masses(masses, sibling_ids)
    return SA[centre_of_mass(view(positions, :, sibling_ids), m),
              centre_of_mass_velocity(view(velocities, :, sibling_ids), m)]
           #m
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

@inline function get_velocities(velocities, masses, sibling::Particle, sibling_ids::SVector{N, Int} where N)
    return SA[velocities[1,sibling.key.i], velocities[2,sibling.key.i], velocities[3,sibling.key.i]]
end

@inline function get_velocities(velocities, masses, total_mass, sibling::T where T <: BinaryIndex, sibling_ids::SVector{N, Int}) where N
    mapreduce(+, sibling_ids) do k
        v = SA[velocities[1,k], velocities[2,k], velocities[3,k]]
        v * masses[k] / total_mass
    end
end



"""
Returns a callback for checking if system has become unbound.
"""
function unbound_callback!(integrator, retcode; max_a_factor=100, G=upreferred(GRAVCONST).val)

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
        M_bin = SA[ustrip(integrator.p.M[binary_ids[1]]), ustrip(integrator.p.M[binary_ids[2]])]
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

        M₁ = ustrip(integrator.p.M[i])
        M₂ = total_mass(integrator.p.M, sibling_ids) .|> ustrip

        com = get_positions(u.x[2], ustrip.(integrator.p.M), M₂, sibling, sibling_ids)
       
        r_rel = position - com

        d = norm(r_rel)

        R_roche = roche_radius(d, M₁/M₂)
        rlof = isless(R_roche, ustrip(integrator.p.R[i]))
        if rlof
            retcode[rcode] = upreferred(1.0u"s")*integrator.t
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
        velocity = SA[u.x[1][1,i], u.x[1][2,i], u.x[1][3,i]]

        # sibling_ids = [part for part in particles if part != i]
        # distances = [norm(position - SA[u.x[2][1,j], u.x[2][2,j], u.x[2][3,j]]) for j in sibling_ids]
        # sibling = sibling_ids[argmin(distances)]

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
        v = SA[u.x[1][1,sibling], u.x[1][2,sibling], u.x[1][3,sibling]]
        
        M₁ = ustrip(integrator.p.M[i])
        M₂ = ustrip(integrator.p.M[sibling])

        r_rel = position - r
        v_rel = velocity - v

        d = norm(r_rel)
        # v² = norm(v_rel)^2

        # a = semi_major_axis(d, v², M₁ + M₂, G)
        # e = eccentricity(r_rel, v_rel, a, M₁ + M₂, G)

        R_roche = roche_radius(d, M₁/M₂)#*(1 - e)

        rlof = R_roche <= ustrip(integrator.p.R[i])
        if rlof
            retcode[rcode] = integrator.t * upreferred(1.0u"s")
        end
    end
    nothing
end

function move_to_com_callback!(integrator)

    com = centre_of_mass(integrator.u.x[2], ustrip(integrator.p.M))
    com_vel = centre_of_mass_velocity(integrator.u.x[1], ustrip(integrator.p.M))

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
    if upreferred(1.0u"s")*integrator.t > 13.8u"Gyr"
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
function democratic_check_callback_binary!(integrator, retcode, system)
    if haskey(retcode, :Democratic_sma) # only need to raise flags once
        return
    end

    pairs = SA[(1,2), (2,3), (1,3)]
    sma = Inf
    bound_pair = (1, 2)
    @inbounds for pair in pairs
        i, j = pair

        ri = SA[integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        rj = SA[integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]

        vi = SA[integrator.u.x[1][1, i], integrator.u.x[1][2, i], integrator.u.x[1][3, i]]
        vj = SA[integrator.u.x[1][1, j], integrator.u.x[1][2, j], integrator.u.x[1][3, j]]

        Mi, Mj = ustrip(integrator.p.M[i]), ustrip(integrator.p.M[j])
        M = Mi + Mj

        r_rel = rj - ri
        v_rel = vj - vi
        
        K = 0.5*(Mi*norm(vi)^2 + Mj*norm(vj)^2)
        U = -upreferred(GRAVCONST).val*Mi*Mj/norm(r_rel)
    
        (K + U) > 0 && continue # if not bound, skip

        d = norm(r_rel)
        v² = norm(v_rel)^2
    
        a = semi_major_axis(d, v², M, upreferred(GRAVCONST).val)
        a < zero(a) && continue

        if a < sma
            sma = a
            bound_pair = pair
        end

    end

    if bound_pair != (1, 2)
        retcode[:Democratic_sma] = (true, integrator.t)
    end


end

"""
Check if any of the binaries have hyperbolic orbits (e > 1).
"""
function democratic_check_callback_hyperbolic(integrator, retcode, system)
    if haskey(retcode, :Democratic_ecc) # only need to raise flags once
        return
    end

    u = integrator.u

    r1 = SA[u.x[2][1,1], u.x[2][2,1], u.x[2][3,1]]
    r2 = SA[u.x[2][1,2], u.x[2][2,2], u.x[2][3,2]]

    v1 = SA[u.x[1][1,1], u.x[1][2,1], u.x[1][3,1]]
    v2 = SA[u.x[1][1,2], u.x[1][2,2], u.x[1][3,2]]

    M1, M2 = ustrip(integrator.p.M[1]), ustrip(integrator.p.M[2])
    M12 = M1 + M2

    r_rel = r2 - r1
    v_rel = v2 - v1

    d = norm(r_rel)
    v² = norm(v_rel)^2

    a = semi_major_axis(d, v², M12, upreferred(GRAVCONST).val)
    
    e = eccentricity(r_rel, v_rel, a, M12, upreferred(GRAVCONST).val)
    if e >= 1
        retcode[:Democratic_ecc] = (true, integrator.t)
    end

end

"""
Check if triple system has ionised (all three stars have become unbound).
"""
function ionization_callback!(integrator, retcode, max_distance, G=upreferred(GRAVCONST).val)

    u = integrator.u
			
    r1 = SA[u.x[2][1,1], u.x[2][2,1], u.x[2][3,1]]
    r2 = SA[u.x[2][1,2], u.x[2][2,2], u.x[2][3,2]]
    r3 = SA[u.x[2][1,3], u.x[2][2,3], u.x[2][3,3]]

    v1 = SA[u.x[1][1,1], u.x[1][2,1], u.x[1][3,1]]
    v2 = SA[u.x[1][1,2], u.x[1][2,2], u.x[1][3,2]]
    v3 = SA[u.x[1][1,3], u.x[1][2,3], u.x[1][3,3]]

    m = SA[ustrip(integrator.p.M[1]), ustrip(integrator.p.M[2]), ustrip(integrator.p.M[3])]

    r12 = r2 .- r1
    r23 = r3 .- r2
    r13 = r3 .- r1
    
    com = Syzygy.centre_of_mass(SA[r1, r2, r3], m)

    d1 = norm(r1 .- com)
    d2 = norm(r2 .- com)
    d3 = norm(r3 .- com)

    distances_from_COM = SA[d1, d2, d3]
    distances_now = SA[norm(r12), norm(23), norm(13)]

    K = 0.5*m .* SA[norm(v1)^2, norm(v2)^2, norm(v3)^2]
    U = -G*m .* SA[m[2]/norm(r12) + m[3]/norm(r13),
                        m[1]/norm(r12) + m[3]/norm(r23),
                        m[1]/norm(r13) + m[2]/norm(r23)]

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
        retcode[:Ionization] = true
        Syzygy.OrdinaryDiffEq.terminate!(integrator)        
    end

end

# function equilibrium_tidal_spin_evolution_callback(u, t, integrator, tidal_potential)

#     k = tidal_potential.k
#     τ = tidal_potential.τ

#     R = integrator.p.R
#     Ms = integrator.p.M
#     # T = R .^ 3 / (tidal_potential.G*)

#     @inbounds for i in eachindex(R)
#         # T = R[j]
#         ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
#         vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]
#         T = R[i]^3/(tidal_potential.G*M[i]*τ)
        
#         # I = moment_of_inertia(...)
#         rg = √(I/(Ms[i]*R^2)) # radius of gyrations
#         for j in eachindex(R)
#             if i != j
#                 rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
#                 vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]
#                 rij, vij = ri - rj, vi - vj
#                 d = norm(rij)
#                 v = norm(vij)
#                 a = semi_major_axis(d, v^2, M, potential.G)
#                 e = eccentricity(rij, vij, a, M, potential.G)
#                 q = Ms[j]/Ms[i]
#                 M = ms[i] + ms[j]
               
#                 n = √(G*M)*(-1/a^3) # mean orbital angular velocity

#                 e² = e^2

#                 f₂ = f₂(e²)
#                 f₅ = f₅(e²)

#                 dΩ_dt = 3k/T*q^2/rg^2*(Rs[i]/a)^6*n/(1 - e²)^6*(
#                         f₂ - √((1 - e)^3)*f₅*Ω/n)
#             end
#         end
    
#     end
# end


# function f₁(e²)
#     1 + 31e²/2 + 255e²^2/8 + 185e²^3/16 + 25e²^4/65
# end

# function f₂(e²)
#     1 + 15e²/2 + 45e²^2/8 + 5e²^3/16
# end

# function f₃(e²)
#     1 + 15e²/4 + 15e²^2/8 + 5e²^3/64
# end

# function f₄(e²)
#     1 + 3e²/2 + e²^2/8
# end

# function f₅(e²)
#     1 + 3e² + 3e²^2/8
# end


function setup_callbacks_old(stopping_conditions, system, p, retcode, G, args; start_time=0)
    isnothing(stopping_conditions) && return stopping_conditions
    isempty(stopping_conditions) && return nothing
    cbs = []

    n = system.n
    for condition in stopping_conditions
        if condition isa String
            if condition == "collision_old"
                condition_collision_old(u, t, integrator) = true
                affect_collision_old!(integrator) = collision_callback_old!(integrator, n, retcode)
                callback_collision = DiscreteCallback(condition_collision_old, affect_collision_old!, save_positions=(false, false))
                push!(cbs, callback_collision)

            elseif condition == "collision"
                condition_collision(u, t, integrator) = true
                affect_collision!(integrator) = collision_callback!(integrator, n, retcode)
                callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
                push!(cbs, callback_collision)

            elseif condition == "escape"
                @assert system.n == 3 "Escape check callback is currently only available for triples."
                function condition_escape(u, t, integrator)
                    (integrator.iter % 100) == 0  # check for escape every 100 iteration
                end
                affect_escape!(integrator) = unbound_callback!(integrator, retcode, G=G)
                callback_escape = DiscreteCallback(condition_escape, affect_escape!, save_positions=(false, false))
                push!(cbs, callback_escape)

            elseif condition == "rlof"
                function condition_rlof(u, t, integrator)
                    (integrator.iter % 100) == 0  # check for rlof every 100 iteration
                end
                rlof_rcodes = [Symbol(:RLOF_, i) for i = 1:n]
                rlof_rcodes = SA[rlof_rcodes...]
                affect_rlof_hier!(integrator) =  rlof_callback_hierarchical!(integrator, retcode, 
                                                                            system.particles, system.binaries, 
                                                                            n, rlof_rcodes)
                affect_rlof_demo!(integrator) =  rlof_callback_democratic!(integrator, retcode, n, rlof_rcodes)
                callback_rlof_demo = DiscreteCallback(condition_rlof, affect_rlof_demo!, save_positions=(false, false))
                callback_rlof_hier = DiscreteCallback(condition_rlof, affect_rlof_hier!, save_positions=(false, false))

                push!(cbs, callback_rlof_demo)
                push!(cbs, callback_rlof_hier)

            elseif condition == "tidal_disruption"
                function condition_td(u, t, integrator)
                    (integrator.iter % 100) == 0  # check for tidal disruption every 100 iteration
                end
                affect_td!(integrator) =  tidal_disruption_callback!(integrator, retcode, system, G)
                callback_td = DiscreteCallback(condition_td, affect_td!, save_positions=(false, false))
                push!(cbs, callback_td)

            elseif condition == "max_cpu_time"
                condition_cpu_time(u, t, integrator) = true
                affect_cpu_time!(integrator) = max_cpu_time_callback!(integrator, retcode, start_time, args[:max_cpu_time])
                callback_cpu_time = DiscreteCallback(condition_cpu_time, affect_cpu_time!, save_positions=(false, false))
                push!(cbs, callback_cpu_time)

            elseif condition == "manifold"
                r0 = [upreferred.(p.position) for p in values(system.particles)] 
                v0 = [upreferred.(p.velocity) for p in values(system.particles)] 
                masses = [upreferred(p.structure.m) for p in values(system.particles)] 
                Einit = total_energy(r0, v0, masses) |> upreferred |> ustrip

                function g(resid, u, p, t)
                    E = total_energy(u.x[2], u.x[1], p.M)
                    resid[1] = Einit - E
                    resid[2] = 0
                end
            
                callback_manifold = ManifoldProjection(g)
                push!(cbs, callback_manifold)

            elseif condition == "com"
                
                affect_com!(integrator) = move_to_com_callback!(integrator)
                condition_com(u, t, integrator) = (integrator.iter % 1000) == 0
                callback_com = DiscreteCallback(condition_com, affect_com!, save_positions=(false, false))
                push!(cbs, callback_com)

            elseif condition == "hubbletime"
                affect_hubble!(integrator) = hubble_time_callback!(integrator, retcode)
                condition_hubble(u, t, integrator) = true
                callback_hubble = DiscreteCallback(condition_hubble, affect_hubble!, save_positions=(false, false))
                push!(cbs, callback_hubble)

            elseif condition == "democratic"
                @assert system.n == 3 "Democratic interaction callback is currently only available for triples."
                affect_democratic1!(integrator) = democratic_check_callback_distance!(integrator, retcode, system)
                affect_democratic2!(integrator) = democratic_check_callback_binary!(integrator, retcode, system)
                affect_democratic3!(integrator) = democratic_check_callback_hyperbolic(integrator, retcode, system)
                condition_democratic(u, t, integrator) = true
                callback_democratic1 = DiscreteCallback(condition_democratic, affect_democratic1!, save_positions=(false, false))
                callback_democratic2 = DiscreteCallback(condition_democratic, affect_democratic2!, save_positions=(false, false))
                callback_democratic3 = DiscreteCallback(condition_democratic, affect_democratic3!, save_positions=(false, false))

                push!(cbs, callback_democratic1)
                push!(cbs, callback_democratic2)
                push!(cbs, callback_democratic3)

            elseif condition == "ionization"
            
                condition_ionization(u, t, integrator) = true
                max_distance =  100*upreferred(system.binaries[2].elements.a).val
                affect_ionization!(integrator) = ionization_callback!(integrator, retcode, max_distance)
                callback_ionization = Syzygy.OrdinaryDiffEq.DiscreteCallback(condition_ionization, affect_ionization!, save_positions=(false, false))
                push!(cbs, callback_ionization)
            else
                continue
            end
        else
            push!(cbs, condition)
        end
    end

    return cbs

end


# function tidal_disruption_callback!(integrator, retcode, system, G=upreferred(GRAVCONST).val)

#     @inbounds for i ∈ 1:system.n
#         stellar_type = round(Int, ustrip(integrator.p.stellar_types[i]))

#         # check if particle is a black hole, supernova, or unknown type
#         if stellar_type == 14 || stellar_type == 15 || stellar_type == 16
#             continue
#         end
#         ri = @SVector [integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
#         for j ∈ i:system.n
#             stellar_type_j = round(Int, ustrip(integrator.p.stellar_types[j]))
#             if stellar_type_j == 14 && i != j
#                 rj = @SVector [integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
#                 d = norm(ri - rj)

#                 tidal_disruption_radius = integrator.p.R[i]*cbrt(integrator.p.M[j]/integrator.p.M[i]) |> ustrip
#                 if (d - ustrip(integrator.p.R[i])) < tidal_disruption_radius
#                     t = integrator.t * upreferred(1.0u"s")
#                     retcode[:TidalDisruption] = (SA[i, j], t)
#                     terminate!(integrator)
#                 end
#             end
#         end
#     end
# end