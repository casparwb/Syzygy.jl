####################################################
# This file contains the setup functions for various
# callbacks used in the simulations, such as collision, ejection, and
# roche-lobe overflows checks.
####################################################

using DiffEqCallbacks: ManifoldProjection


function setup_callbacks(stopping_conditions, system, p, retcode, G, args; start_time=0)
    isnothing(stopping_conditions) && return stopping_conditions
    isempty(stopping_conditions) && return nothing
    cbs = []

    for condition in stopping_conditions
        if condition isa String
            if condition == "collision"
                n = system.n
                condition_collision(u, t, integrator) = true
                affect_collision!(integrator) = collision_callback!(integrator, n, retcode)
                callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
                push!(cbs, callback_collision)

            elseif condition == "escape"
                function condition_escape(u, t, integrator)
                    (integrator.iter % 100) == 0  # check for escape every 100 iteration
                end
                affect_escape!(integrator) = unbound_callback!(integrator, retcode, system, G=G)
                callback_escape = DiscreteCallback(condition_escape, affect_escape!, save_positions=(false, false))
                push!(cbs, callback_escape)

            elseif condition == "rlof"
                function condition_rlof(u, t, integrator)
                    (integrator.iter % 100) == 0  # check for rlof every 100 iteration
                end
                affect_rlof!(integrator) =  rlof_callback!(integrator, retcode, system, G)
                callback_rlof = DiscreteCallback(condition_rlof, affect_rlof!, save_positions=(false, false))
                push!(cbs, callback_rlof)

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
                n = system.n
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
                affect_democratic1!(integrator) = democratic_check_callback1!(integrator, retcode, system)
                affect_democratic2!(integrator) = democratic_check_callback2!(integrator, retcode, system)
                # affect_democratic3!(integrator) = democratic_check_callback3!(integrator, retcode, system)
                condition_democratic(u, t, integrator) = true
                callback_democratic1 = DiscreteCallback(condition_democratic, affect_democratic1!, save_positions=(false, false))
                callback_democratic2 = DiscreteCallback(condition_democratic, affect_democratic2!, save_positions=(false, false))
                # callback_democratic3 = DiscreteCallback(condition_democratic, affect_democratic3!, save_positions=(false, false))

                push!(cbs, callback_democratic1)
                push!(cbs, callback_democratic2)
                # push!(cbs, callback_democratic3)
            
            else
                continue
            end
        else
            push!(cbs, condition)
        end
    end

    return cbs

end

"""

Returns a callback for checking if collision has occured in system.
"""
function collision_callback!(integrator, n, retcode)
    # k = 1
    @inbounds for i ∈ 1:n
        ri = @SVector [integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        for j ∈ i:n
            if i != j
                rj = @SVector [integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
                d = norm(ri - rj)

                if (d - integrator.p.R[j]) < integrator.p.R[i]
                    t = u"kyr"(integrator.t * u"s")
                    retcode[:Collision] = ([i, j], t)
                    terminate!(integrator)
                end
                # k += 1
            end
        end
    end
end



function get_state_vectors(positions, velocities, masses, particle::Particle, sibling_ids)
    return SA[positions[1,particle.key.i], positions[2,particle.key.i], positions[3,particle.key.i]],
           SA[velocities[1,particle.key.i], velocities[2,particle.key.i], velocities[3,particle.key.i]]
    #view(positions, :, particle.key.i)#, view(velocities, :, particle.key.i), view(masses, particle.key.i)
end

function get_state_vectors(positions, velocities, masses, binary::Binary, sibling_ids)
    m = masses[sibling_ids]
    return SA[centre_of_mass(view(positions, :, sibling_ids), m),
              centre_of_mass_velocity(view(velocities, :, sibling_ids), m)]
           #m
end

function get_positions(positions, masses, particle::Particle, sibling_ids)
    return SA[positions[1,particle.key.i], positions[2,particle.key.i], positions[3,particle.key.i]]
end

function get_positions(positions, masses, binary::Binary, sibling_ids)
    m = masses[sibling_ids]
    return centre_of_mass(view(positions, :, sibling_ids), m)
end

function get_velocities(velocities, masses, particle::Particle, sibling_ids)
    return SA[velocities[1,particle.key.i], velocities[2,particle.key.i], velocities[3,particle.key.i]]
end

function get_velocities(velocities, masses, binary::Binary, sibling_ids)
    m = masses[sibling_ids]
    return centre_of_mass_velocity(view(velocities, :, sibling_ids), m)
end


"""

Returns a callback for checking if system has become unbound.
"""
function unbound_callback!(integrator, retcode, system; max_a=100, G=upreferred(𝒢).val)
    # positions = u.x[2]
    # velocities = u.x[1]
    # criteria = @MArray zeros(Bool, 4, system.n)
    u = integrator.u
    @inbounds for i ∈ 1:system.n#particle ∈ system.particles
        particle = system.particles[i]
        
        position = SA[u.x[2][1,i], u.x[2][2,i], u.x[2][3,i]]
        velocity = SA[u.x[1][1,i], u.x[1][2,i], u.x[1][3,i]]

        sibling = particle.sibling
        sibling_ids = get_particle_ids(system[sibling])

        # r = get_state_vectors(u.x[2], u.x[1], integrator.p.M, system[sibling], sibling_ids)[1]
        r = get_positions(u.x[2], integrator.p.M, system[sibling], sibling_ids)
        r_rel = position  - r#, velocity  - v

        d = norm(r_rel)

        new_position = position + velocity .* integrator.dt

        criteria_1 = d > (max_a * integrator.p.a[i])    # Body is a certain distance from sibling's centre of mass
        criteria_2 = norm(new_position - r) > d       # Body is moving away from centre of mass
        
        escape = criteria_1 && criteria_2
        if escape
            escapee = i
            T = kinetic_energy(velocity, integrator.p.M[i]) 
            
            U = -integrator.p.M[i]*integrator.p.M[sibling_ids[1]]/(norm(position - u.x[2][:,sibling_ids[1]]))
            @inbounds for j in sibling_ids[2:end]
                U -= integrator.p.M[i]*integrator.p.M[j]/(norm(position - u.x[2][:,j]))
            end
            U *= G
            
            Etot = T + U
            if Etot > zero(Etot)
                retcode[:Escape] = escapee
                terminate!(integrator)
            elseif d >= upreferred(1.0u"pc").val      # Body is 1 parsec away from its sibling
                retcode[:Drifter] = escapee
                terminate!(integrator)
            end
        end

    end

end

do_nothing!(integrator) = nothing

function rlof_callback!(integrator, retcode, system, G=upreferred(𝒢).val)

    u = integrator.u
    @inbounds for i ∈ 1:system.n#particle ∈ system.particles
        particle = system.particles[i]
        # particle = system.particles[i]
        
        # position, velocity = u.x[2][:,i], u.x[1][:,i]

        position = SA[u.x[2][1,i], u.x[2][2,i], u.x[2][3,i]]
        velocity = SA[u.x[1][1,i], u.x[1][2,i], u.x[1][3,i]]

        sibling = particle.sibling
        siblings = system[sibling]
        sibling_ids = get_particle_ids(siblings)
        
        r, v = get_state_vectors(u.x[2], u.x[1], integrator.p.M, siblings, sibling_ids)

        # m = integrator.p.M[sibling_ids]
        M₁ = integrator.p.M[i]
        M₂ = zero(M₁)

        for k in integrator.p.M[sibling_ids]
            M₂ += k
        end
        # M₂ = sum(m)
        r_rel = position - r
        v_rel = velocity - v

        d = norm(r_rel)
        v² = norm(v_rel)^2

        a = semi_major_axis(d, v², M₁ + M₂, G)
        e = eccentricity(r_rel, v_rel, a, M₁ + M₂, G)

        R_roche = roche_radius(a, M₁/M₂)*(1 - e)

        rlof = R_roche <= integrator.p.R[i]
        if rlof
            rcode = Symbol("RLOF_$i")
            if !(rcode in keys(retcode))
                retcode[rcode] = (upreferred(1.0u"s")*integrator.t)
            end
            # outcome(integrator)
        end
    end
end

function tidal_disruption_callback!(integrator, retcode, system, G=upreferred(𝒢).val)

    @inbounds for i ∈ 1:system.n
        stellar_type = system.particles[i].structure.type.index

        # check if particle is a black hole, supernova, or unknown type
        if stellar_type == 14 || stellar_type == 15 || stellar_type == 16
            continue
        end
        ri = @SVector [integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        for j ∈ i:system.n
            if i != j && system.particles[j].structure.type isa BlackHole
                rj = @SVector [integrator.u.x[2][1, j], integrator.u.x[2][2, j], integrator.u.x[2][3, j]]
                d = norm(ri - rj)

                tidal_disruption_radius = integrator.p.R[i]*cbrt(integrator.p.M[j]/integrator.p.M[i])
                if (d - integrator.p.R[i]) < tidal_disruption_radius
                    t = u"kyr"(integrator.t * u"s")
                    retcode[:TidalDisruption] = ([i, j], t)
                    terminate!(integrator)
                end
            end
        end
    end
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
    if upreferred(1.0u"s")*integrator.t > 13.8u"Gyr"
        retcode[:HubbleTime] = true
        terminate!(integrator)
    end
end


"""
Check if the pair with the smallest distance is no longer the initial inner binary.
"""
function democratic_check_callback1!(integrator, retcode, system)

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
Check if the original inner binary is no longer the one with the smallest semi-major axis,
or if any of the binaries have hyperbolic orbits (e > 1).
"""
function democratic_check_callback2!(integrator, retcode, system)
    if haskey(retcode, :Democratic_sma) && haskey(retcode, :Democratic_ecc) # only need to raise flags once
        return
    end

    u = integrator.u

    r1 = SA[u.x[2][1,1], u.x[2][2,1], u.x[2][3,1]]
    r2 = SA[u.x[2][1,2], u.x[2][2,2], u.x[2][3,2]]
    r3 = SA[u.x[2][1,3], u.x[2][2,3], u.x[2][3,3]]

    v1 = SA[u.x[1][1,1], u.x[1][2,1], u.x[1][3,1]]
    v2 = SA[u.x[1][1,2], u.x[1][2,2], u.x[1][3,2]]
    v3 = SA[u.x[1][1,3], u.x[1][2,3], u.x[1][3,3]]

    M1, M2, M3 = integrator.p.M[1], integrator.p.M[2], integrator.p.M[3]
    M12 = M1 + M2
    M123 = M12 + M3

    inner_com = centre_of_mass(SA[r1, r2], SA[M1, M2])
    inner_com_vel = centre_of_mass_velocity(SA[v1, v2], SA[M1, M2])

    r_3_inner_rel = r3 - inner_com
    v_3_inner_rel = v3 - inner_com_vel

    r_12_rel = r2 - r1
    v_12_rel = v2 - v1

    d_3_inner = norm(r_3_inner_rel)
    v²_3_inner = norm(v_3_inner_rel)^2

    d_12 = norm(v_12_rel)
    v²_12_rel = norm(v_12_rel)^2

    a_12 = semi_major_axis(d_12, v²_12_rel, M12, upreferred(𝒢.val))
    a_3_inner = semi_major_axis(d_3_inner, v²_3_inner, M123, upreferred(𝒢).val)
    
    e_in = eccentricity(r_12_rel, v_12_rel, a_12, M12, upreferred(𝒢).val)
    e_out = eccentricity(r_3_inner_rel, v_3_inner_rel, a_3_inner, M123, upreferred(𝒢).val)

    if a_12 > a_3_inner
        if !(haskey(retcode, :Democratic_sma))
            retcode[:Democratic_sma] = (true, integrator.t)
        end
    elseif e_in >= 1 || e_out >= 1
        if !(haskey(retcode, :Democratic_ecc))
            retcode[:Democratic_ecc] = (true, integrator.t)
        end
    end

end


"""
TO DO
"""
function supernova_kick_callback(u, t, integrator, t_sn)


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
