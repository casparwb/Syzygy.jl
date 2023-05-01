####################################################
# This file contains the setup functions for various
# callbacks used in the simulations, such as collision, ejection, and
# roche-lobe overflows checks.
####################################################

using DiffEqCallbacks: ManifoldProjection


function setup_callbacks(stopping_conditions, nbody, p, retcode, G, args; start_time=0)
    isnothing(stopping_conditions) && return stopping_conditions
    isempty(stopping_conditions) && return nothing
    cbs = []

    pin = (minimum([b.elements.P for b in values(nbody.binaries)]) |> upreferred).val
    for condition in stopping_conditions
        if condition isa String
            if condition == "collision"
                n = nbody.n
                condition_collision(u, t, integrator) = true
                affect_collision!(integrator) = collision_callback!(integrator, n, retcode)
                callback_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(false, false))
                push!(cbs, callback_collision)
            elseif condition == "escape"
                function condition_escape(u, t, integrator)
                    (integrator.iter % 1000) == 0  # check for escape every 1000 iteration
                end
                affect_escape!(integrator) = unbound_callback!(integrator, retcode, nbody, G=G)
                callback_escape = DiscreteCallback(condition_escape, affect_escape!, save_positions=(false, false))
                push!(cbs, callback_escape)
            elseif condition == "stellar"
                condition_evolve(u, t, integrator) = evolve_callback(u, t, integrator, p)
                affect_s!(integrator) = nothing
                callback_evolve = DiscreteCallback(condition_evolve, affect_s!)
                push!(cbs, callback_evolve)
            elseif condition == "rlof"
                function condition_rlof(u, t, integrator)
                    (integrator.iter % 1000) == 0  # check for escape every 1000 iteration
                end
                affect_rlof!(integrator) =  rlof_callback!(integrator, retcode, nbody, G)
                callback_rlof = DiscreteCallback(condition_rlof, affect_rlof!, save_positions=(false, false))
                push!(cbs, callback_rlof)
            elseif condition == "tidal_disruption"
                function condition_td(u, t, integrator)
                    (integrator.iter % 1000) == 0  # check for escape every 1000 iteration
                end
                affect_td!(integrator) =  tidal_disruption_callback!(integrator, retcode, nbody, G)
                callback_td = DiscreteCallback(condition_td, affect_td!, save_positions=(false, false))
                push!(cbs, callback_td)

            elseif condition == "max_cpu_time"
                condition_cpu_time(u, t, integrator) = true
                affect_cpu_time!(integrator) = max_cpu_time_callback!(integrator, retcode, start_time, args[:max_cpu_time])
                callback_cpu_time = DiscreteCallback(condition_cpu_time, affect_cpu_time!, save_positions=(false, false))
                push!(cbs, callback_cpu_time)
            elseif condition == "manifold"
                Einit = total_energy(u"m".(reduce(hcat, nbody.r₀))   |> ustrip |> Matrix, 
                                    u"m/s".(reduce(hcat, nbody.v₀)) |> ustrip |> Matrix, 
                                    u"kg".(nbody.m) |> ustrip)
                function g(resid, u, p, t)
                    E = total_energy(u.x[2], u.x[1], p.M)
                    resid[1] = Einit - E
                    resid[2] = 0
                end
            
                callback_manifold = ManifoldProjection(g)
                # callback_manifold = DiscreteCallback(condition_manifold, affect_s!)
                push!(cbs, callback_manifold)

            elseif condition == "com"
                n = nbody.n
                affect!(integrator) = move_to_com_callback!(integrator, n)
                condition_com(u, t, integrator) = true
                callback_com = DiscreteCallback(condition_com, affect!)
                push!(cbs, callback_com)
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
                # d = norm(integrator.u.x[2][:,i] - integrator.u.x[2][:,j])
                # d = norm(view(integrator.u.x[2], :, i) - view(integrator.u.x[2], :,j))
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

        r = get_state_vectors(u.x[2], u.x[1], integrator.p.M, system[sibling], sibling_ids)[1]
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
                # return true
                # @info "Escape"
                terminate!(integrator)
            elseif d >= upreferred(1.0u"pc").val      # Body is 1 parsec away from its sibling
                retcode[:Drifter] = escapee
                # return true
                # @info "Drift"
                terminate!(integrator)
            end
        end

    end

end

function evolve_callback(u, t, integrator, p)
    integrator.p[2] = p.R(t)
    integrator.p[3] = p.m(t)
    false
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
        if any(system.particles[i].structure.type .== (14, 15, 16))
            continue
        end
        ri = @SVector [integrator.u.x[2][1, i], integrator.u.x[2][2, i], integrator.u.x[2][3, i]]
        for j ∈ i:system.n
            if i != j && system.particles[j].structure.type == 14
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

function max_cpu_time_callback!(integrator, retcode, start_time, max_cpu_time)

    if (time() - start_time) >= max_cpu_time
        retcode[:MaxCPUTime] = true
        terminate!(integrator)
    end

end

function supernova_kick_callback(u, t, integrator, t_sn)


end

function equilibrium_tidal_spin_evolution_callback(u, t, integrator, tidal_potential)

    k = tidal_potential.k
    τ = tidal_potential.τ

    R = integrator.p.R
    Ms = integrator.p.M
    # T = R .^ 3 / (tidal_potential.G*)

    @inbounds for i in eachindex(R)
        # T = R[j]
        ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
        vi = @SVector [vs[1, i], vs[2, i], vs[3, i]]
        T = R[i]^3/(tidal_potential.G*M[i]*τ)
        
        # I = moment_of_inertia(...)
        rg = √(I/(Ms[i]*R^2)) # radius of gyrations
        for j in eachindex(R)
            if i != j
                rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
                vj = @SVector [vs[1, j], vs[2, j], vs[3, j]]
                rij, vij = ri - rj, vi - vj
                d = norm(rij)
                v = norm(vij)
                a = semi_major_axis(d, v^2, M, potential.G)
                e = eccentricity(rij, vij, a, M, potential.G)
                q = Ms[j]/Ms[i]
                M = ms[i] + ms[j]
               
                n = √(G*M)*(-1/a^3) # mean orbital angular velocity

                e² = e^2

                f₂ = f₂(e²)
                f₅ = f₅(e²)

                dΩ_dt = 3k/T*q^2/rg^2*(Rs[i]/a)^6*n/(1 - e²)^6*(
                        f₂ - √((1 - e)^3)*f₅*Ω/n)
            end
        end
    
    end
end


function f₁(e²)
    1 + 31e²/2 + 255e²^2/8 + 185e²^3/16 + 25e²^4/65
end

function f₂(e²)
    1 + 15e²/2 + 45e²^2/8 + 5e²^3/16
end

function f₃(e²)
    1 + 15e²/4 + 15e²^2/8 + 5e²^3/64
end

function f₄(e²)
    1 + 3e²/2 + e²^2/8
end

function f₅(e²)
    1 + 3e² + 3e²^2/8
end
