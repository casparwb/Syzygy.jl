using LinearAlgebra, StaticArrays, AxisArrays   


function Base.show(io::IO, sol::FewBodySolution)
    print(io, "n datapoints = ")
    show(io, length(sol.t))
    println(io)
    print(io, "Time span: ")
    timespan = u"kyr".([sol.t[1], sol.t[end]])
    show(io, (timespan[1], timespan[2]))

    if !ismissing(sol.ode_system)
        println()
        print(io, "ODE parameters: \n")
        for (k, v) in sol.ode_system
            # show(io, (k, v))
            @printf(io, "   %-16s %s", "$k", "$v")
            println()
        end
    end

end

function setup_simulation_solution(n_steps, n_bodies, n_binaries, masses)

    elements = OrbitalElements[]

    radius_matrix = Matrix{typeof((1.0u"m"))}(undef, n_bodies, n_steps)
    mass_matrix = Matrix{typeof(1.0u"kg")}(undef, n_bodies, n_steps)
    type_matrix = Matrix{Int8}(undef, n_bodies, n_steps)
    lum_matrix = Matrix{typeof(1.0u"Lsun")}(undef, n_bodies, n_steps)
    spin_matrix = Matrix{typeof(1.0u"1/s")}(undef, n_bodies, n_steps)

    structure = StellarStructure(type_matrix, mass_matrix, 
                                 radius_matrix, spin_matrix, lum_matrix,
                                 similar(radius_matrix), similar(mass_matrix), 
                                 similar(radius_matrix), similar(mass_matrix))

    quantities = PhysicalQuantities(
                                    Array{typeof(1.0u"m^2/s"), 3}(undef, 3, n_binaries, n_steps),
                                    Vector{typeof(1.0u"J")}(undef, n_steps),
                                    Vector{typeof(1.0u"J")}(undef, n_steps),
                                    Vector{typeof(1.0u"J")}(undef, n_steps)
                                    )                           
    let
        sma_vec    = Vector{typeof(1.0u"m")}(undef, n_steps)
        period_vec = similar(sma_vec, typeof(1.0u"d"))
        ecc_vec    = similar(sma_vec, typeof(1.0u"1"))
        ω_vec      = similar(sma_vec, typeof(1.0u"°"))
        inc_vec    = similar(ω_vec)
        Ω_vec      = similar(inc_vec)

        for b ∈ 1:n_binaries
            binary_elements = OrbitalElements(similar(sma_vec),
                                              similar(period_vec),
                                              similar(ecc_vec),
                                              similar(ω_vec),
                                              similar(inc_vec),
                                              similar(Ω_vec),
                                              similar(Ω_vec))
            
            
            push!(elements, binary_elements)
        end
    end

    return elements, structure, quantities
end


function add_orbital_elements!(elements, new_elements, binary_key, idx)

    # Inner binary orbital elements
    a, P, e, i, Ω, ω, ν = new_elements

    elements[binary_key].a[idx] = a
    elements[binary_key].P[idx] = P
    elements[binary_key].e[idx] = e
    elements[binary_key].ω[idx] = ω
    elements[binary_key].i[idx] = i
    elements[binary_key].Ω[idx] = Ω
    elements[binary_key].ν[idx] = ν

end

function add_orbital_quantities!(quantities, n_bodies, pos, vel, mass, idx)
    quantities.K[idx] = kinetic_energy(vel, mass)
    quantities.U[idx] = potential_energy(pos, mass)
    quantities.E[idx] = quantities.K[idx]  + quantities.U[idx] 
    # @inbounds for i ∈ 1:n_bodies
    #     quantities.h[:,i,idx] .= angular_momentum(pos[:,i], vel[:,i])
    # end

end

function analyse_simulation(result::SimulationResult)

    # result = sim.result
    system = result.simulation.ic

    n_bodies = system.n
    binaries = system.binaries
    binary_keys = keys(binaries)
    n_binaries = length(binaries)
    time = result.solution.t .* u"s"
    n_steps = length(time)

    masses = result.simulation.params.M .* upreferred(1.0u"kg") 
    radii = result.simulation.params.R  .* upreferred(1.0u"m")
    spins = result.simulation.params.S  .* upreferred(1.0u"1/s")

    luminosities = if :L in propertynames(result.simulation.params)
                       result.simulation.params.L  .* upreferred(1.0u"Lsun")
                   else
                       [system.particles[i].structure.L for i ∈ 1:system.n]
                   end

    if masses[1] isa Number
        mass_vec = Matrix{typeof(upreferred(1.0u"kg"))}(undef, n_bodies, n_steps)
        @inbounds for i = 1:n_bodies
            mass_vec[i,:] .= masses[i]
        end
        masses = mass_vec
    end

    if radii[1] isa Number
        rad_vec = Matrix{typeof(upreferred(1.0u"m"))}(undef, n_bodies, n_steps)
        @inbounds for i = 1:n_bodies
            rad_vec[i,:] .= radii[i]
        end
        radii = rad_vec
    end

    if spins[1] isa Number
        spin_vec = Matrix{typeof(upreferred(1.0u"1/s"))}(undef, n_bodies, n_steps)
        @inbounds for i = 1:n_bodies
            spin_vec[i,:] .= spins[i]
        end
        spins = spin_vec
    end

    if luminosities[1] isa Number
        lum_vec = Matrix{typeof(upreferred(1.0u"Lsun"))}(undef, n_bodies, n_steps)
        @inbounds for i = 1:n_bodies
            lum_vec[i,:] .= luminosities[i]
        end
        luminosities = lum_vec
    end

    elements, structure, quantities = setup_simulation_solution(n_steps, n_bodies, n_binaries, masses)
    @inbounds for i in eachindex(time)
        structure.m[:,i] .= masses[:,i]
        structure.R[:,i] .= radii[:,i]
        structure.L[:,i] .= luminosities[:,i]
        
        structure.S[:,i] .= spins[:,i]
        structure.type[:,i] .= [system.particles[i].structure.type.index for i in 1:n_bodies]

        structure.R_core[:,i] = [system.particles[i].structure.R_core for i in 1:n_bodies]
        structure.m_core[:,i] = [system.particles[i].structure.m_core for i in 1:n_bodies]
        structure.R_env[:,i] = [system.particles[i].structure.R_env for i in 1:n_bodies]
        structure.m_env[:,i] = [system.particles[i].structure.m_env for i in 1:n_bodies]

    end

    r = Array{typeof(upreferred(1.0u"m")), 3}(undef, 3, n_bodies, n_steps)
    v = similar(r, typeof(1.0u"m/s"))

    r = AxisArray(r; dim=1:3, particle=1:n_bodies, time=time)
    v = AxisArray(v; dim=1:3, particle=1:n_bodies, time=time)
    for idx in eachindex(time)

        pos = result.solution.u[idx].x[2] .* upreferred(1.0u"m")
        vel = result.solution.u[idx].x[1] .* upreferred(1.0u"m/s")

        r[:, :, idx] .= pos
        v[:, :, idx] .= vel
        
        @inbounds for i ∈ binary_keys

            bin = binaries[i]
            children = bin.children
            
            if all(x -> x isa Particle, children) # Two particles
                particles = [p.key.i for p in children] 

                r_rel = pos[:,particles[2]] .- pos[:,particles[1]]# |> vec
                v_rel = vel[:,particles[2]] .- vel[:,particles[1]]# |> vec
                h = angular_momentum(r_rel, v_rel)
                mass = masses[particles[1], i] + masses[particles[2], i]
                
                new_elements = binary_orbital_elements(r_rel, v_rel, mass)
                add_orbital_elements!(elements, new_elements, i, idx)
            elseif any(x -> x isa Particle, children) # One binary and one Particle
                particle = findall(x -> x isa Particle, children) |> first
                binary = findall(x -> x isa Binary, children) |> first
                
                binary = children[binary]
                particle = children[particle]

                p_idx = particle.key.i

                sibling_keys = sort(binary.nested_children)
                m = masses[sibling_keys,i]
                
                binary_com = centre_of_mass(pos[:, sibling_keys], m) 
                binary_com_vel = centre_of_mass_velocity(vel[:, sibling_keys], m) 
                
                r_rel = pos[:, p_idx] .- binary_com
                v_rel = vel[:, p_idx] .- binary_com_vel
                h = angular_momentum(r_rel, v_rel)

                mass = sum(m) + particle.mass

                new_elements = binary_orbital_elements(r_rel, v_rel, mass)
                add_orbital_elements!(elements, new_elements, i, idx)
            else # Two binaries
                binary_1 = children[1]
                binary_2 = children[2]

                # binary_1_particles = get_particles_recursive(children[1])
                binary_1_keys = binary_1.nested_children#[p.key.i for p in binary_1_particles]

                binary_1_com = centre_of_mass(view(pos, :, binary_1_keys), masses[binary_1_keys,i]) 
                binary_1_com_vel = centre_of_mass_velocity(view(vel, :, binary_1_keys), masses[binary_1_keys,i]) 

                # binary_2_particles = get_particles_recursive(children[1])
                binary_2_keys = binary_2.nested_children#[p.key.i for p in binary_2_particles]

                binary_2_com = centre_of_mass(view(pos, :, binary_2_keys), masses[binary_2_keys,i]) 
                binary_2_com_vel = centre_of_mass_velocity(view(vel, :, binary_2_keys), masses[binary_2_keys,i]) 

                r_rel = binary_2_com .- binary_1_com 
                v_rel = binary_2_com_vel .- binary_1_com_vel
                h = angular_momentum(r_rel, v_rel)

                mass = sum(children[1].masses) + sum(children[2].masses)
                new_elements = binary_orbital_elements(r_rel, v_rel, mass)

                add_orbital_elements!(elements, new_elements, i, idx)
            end
        end

        add_orbital_quantities!(quantities, n_bodies, pos, vel, structure.m[:,idx], idx)

    end

    ode_solution = Dict(:potential => result.simulation.system.potential |> values,
                        :retcodes => result.retcode,
                        :timesteps => result.solution.destats.naccept,
                        :func_1_evals => result.solution.destats.nf,
                        :func_2_evals => result.solution.destats.nf2,
                        :runtime => result.runtime)
    attributes = merge(ode_solution, result.simulation.args, result.simulation.diffeq_args)
    # FewBodySolution(system, time, r, v, elements, structure, quantities, system)
    FewBodySolution(system, time, r, v, elements, structure, quantities, attributes, result.simulation.params)
end



"""
    binary_orbital_elements(pos, vel, mass)

Calculate orbital elements of binary with relative position `pos`, relative velocities
`vel`, and total mass `mass`.
"""
function binary_orbital_elements(pos, vel, M)

    # Relative velocities
    v²   = norm(vel) .^ 2 

    # Relative distances
    d  = norm(pos)

    orbital_elements_from_kinematics(pos, d, vel, v², M)
end

function orbital_elements_from_kinematics(r, d, v, v², M, G=𝒢)

    # M = sum(masses)

    # Semi major axes
    a =  semi_major_axis(d, v², M)

    # Period 
    P = orbital_period(a, M, G)

    # Eccentricities
    e = eccentricity(r, v, a, M)

    h = angular_momentum(r, v)

    # Inclination
    i = inclination(h)

    # Longitude of ascending node
    Ω = longitude_of_ascending_node(h)

    # Argument of periapsis
    ω = argument_of_periapsis(r, v, h, M)

    # True anomaly
    ν = true_anomaly(r, v, h, M)

    return a, P, e, i, Ω, ω, ν
end

function multibodysystem(sol::FewBodySolution, time)

    time_index = argmin(abs.(time .- sol.t))
    # positions = sol.r[time=time]
    # velocities = sol.v[time=time]
    masses = sol.structure.m[:,time_index]

    n_bodies = length(masses)
    multibodysystem(masses, R = sol.structure.R[:,time_index], 
                            L = sol.structure.L[:,time_index],
                            S = sol.structure.S[:,time_index],
                            R_core = sol.structure.R_core[:,time_index],
                            m_core = sol.structure.m_core[:,time_index],
                            R_env = sol.structure.R_env[:,time_index],
                            m_env = sol.structure.m_env[:,time_index],
                            types = Int.(sol.structure.type[:,time_index]),
                            a = [e.a[time_index] for e in sol.elements],
                            e = [e.e[time_index] for e in sol.elements], 
                            ω = [e.ω[time_index] for e in sol.elements], 
                            i = [e.i[time_index] for e in sol.elements], 
                            Ω = [e.Ω[time_index] for e in sol.elements], 
                            ν = [e.ν[time_index] for e in sol.elements], 
                            hierarchy = sol.initial_conditions.hierarchy,
                            t0 = time
                            )
end

"""

add callback for adding supernova kick

"""

# begin
#     res = simulate(m, dtmin=1.0e-9, npoints=10000);
#     time =res.solution.t

#     n_bodiess = m.n
#     n_steps = length(time)
#     r = Array{typeof(1.0u"m"), 3}(undef, 3, n_bodies, n_steps)
#     v = similar(r, typeof(1.0u"m/s"))
#     for idx in eachindex(time)
           
#                    pos = res.solution.u[idx].x[2] .* u"m"
#                    vel = res.solution.u[idx].x[1] .* u"m/s"
           
#                    r[:, :, idx] .= pos
#                    v[:, :, idx] .= vel
#                    end
#     plt = plot(aspect_ratio=1);
#     for i = 1:m.n
#            plot!(plt, r[1,i,:], r[2,i,:], label="$i")
#     end
#     plt
# end