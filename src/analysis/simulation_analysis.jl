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

function setup_simulation_solution(n_steps, n_bodies, n_binaries)

    elements = OrbitalElements[]

    radius_matrix = Matrix{typeof((1.0u"m"))}(undef, n_bodies, 2)
    mass_matrix = Matrix{typeof(1.0u"kg")}(undef, n_bodies, 2)
    type_matrix = Matrix{Int8}(undef, n_bodies, 2)
    lum_matrix = Matrix{typeof(1.0u"Lsun")}(undef, n_bodies, 2)
    spin_matrix = Matrix{typeof(1.0u"1/s")}(undef, n_bodies, 2)

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

    initial_stellar_parameters = begin
        tmp = Dict()
        for param in propertynames(system.particles[1].structure)
            tmp[param] = [getproperty(system.particles[i].structure, param) for i = 1:n_bodies]
        end

        tmp[:type] = [t.index for t in tmp[:type]] 
        tmp
    end

    final_stellar_parameters = begin
        tmp = Dict()
        for param in propertynames(result.ode_params)
            tmp[param] = convert.(Unitful.Quantity, getproperty(result.ode_params, param))
        end

        tmp[:type] = pop!(tmp, :stellar_types)
        tmp[:m] = pop!(tmp, :M)
        tmp
    end

    elements, structure, quantities = setup_simulation_solution(n_steps, n_bodies, n_binaries)

    for param in keys(initial_stellar_parameters)
        if !(param ∈ keys(final_stellar_parameters))
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 1)
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 2)
        else
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 1)
            setindex!(getproperty(structure, param), final_stellar_parameters[param], 1:n_bodies, 2)
        end
    end

    r = Array{typeof(upreferred(1.0u"m")), 3}(undef, 3, n_bodies, n_steps)
    v = similar(r, typeof(1.0u"m/s"))

    r = AxisArray(r; dim=1:3, particle=1:n_bodies, time=time)
    v = AxisArray(v; dim=1:3, particle=1:n_bodies, time=time)
    for idx in eachindex(time)

        str_idx = min(idx, 2)

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
                mass = structure.m[particles[1], str_idx] + structure.m[particles[2], str_idx]
                
                new_elements = binary_orbital_elements(r_rel, v_rel, mass)
                add_orbital_elements!(elements, new_elements, i, idx)
            elseif any(x -> x isa Particle, children) # One binary and one Particle
                particle = findall(x -> x isa Particle, children) |> first
                binary = findall(x -> x isa Binary, children) |> first
                
                binary = children[binary]
                particle = children[particle]

                p_idx = particle.key.i

                sibling_keys = sort(binary.nested_children)
                m = structure.m[sibling_keys,str_idx]

                binary_com = centre_of_mass(pos[:, sibling_keys], m) 
                binary_com_vel = centre_of_mass_velocity(vel[:, sibling_keys], m) 
                
                r_rel = pos[:, p_idx] .- binary_com
                v_rel = vel[:, p_idx] .- binary_com_vel
                h = angular_momentum(r_rel, v_rel)

                mass = sum(m) + structure.m[p_idx, str_idx]

                new_elements = binary_orbital_elements(r_rel, v_rel, mass)
                add_orbital_elements!(elements, new_elements, i, idx)
            else # Two binaries
                binary_1 = children[1]
                binary_2 = children[2]

                # binary_1_particles = get_particles_recursive(children[1])
                binary_1_keys = binary_1.nested_children#[p.key.i for p in binary_1_particles]

                binary_1_com = centre_of_mass(view(pos, :, binary_1_keys), structure.m[binary_1_keys,str_idx]) 
                binary_1_com_vel = centre_of_mass_velocity(view(vel, :, binary_1_keys), structure.m[binary_1_keys,str_idx]) 

                # binary_2_particles = get_particles_recursive(children[1])
                binary_2_keys = binary_2.nested_children#[p.key.i for p in binary_2_particles]

                binary_2_com = centre_of_mass(view(pos, :, binary_2_keys), structure.m[binary_2_keys,str_idx]) 
                binary_2_com_vel = centre_of_mass_velocity(view(vel, :, binary_2_keys), structure.m[binary_2_keys,str_idx]) 

                r_rel = binary_2_com .- binary_1_com 
                v_rel = binary_2_com_vel .- binary_1_com_vel
                h = angular_momentum(r_rel, v_rel)

                # mass = sum(children[1].masses) + sum(children[2].masses)
                mass = sum(structure.m[binary_1_keys,str_idx]) + sum(structure.m[binary_2_keys,str_idx])
                new_elements = binary_orbital_elements(r_rel, v_rel, mass)

                add_orbital_elements!(elements, new_elements, i, idx)
            end
        end

        add_orbital_quantities!(quantities, n_bodies, pos, vel, structure.m[:,str_idx], idx)

    end

    ode_solution = Dict(:potential => result.simulation.system.potential |> values,
                        :retcodes => result.retcode,
                        :timesteps => result.solution.destats.naccept,
                        :func_1_evals => result.solution.destats.nf,
                        :func_2_evals => result.solution.destats.nf2,
                        :runtime => result.runtime)
    attributes = merge(ode_solution, result.simulation.args, result.simulation.diffeq_args)
    # FewBodySolution(system, time, r, v, elements, structure, quantities, system)
    FewBodySolution(system, time, r, v, elements, structure, 
                    quantities, attributes, result.ode_params)
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

    a =  semi_major_axis(d, v², M)

    P = orbital_period(a, M, G)

    e = eccentricity(r, v, a, M)

    h = angular_momentum(r, v)

    i = inclination(h)

    Ω = longitude_of_ascending_node(h)

    ω = argument_of_periapsis(r, v, h, M)

    ν = true_anomaly(r, v, h, M)

    return a, P, e, i, Ω, ω, ν
end

function multibodysystem(sol::FewBodySolution, time)

    time_index = argmin(abs.(time .- sol.t))
    # positions = sol.r[time=time]
    # velocities = sol.v[time=time]
    quant_time_index = size(sol.structure.m, 2) == length(sol.t) ? time_index : 2
    masses = sol.structure.m[:,quant_time_index]

    n_bodies = length(masses)
    # structure_args = Dict(prop => getproperty(sol.structure, prop)[:,time_index] for prop in fieldnames(StellarStructure))
    # structure_args[:types] = Int.(structure_args[:types])
    # orbital_args = Dict(prop => getindex.(getproperty(sol.elements, prop), time_index) for prop in fieldnames(OrbitalElements))
    multibodysystem(masses, R = sol.structure.R[:,quant_time_index], 
                            L = sol.structure.L[:,quant_time_index],
                            S = sol.structure.S[:,quant_time_index],
                            R_core = sol.structure.R_core[:,quant_time_index],
                            m_core = sol.structure.m_core[:,quant_time_index],
                            R_env = sol.structure.R_env[:,quant_time_index],
                            m_env = sol.structure.m_env[:,quant_time_index],
                            types = Int.(sol.structure.type[:,quant_time_index]),
                            a = [e.a[quant_time_index] for e in sol.elements],
                            e = [e.e[quant_time_index] for e in sol.elements], 
                            ω = [e.ω[quant_time_index] for e in sol.elements], 
                            i = [e.i[quant_time_index] for e in sol.elements], 
                            Ω = [e.Ω[quant_time_index] for e in sol.elements], 
                            ν = [e.ν[quant_time_index] for e in sol.elements], 
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