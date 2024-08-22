using LinearAlgebra: norm
using StaticArrays, AxisArrays   




"""
    get_stellar_structure(result::SimulationResult)

Return a `StellarStructure` type containing the initial and final stellar structure properties of the components. If stellar
evolution was not included, the initial and final values will be the same.
"""
function get_stellar_structure(result)
    system = result.simulation.ic
    n_bodies = system.n

    initial_stellar_parameters = begin
        tmp = Dict()
        for param in propertynames(system.particles[1].structure)
            tmp[param] = [getproperty(system.particles[i].structure, param) for i = 1:n_bodies]
        end

        tmp[:stellar_type] = [t.index for t in tmp[:stellar_type]] 
        tmp
    end

    final_stellar_parameters = begin
        tmp = Dict()
        for param in propertynames(result.ode_params)
            # tmp[param] = convert.(Unitful.Quantity, getproperty(result.ode_params, param))
            tmp[param] =  getproperty(result.ode_params, param)
        end

        tmp[:stellar_type] = pop!(tmp, :stellar_types)
        tmp[:m] = pop!(tmp, :M)
        tmp
    end

    for (k, v) in initial_stellar_parameters
        !(k in keys(final_stellar_parameters)) && continue
        final_stellar_parameters[k] = final_stellar_parameters[k] .* unit(v[1])
    end

    # elements, structure, quantities = setup_simulation_solution(n_steps, n_bodies, n_binaries)

    radius_matrix = Matrix{typeof((1.0u"m"))}(undef, n_bodies, 2)
    mass_matrix = Matrix{typeof(1.0u"kg")}(undef, n_bodies, 2)
    stellar_type_matrix = Matrix{Int8}(undef, n_bodies, 2)
    lum_matrix = Matrix{typeof(1.0u"Lsun")}(undef, n_bodies, 2)
    spin_matrix = nothing#Matrix{typeof(1.0u"1/s")}(undef, n_bodies, 2)

    structure = StellarStructure(stellar_type_matrix, mass_matrix, 
                                radius_matrix, spin_matrix, lum_matrix,
                                similar(radius_matrix), similar(mass_matrix), 
                                similar(radius_matrix), similar(mass_matrix))

    for param in keys(initial_stellar_parameters)
        param == :S && continue
        if !(param ∈ keys(final_stellar_parameters))
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 1)
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 2)
        else
            setindex!(getproperty(structure, param), initial_stellar_parameters[param], 1:n_bodies, 1)
            setindex!(getproperty(structure, param), final_stellar_parameters[param], 1:n_bodies, 2)
        end
    end

    return structure
end

"""
    to_solution(result::SimulationResult)

Convert a `SimulationResult` to a `MultiBodySolution` type with easy access to 
state vectors and other quantities.
"""	
function to_solution(result::SimulationResult)

    system = result.simulation.ic
    n_bodies = system.n

    time = result.solution.t .* upreferred(u"s")
    
    n_steps = length(time)

    S_unit = unit(system.particles.S[1][1])

    r = Array{typeof(upreferred(1.0u"m")), 3}(undef, 3, n_bodies, n_steps)
    v = similar(r, typeof(upreferred(1.0u"m/s")))
    S = similar(r, typeof(upreferred(1.0*S_unit)))
    Sv = similar(r, typeof(upreferred(1.0*S_unit/upreferred(u"yr"))))


    r = AxisArray(r; dim=1:3, particle=1:n_bodies, time=time)
    v = AxisArray(v; dim=1:3, particle=1:n_bodies, time=time)
    S = AxisArray(S; dim=1:3, particle=1:n_bodies, time=time)
    Sv = AxisArray(Sv; dim=1:3, particle=1:n_bodies, time=time)

    
    for idx in eachindex(time)
        pos = result.solution.u[idx].x[2][1:3,:] .* upreferred(u"m")
        vel = result.solution.u[idx].x[1][1:3,:] .* upreferred(u"m/s")

        spin = result.solution.u[idx].x[2][4:6,:] .* upreferred(S_unit)
        spin_vel = result.solution.u[idx].x[1][4:6,:] .* upreferred(S_unit/upreferred(u"yr"))
        
        r[:, :, idx] .= pos
        v[:, :, idx] .= vel

        S[:, :, idx] .= spin
        Sv[:, :, idx] .= spin_vel
    end

    ode_solution = Dict(:potential => result.simulation.potential |> values,
                        :retcodes => result.retcode,
                        :timesteps => result.solution.stats.naccept,
                        :func_1_evals => result.solution.stats.nf,
                        :func_2_evals => result.solution.stats.nf2,
                        :runtime => result.runtime)
    attributes = merge(ode_solution, result.simulation.args, result.simulation.diffeq_args)

    stellar_structure = get_stellar_structure(result)

    return MultiBodySolution(system, time, r, v, S, Sv, stellar_structure, attributes, result.ode_params)

end


# """
# TO DO: allow elements to be supplied OR calculate orbital elements given positions
# and velocities. How to generalise to higher-order multiples?
# """
# function multibodysystem(sol::MultiBodySolution, time)

#     time_index = argmin(abs.(time .- sol.t))
#     if !isnothing(structure)
#         quant_time_index = size(structure.m, 2) == length(sol.t) ? time_index : 2
#         masses = structure.m[:,quant_time_index]

#         R = structure.R[:,quant_time_index], 
#         L = structure.L[:,quant_time_index],
#         S = structure.S[:,quant_time_index],
#         R_core = structure.R_core[:,quant_time_index],
#         m_core = structure.m_core[:,quant_time_index],
#         R_env = structure.R_env[:,quant_time_index],
#         m_env = structure.m_env[:,quant_time_index],
#         stellar_types = structure.type[:,quant_time_index]
#     end

#     n_bodies = length(sol.ic.masses)
#     multibodysystem(masses, R = structure.R[:,quant_time_index], 
#                             L = structure.L[:,quant_time_index],
#                             S = structure.S[:,quant_time_index],
#                             R_core = structure.R_core[:,quant_time_index],
#                             m_core = structure.m_core[:,quant_time_index],
#                             R_env = structure.R_env[:,quant_time_index],
#                             m_env = structure.m_env[:,quant_time_index],
#                             stellar_types = Int.(structure.type[:,quant_time_index]),
#                             a = [e.a[quant_time_index] for e in sol.elements],
#                             e = [e.e[quant_time_index] for e in sol.elements], 
#                             ω = [e.ω[quant_time_index] for e in sol.elements], 
#                             i = [e.i[quant_time_index] for e in sol.elements], 
#                             Ω = [e.Ω[quant_time_index] for e in sol.elements], 
#                             ν = [e.ν[quant_time_index] for e in sol.elements], 
#                             hierarchy = sol.initial_conditions.hierarchy,
#                             t0 = time
#                             )
# end


