using LinearAlgebra, StaticArrays, AxisArrays   



function Base.show(io::IO, sol::MultiBodySolution)
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

        tmp[:type] = [t.index for t in tmp[:type]] 
        tmp
    end

    final_stellar_parameters = begin
        tmp = Dict()
        for param in propertynames(result.ode_params)
            # tmp[param] = convert.(Unitful.Quantity, getproperty(result.ode_params, param))
            tmp[param] =  getproperty(result.ode_params, param)
        end

        tmp[:type] = pop!(tmp, :stellar_types)
        tmp[:m] = pop!(tmp, :M)
        tmp
    end

    # elements, structure, quantities = setup_simulation_solution(n_steps, n_bodies, n_binaries)

    radius_matrix = Matrix{typeof((1.0u"m"))}(undef, n_bodies, 2)
    mass_matrix = Matrix{typeof(1.0u"kg")}(undef, n_bodies, 2)
    stellar_type_matrix = Matrix{Int8}(undef, n_bodies, 2)
    lum_matrix = Matrix{typeof(1.0u"Lsun")}(undef, n_bodies, 2)
    spin_matrix = Matrix{typeof(1.0u"1/s")}(undef, n_bodies, 2)

    structure = StellarStructure(stellar_type_matrix, mass_matrix, 
                                radius_matrix, spin_matrix, lum_matrix,
                                similar(radius_matrix), similar(mass_matrix), 
                                similar(radius_matrix), similar(mass_matrix))

    for param in keys(initial_stellar_parameters)
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

    r = Array{typeof(upreferred(1.0u"m")), 3}(undef, 3, n_bodies, n_steps)
    v = similar(r, typeof(1.0u"m/s"))

    r = AxisArray(r; dim=1:3, particle=1:n_bodies, time=time)
    v = AxisArray(v; dim=1:3, particle=1:n_bodies, time=time)

    for idx in eachindex(time)
        pos = result.solution.u[idx].x[2] .* upreferred(1.0u"m")
        vel = result.solution.u[idx].x[1] .* upreferred(1.0u"m/s")

        r[:, :, idx] .= pos
        v[:, :, idx] .= vel
    end

    ode_solution = Dict(:potential => result.simulation.system.potential |> values,
                        :retcodes => result.retcode,
                        :timesteps => result.solution.destats.naccept,
                        :func_1_evals => result.solution.destats.nf,
                        :func_2_evals => result.solution.destats.nf2,
                        :runtime => result.runtime)
    attributes = merge(ode_solution, result.simulation.args, result.simulation.diffeq_args)

    stellar_structure = get_stellar_structure(result)

    return MultiBodySolution(system, time, r, v, stellar_structure, attributes, result.ode_params)

end


"""
TO DO: allow elements to be supplied OR calculate orbital elements given positions
and velocities. How to generalise to higher-order multiples?
"""
function multibodysystem(sol::MultiBodySolution, time; elements=nothing)

    time_index = argmin(abs.(time .- sol.t))
    if !isnothing(structure)
        quant_time_index = size(structure.m, 2) == length(sol.t) ? time_index : 2
        masses = structure.m[:,quant_time_index]

        R = structure.R[:,quant_time_index], 
        L = structure.L[:,quant_time_index],
        S = structure.S[:,quant_time_index],
        R_core = structure.R_core[:,quant_time_index],
        m_core = structure.m_core[:,quant_time_index],
        R_env = structure.R_env[:,quant_time_index],
        m_env = structure.m_env[:,quant_time_index],
        stellar_types = structure.type[:,quant_time_index]
    end

    n_bodies = length(sol.ic.masses)
    # multibodysystem(masses, R = structure.R[:,quant_time_index], 
    #                         L = structure.L[:,quant_time_index],
    #                         S = structure.S[:,quant_time_index],
    #                         R_core = structure.R_core[:,quant_time_index],
    #                         m_core = structure.m_core[:,quant_time_index],
    #                         R_env = structure.R_env[:,quant_time_index],
    #                         m_env = structure.m_env[:,quant_time_index],
    #                         stellar_types = Int.(structure.type[:,quant_time_index]),
    #                         a = [e.a[quant_time_index] for e in sol.elements],
    #                         e = [e.e[quant_time_index] for e in sol.elements], 
    #                         ω = [e.ω[quant_time_index] for e in sol.elements], 
    #                         i = [e.i[quant_time_index] for e in sol.elements], 
    #                         Ω = [e.Ω[quant_time_index] for e in sol.elements], 
    #                         ν = [e.ν[quant_time_index] for e in sol.elements], 
    #                         hierarchy = sol.initial_conditions.hierarchy,
    #                         t0 = time
    #                         )
end


