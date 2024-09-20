function Base.show(io::IO, system::T where T <: MultiBodySystem)
    print(io, "\nN binaries: ")
    show(io, size(system.hierarchy, 1)-1)
    println(io)

    print(io, "N levels: ")
    show(io, length(system.levels))
    println(io)

    print(io, "N particles: ")
    show(io, system.n)
    println(io)

    print(io, "Hierarchy: ")
    show(io, system.hierarchy)
    print(io, "\n\n")

    show(io, system.root)

end

function Base.show(io::IO, binary::T where T <: Binary)

    indent = get(io, :indent, 2 + binary.level*2)

    # print(io, " "^indent, "Binary key: ")
    printstyled(io, " "^indent, "Binary ", color=:blue)
    print(io,  "key: ")

    show(io, binary.key)
    println(io)

    # print(io, " "^indent, "Binary level: ")
    printstyled(io, " "^indent, "Binary ", color=:blue)
    print(io,  "level: ")

    show(io, binary.level)
    println(io)

    # print(io, " "^indent, "Binary parent: ")
    printstyled(io, " "^indent, "Binary ", color=:blue)
    print(io,  "parent: ")

    show(io, binary.parent)
    println(io)

    # println(io, " "^indent, "Binary elements: ")
    printstyled(io, " "^indent, "Binary ", color=:blue)
    printstyled(io, "elements:\n")
    show(IOContext(io, :indent => indent+2), binary.elements)
    println(io)


    printstyled(io, " "^indent, "Children:\n\n", color=:green)
    for child in binary.children
        show(IOContext(io, :indent => indent+4), child)
        println()
    end

end

function Base.show(io::IO, particle::Particle)
    indent = get(io, :indent, 2)


    # print(io,  " "^indent, "Particle key: ")
    printstyled(io,  " "^indent, "Particle ", color=:yellow)
    print(io, "key: ")
    print(io, particle.key)
    println(io)

    # print(io,  " "^indent, "Particle parent: ")
    printstyled(io,  " "^indent, "Particle ", color=:yellow)
    print(io, "parent: ")
    print(io, particle.parent)
    println(io)

    # print(io,  " "^indent, "Particle mass: ")
    printstyled(io,  " "^indent, "Particle ", color=:yellow)
    print(io, "mass: ")
    print(io, particle.mass)
    println(io)


    # print(io,  " "^indent, "Particle type: ")
    printstyled(io,  " "^indent, "Particle ", color=:yellow)
    print(io, "type: ")
    print(io, stellar_type_index[particle.structure.stellar_type.index])
    println(io)

    # println(io,  " "^indent, "Particle structure: ")
    printstyled(io,  " "^indent, "Particle ", color=:yellow)
    print(io, "structure: ")
    show(io, particle.structure)

end

function Base.show(io::IO, elements::OrbitalElements)#{T, T, T, T, T, T, T}) where T <: Number
    if elements.a isa AbstractArray
        # show(io, elements)
        return
    end
    indent = get(io, :indent, 2)
    els = propertynames(elements)
    print(io, " "^indent, "|")

    for el in els
        print(io, " $el: ")
        prop = getproperty(elements, el)
        prop = prop isa Quantity ? round(prop.val, digits=2)*unit(prop) : round(prop, digits=2)
        # prop = round(prop.val, digits=2)*unit(prop)
        show(io, prop)
        print(io, " | ")
    end 
    
    println(io)
end

# function Base.show(io::IO, elements::OrbitalElements{T, T, T, T, T, T, T}) where T <: Number

function Base.show(io::IO, structure::StellarStructure)
    # if structure.R isa AbstractArray
    #     # show(io, elements)
    #     return
    # end
    indent = get(io, :indent, 0) + 2

    str_names = propertynames(structure)
    print(io, " "^indent, "|")
    for str in str_names
        # str === :type && continue
        print(io, " $str: ")
        prop = getproperty(structure, str)
        if prop isa AbstractArray
            # show(io, prop)
            print("["*join(prop,",")*"]")
        elseif prop isa StellarType
            show(io, prop.index)
        else
            prop = prop isa Quantity ? round(prop.val, digits=2)*unit(prop) : round(prop, digits=2)
            show(io, prop)
        end
        print(io, " | ")
        # println(io)   
    end
    println(io)
end

function Base.show(io::IO, sim::MultiBodySimulation)

    println(io, "\nSimulation setup\n-------------------------------")

    print(io, "\nTime span: ")
    timespan = u"kyr".(sim.tspan .* u"s")
    print(io, timespan)

    println("\n")

    println(io, "Arguments ")
    for (arg, val) in sim.args
        arg == :potential && continue 
        @printf(io, "   %-16s %s", "$arg", "$val")
        # @printf(io, "%7s ", "$val")
        println(io)
    end

    println(io)

    print(io, "DiffEq arguments: \n")
    for (arg, val) in sim.diffeq_args
        @printf(io, "   %-16s %s", "$arg", "$val")
        println(io)
    end

    print(io, "\nPotential(s): \n")
    if sim.args[:potential] isa Vector
        for pot in sim.args[:potential]
            pot = nameof(typeof(pot))
            @printf(io, "   %-16s", "$pot")
            println(io)
        end
    else
        pot = nameof(typeof(sim.args[:potential]))
        @printf(io, "   %-16s", "$pot")
    end

    println(io)
end

function Base.show(io::IO, sim::SimulationResult)

    println(io, "Simulation result:\n====================================\n")
    println(io, "Retcodes: ")
    for (k, v) in sim.retcode
        println(io, "   $k", " ", v)
    end

    println(io, "\nRuntime: $(sim.runtime)\n")

    println(io, "Number of datapoints: $(length(sim.solution.t))\n")
    println(io, "ODE Retcode: $(sim.solution.retcode)\n")
    println(io, "$(sim.solution.destats)")
end


function Base.show(io::IO, params::DefaultSimulationParams)
    print(nameof(typeof(params)))
    print(":")
    println()
    for prop in propertynames(params)
        val = getproperty(params, prop) 
        un = unit(val[1])
        val = ustrip(val)

        @printf(io, "   %-16s %s %s", "$prop", "$val", "$un")
        println(io)
    end
end

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