using LinearAlgebra


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
    print(io, stellar_type_index[particle.structure.type.index])
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

function Base.show(io::IO, structure::StellarStructure2)
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
            show(io, prop)
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

# function Base.getproperty(system::MultiBodySystem, prop)

# end

function Base.getindex(system::MultiBodySystem, idx::ParticleIndex)
    # get_particle(system, idx.i)
    system.particles[idx.i]
end

function Base.getindex(system::MultiBodySystem, idx::BinaryIndex)
    # get_binary(system, idx.i)
    system.binaries[idx.i]
end


function elements_to_cartesian(hierarchy, elements, masses)
    r, v = keplers_problem(hierarchy, masses, elements)
    return r, v
end

"""
    multibodysystem(masses::Vector{<:Quantity}; <kwargs>)

Compute and return a multibody system with given structures, orbital elements, and
specified hierarchy.

Structure and element arguments can be given as vectors of length # of bodies or # of binaries
respectively, or as numbers, in which case the same value will be used for each body and/or binary.

...
# Particle structure arguments
- `R = 1.0u"Rsun"`: radius of each particle.
- `S = 0.0u"1/yr"`: spin magnitude of each particle. If given as a negative number, the spin will be calculated using [`stellar_spin`](@ref).
- `L = 1.0u"Lsun"`: luminosity of each particle.
- `types = 1`: type of each particle. See [`FewBodySimulator.stellar_type_index`](@ref).

# Binary elements arguments
- `a = 1.0u"AU"`: semi-major axis.
- `e = 0.1`: eccentricity.
- `ω = 0.0u"rad"`: argument of periapsis.
- `i = 0.0u"rad"`: (mutual) inclination. The inclination of the first binary is with respect to the xy-plane.
- `Ω = 0.0u"rad"`: longitude of the ascending node.
- `ν = (π)u"rad"`: true anomaly, with `ν = 0u"rad"` corresponding to periapsis.

# Other arguments
- `hierarchy`: specification of the hierarchy structure. First element is total number of bodies, followed 
               by number of binaries on each level. If given a number, system is assumed
               to be hierarchichal. 
- `t0 = 0.0u"s"`: time of the system.
- `verbose = false`
...

# Examples
```jldoctest
julia> binary = multibodysystem([1.0, 1.0]u"Msun");
julia> binary = multibodysystem([1.0, 1.3]u"Msun", R=[1.0, 0.13]u"Rsun", a=1.0u"AU", e=0.4);
julia> triple = multibodysystem([1.0, 1.0, 1.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.1, 0.7], i=[π/2, 0.0]u"rad");
julia> quadruple = multibodysystem([1.0, 1.0, 1.0, 1.0]u"Msun");
julia> quadruple = multibodysystem([1.0, 1.0, 1.0, 1.0]u"Msun", hierarchy=[4, 2, 1]);

```
"""
function multibodysystem(masses::Vector{<:Quantity}; R      = 1.0u"Rsun", 
                                                     S      = 0.0u"1/yr", 
                                                     L      = 1.0u"Lsun", 
                                                     R_core = 0.0u"Rsun",
                                                     m_core = 0.0u"Msun",
                                                     R_env  = 0.0u"Rsun",
                                                     m_env  = 0.0u"Msun",
                                                     types  = 1, 
                                                     a      = 1.0u"AU", 
                                                     e      = 0.1, 
                                                     ω      = 0.0u"rad", 
                                                     i      = 0.0u"rad", 
                                                     Ω      = 0.0u"rad", 
                                                     ν      = (π)u"rad", 
                                                     t0::Quantity = 0.0u"s", 
                                                     verbose::Bool = false, 
                                                     hierarchy = [length(masses), repeat([1], length(masses)-1)...], 
                                                     extras::Dict=Dict())
    n_bodies = length(masses)
    n_bins = n_bodies - 1

    R = ifelse(R isa Number, repeat([R], n_bodies), R)
    S = ifelse(S isa Number, repeat([S], n_bodies), S)
    L = ifelse(L isa Number, repeat([L], n_bodies), L)
    R_core = ifelse(R_core isa Number, repeat([R_core], n_bodies), R_core)
    m_core = ifelse(m_core isa Number, repeat([m_core], n_bodies), m_core)
    R_env  = ifelse(R_env  isa Number, repeat([R_env ], n_bodies), R_env )
    m_env  = ifelse(m_env  isa Number, repeat([m_env ], n_bodies), m_env )

    a = ifelse(a isa Number, repeat([a], n_bins), a)
    e = ifelse(e isa Number, repeat([e], n_bins), e)
    ω = ifelse(ω isa Number, repeat([ω], n_bins), ω)
    i = ifelse(i isa Number, repeat([i], n_bins), i)
    Ω = ifelse(Ω isa Number, repeat([Ω], n_bins), Ω)
    ν = ifelse(ν isa Number, repeat([ν], n_bins), ν)

    if any(x -> x < zero(x), S)
        idx = findall(x -> x < zero(x), S)
        S[idx] .= stellar_spin.(masses[idx], R[idx])
    end

    types = ifelse(types isa Number, repeat([types], n_bodies), types)
    multibodysystem(masses, R, S, L, types, R_core, m_core, R_env, m_env, a, e, ω, i, Ω, ν, hierarchy, t0, verbose=verbose, extras=extras)
end


function multibodysystem(masses::Vector{<:Quantity}, 
                         R::Vector{<:Quantity}, S::Vector{<:Quantity}, L::Vector{<:Quantity}, types::Vector{Int}, 
                         R_core::Vector{<:Quantity}, m_core::Vector{<:Quantity}, R_env::Vector{<:Quantity}, m_env::Vector{<:Quantity},
                         a::Vector{<:Quantity}, e::Vector{<:Real},     ω::Vector{<:Quantity}, 
                         i::Vector{<:Quantity}, Ω::Vector{<:Quantity}, ν::Vector{<:Quantity}, 
                         hierarchy::Vector{Int}, t0::Quantity; 
                         verbose::Bool = false, extras::Dict = Dict())

    n_bins = sum(hierarchy[2:end])
    n_particles = length(masses)
    
    @assert length(masses) == length(S) == length(L) == length(types) "Must give structural property for each particles."
    @assert n_bins == n_particles - 1 "Number of binary elements must equal N particles - 1."

    elements = OrbitalElements[]
    for idx = 1:n_bins
        push!(elements, OrbitalElements(a[idx], 0.0u"d", e[idx], ω[idx], i[idx], Ω[idx], ν[idx]))
    end

    structures = StellarStructure2[]
    for idx = 1:n_particles
        stellar_type = stellar_type_from_index(types[idx])
        # push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx]))
        push!(structures, StellarStructure2(stellar_type, masses[idx], R[idx], S[idx], L[idx],
                                            R_core[idx], m_core[idx], R_env[idx], m_env[idx]))
    end


    MultiBodySystem(masses, hierarchy, elements, structures, t0, verbose=verbose)
end



function MultiBodySystem(masses::Vector, hierarchy::Vector, 
                         elements::AbstractVector{T} where T <: OrbitalElements,
                         structures::AbstractVector{T} where T <: StellarStructure2,
                         t = 0.0u"yr"; verbose=false)

    hierarchy_mat = hierarchy_matrix(hierarchy)
    positions, velocities = keplers_problem(hierarchy_mat, masses, elements)
    A_matrix = mass_ratio_matrix(hierarchy_mat, masses)
    Ainv = inv(A_matrix)

    positions = permutedims(Ainv*positions)
    velocities = permutedims(Ainv*velocities)

    n_bins = sum(hierarchy[2:end])
    n_bodies = length(masses)

    level = length(hierarchy[2:end]) - 1
    levels = Int[level]
    added_ids = Int[]
    parent_key = if iszero(level)
                    -1
                 else
                    hierarchy[2] + 1
                 end
    first_binary = let idx = hierarchy_mat[1,:] .|> abs .|> Bool
                        verbose && println("Initial binary on level $level")
                        append!(added_ids, findall(idx))
                        binary_elements = elements[1]# + 1 + bin]
                        μ = sum(masses[idx])
                        P = orbital_period(binary_elements.a, μ) |> u"d"
                        binary_elements = let e = binary_elements
                            OrbitalElements(e.a, P, e.e, e.ω, e.i, e.Ω, e.ν)
                        end

                        # sibling_key = if 

                        create_binary(positions[:,idx], velocities[:,idx], 
                                      masses[idx], structures[idx], binary_elements, 
                                      [1, 2], 1, level, parent_key) 
                    end

    i = 2
    # bin = -1
    binaries = Binary[first_binary]
    tot_bodies = 2
    binary_key = 2

    while i < n_bodies
        binary_elements = elements[i]# + 1 + bin]
        idx = hierarchy_mat[i,:] .|> abs .|> Bool
        μ = sum(masses[idx])
        P = orbital_period(binary_elements.a, μ) |> u"d"
        binary_elements = let e = binary_elements
            OrbitalElements(e.a, P, e.e, e.ω, e.i, e.Ω, e.ν)
        end


        # new binary
        if iszero(first(idx))# && sum(idx) == 2
            verbose && println("New binary on level $level")
            parent_key = hierarchy[level+1] + 1

            ids = [tot_bodies+1, tot_bodies+2]
            binary = create_binary(positions[:,idx], velocities[:,idx], 
                                   masses[idx], structures[idx], binary_elements, ids, binary_key, level, parent_key)
            # @show binary
            push!(binaries, binary)
            append!(added_ids, findall(idx))
            tot_bodies += 2
            binary_key += 1
        else 
            if i == n_bodies-1  # last level
                if tot_bodies == n_bodies # two binaries
                    level -= 1
                    verbose && println("Double binary on level $level")
                    child1 = binaries[end-1]
                    parent_key = iszero(level) ? -1 : child1.parent + 1#iszero(level) ? -1 : hierarchy[level] + 1  
                    child2 = binaries[end]
                    binary = create_binary(child1, child2, binary_elements, level, binary_key, parent_key)
                    push!(binaries, binary)
                    binary_key += 1
                    push!(levels, level)
                else  # one binary and one particle
                    level -= 1
                    verbose && println("Binary-particle binary on level $level")
                    child1 = last(binaries)
                    parent_key = iszero(level) ? -1 : child1.parent + 1#hierarchy[level] + 1
                    child_2_idx = last(findall(idx))
                    child2 = Particle(ParticleIndex(tot_bodies+1), BinaryIndex(binary_key), child1.key, masses[child_2_idx],
                                      positions[:,child_2_idx], 
                                      velocities[:, child_2_idx], 
                                      structures[child_2_idx], Dict())
    
                    binary = create_binary(child1, child2, binary_elements, level, binary_key, parent_key)
                    push!(binaries, binary)
                    tot_bodies += 1
                    binary_key += 1
                    push!(levels, level)
                    append!(added_ids, findall(idx)[end])
                end
            else 
                if all(x -> x in added_ids, findall(idx)) # connect double binary
                    # level -= 1
                    verbose && println("Double binary on level $level")
                    child1 = binaries[end-1]
                    parent_key = iszero(level) ? -1 : child1.parent + 1#iszero(level) ? -1 : hierarchy[level] + 1  
                    child2 = binaries[end]
                    binary = create_binary(child1, child2, binary_elements, level, binary_key, parent_key, binary_key)
                    push!(binaries, binary)
                    push!(levels, level)
                    binary_key += 1
                    level -= 1
                else # particle-binary binary
                    level -= 1
                    verbose && println("Binary-particle binary on level $level")
                    child1 = last(binaries)
                    parent_key = iszero(level) ? -1 : child1.parent.i + 1#hierarchy[level] + 1
                    child_2_idx = last(findall(idx))
                    child2 = Particle(ParticleIndex(tot_bodies+1), BinaryIndex(binary_key), child1.key, masses[child_2_idx],
                                      positions[:,child_2_idx], 
                                      velocities[:, child_2_idx], 
                                      structures[child_2_idx], Dict())
    
                    binary = create_binary(child1, child2, binary_elements, level, binary_key, parent_key)
                    push!(binaries, binary)
                    tot_bodies += 1
                    binary_key += 1
                    push!(levels, level)
                    append!(added_ids, findall(idx)[end])
                end
            end
        end

        i += 1
    end
    # binaries = SA[reverse(binaries)...]
    binaries = Dict{Int, Binary}(b.key.i => b for b in binaries)
    root_bin = [bin for bin in values(binaries) if bin.level == 0][1]


    bodies = get_particles(root_bin)
    bodies = Dict{Int, Particle}(p.key.i => p for p in bodies)
    # bodies = SA[bodies[sortperm([b.key.i for b in bodies])]...]

    levels = SA[sort(unique(levels))...]

    MultiBodySystem(n_bodies, t, bodies, binaries, levels, root_bin, hierarchy, nothing)
end

function create_binary(positions, velocities, masses, structures, elements, 
                      particle_keys, binary_key, level, parent_key, sibling_key=nothing)
    children = Particle[]

    for i ∈ eachindex(masses)
        sibling_key = findall(particle_keys .!= particle_keys[i])[1]
        child = Particle(ParticleIndex(particle_keys[i]), 
                         BinaryIndex(binary_key), 
                         ParticleIndex(sibling_key),
                         masses[i], SA[positions[:,i]...], 
                         SA[velocities[:,i]...], 
                         structures[i], Dict())
        push!(children, child)
    end

    com_position = centre_of_mass(positions, masses)
    com_velocity = centre_of_mass_velocity(velocities, masses)
    
    Binary(BinaryIndex(binary_key), level, BinaryIndex(parent_key), sibling_key,
           SA[children...], SA[[c.key.i for c in children]...], com_position, com_velocity, masses, elements)
end

function create_binary(child1::Binary, child2::Particle, elements, level, binary_key, parent_key, sibling_key=nothing)
    bodies = [get_particles(child1)..., child2]
    masses = [body.mass for body in bodies]
    positions = [body.position for body in bodies]
    velocities = [body.velocity for body in bodies]

    com_position = centre_of_mass(positions, masses)
    com_velocity = centre_of_mass_velocity(velocities, masses)
    
    binary_masses = [sum(masses[1:end-1]), masses[end]]
    children = SA[child1, child2]

    nested_children = SA[child2.key.i, get_particle_ids(child1)...]
    Binary(BinaryIndex(binary_key), level, BinaryIndex(parent_key), sibling_key, 
           children, nested_children, com_position, com_velocity, binary_masses, elements)
end

function create_binary(child1::Binary, child2::Binary, elements, level, binary_key, parent_key, sibling_key=nothing)
    children = SA[child1, child2]
    bodies = [get_particles(child1)..., get_particles(child2)...]
    masses = [body.mass for body in bodies]
    positions = [body.position for body in bodies]
    velocities = [body.velocity for body in bodies]
    com_position = centre_of_mass(positions, masses)
    com_velocity = centre_of_mass_velocity(velocities, masses)
    
    child1_bodies = length(child1.children)
    binary_masses = [sum(masses[1:child1_bodies]), sum(masses[child1_bodies+1:end])]
    nested_children = SA[get_particle_ids(child1)..., get_particle_ids(child2)...]
    Binary(BinaryIndex(binary_key), level, BinaryIndex(parent_key), sibling_key, 
           children, nested_children, com_position, com_velocity, binary_masses, elements)
end



get_particle_ids(binary::Binary) = binary.nested_children
get_particle_ids(particle::Particle) = SA[particle.key.i]

function get_particles(binary::Binary, particles=Particle[])
    children = binary.children
    # particles = Particle[]
    for child in children
        if child isa Particle
            push!(particles, child)
        else
            get_particles(child, particles)
        end
    end

    return particles
end


function get_particles(system::MultiBodySystem)
    bins = values(system.binaries)
    particles = Particle[]
    for bin in bins
        append!(particles, get_particles(bin))
    end

    return particles[sortperm([p.key.i for p in particles])]
end

function get_binaries(binary::Binary)
    children = binary.children
    binaries = Binary[]
    for child in children
        if child isa Binary
            push!(binaries, child)
            append!(binaries, get_binaries(child))
        else
            continue
        end
    end

    return binaries
end

function get_binaries(system::MultiBodySystem)
    root = system.root
    binaries = Binary[root]
    append!(binaries, get_binaries(root))

    # levels = [bin.level for bin in binaries]
    keys = [bin.key.i for bin in binaries]
    binaries = binaries[sortperm(keys)]
    return binaries
end

# match_key(x) = x.ke

function get_binary(system::MultiBodySystem, key::Int)
    system.binaries[findfirst(x -> x.key.i == key, system.binaries)]
end

function get_binary(binary::Binary, key::Int)

    @inbounds for child in binary.children
        if child isa Binary && child.key.i == key
            return child
        end
    end

    @info "Binary with key $key not found."
    return nothing
end

function get_particle(system::MultiBodySystem, key::Int)

    # match_key(x::Particle)::Bool = x.key.i == key
    system.particles[findfirst(x -> x.key.i == key, system.particles)]
    # @inbounds for body in system.particles
    #     if body.key.i == key
    #         return body
    #     end
    # end

    # @info "Particle with key $key not found."
    # return nothing
end

function get_particle(binary::Binary, key::Int)
    # findfirst(x -> x.key.i == key, binary.particles)
    # local child::Particle
    @inbounds for child in binary.children
        if child isa Particle && child.key.i == key
            return child
        end
    end

    # @info "Particle with key $key not found."
    # return nothing
    throw(BoundsError)
end

function getparent(system::MultiBodySystem, object::Union{Binary, Particle})
    parent_key = object.parent
    if parent_key == -1
        tpe = nameof(typeof(object))
        @info "$tpe is root."
    end

    for (level, bin) in system.binaries
        # @show level, bin.key
        if bin.key == parent_key
            return bin
        end
    end

end

"""
    multibodysystem(system::MultiBodySystem; new_params...)

Remake a given system with new parameters. 
"""
function multibodysystem(system::MultiBodySystem; new_params...)

    possible_kwargs = [:R, :S, :L, :types, :a, :e, :ω, :i, :Ω, :ν, 
                       :hierarchy, :t0, :verbose, :extras]

    new_kwargs = Dict(new_params)
    unchanched_kwargs = [kw for kw in possible_kwargs if (!in(kw, keys(new_kwargs)))]


    particle_keys = collect(keys(system.particles)) |> sort
    binary_keys = collect(keys(system.binaries)) |> sort
    all_args = Dict{Symbol, Any}()
    for arg in unchanched_kwargs
        # @show arg
        if any(arg .=== (:R, :S, :L))
            quantity = [getproperty(system.particles[p].structure, arg) for p in particle_keys] |> Vector
            all_args[arg] = quantity
        elseif any(arg .=== (:a, :e, :ω, :i, :Ω, :ν))
            element = [getproperty(system.binaries[p].elements, arg) for p in binary_keys] |> Vector
            all_args[arg] = element
        elseif arg === :hierarchy
            all_args[arg] = system.hierarchy
        elseif arg === :t0
            all_args[arg] = system.time
        else
            continue
        end
    end


    masses = pop!(new_kwargs, :masses, [system.particles[p].mass for p in particle_keys]) |> Vector

    args = merge(new_kwargs, all_args)

    multibodysystem(masses; args...)
end


"""
m = FewBodySimulator.multibodysystem([1.0, 1.0, 1.0, 1.0]u"Msun", a=[0.1, 0.1, 5.0]u"AU", e=[0.1, 0.1, 0.1], ω=[0.0, 0.0, 0.0]u"rad", 
                                      i = [0.0, 0.0, 0.0]u"rad", Ω=[0.0, 0.0, 0.0]u"rad", 
                                      ν=[0.0, 0.0, 0.0]u"rad", hierarchy=[4, 2, 1]);


m = FewBodySimulator.multibodysystem([1.0, 1.0]u"Msun", a=[0.1]u"AU", e=[0.1], 
                                ω=[0.0]u"rad", i=[0.0]u"rad", Ω=[0.0]u"rad", 
                                ν=[0.0]u"rad", hierarchy=[2,1])
"""