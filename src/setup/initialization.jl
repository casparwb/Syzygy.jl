using LinearAlgebra: norm



"""
    multibodysystem(masses::Vector{<:Quantity}; <kwargs>)

Compute and return a multibody system with given structures, orbital elements, and
specified hierarchy.

Structure and element arguments can be given as vectors of length # of bodies or # of binaries
respectively, or as numbers, in which case the same value will be used for each body and/or binary.

...
# Binary elements arguments
- `a = 1.0u"AU"`: semi-major axis.
- `e = 0.1`: eccentricity.
- `ω = 0.0u"rad"`: argument of periapsis.
- `i = 0.0u"rad"`: (mutual) inclination. The inclination of the first binary is with respect to the xy-plane.
- `Ω = 0.0u"rad"`: longitude of the ascending node.
- `ν = (π)u"rad"`: true anomaly, with `ν = 0u"rad"` corresponding to periapsis.

# Particle structure arguments
- `R = 1.0u"Rsun"`: radius of each particle.
- `S = 0.0u"1/yr"`: spin magnitude of each particle. If given as a negative number, the spin will be calculated using [`stellar_spin`](@ref).
- `L = 1.0u"Lsun"`: luminosity of each particle.
- `stellar_types = 1`: type of each particle. See [`Syzygy.stellar_type_index`](@ref).
- `R_core = 0.0u"Rsun"`: stellar core radius. Only used if tidal potential is included.
- `m_core = 0.0u"Msun"`: stellar core mass. Only used if tidal potential is included.
- `R_env = 0.0u"Rsun"`: stellar envelope radius. Only used if tidal potential is included.
- `m_env = 0.0u"Msun"`: stellar envelope mass. Only used if tidal potential is included.

# Other arguments
- `hierarchy`: specification of the hierarchy structure. First element is total number of bodies, followed 
               by number of binaries on each level. If given a number, system is assumed
               to be hierarchichal. 
- `time = 0.0u"s"`: time of the system.
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
function multibodysystem(masses::Vector{<:Unitful.Mass}; R  = 1.0u"Rsun", 
                                                         S      = 0.0u"1/yr", 
                                                         L      = 1.0u"Lsun", 
                                                         R_core = unit(R[1])*(0.0),
                                                         m_core = unit(masses[1])*(0.0),
                                                         R_env  = unit(R[1])*(0.0),
                                                         m_env  = unit(masses[1])*(0.0),
                                                         stellar_types  = 1, 
                                                         a      = 1.0u"AU", 
                                                         e      = 0.1, 
                                                         ω      = 0.0u"rad", 
                                                         i      = 0.0u"rad", 
                                                         Ω      = 0.0u"rad", 
                                                         ν      = (π)u"rad", 
                                                         time::Quantity = 0.0u"s", 
                                                         verbose::Bool = false, 
                                                         hierarchy = [length(masses), repeat([1], length(masses)-1)...])
    n_bodies = length(masses)
    n_bins = n_bodies - 1

    R = ifelse(R isa Number, repeat([R], n_bodies), R)
    L = ifelse(L isa Number, repeat([L], n_bodies), L)
    R_core = ifelse(R_core isa Number, repeat([R_core], n_bodies), R_core)
    m_core = ifelse(m_core isa Number, repeat([m_core], n_bodies), m_core)
    R_env  = ifelse(R_env  isa Number, repeat([R_env ], n_bodies), R_env )
    m_env  = ifelse(m_env  isa Number, repeat([m_env ], n_bodies), m_env )

    S = if S isa Number 
        Sv = [[S, zero(S), zero(S)] for i in 1:n_bodies]
    elseif S isa Vector{<:Number}
        Sv = [[ss, zero(ss), zero(ss)] for ss in S]
    else
        S
    end

    a = ifelse(a isa Number, repeat([a], n_bins), a)
    e = ifelse(e isa Number, repeat([e], n_bins), e)
    ω = ifelse(ω isa Number, repeat([ω], n_bins), ω)
    i = ifelse(i isa Number, repeat([i], n_bins), i)
    Ω = ifelse(Ω isa Number, repeat([Ω], n_bins), Ω)
    ν = ifelse(ν isa Number, repeat([ν], n_bins), ν)

    # if any(x -> x < zero(x), S)
    #     idx = findall(x -> x < zero(x), S)
    #     S[idx] .= stellar_spin.(masses[idx], R[idx])
    # end

    stellar_types = ifelse(stellar_types isa Number, repeat([stellar_types], n_bodies), stellar_types)
    multibodysystem(masses, R, S, L, stellar_types, R_core, m_core, R_env, m_env, a, e, ω, i, Ω, ν, hierarchy, time, verbose=verbose)
end


function multibodysystem(masses::Vector{<:Unitful.Mass}, 
                         R::Vector{<:Unitful.Length}, 
                         S::Vector{<:Union{<:Vector{<:Quantity}, <:Quantity}}, 
                         L::Vector{<:Unitful.Power}, 
                         stellar_types::Vector{Int}, 
                         R_core::Vector{<:Unitful.Length}, 
                         m_core::Vector{<:Unitful.Mass}, 
                         R_env::Vector{<:Unitful.Length}, 
                         m_env::Vector{<:Unitful.Mass},
                         a::Vector{<:Unitful.Length}, e::Vector{<:Real}, ω::Vector{<:Quantity}, 
                         i::Vector{<:Quantity}, Ω::Vector{<:Quantity}, ν::Vector{<:Quantity}, 
                         hierarchy::Vector{Int}, time::Quantity; 
                         verbose::Bool = false)

    n_bins = sum(hierarchy[2:end])
    n_particles = length(masses)
    
    @assert length(masses) == length(S) == length(L) == length(stellar_types) "Must give structural property for each particles."
    @assert n_bins == n_particles - 1 "Number of binary elements must equal N particles - 1."

    elements = OrbitalElements[]
    for idx = 1:n_bins
        push!(elements, OrbitalElements(a[idx], 0.0u"d", e[idx], ω[idx], i[idx], Ω[idx], ν[idx]))
    end

    # if unit(masses)[1] != unit()
    structures = StellarStructure[]
    for idx ∈ 1:n_particles
        stellar_type = stellar_type_from_index(stellar_types[idx])
        # push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx]))
        push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx],
                                            R_core[idx], m_core[idx], R_env[idx], m_env[idx]))
    end


    multibodysystem(masses, hierarchy, elements, structures, time, verbose=verbose)
end

# function multibodysystem(masses::Vector{<:Unitful.Mass}, 
#                          R::Vector{<:Unitful.Length}, 
#                          S::Vector{Vector{<:Quantity}}, 
#                          L::Vector{<:Quantity}, 
#                          stellar_types::Vector{Int}, 
#                          R_core::Vector{<:Unitful.Length}, 
#                          m_core::Vector{<:Unitful.Mass}, 
#                          R_env::Vector{<:Unitful.Length}, 
#                          m_env::Vector{<:Unitful.Mass},
#                          a::Vector{<:Unitful.Length}, e::Vector{<:Real}, ω::Vector{<:Quantity}, 
#                          i::Vector{<:Quantity}, Ω::Vector{<:Quantity}, ν::Vector{<:Quantity}, 
#                          hierarchy::Vector{Int}, time::Quantity; 
#                          verbose::Bool = false)

#     n_bins = sum(hierarchy[2:end])
#     n_particles = length(masses)
    
#     @assert length(masses) == length(S) == length(L) == length(stellar_types) "Must give structural property for each particles."
#     @assert n_bins == n_particles - 1 "Number of binary elements must equal N particles - 1."

#     elements = OrbitalElements[]
#     for idx = 1:n_bins
#         push!(elements, OrbitalElements(a[idx], 0.0u"d", e[idx], ω[idx], i[idx], Ω[idx], ν[idx]))
#     end

#     # if unit(masses)[1] != unit()

#     structures = StellarStructure[]
#     for idx = 1:n_particles
#         stellar_type = stellar_type_from_index(stellar_types[idx])
#         # push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx]))
#         push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx],
#                                             R_core[idx], m_core[idx], R_env[idx], m_env[idx]))
#     end

#     multibodysystem(masses, hierarchy, elements, structures, time, verbose=verbose)
# end


function multibodysystem(masses::Vector, hierarchy::Vector, 
                         elements::AbstractVector{T} where T <: OrbitalElements,
                         structures::AbstractVector{T} where T <: StellarStructure,
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
                                      structures[child_2_idx])
    
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
                                      structures[child_2_idx])
    
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

    pairs = Tuple{Int, Int}[]
    for i in 1:n_bodies
        for j = (i+1):n_bodies
            if i != j
                push!(pairs, (i, j))
            end
        end
    end

    MultiBodySystem(n_bodies, t, bodies, pairs, binaries, levels, root_bin, hierarchy)
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
                         structures[i])
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

function get_particle_ids(binary::Binary)
    binary.nested_children
end
function get_particle_ids(particle::Particle)
    SA[particle.key.i]
end

function get_particles(binary::Binary, particles=Particle[])
    children = binary.children
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

function get_particle(system::MultiBodySystem, key::Int)
    system.particles[key]
end

function get_particle(binary::Binary, key::Int)

    @inbounds for child in binary.children
        if child isa Particle && child.key.i == key
            return child
        end
    end

    throw(BoundsError)
end

"""
    multibodysystem(masses, positions, velocities; kwargs...)

Create a new MultiBodySystem from masses, positions, and velocities. The resulting object
will not have a specific hierarchy, and therefore the binary fields should not be accessed.


# Examples
```jldoctest
julia> masses = rand(3)u"kg"
julia> positions = [rand(3)u"m" for i = 1:3]
julia> velocities = [rand(3)u"m/s" for i = 1:3]
julia> multibodysystem(masses, positions, velocities);
julia> multibodysystem(masses, positions, velocities, R=ones(3)u"m") # you an also supply any of the structure arguments
```

"""
function multibodysystem(masses, positions, velocities;
                        R      = 1.0u"Rsun", 
                        S      = 0.0u"1/yr", 
                        L      = 1.0u"Lsun", 
                        R_core = zero(R[1]),
                        m_core = zero(masses[1]),
                        R_env  = zero(R[1]),
                        m_env  = zero(masses[1]),
                        stellar_types  = 1, 
                        time::Quantity = 0.0u"s", 
                        verbose::Bool = false)

    n_bodies = length(masses)

    particle_keys = 1:n_bodies

    R = ifelse(R isa Number, repeat([R], n_bodies), R)
    S = ifelse(S isa Number, repeat([S], n_bodies), S)
    L = ifelse(L isa Number, repeat([L], n_bodies), L)
    R_core = ifelse(R_core isa Number, repeat([R_core], n_bodies), R_core)
    m_core = ifelse(m_core isa Number, repeat([m_core], n_bodies), m_core)
    R_env  = ifelse(R_env  isa Number, repeat([R_env ], n_bodies), R_env )
    m_env  = ifelse(m_env  isa Number, repeat([m_env ], n_bodies), m_env )
    stellar_types = ifelse(stellar_types isa Number, repeat([stellar_types], n_bodies), stellar_types)

    S = if S isa Number 
        Sv = [[S, zero(S), zero(S)] for i in 1:n_bodies]
    elseif S isa Vector{<:Number}
        Sv = [[ss, zero(ss), zero(ss)] for ss in S]
    else
        S
    end

    # particle_structures = StellarStructure[]
    # for idx in eachindex(masses)
    #     structure = StellarStructure(stellar_types[idx], masses[idx], R[idx], S[idx], L[idx],
    #                                 R_core[idx], m_core[idx], R_env[idx], m_env[idx])
    #     push!(particle_structures, structure)
    # end

    particles = Particle[]
    for idx in eachindex(masses)
        stellar_type = stellar_type_from_index(stellar_types[idx])
        structure = StellarStructure(stellar_type, 
                                     masses[idx], 
                                     R[idx], S[idx], L[idx],
                                     R_core[idx], m_core[idx], 
                                     R_env[idx], m_env[idx])

        p = Particle(ParticleIndex(particle_keys[idx]), 
                     BinaryIndex(-1), -1, 
                     masses[idx], 
                     positions[idx], 
                     velocities[idx], 
                     structure)
        push!(particles, p)
    end

    pairs = Tuple{Int, Int}[]
    for i in 1:n_bodies
        for j = (i+1):n_bodies
            if i != j
                push!(pairs, (i, j))
            end
        end
    end

    # binaries = nothing
    # levels = SA[1]
    # root = Binary(BinaryIndex(1), 1, BinaryIndex(-1), 1, SA[1, 2], SA[particle_keys...], 
    #               centre_of_mass(positions, masses), 
    #               centre_of_mass_velocity(velocities, masses), 
    #               masses, OrbitalElements())
    # hierarchy = [n_bodies, repeat([1], n_bodies-1)...]

    # bodies = Dict{Int, Particle}(i => p for (i, p) in enumerate(particles))
    # binaries = Dict(1 => root)

    particles = Dict{Int, Particle}(i => p for (i, p) in enumerate(particles))

    NonHierarchicalSystem(n_bodies, time, particles, pairs)
end


"""
    multibodysystem(system::MultiBodySystem; new_params...)

Remake a given system with new parameters. 
"""
function multibodysystem(system::MultiBodySystem; new_params...)

    possible_kwargs = [:R, :S, :L, :stellar_types, :a, :e, :ω, :i, :Ω, :ν, 
                       :hierarchy, :time, :verbose]

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
        elseif arg === :time
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
    multibodysystem(sol::MultiBodySolution, time)

Remake a system at a specific point in time from a solution object.

# Examples
```jldoctest
julia> triple = multibodysystem([3, 2, 1]u"Msun")
julia> result = simulate(triple, t_sim=10u"yr", save_everystep=true)
julia> solution = to_solution(result)
julia> triple_new = multibodystem(solution, 5.0u"yr")
```
"""
function multibodysystem(sol::MultiBodySolution, time)

    time_index = argmin(abs.(time .- sol.t))

    positions = sol.r[time=time]
    velocities = sol.v[time=time]

    structure = sol.structure
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

    return multibodysystem(masses, positions, velocities, R=R, L=L, S=S, 
                           R_core=R_core, m_core=m_vore, R_env=R_env, m_env=m_env, 
                           stellar_types=stellar_types, time=time)

end

