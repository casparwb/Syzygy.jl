using LinearAlgebra: norm

const multibodysystem_parameter_aliases = Dict(:radii => :R,  
                                               :spins => :S,  
                                               :luminositites => :L, 
                                               :core_radii => :R_core,  
                                               :core_masses => :m_core,
                                               :envelope_radii => :R_env,  
                                               :envelope_masses => :m_env,
                                               :sma => :a,  
                                               :semi_major_axis => :a,  
                                               :semi_major_axes => :a,
                                               :semimajor_axis => :a,
                                               :semimajor_axes => :a,
                                               :eccentricity => :e,  
                                               :eccentricities => :e,  
                                               :argument_of_periapsis => :ω,  
                                               :argument_of_pericenter => :ω,
                                               :argument_of_pericentre => :ω,  
                                               :aop => :ω,
                                               :inclination => :i,
                                               :longitude_of_ascending_node => :Ω,  
                                               :loan => :Ω,
                                               :true_anomaly => :ν,  
                                               :true_anomalies => :ν,  
                                               :anomaly => :ν,  
                                               :anomalies => :ν,
                                               :R => :R,
                                               :S => :S,
                                               :L => :L,
                                               :R_core => :R_core,
                                               :m_core => :m_core,
                                               :R_env => :R_env,
                                               :m_env => :m_env,
                                               :a => :a,
                                               :ω => :ω,
                                               :i => :i,
                                               :Ω => :Ω,
                                               :ν => :ν,
                                               :e => :e,
                                               :stellar_types => :stellar_types)

function parse_multibodysystem_args(keyword_arguments)
    keyword_arguments = Dict(keyword_arguments)
    default_param_values = Dict{Symbol, Any}(:R      => 1.0u"Rsun", 
                                             :S      => 0.0u"1/yr", 
                                             :L      => 1.0u"Lsun", 
                                             :R_core => 0.0u"Rsun",
                                             :m_core => 0.0u"Msun",
                                             :R_env  => 0.0u"Rsun",
                                             :m_env  => 0.0u"Msun",
                                             :a      => 1.0u"AU", 
                                             :ω      => 0.0u"rad", 
                                             :i      => 0.0u"rad", 
                                             :Ω      => 0.0u"rad", 
                                             :ν      => (π)u"rad",
                                             :e      => 0.1, 
                                             :stellar_types  => 1)

    args_out = copy(default_param_values)

    for (k, v) in keyword_arguments
        param_key = get(multibodysystem_parameter_aliases, k, nothing)
        isnothing(param_key) && throw(ArgumentError("Parameter $k is not valid. See `keys(Syzygy.multibodysystem_parameter_aliases)` for all possible parameters."))
        args_out[param_key] = v
    end

    return args_out
end

function check_units(masses, radii, spins, luminosities, 
                     core_radii, core_masses, 
                     envelope_radii, envelope_masses, 
                     semi_major_axes, argument_of_periapses, inclinations, 
                     longitude_of_ascending_nodes, 
                     true_anomalies)

    𝐋, 𝐌, 𝐓 = Unitful.𝐋, Unitful.𝐌, Unitful.𝐓
    Length, Mass, Time = 𝐋, 𝐌, 𝐓
    Power = Length^2*Mass/Time^3

    @assert all(x -> dimension(x) == Mass, masses) "Unit of masses is not correct. Got $(unit.(masses)), expected mass"
    @assert all(x -> dimension(x) == Length, radii) "Unit of radii is not correct. Got $(unit.(radii)), expected length"
    # @assert all(x -> dimension(x) == 1/Time, spins) "Unit of spins is not correct. Got $(unit.(spins)), expected mass"
    @assert all(x -> dimension(x) == Power, luminosities) "Unit of luminosities is not correct. Got $(unit.(luminosities)), expected power"
    @assert all(x -> dimension(x) == Length, core_radii) "Unit of core radii is not correct. Got $(unit.(core_radii)), expected length"
    @assert all(x -> dimension(x) == Mass, core_masses) "Unit of core masses is not correct. Got $(unit.(core_masses)), expected mass"
    @assert all(x -> dimension(x) == Length, envelope_radii) "Unit of envelope radii is not correct. Got $(unit.(envelope_radii)), expected length"
    @assert all(x -> dimension(x) == Mass, envelope_masses) "Unit of envelope masses is not correct. Got $(unit.(envelope_masses)), expected mass"
    @assert all(x -> dimension(x) == Length, semi_major_axes) "Unit of semi-major axes is not correct. Got $(unit.(semi_major_axes)), expected length"

    @assert all(x -> (unit(x) == u"rad" || unit(x) == u"°"), argument_of_periapses) "Unit of argument of periapses is not correct. Got $(unit.(argument_of_periapses)), expected angle (rad or °)"
    @assert all(x -> (unit(x) == u"rad" || unit(x) == u"°"), inclinations) "Unit of inclinations is not correct. Got $(unit.(inclinations)), expected angle (rad or °)"
    @assert all(x -> (unit(x) == u"rad" || unit(x) == u"°"), longitude_of_ascending_nodes) "Unit of longitude_of_ascending_nodes is not correct. Got $(unit.(longitude_of_ascending_nodes)), expected angle (rad or °)"
    @assert all(x -> (unit(x) == u"rad" || unit(x) == u"°"), true_anomalies) "Unit of true_anomalies is not correct. Got $(unit.(true_anomalies)), expected angle (rad or °)"

end

""" 
The following macro and function are taken from https://github.com/mauro3/UnPack.jl. All credits go 
to the authors of that package.
"""
@inline unpack(x::AbstractDict{Symbol}, ::Val{k}) where {k} = x[k]

macro unpack(args)
    items, dict = args.args
    items = isa(items, Symbol) ? [items] : items.args
    dict_instance = gensym()
    kd = [:( $key = $unpack($dict_instance, Val{$(Expr(:quote, key))}()) ) for key in items]
    kdblock = Expr(:block, kd...)
    expr = quote
        local $dict_instance = $dict 
        $kdblock
        $dict_instance 
    end
    esc(expr)
end

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
function multibodysystem(masses::Vector{<:Unitful.Mass};
                         time::Quantity = 0.0u"s", 
                         verbose::Bool = false, 
                         hierarchy = [length(masses), repeat([1], length(masses)-1)...],
                         system_params...)

    system_params = parse_multibodysystem_args(system_params)
    n_bodies = length(masses)
    n_bins = n_bodies - 1

    @unpack R, S, L, stellar_types, R_core, m_core, R_env, m_env, a, e, ω, i, Ω, ν = system_params
    
    R = ifelse(R isa Number, repeat([R], n_bodies), R)
    L = ifelse(L isa Number, repeat([L], n_bodies), L)
    R_core = ifelse(R_core isa Number, repeat([R_core], n_bodies), R_core)
    m_core = ifelse(m_core isa Number, repeat([m_core], n_bodies), m_core)
    R_env  = ifelse(R_env  isa Number, repeat([R_env ], n_bodies), R_env )
    m_env  = ifelse(m_env  isa Number, repeat([m_env ], n_bodies), m_env )
    stellar_types = ifelse(stellar_types isa Number, repeat([stellar_types], n_bodies), stellar_types)


    S = if S isa Number 
        [[S, zero(S), zero(S)] for i in 1:n_bodies]
    elseif S isa Vector{<:Number}
       [[ss, zero(ss), zero(ss)] for ss in S]
    else
        S
    end

    a = ifelse(a isa Number, repeat([a], n_bins), a)
    e = ifelse(e isa Number, repeat([e], n_bins), e)
    ω = ifelse(ω isa Number, repeat([ω], n_bins), ω)
    i = ifelse(i isa Number, repeat([i], n_bins), i)
    Ω = ifelse(Ω isa Number, repeat([Ω], n_bins), Ω)
    ν = ifelse(ν isa Number, repeat([ν], n_bins), ν)

    check_units(masses, R, S, L, 
                R_core, m_core, 
                R_env, m_env, 
                a, ω, i, Ω, ν)

    hierarchical_multibodysystem_preprocess(masses, R, S, L, 
                                            stellar_types, 
                                            R_core, m_core, 
                                            R_env, m_env, 
                                            a, e, ω, i, Ω, ν, 
                                            hierarchy, time, 
                                            verbose=verbose)
end


function hierarchical_multibodysystem_preprocess(masses::Vector{<:Unitful.Mass}, 
                                                 R::Vector{<:Unitful.Length}, 
                                                 S::Vector{<:Union{<:AbstractVector{<:Quantity}, <:Quantity}}, 
                                                 L::Vector{<:Unitful.Power}, 
                                                 stellar_types::Vector{Int}, 
                                                 R_core::Vector{<:Unitful.Length}, 
                                                 m_core::Vector{<:Unitful.Mass}, 
                                                 R_env::Vector{<:Unitful.Length}, 
                                                 m_env::Vector{<:Unitful.Mass},
                                                 a::Vector{<:Unitful.Length}, e::Vector{<:Real}, ω::Vector{<:Quantity}, 
                                                 i::Vector{<:Quantity}, Ω::Vector{<:Quantity}, ν::Vector{<:Quantity}, 
                                                 hierarchy::Vector{Int}, time::Unitful.Time; 
                                                 verbose::Bool = false)

    n_bins = sum(hierarchy[2:end])
    n_particles = length(masses)
    
    @assert length(masses) == length(S) == length(L) == length(stellar_types) "Must give structural property for each particles."
    @assert n_bins == n_particles - 1 "Number of binary elements must equal N particles - 1."

    elements = OrbitalElements[]
    for idx = 1:n_bins
        push!(elements, OrbitalElements(a[idx], 0.0u"d", e[idx], ω[idx], i[idx], Ω[idx], ν[idx]))
    end

    structures = StellarStructure[]
    for idx ∈ 1:n_particles
        stellar_type = stellar_type_from_index(stellar_types[idx])
        push!(structures, StellarStructure(stellar_type, masses[idx], R[idx], S[idx], L[idx],
                                            R_core[idx], m_core[idx], R_env[idx], m_env[idx]))
    end


    hierarchical_multibodysystem(masses, hierarchy, elements, structures, time, verbose=verbose)
end


function hierarchical_multibodysystem(masses::Vector, hierarchy::Vector, 
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
                    child2 = Particle(ParticleIndex(tot_bodies+1), BinaryIndex(binary_key), child1.key,
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
                    parent_key = iszero(level) ? -1 : child1.parent.i + 1#iszero(level) ? -1 : hierarchy[level] + 1  
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
                    child2 = Particle(ParticleIndex(tot_bodies+1), BinaryIndex(binary_key), child1.key,
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

    HierarchicalMultiple(n_bodies, t, bodies, pairs, binaries, levels, root_bin, hierarchy)
end

function create_binary(positions, velocities, masses, structures, elements, 
                      particle_keys, binary_key, level, parent_key, sibling_key=nothing)
    children = Particle[]

    for i ∈ eachindex(masses)
        sibling_key = findall(particle_keys .!= particle_keys[i])[1]
        child = Particle(ParticleIndex(particle_keys[i]), 
                         BinaryIndex(binary_key), 
                         ParticleIndex(sibling_key), SA[positions[:,i]...], 
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

function get_particles(system::HierarchicalMultiple)
    bins = values(system.binaries)
    particles = Particle[]
    for bin in bins
        append!(particles, get_particles(bin))
    end

    return particles[sortperm([p.key.i for p in particles])]
end

function get_particle(system::HierarchicalMultiple, key::Int)
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

function parse_stellar_structure_args(keyword_arguments)
    keyword_arguments = Dict(keyword_arguments)
    default_param_values = Dict{Symbol, Any}(:R      => 1.0u"Rsun", 
                                             :S      => 0.0u"1/yr", 
                                             :L      => 1.0u"Lsun", 
                                             :R_core => 0.0u"Rsun",
                                             :m_core => 0.0u"Msun",
                                             :R_env  => 0.0u"Rsun",
                                             :m_env  => 0.0u"Msun",
                                             :stellar_types  => 1)

    args_out = copy(default_param_values)

    for (k, v) in keyword_arguments
        param_key = get(multibodysystem_parameter_aliases, k, nothing)
        isnothing(param_key) && throw(ArgumentError("Parameter $k is not valid. See `keys(Syzygy.multibodysystem_parameter_aliases)` for all possible parameters."))
        args_out[param_key] = v
    end

    return args_out
end

"""
    multibodysystem(masses, positions, velocities; kwargs...)

Create a NonHierarchicalSystem from masses, positions, and velocities. 

# Examples
```jldoctest
julia> masses = rand(3)u"kg"
julia> positions = [rand(3)u"m" for i = 1:3]
julia> velocities = [rand(3)u"m/s" for i = 1:3]
julia> multibodysystem(masses, positions, velocities);
julia> multibodysystem(masses, positions, velocities, radii=ones(3)u"m") # you can also supply any of the stellar structure arguments
```

"""
function multibodysystem(masses, positions, velocities;
                         time::Quantity = 0.0u"s", 
                         verbose::Bool = false,
                         stellar_structure_params...)

    stellar_structure_params = parse_stellar_structure_args(stellar_structure_params)
    @unpack R, S, L, R_core, m_core, R_env, m_env, stellar_types = stellar_structure_params

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
        [[S, zero(S), zero(S)] for i in 1:n_bodies]
    elseif S isa Vector{<:Number}
        [[ss, zero(ss), zero(ss)] for ss in S]
    else
        S
    end

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

    particles = Dict{Int, Particle}(i => p for (i, p) in enumerate(particles))

    NonHierarchicalSystem(n_bodies, time, particles, pairs)
end


"""
    multibodysystem(system::HierarchicalMultiple; new_params...)

Remake a given system with new parameters. 
"""
function multibodysystem(system::HierarchicalMultiple; new_params...)

    possible_kwargs = [:R, :S, :L, :R_core, :m_core, :R_env, :m_env, 
                       :stellar_types, :a, :e, :ω, :i, :Ω, :ν, 
                       :hierarchy, :time, :verbose]

    new_kwargs = Dict(new_params)
    unchanched_kwargs = [kw for kw in possible_kwargs if (!in(kw, keys(new_kwargs)))]


    particle_keys = collect(keys(system.particles)) |> sort
    binary_keys = collect(keys(system.binaries)) |> sort
    all_args = Dict{Symbol, Any}()
    for arg in unchanched_kwargs
        # @show arg
        if any(arg .=== (:R, :S, :L, :R_core, :m_core, :R_env, :m_env))
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
                           R_core=R_core, m_core=m_core, R_env=R_env, m_env=m_env, 
                           stellar_types=stellar_types, time=time)

end

function initialize_from_file(filepath)

    if !(endswith(filepath, "jld2"))
        @error "Only JLD2 files are currently supported."
    end

    data = JLD2.load(filepath)
    datakeys = collect(keys(data))

    if !(masses in datakeys)
        @error "Datafile must contain masses"
    end

    data = if eltype(datakeys) == String
        data_new = Dict{Symbol, Any}()
        for (k, v) in data
            data_new[Symbol(k)] = v
        end
        data_new
    else
        data
    end

    masses = try
        pop!(data, :masses)
    catch e
        pop!(data, :m)
    end

    if "positions" in datakeys
        positions, velocities = pop!(data, "positions"), pop!(data, "velocities")
        return multibodysystem(masses, positions, velocities; data...)
    else
        return multibodysystem(masses; data...)
    end
end

function initialize_from_dict(d::Dict)

    d = copy(d)
    masses = if haskey(d, :masses)
        pop!(d, :masses)
    else
        pop!(d, :m)
    end

    if "positions" in keys(d)
        positions, velocities = pop!(d, "positions"), pop!(d, "velocities")
        return multibodysystem(masses, positions, velocities; d...)
    else
        return multibodysystem(masses; d...)
    end
end