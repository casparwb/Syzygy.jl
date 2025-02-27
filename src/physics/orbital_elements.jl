
@doc """

    semi_major_axis(d, v², M)

Return the semi-major axis of a binary system with total mass 'M'.
Here, 'd' denotes the relative distance, while 'v²' is the relative
velocity magnitude squared.

```math
    a = \\frac{GM d}{2GM - d v²} 
```
"""
function semi_major_axis(d, v², M, G)
    GM = G*M
    GM*d/(2GM - d*v²)
end

semi_major_axis(d, v², M) = semi_major_axis(d, v², M, GRAVCONST)

"""
    orbital_period(a, M, G)

Get the orbital period of binary with semi-major axis `a` and total mass `M`.
"""
function orbital_period(a, M, G)
    a < zero(a) && @warn "Given semi-major axis is negative: " a
    2π*√(abs(a)^3/(G*M))
end

orbital_period(a, M::Unitful.Mass) = orbital_period(a, M, GRAVCONST)

"""
    eccentricity_vector(r, v, d, μ)

Eccentricity vector of binary orbit with relative position vector `r`, separation `d`, relative
velocity `v`, and standard gravitational parameter `μ`, where `μ ≡ G(m₁ + m₂)`.

```math
    \\vec{e} = \\frac{\\vec{v} \\times (\\vec{r} \\times \\vec{v})}{\\mu} - \\frac{\\vec{r}}{d} 
```
"""
function eccentricity_vector(r, v, d, μ)
    (v × (r × v))/μ - r/d
end

"""
    eccentricity_vector(r, v, d, m::AbstractVector, G=GRAVCONST)
"""
function eccentricity_vector(r, v, d, masses::AbstractVector, G=GRAVCONST)
    μ = G*(masses[1] + masses[2])
    eccentricity_vector(r, v, d, μ)
end

# function eccentricity_vector(r, v, d, M::Number, G=GRAVCONST)
#     μ = G*M
#     eccentricity_vector(r, v, d, μ)
# end

function eccentricity(r, v, d, m::AbstractVector, G=GRAVCONST)
    μ = G*(m[1] + m[2])
    return norm(eccentricity_vector(r, v, d, μ))
end

function eccentricity(r, v, d, M::Number, G=GRAVCONST)
    μ = G*M
    return norm(eccentricity_vector(r, v, d, μ))
end

"""
    eccentricity(r, v, a, M)

Eccentricity of binary orbit with total mass `M`, semi-major axis `a`,
relative positions `r`, and relative velocities `v`.

``e = \\sqrt{1 - \\frac{|r×v|²}{GMa}}`` 
"""
function eccentricity_old(r, v, a, M, G)
    e² = 1 - sum(abs2, r × v)/(G*M*a)
    e² < 0 && return 1e-4
    return sqrt(e²)
end

eccentricity_old(r, v, a, M) = eccentricity(r, v, a, M, GRAVCONST)

function angular_momentum(r, v)
    r×v
end

inclination(h) = acos(h[end]/norm(h))
    
"""
    mutual_inclination(h₁, h₂)

Return the mutual inclination angle between two binaries with orbital angular momenta h₁ and h₂.
"""
function mutual_inclination(h₁, h₂)
    return acos(dot(h₁, h₂)/(norm(h₁)*norm(h₂)))
end


"""
    longitude_of_ascending_node(h)

Return the longitude of ascending node of a body in an orbit with angular momentum `h`.
"""
function longitude_of_ascending_node(h)
    n = SA[-h[2], h[1], zero(h[1])]
    Ω = acos(n[1]/norm(n))
    isnan(Ω) && return 0.0u"rad"
    n[2] >= n[3] && return (Ω)u"rad"
    return (2π - Ω)u"rad"
end

"""
    argument_of_periapsis(r, v, h, m, G)

Return the argument of periapsis of a body in an orbit with relative position `r`,
velocity `r`, angular momentum `h` and mass `m`.
"""
function argument_of_periapsis(r, v, h, m, G)
    μ = G*m
    n = SA[-h[2], h[1], zero(h[1])]
    # n = 
    e = (v × h)/μ .- r/norm(r)

    ω = acos(dot(e, n)/(norm(n)*norm(e)))
    if isnan(ω) 
        ω = atan(e[2], e[1])
        return ifelse(h[3] < zero(h[3]), 2π - ω, ω)u"rad"
    end

    return ifelse(e[3] < zero(e[3]), 2π - ω, ω)u"rad"
end

argument_of_periapsis(r, v, h, M) = argument_of_periapsis(r, v, h, M, GRAVCONST)

"""
    true_anomaly(r, v, h, M, G)

Return the true anomaly of a body in an orbit with relative position `r`,
velocity `r`, angular momentum `h` and mass `m`.
"""
function true_anomaly(r, v, h, m, G)
    μ = G*m
    r_norm = norm(r)
    e = (v × h)/μ .- r/r_norm
    
    vᵣ = dot(r/r_norm, v)
    arg = dot(e, r)/(norm(e)*r_norm)
    ν = acos(arg)

    # ν = acos(min(1.0, abs(arg)))
    # rdotv = dot(r, v)
    # println(ν, " ", rdotv)
    # ν = ifelse(rdotv < zero(rdotv), 2π - ν, ν)
    # println(2π - ν, " ", vᵣ)
    vᵣ < zero(vᵣ) && return 2π - ν
    return ν
    # return ifelse(vᵣ >= zero(vᵣ), ν, 2π - ν)
end

true_anomaly(r, v, h, M) = true_anomaly(r, v, h, M, GRAVCONST)


"""
    elements_from_cartesian(positions, velocities, masses)

Get orbital elements for given state vectors by checking which bodies are bound to each other
and calculating their keplerian elements. WORK IN PROGRESS 
"""
function elements_from_cartesian(positions, velocities, masses)

    n_bodies = length(masses)
    bodies = 1:n_bodies
    pairs = SVector{2, Int}[]
    for i in bodies
        for j in (i+1):n_bodies
            if i != j
                push!(pairs, SA[i, j])
            end
        end
    end

    bound_pairs = SVector{2, Int}[]
    for pair in pairs
        if is_binary(positions[pair], velocities[pair], masses[pair])
            push!(bound_pairs, pair)
        end
    end

    unbound_bodies = Int[]
    for body in bodies
        in_binary = false
        for pair in bound_pairs
            if body in pair
                in_binary = true
                break
            end
        end

        if !in_binary
            push!(unbound_bodies, body)
        end
    end
            

    @info "Bound pairs: " bound_pairs
    @info "Unbound bodies: " unbound_bodies

    # return pairs
end

"""
    binary_elements(positions, velocities, masses)

Calculate binary properties of bodies with given positions and velocities.
Assumes the bodies are gravitationally bound.
"""
function binary_elements(positions, velocities, masses)

    r1 = positions[1]
    r2 = positions[2]

    v1 = velocities[1]
    v2 = velocities[2]
    
    M1, M2 = masses
    
    r_rel = r2 - r1
    v_rel = v2 - v1
    h = angular_momentum(r_rel, v_rel)

    d = norm(r_rel)
    v = norm(v_rel)
    M = M1 + M2

    a = semi_major_axis(d, v^2, M, GRAVCONST) |> u"Rsun"
    e = eccentricity(r_rel, v_rel, d, M, GRAVCONST) 

    P = orbital_period(a, M, GRAVCONST) |> u"d"
    h = angular_momentum(r_rel, v_rel)
    ω = argument_of_periapsis(r_rel, v_rel, h, M, GRAVCONST) 
    i = inclination(h)
    Ω = longitude_of_ascending_node(h)
    ν = true_anomaly(r_rel, v_rel, h, M, GRAVCONST)
    # ν = true_anomaly(r1, v1, h1, M1, GRAVCONST)

    OrbitalElements(a, P, e, ω, i, Ω, ν)
end

"""
    is_binary(positions, velocities, masses)

Check if two bodies with the given positions, velocities, and masses
are gravitationally bound and form a binary.
"""
function is_binary(positions::AbstractVector, velocities::AbstractVector, masses::AbstractVector)

    r1 = positions[1]
    r2 = positions[2]

    v1 = velocities[1]
    v2 = velocities[2]

    M1, M2 = masses

    r_rel = r2 - r1

    v₁² = norm(v1)^2
    v₂² = norm(v2)^2

    K = 0.5*(M1*v₁² + M2*v₂²)
    U = -GRAVCONST*M1*M2/norm(r_rel)

    (K + U) < zero(U)
end

function is_binary(positions::AbstractMatrix, velocities::AbstractMatrix, masses)
    return is_binary(eachcol(positions), eachcol(velocities), masses)
end