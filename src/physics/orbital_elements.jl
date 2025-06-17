
@doc """

    semi_major_axis(d, v², M)

Return the semi-major axis of a binary system with total mass 'M'.
Here, 'd' denotes the relative separation, while 'v²' is the relative
velocity magnitude squared.

```math
    a = \\frac{GM d}{2GM - d v²} 
```
"""
function semi_major_axis(d, v², M, G)
    GM = G*M
    GM*d/(2GM - d*v²)
end

semi_major_axis(d::Quantity, v², M) = semi_major_axis(d, v², M, GRAVCONST)
semi_major_axis(d::Real, v², M) = semi_major_axis(d, v², M, UNITLESS_G)

"""
    orbital_period(a, M)

Get the orbital period of binary with semi-major axis `a` and total mass `M`.
"""
function orbital_period(a, M, G)
    if a < zero(a) 
        @warn "Given semi-major axis is negative: " a
        return Inf
    end
    2π*√(abs(a)^3/(G*M))
end

orbital_period(a, M::Unitful.Mass) = orbital_period(a, M, GRAVCONST)
orbital_period(a, M::Real) = orbital_period(a, M, UNITLESS_G)

"""
    eccentricity_vector(r, v, d, M)

Eccentricity vector of binary orbit with relative position vector `r`, separation `d`, relative
velocity `v`, and total mass `M`.

```math
    \\vec{e} = \\frac{\\vec{v} \\times (\\vec{r} \\times \\vec{v})}{\\mu} - \\frac{\\vec{r}}{d} 
```
"""
function eccentricity_vector(r, v, d, M, G)
    μ = G*M
    (v × (r × v))/μ - r/d
end

"""
    eccentricity_vector(r, v, d, m::AbstractVector, G=GRAVCONST)
"""
# function eccentricity_vector(r, v, d, masses::AbstractVector, G)
#     μ = G*(masses[1] + masses[2])
#     eccentricity_vector(r, v, d, μ)
# end

eccentricity_vector(r, v, d::Unitful.Length, M) = eccentricity_vector(r, v, d, M, GRAVCONST)
eccentricity_vector(r, v, d::Real, M) = eccentricity_vector(r, v, d, M, UNITLESS_G)


# function eccentricity(r, v, d, m::AbstractVector, G)
#     μ = G*(m[1] + m[2])
#     return norm(eccentricity_vector(r, v, d, μ))
# end
"""
    eccentricity(r, v, d, M)

Eccentricity of binary orbit with relative position vector `r`, relative
velocity `v`,  separation `d`, and total mass `M`, where `M ≡ m₁ + m₂`.
The eccentricity is calculated as the norm of the eccentricity vector

```math
    e = |\\vec{e}| = \\frac{\\vec{v} \\times (\\vec{r} \\times \\vec{v})}{\\mu} - \\frac{\\vec{r}}{d} 
```
"""
function eccentricity(r, v, d, M, G)
    return norm(eccentricity_vector(r, v, d, M, G))
end

eccentricity(r, v, d::Unitful.Length, M) = eccentricity(r, v, d, M, GRAVCONST)
eccentricity(r, v, d::Real, M) = eccentricity(r, v, d, M, UNITLESS_G)



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
    isnan(Ω) && return 0.0
    n[2] >= n[3] && return (Ω)
    return (2π - Ω)
end

"""
    argument_of_periapsis(r, v, h, m, G)

Return the argument of periapsis of a body in an orbit with relative position `r`,
velocity `v`, angular momentum `h` and mass `m`.
"""
function argument_of_periapsis(r, v, h, m, G)
    μ = G*m
    n = SA[-h[2], h[1], zero(h[1])]
    e = (v × h)/μ .- r/norm(r)

    ω = acos(dot(e, n)/(norm(n)*norm(e)))
    if isnan(ω) 
        ω = atan(e[2], e[1])
        return ifelse(h[3] < zero(h[3]), 2π - ω, ω)
    end

    return ifelse(e[3] < zero(e[3]), 2π - ω, ω)
end

argument_of_periapsis(r, v, h, M::Quantity) = argument_of_periapsis(r, v, h, M, GRAVCONST)*u"rad"
argument_of_periapsis(r, v, h, M::Real) = argument_of_periapsis(r, v, h, M, UNITLESS_G)

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

true_anomaly(r, v, h, M::Quantity) = true_anomaly(r, v, h, M, GRAVCONST)
true_anomaly(r, v, h, M::Real) = true_anomaly(r, v, h, M, UNITLESS_G)


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
Assumes the bodies are gravitationally bound. Returns an instance of [`OrbitalElements`](@ref).
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

    a = semi_major_axis(d, v^2, M)
    e = eccentricity(r_rel, v_rel, d, M) 

    P = orbital_period(a, M)
    h = angular_momentum(r_rel, v_rel)
    ω = argument_of_periapsis(r_rel, v_rel, h, M) 
    i = inclination(h)
    Ω = longitude_of_ascending_node(h)
    ν = true_anomaly(r_rel, v_rel, h, M)

    OrbitalElements(a, P, e, ω, i, Ω, ν)
end

"""
    is_binary(r1, r2, v1, v2, m1, m1)

Check if two bodies with the given positions (r1, r2), velocities (v1, v2), and masses (m1, m2)
are gravitationally bound and form a binary.
"""
function is_binary(r1, r2, v1, v2, m1, m2, G)
    r_rel = r2 - r1

    v₁² = norm(v1)^2
    v₂² = norm(v2)^2

    K = 0.5*(m1*v₁² + m2*v₂²)
    U = -G*m1*m2/norm(r_rel)

    (K + U) < zero(U)
end

is_binary(r1, r2, v1, v2, m1, m2::Quantity) = is_binary(r1, r2, v1, v2, m1, m2, GRAVCONST)
is_binary(r1, r2, v1, v2, m1, m2::Real) = is_binary(r1, r2, v1, v2, m1, m2, UNITLESS_G)


function is_binary(positions::AbstractVector, velocities::AbstractVector, masses::AbstractVector)

    r1 = positions[1]
    r2 = positions[2]

    v1 = velocities[1]
    v2 = velocities[2]

    m1, m2 = masses

    is_binary(r1, r2, v1, v2, m1, m2)

end

function is_binary(positions::AbstractMatrix, velocities::AbstractMatrix, masses)
    return is_binary(eachcol(positions), eachcol(velocities), masses)
end