
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
    G*M*d/(2G*M - d*v²)
end

semi_major_axis(d, v², M) = semi_major_axis(d, v², M, 𝒢)

# function semi_major_axis(P, M, G)
#     return cbrt(P^2/(4π^2)*(G*M))
# end

function orbital_period(a, M, G)
    2π*√(abs(a)^3/(G*M))
end

orbital_period(a, M) = orbital_period(a, M, 𝒢)

function eccentricity_vector(r, d, v, α::Number, G)
    # α = G*(m[1] + m[2])
    e = (v × (r × v))/α - r/d
end


function eccentricity_vector(r, d, v, m, G)
    α = G*(m[1] + m[2])
    eccentricity_vector(r, d, v, α, G)
end


"""
    eccentricity(r, v, a, M)

Eccentricity of binary orbit with total mass `M`, semi-major axis `a`,
relative positions `r`, and relative velocities `v`.

``e = \\sqrt{1 - \\frac{|r×v|²}{GMa}}`` 
"""
function eccentricity(r, v, a, M, G)
    e² = 1 - sum(abs2, r × v)/(G*M*a)
    e² < 0 && return 1e-4
    return sqrt(e²)
end

eccentricity(r, v, a, M) = eccentricity(r, v, a, M, 𝒢)

angular_momentum(r, v) = r×v

inclination(h) = acos(h[end]/norm(h))
    
function mutual_inclination(h₁, h₂)
    return acos(dot(h₁, h₂)/(norm(h₁)*norm(h₂)))
end

function longitude_of_ascending_node(h)
    n = SA[-h[2], h[1], zero(h[1])]
    Ω = acos(n[1]/norm(n))
    isnan(Ω) && return 0.0
    n[2] >= n[3] && return Ω
    return 2π - Ω
end

function argument_of_periapsis(r, v, h, M, G)
    μ = G*M
    n = SA[-h[2], h[1], zero(h[1])]
    e = (v × h)/μ .- r/norm(r)
    ω = acos(min(1.0, dot(n, e)/(norm(n)*norm(e))))
    isnan(ω) && return atan(e[2], e[1])
    e[3] < zero(e[3]) && return 2π - ω
    return ω
end

argument_of_periapsis(r, v, h, M) = argument_of_periapsis(r, v, h, M, 𝒢)

function true_anomaly(r, v, h, M, G)
    μ = G*M
    e = (v × h)/μ .- r/norm(r)
    arg = dot(e, r)/(norm(e)*norm(r))
    # if dot(r, v) < 0
    ν = acos(min(1.0, abs(arg)))
    rdotv = dot(r, v)
    ν = ifelse(rdotv < zero(rdotv), 2π - ν, ν)
end

true_anomaly(r, v, h, M) = true_anomaly(r, v, h, M, 𝒢)
