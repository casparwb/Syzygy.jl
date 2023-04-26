""" 
Functions for calculating cartesian positions and velocities 
from orbital elements of a binary.
"""

# using NonlinearSolve

# """
#     keplers_equation(e, M)

# Solve Kepler's equation for the eccentric anomaly E given
# eccentricity e and mean anomaly M.
# """
# function eccentric_anomaly(e, M)

#     function f(dE, E, p)
#         dE .= E[1] - e*sin(E[1]) - M
#     end

#     E₀ = [π/2]

#     probN = NonlinearProblem(f, E₀)
#     sol = solve(probN, NewtonRaphson(), reltol=1e-10)

#     return first(sol.u)
# end


"""
    keplerian_to_cartesian(elements::OrbitalElements, mass::Quantity; G=𝒢)

Convert Keplerian orbital elements to cartesian position and velocity following the formalism of 
[](https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf). 
Motion is by default assumed to be in the xy-plane. 
"""
function keplerian_to_cartesian(elements::OrbitalElements, mass; G=𝒢)

    # semi-major axis
    # time
    # eccentricity
    # argument of periastron
    # inclination
    # longitude of ascending node
    # true anomaly
    # els = (:a, :P, :e, :ω, :i, :Ω, :ν)
    a, e, ω, i, Ω, ν = let e = elements
        e.a, e.e, e.ω, e.i, e.Ω, e.ν
    end

    # P = orbital_period(a, mass, G)
    μ = G*mass

    sqrt_1_min_e² = sqrt(1.0 - e^2)

    E = atan((sqrt_1_min_e²*sin(ν))/(e + cos(ν)))
    # M = atan(sqrt_1_min_e²*sin(ν), -e - cos(ν)) + π - e*(sqrt_1_min_e²*sin(ν))/(1 + e*cos(ν))
    # E = eccentric_anomaly(e, M)

    # ecosE = 1 + e*cos(E)
    # Edot = n*ecosE^2/sqrt_1_min_e²^3

    # r = a*(1 - e^2)/ecosE
    r = a*(1 - e^2)/(1 + e*cos(ν))
    h = √(μ*a*(1 - e^2))
    
    # rdot = a*n/sqrt_1_min_e²*e*sin(E)

    o = zeros(typeof(1.0u"m"  ), 3)
    o_dot = zeros(typeof(1.0u"m/s"  ), 3)
    x = zeros(typeof(1.0u"m"  ), 3)
    v = zeros(typeof(1.0u"m/s"), 3)

    cosΩ, sinΩ = cos(Ω), sin(Ω)
    # cosΩ, sinΩ = if Ω ≈ π
    #                  cos(Ω), sin(Ω)
    #              else
    #                  -1.0, 0.0
    #              end

    cosω = cos(ω)
    sinω = sin(ω)

    cosω_plus_ν = cos(ω + ν)
    sinω_plus_ν = sin(ω + ν)
    # @show E
    cosi, sini = cos(i), sin(i)
    # cosi, sini = if i ≈ π/2
    #                 #  cos(i), sin(i)
    #                 0.0, 1.0
    #              else
    #                  0.0, 1.0
    #              end


    # Positions 
    o[1] = r*cos(ν)
    o[2] = r*sin(ν)
    o[3] = 0.0u"m"

    f = √(μ*a)/r
    o_dot[1] = -f*sin(E)
    o_dot[2] = f*sqrt_1_min_e²*cos(E)
    o_dot[3] = 0.0u"m/s"

    x[1] = o[1]*(cosω*cosΩ - sinω*cosi*sinΩ) - o[2]*(sinω*cosΩ + cosω*cosi*sinΩ)
    x[2] = o[1]*(cosω*sinΩ + sinω*cosi*cosΩ) + o[2]*(cosω*cosi*cosΩ - sinω*sinΩ)
    x[3] = o[1]*sinω*sini + o[2]*cosω*sini

    # x[1] = r*(cosΩ*cosω_plus_ν - sinΩ*sinω_plus_ν*cosi)
    # x[2] = r*(sinΩ*cosω_plus_ν + cosΩ*sinω_plus_ν*cosi)
    # x[3] = r*sinω_plus_ν*sini
    # # rdot_r = rdot/r

    # # Velocities
    # # rEdot = r*Edot
    # p = a*(1 - e^2)
    # rp = r*p
    # he = h*e

    v[1] = o_dot[1]*(cosω*cosΩ - sinω*cosi*sinΩ) - o_dot[2]*(sinω*cosΩ + cosω*cosi*sinΩ)
    v[2] = o_dot[1]*(cosω*sinΩ + sinω*cosi*cosΩ) + o_dot[2]*(cosω*cosi*cosΩ - sinω*sinΩ)
    v[3] = o_dot[1]*sinω*sini + o_dot[2]*cosω*sini

    # v[1] = x[1]*he/rp*sin(ν) - h/r*(cosΩ*sinω_plus_ν + sinΩ*cosω_plus_ν*cosi)
    # v[2] = x[2]*he/rp*sin(ν) - h/r*(sinΩ*sinω_plus_ν + cosΩ*cosω_plus_ν*cosi)
    # v[3] = x[3]*he/rp*sin(ν) - h/r*sini*cosω_plus_ν
    return x,v
end


"""
    keplers_problem(hierarchy::AbstractMatrix, masses::AbstractVector, elements::AbstractVector{OrbitalElements})

Solve Kepler's problem for a system with a given hierarchy. Returns two 
(n×3) matrices of positions and velocities for each of the 'n' bodies in the system.
"""
function keplers_problem(hierarchy, masses, elements)

    n = length(masses)
    positions = zeros(typeof(1.0u"m"), n, 3)
    velocities = zeros(typeof(1.0u"m/s"), n, 3)

    # positions = zeros(typeof(1.0u"m"), 3, n)
    # velocities = zeros(typeof(1.0u"m/s"), 3, n)

    i = 1
    bin = -1

    # elements = reverse(elements)
    while i < n
        idx = hierarchy[i,:] .|> abs .|> Bool
        μ = sum(masses[idx])
        # @show masses[i], μ

        if iszero(first(idx))
            bin =+ 1
        end
        # @show elements
        binary_elements = elements[i]
        # @show 
        pos, vel = keplerian_to_cartesian(binary_elements, μ)
        positions[i,:] .= pos
        velocities[i,:] .= vel

        # positions[:,i] .= pos
        # velocities[:,i] .= vel

        if bin > 0
            bin -= 2
        elseif bin < -1
            bin = 0
        end

        i += 1
    end

    return positions, velocities
end


"""
    mass_ratio_matrix(hierarchy::AbstractMatrix, masses::AbstractVector)

Compute the mass ratio matrix A as defined in [Hamers and Portegies Zwart 2016](https://doi.org/10.1093/mnras/stw784)
given a hierarchy matrix and mass vector.
"""
function mass_ratio_matrix(hierarchy::AbstractMatrix, masses::AbstractVector)

    A = zeros(Float64, size(hierarchy))
    n = size(A, 1)

    for j = 1:n
        for i = 1:n
            m = zero(typeof(1.0u"kg"))
            for l = 1:n
                m += masses[l]*ifelse(hierarchy[i, j] == hierarchy[i, l], 1.0, 0.0)
            end
            A[i, j] = (hierarchy[i, j]*masses[j])/m

        end
    end

    return A
end