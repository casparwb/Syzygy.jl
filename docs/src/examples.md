
## Quadruples

For a quadruple stellar system, there are two possible stable configurations: 3+1 or 2+2, where the former defines a hierarchical triple with a more distant quaternary companion, and the latter is a system with two binaries orbiting a common centre of mass. `Syzygy` let's you set up both kinds of systems using the `hierarchy` keyword argument for [`multibodysystem`](@ref). First, we'll set up the masses of the stars.

```@example quad
using Syzygy, Plots # hide
masses = ones(4)u"Msun"
```

### Hierarchical quadruple

A hierarchical configuration is the default way of instantiating a multiple system, so we don't need to set any value for the `hierarchy` parameter. Let's define the orbital properties. For this example, we'll define two binaries with semi-major axes $1$ and $0.5$ AU, with a semi-major axis of $10$ AU separating the centre of mass of the two binaries.

```@example quad

a_bin_1 = 0.1u"AU"
a_bin_2 = 1.0u"AU"
a_bin_3 = 6.0u"AU"

e_bin_1 = 0.5
e_bin_2 = 0.0
e_bin_3 = 0.4

quad = multibodysystem(masses, a=[a_bin_1, a_bin_2, a_bin_3], 
                               e=[e_bin_1, e_bin_2, e_bin_3])
res = simulate(quad, t_sim=1)
sol = to_solution(res)
orbitplot(sol, dims=[1, 2])
```



### 2+2 quadruple

For this example, we'll define two binaries with semi-major axes $1$ and $0.5$ AU, with a semi-major axis of $10$ AU separating the centre of mass of the two binaries. So set up this binary, we have to set the `hierachy` argument such that we get the correct configuration.

```@example quad

a_bin_1 = 1.0u"AU"
a_bin_2 = 0.5u"AU"
a_bin1_bin2 = 10.0u"AU"

e_bin_1 = 0.5
e_bin_2 = 0.0
e_bin1_bin2 = 0.4

hierarchy = [4, 2, 1]

quad = multibodysystem(masses, a=[a_bin_1, a_bin_2, a_bin1_bin2], 
                               e=[e_bin_1, e_bin_2, e_bin1_bin2], 
                               hierarchy=hierarchy)
res = simulate(quad, t_sim=1)
sol = to_solution(res)
orbitplot(sol, dims=[1, 2], label=false)
```

## The solar system

In this example we show how the code can be used to simulate plantary systems, using our own solar system as the test case. First, we'll define the masses and orbital properties of all the planets.

```@example solar_system
using Syzygy, Plots # hide
mass_sun = 1u"Msun"
# mass_planets = [0.330, 4.87, 5.97, 0.642, 1898, 568, 86.8, 102]*1e24u"kg"
mass_planets = [0.0553, 0.8154, 0.9996, 0.1075, 
                317.8075, 95.1078, 14.5341, 17.0792]u"Mearth"
masses = [mass_sun, mass_planets...]

semi_major_axes = [0.3871, 0.7233, 1.0, 1.5273, 
                   5.2028, 9.5388, 19.1914, 30.0611]u"AU"

eccentricites = [0.206, 0.007, 0.017, 0.093, 
                   0.048, 0.056, 0.046, 0.010]

inclinations = [7.0, 3.39, 0.0, 1.85, 
                   1.31, 2.49, 0.77, 1.77]u"°"

planet_names = ["Sun" "Mercury" "Venus" "Earth" "Mars" "Jupiter" "Saturn" "Uranus" "Neptune"]
nothing # hide
```

Then, we'll create an instance of `multibodysystem`, and simulate the system for one neptune orbit.

```@example solar_system
solar_system = multibodysystem(masses, sma=semi_major_axes, e=eccentricites, i=inclinations)

orbital_periods = [bin.elements.P for (_, bin) in sort(solar_system.binaries)]
sol = to_solution(simulate(solar_system, t_sim=orbital_periods[end])) ;

orbitplot(sol, dims=[1, 2], label=planet_names)
```

Let's just simulate the inner solar system

```@example solar_system
inner_solar_system = multibodysystem(masses[1:4], sma=semi_major_axes[1:4], 
                                                  e=eccentricites[1:4], 
                                                  i=inclinations[1:4])

orbital_periods = [bin.elements.P for (_, bin) in sort(inner_solar_system.binaries)]
sol = to_solution(simulate(inner_solar_system, t_sim=orbital_periods[end])); # 1 orbital period of mars

orbitplot(sol, dims=[1, 2], label=reshape(planet_names[1:4], 1, :))
```

## Simulating a precessing binary black-hole system

For this example, we're going to simulate a binary black-hole system undergoing orbital precession due to general relativistic effects. 

```@example bbh
using Syzygy, Plots # hide
masses = [27, 25]u"Msun"
radii = Syzygy.gravitational_radius.(masses)
stellar_types = [14, 14]

separation = 3.0u"Rsun"
ecc = 0.9
bbh = multibodysystem(masses, semi_major_axis=separation, e=ecc, radii=radii, stellar_types=stellar_types)
nothing # hide
```

```@example bbh
# set up callbacks and potentials

callbacks = [CollisionCB(1, 50)] # we'll reduce the separation for checking for collision
potentials = [PureGravitationalPotential(), PNPotential()]

# the default algorithm (DPRKN8) does not support velocity-dependent 
# acceleration functions, so we have to use a different one.
alg = Syzygy.FineRKN5() 
P = bbh.binaries[1].elements.P
T = 100P # simulate for 100 orbital periods
res = simulate(bbh, t_sim=T,  
                    alg=alg, callbacks=callbacks,
                    potential=potentials)
sol = to_solution(res);

orbitplot(sol, dims=[1, 2])
```

## High precision simulations

In this example we'll simulate a chaotic, unstable triple system using higher-precision numerics. For this, we'll use the `Double64` number type from the [DoubleFloats.jl](https://github.com/JuliaMath/DoubleFloats.jl) package. As we can see, this type has twice as many bits in the significand compared to the standard `Float64`.

```@example highprec
using Syzygy

println("Float64 precision: ", precision(Float64), " \nDouble64 precision: ", precision(Syzygy.Double64))
```

To specify this datatype for the simulation, pass the keyword argument `precision=:Double64` to [`simulate`](@ref).

```@example highprec

masses = [5, 1, 10]u"Msun"
a_in = 1.0u"AU"
a_out = 3.0u"AU"
e_in = 0.7
e_out = 0.5

triple = multibodysystem(masses, a=[a_in, a_out], e=[e_in, e_out])

res = simulate(triple, t_sim=1, precision=:Double64, abstol=1e-20, reltol=1e-20, save_everystep=false)
E_init = Syzygy.total_energy(res, res.solution.t[1])
E_final = Syzygy.total_energy(res, res.solution.t[end])

E_init/E_final - 1
```

## Adding single stellar evolution

In this example we'll show how you can take advantage of the extensibility of `Syzygy` to include stellar evolution, allowing the properties of the stellar bodies change during the simulation. There are a few things we need to set up in order to get this working. The checklist is:

    1. A struct whose fields are mutable vectors of stellar parameters.
    2. A function for instantiating the parameter struct.
    3. A function for setting up the callback where the stellar parameters are evolved.

We'll first do (1) and (2). The struct needs to have  _at least_ the `masses`, `radii`, `stellar_types`, and `stellar_type_numbers` fields. 

```julia

struct StellarEvolutionParams <: Syzygy.SimulationParams
    radii::MVector{3, Float64} 
    masses::MVector{3, Float64} 
    stellar_types::MVector{3, Syzygy.StellarType} 
    stellar_type_numbers::MVector{3, Int}
end

```

We now need to add a method to [`Syzygy.setup_params`](@ref) using the type signature of our new struct to make sure that the solver will use these params.

```julia

function Syzygy.setup_params(::Type{<:StellarEvolutionParams}, system, datatype=Float64)
    particles = system.particles

    masses  = datatype[]
    radii  = datatype[]
    stellar_type_nums = Int[]
    stellar_types = Syzygy.StellarType[]

    particle_keys = keys(particles) |> collect |> sort

    for i in particle_keys
        p = particles[i]
        
        mass         = ustrip(unit_mass, p.structure.m) 
        radius       = ustrip(unit_radius, p.structure.R) 
        stellar_type = p.structure.stellar_type 

        push!(masses, mass)
        push!(radii, radius)
        push!(stellar_types, stellar_type)
        push!(stellar_type_nums, stellar_type.number)
    end

    radii = MVector(radii...)
    masses = MVector(masses...)
    stellar_types = MVector(stellar_types...)
    stellar_type_nums = MVector(stellar_type_nums...)

    ode_params = StellarEvolutionParams(radii, masses, 
                                        stellar_types,
                                        stellar_type_nums)

    return ode_params
end 

```

Next step is to set up the callback that actually mutates the stellar structure parameters. For this example, we'll assume we have pre-evolved these stars and saved the mass, radius, stellar type evolution to a file. We can then read this file, construct and interpolator, and use this in the callback to set the stellar properties at each time step. Similar to the previous steps, we first need to define a struct, followed by overloading the function `Syzygy.get_callback` using the type signature of this struct. For this example, we'll use [DataInterpolations.jl](https://docs.sciml.ai/DataInterpolations/stable/) to do the interpolation.

```julia
using DataInterpolations, DiffEqCallbacks

struct StellarEvolutionCallback :> Syzygy.AbstractSyzygyCallback

function Syzygy.get_callback(cb::StellarEvolutionCallback, system, retcode, args)

    pre_calculated_evolution = JLD2.load("evolution.jld2")
    times = pre_calculated_evolution["time"]
    masses = pre_calculated_evolution["mass"] # a 3xN matrix
    radii = pre_calculated_evolution["radius"]
    stellar_type_nums = pre_calculated_evolution["stellar_type"]

    itp_mass = LinearInterpolation(masses, times)
    itp_radius = LinearInterpolation(radii, times)
    itp_stellar_type = ConstantInterpolation(stellar_type_nums, times)

    function do_stellar_evolution!(u, t, integrator)
        integrator.p.masses .= SVector{3, Float64}(itp_mass(t))
        integrator.p.radii .= SVector{3, Float64}(itp_radius(t))
        integrator.p.stellar_type_numbers .= SVector{3, Float64}(itp_stellar_type(t))
        integrator.p.stellar_types .= Syzygy.stellar_type_from_index.(integrator.p.stellar_type_numbers)
    end

    return FunctionCallingCallback(do_stellar_evolution)
end
```

Now, we can put everything together and simulate the system

```julia

masses = [2, 1, 3]u"Msun"
radii = [2, 1, 3]u"Rsun"
stellar_types = [1, 1, 1]

a_in = 1.0u"AU"
a_out = 10.0u"AU"

triple = multibodysystem(masses, a=[a_in, a_out], radii=radii, stellar_types=stellar_types)

callbacks = [StellarEvolutionCallback()]
params = StellarEvolutionParams

res = simulate(triple, t_sim=10u"Myr", 
                       save_everystep=false, 
                       callbacks=callbacks, 
                       params=params)

```