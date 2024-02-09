# Syzygy.jl

> In astronomy, a syzygy (/ˈsɪzədʒi/ SIZ-ə-jee; from Ancient Greek συζυγία (suzugía) 'union, yoke') is a roughly straight-line configuration of three or more celestial bodies in a gravitational system.

`Syzygy.jl` is a [Julia](https://julialang.org/) code for simulating and visualizing the dynamics of hierarchical multistar systems. This package is being developed by someone with very little experience in software development, and is not meant for use by the general public. 

The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations. This makes the package highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances. 

# Usage

## System initialization

A simulation in `Syzygy.jl` begins by setting up the system you want to simulate. This is done by the `multibodysystem` function, which takes in the structural arguments of the bodies in the system - masses, radii, stellar types, luminositites, etc... -, and the orbital parameters of the binaries - semi-major axes, eccentricities, etc.. . The masses are set as the first positional argument, while all other parameters are set using keyword arguments. `Syzygy.jl` uses units by utilizing [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro.jl/stable/), which of which are re-exported upon loading `Syzygy.jl`. Arguments that are not unitless need to be defined with a unit when initializing the system. The system is set up following [Hamers & Portegies Zwart 2016](https://doi.org/10.1093/mnras/stw784), with labeling being done as shown in this figure (Evans 1968):


<img src="https://github.com/casparwb/Syzygy.jl/assets/42641901/971bfd0d-d206-4912-963c-5edd2eeee186" width="500" />

> The code for setting up the hierarchy (converting the orbital elements into state vectors) is either heavily inspired by, or directly taken from [NbodyGradient.jl](https://github.com/ericagol/NbodyGradient.jl). All credits go to the authors.

### Examples


```julia

binary = multibodysystem([1.0, 1.0]u"Msun", a=1.0u"AU", e=0.4) # set up a binary system with two 1 solar-mass stars, in an orbit with 1 semi-major axis of 1 AU and an eccentricity of 0.4
triple = multibodysystem([2.0, 1.0, 3.0]u"Msun", a=[0.1, 0.5]u"Rsun", e=[0.1, 0.4], i=[π/2, 0.0]u"rad") # hierarchical triple
quadruple = multibodysystem([1.0, 1.0, 1.0, 1.0]u"Msun", a=[0.1, 0.5, 10.0]u"Rsun", e=[0.1, 0.4, 0.2], i=[90.0, 45.0, 0.0]u"degree", hierarchy=[4, 2, 1]) # 2+2 quadruple

# set up a binary black hole system
bh1_mass = 9.62u"Msun"
bh2_mass = 8.4u"Msun"

bh1_radius = 2*𝒢*bh1_mass/c² # the gravitational constant G can be accessed using \scrG in the REPL. 
bh2_radius = 2*𝒢*bh2_mass/c²

binary_blackholes = multibodysystem([bh1_mass, bh2_mass], a=15.3u"Rsun", R=[bh1_radius, bh2_radius], type=[14, 14]) 
```

If you wanted to set up the system in the above figure, a quintuple system, you would do

```julia
masses = ones(5)u"Msun"
r1 = 1.0u"Rsun"
r2 = 5.0u"Rsun"
r3 = 1.0u"Rsun"
r4 = 10.0u"Rsun"
hierarchy = [5, 1, 2, 1] # 1 binary on level 2, 2 on 1, and 1 on 2

quintuple = multibodysystem(masses, a=[r1, r2, r3, r4], hierarchy=hierarchy)
```

The components of the system can be accessed via `system.binaries[index]`, with `index` being an integer from 1 to number of binaries, and `system.particles[index]`, for the individual particles. A binary or particle are instances of the `Binary` or `Particle` type, both of which have a number of fields containing information about their state and internal structure.



## Simulation setup

Once the system has been initialized, you can set up the parameters of the n-body simulation using the `simulation` function, which takes in the system as the first and only positional argument, followed by keyword arguments for all the other properties. `simulation` accepts almost any argument supported by the common solver options of the `solve` function from [DifferentialEquations.jl](https://diffeq.sciml.ai/). See the [Docs](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options) for a full overview. Note that solvers are not re-exported, thus they can be specified with the `Syzygy` prefix, like `Syzygy.DPRKN6()`.

### Important keyword arguments
- `t_sim`: total simulation time, which can either be a `Number`, in which case the simulation time will be `t_sim` multiplied by the period of the outermost binary, or it can be a `Quantity`, i.e., a number with a (time) unit, in which case the simulation time is simply `t_sim`.
- `npoints`: number of snapshots to save. With this parameter, a snapshot will be taken every `(t_end - t_start)/npoints`, and will therefore be uniformly distributed in time. You can also use the `saveat` argument from `DiffEq`, which can take in an array of specific times, or other types.
- `alg`: the solver to use. Default is the 8-th order Runge-Kyutta-Nyström method as described in the DiffEq 2nd order ODE solver [docs](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dynamical_solve/). This is an adaptive solver, thus the only solver arguments are the absolute and relative error tolerances, which can be set using `abs_tol` and `rel_tol` respectively. By default they are both set to  `1e-10`.
- `dt`: time step. Only used if the solver is not adaptive, such as `VerletLeapfrog`, `McAte5`, `Yoshida6`, etc.. Can either be a `Number` or a `Quantity`, where a `Number` would specify a timestep equal to `dt` multiplied by the innermost orbital period.
- `showprogress`: whether to display a progress tracker, showing the current absolute system time and its percentage of `t_sim`. Will slow down the simulation.
- `verbose`: whether to output information about the setup and the final outcome, including total runtime and energy loss.
- `callbacks`: see section on callbacks.
- `potentials`: see section on potentials.

### Examples
```julia
sim = simulation(triple, t_sim=10, npoints=10_000) # simulate for 10 outer orbits, and save 10 000 snapshots
sim = simulation(triple, t_sim=10u"kyr", save_everystep=false) # simulate for 10 000 years, and only save the initial and final state of the system
sim = simulation(triple, t_sim=1, alg=Syzygy.McAte5(), dt=1e-5) # use a timestep of 1e-5 * P_in
```

## Callbacks

The package contains several pre-defined callbacks that can be specified when setting up the simulation. They are mostly stopping conditions, but can also be flags or other code injections. They can be specified using the `callbacks` argument, and should be given either as a `String`, which specifies the pre-defined callbacks, or it can be a callback type as described in [Event Handling and Callback Functions](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/#The-Callback-Types) if the user creates their own (see Advanced usage for more details). The most important pre-defined callbacks are

- `collision`: a stopping condition that is invoked when two particles have overlapping radii. If the systems contains compact objects, the [tidal disruption radius](https://en.wikipedia.org/wiki/Tidal_disruption_event#Tidal-disruption_radius) or the [gravitational radius](https://en.wikipedia.org/wiki/Schwarzschild_radius) is used to check for collisions between a CO and a Star, or two COs.
- `escape`: a stopping condition for checking whether an object has been ejected from the system. The callback has three checks that need to be passed in order for the condition to be invoked. See Standish & Myles 1971 for more details. Currently this callback is only possible to use for a triple system.

To get the full list of callbacks, you can call `Syzygy.callbacks()`.

## Potentials

`Syzygy.jl` has support for including/varying the potentials and corresponding acceleration functions in a simulation.  In this package, you can specify the acceleration functions to include by setting the `potentials` keyword argument when calling `simulate` or `simulation`. This argument accepts a `Vector{FewBodyPotential}` where `FewBodyPotential` is a supertype of all the potentials defined in the package, and for each of which there exists a unique acceleration function. The total acceleration of a particle at each time step is thus a sum of all the acceleration functions defined by the `potential` vector. Internally, this looks something like

```julia
function acceleration_function(potential::PureGravitationalPotential)
    return (u, p, t) -> newtonian_gravity(u, p, t, potential)
end

function acceleration_function(potential::PostNewtonianPotential)
    return (u, p, t) -> post_newtonian_acceleration(u, p, t, potential)
end
```

Currently, there are 4 possible potentials, with post-newtonian acceleration being being implemented. These are

- `PureGravitationalAcceleration`: Newtonian gravitational acceleration
- `DynamicalTidalAcceleration`: energy dissipation from dynamical tidal, following the prescription of [Samsing et al. 2018](http://arxiv.org/abs/1803.08215).
- `EquilbriumTidalAcceleration`: tidal drag force from equilibrium tides, as described by [Hurley at al. 2002](http://arxiv.org/abs/astro-ph/0201220)
- `StaticEquilibriumTidalAcceleration`: same as above, but with the assumption that masses, radii, and other structural parameters of the particles do not change during the simulation. This is more efficient as it only calculated certain quantities once, and re-uses them throughout the simulation.

Each of these accepts a specific set of arguments related to the acceleration functions. To get an overview of these, enter the `help`-mode in the REPL with `?` and type any of the above names. 

> The numerical implementation of potentials and acceleration functions in this package is again heavily inspired by or completely taken from [NbodySimulator.jl](https://github.com/SciML/NBodySimulator.jl), so all credit go the authors.
### Example

```julia
grav_pot = PureGravitationalAcceleration()
tidal_pot = DynamicalTidalAcceleration(;kwargs...)
res = simulate(triple, potentials=[grav_pot, tidal_pot]) 
```

## Running a simulation

Once the simulation has been instantiated, it can be run using `simulate(sim)`. You can also skip the previous step and directly call `simulate` with the exact same arguments as for `simulation`.

```julia
res = simulate(triple, t_sim=10, npoints=1000)
```


## Analysis and visualization

Once a simulation has finished, it can be analyzed and converted to a `FewBodySolution` type, which allows for easy access to the state vectors and other parameters. This is done using the `analyse_simulation` function, which only takes in the resulting output from `simulate` as its argument. In the `FewBodySolution` type, the state vectors are stored in 3-dimensional `LArray` from [LabelledArrays.jl](https://github.com/SciML/LabelledArrays.jl), with the axes being (dimension, particle, time). The time evolution of the orbital elemenst of each binary can also be accessed from the `elements`-field, which is a `Vector{OrbitalElements{`, with the index corresponding to the original binary indices.

```julia
sol = analyse_simulation(res)
N = length(sol.t)
r1 = sol.r[particle=1] # a (3, N) array
v1 = sol.v[particle=1]

r1 == sol.r[:,1,:]

inner_sma = sol.elements[1].a
```

`Syzygy.jl` contains a few plot recipes for [Plots.jl](https://docs.juliaplots.org/latest/tutorial/), which can be used to create quick and simple visualizations of the orbits. To use this, simply import `Plots.jl`, and use the function `orbitplot`, which takes in a `FewBodySolution` as its first argument, followed by keyword arguments:

- `dims`: which dimensions to plot. Should be a `Vector{Int}`, with any combination of `(1, 2, 3)` (x, y, z).
- `bodies`: which bodies (particles) to plot. Should be a `Vector{Int}`, with any combination of `1:N_particles`.
- `tspan`: the temporal selection to plot. Should be a `Tuple{Quantity, Quantity}`, with the start- and end-points.

```julia
using Plots
orbitplot(sol, bodies=[1, 2], dims=[1, 2]) # plot only particles 1 and 2 in the x-y plane
```

# Advanced Usage

`Syzygy.jl` is designed to be highly composable and flexible, and not only includes additional acceleration functions for including other forces, but also allows the users to define their own potentials and callbacks. 

## Defining your own callbacks

To define you own callback, simply set it up following the the [DiffEq docs](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions), and include it in the `callbacks` vector when initializing a simulation. 


## Adding another potential

To add a new potential and acceleration function, you first need to define the acceleration function itself, which should be an in-place function that calculates the force for one particle and adds it to an acceleration vector (dv). Example:

```julia
function my_acceleration_function!(dv, v, r, i, n, p)
    accel = @SVector zeros(3) # for storing the x,y,z acceleration on particle i

    for k = 1:n
        if k != i
            force = some_interaction(i, k)
            accel .+= force
        end
    end
    dv .+= accel

end
```

You then need to define a potential type, which must be a subtype of the `Syzygy.FewBodyPotential` supertype. 

```julia
struct MyPotential{T1, T2} <: Syzygy.FewBodyPotential
    parameter_1::T1
    parameter_2::T2
end
```

Finally, you add a method to the `Syzygy.get_accelerating_function` such that is returns a wrapper around your acceleration function that can be used by DiffEq:

```julia
function Syzygy.get_accelerating_function(parameters::MyPotential, n)
    (dv, u, v, p, t, i) -> my_acceleration_function!(dv, u, parameters, i, n)
end
```

To use this potential in a simulation, create an instance of your new `MyPotential` type, and include it in the `potentials` vector as a keyword argument when initializing the simulation. 





