# Syzygy.jl

> In astronomy, a syzygy (/ˈsɪzədʒi/ SIZ-ə-jee; from Ancient Greek συζυγία (suzugía) 'union, yoke') is a roughly straight-line configuration of three or more celestial bodies in a gravitational system.

`Syzygy.jl` is a [Julia](https://julialang.org/) code for simulating and visualizing the dynamics of hierarchical multistar systems. This package is being developed by someone with very little experience in software development, and is not meant for use by the general public. 

The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations. This makes the package highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances. 

# Usage

## System initialization

A simulation in `Syzygy.jl` begins by setting up the system you want to simulate. This is done by the `multibodysystem` function, which takes in the structural arguments of the bodies in the system - masses, radii, stellar types, luminositites, etc... -, and the orbital parameters of the binaries - semi-major axes, eccentricities, etc.. . The masses are set as the first positional argument, while all other parameters are set using keyword arguments. `Syzygy.jl` uses units by utilizing [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro.jl/stable/), which of which are re-exported upon loading `Syzygy.jl`. Arguments that are not unitless need to be defined with a unit when initializing the system. The system is set up following [Hamers & Portegies Zwart 2016](https://doi.org/10.1093/mnras/stw784), with labeling being done as shown in this figure (Evans 1968):


<img src="https://github.com/casparwb/Syzygy.jl/assets/42641901/971bfd0d-d206-4912-963c-5edd2eeee186" width="500" />

The code for setting up the hierarchy (converting the orbital elements into state vectors) is either heavily inspired by, or directly taken from [NbodyGradient.jl](https://github.com/ericagol/NbodyGradient.jl). All credits go to the authors.

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

The components of the system can be accessed via `system.binaries[index]`, with `index` being an integer from 1 to number of binaries, and `system.particles[index]`, for the individual particles. A binary or particle are instances of the `Binary` or `Particle` type, both of which have a number of fields containing information about their constituents. For a `Particle`, you can access 

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

`Syzygy.jl` is designed to be highly composable and flexible, and allows not only includes additional acceleration functions for including other forces, but also allows the users to define their own potentials and callbacks.
