# Simulating systems

## Simulation setup

Once the system has been initialized, you can set up the parameters of the n-body simulation using the `simulation` function, which takes in the system as the first and only positional argument, followed by keyword arguments for all the other properties. `simulation` accepts almost any argument supported by the common solver options of the `solve` function from [DifferentialEquations.jl](https://diffeq.sciml.ai/). See the [Docs](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options) for a full overview. Note that solvers are not re-exported, thus they can be specified with the `Syzygy` prefix, like `Syzygy.DPRKN6()`.

````@docs
Syzygy.simulation
````

<!-- ### Important keyword arguments
- `t_sim`: total simulation time, which can either be a `Number`, in which case the simulation time will be `t_sim` multiplied by the period of the outermost binary, or it can be a `Quantity`, i.e., a number with a (time) unit, in which case the simulation time is simply `t_sim`.
- `npoints`: See section on Saving controls.
- `saveat`: See section on Saving controls.
- `save_every`: See section on Saving controls.
- `alg`: the solver to use. Default is the 8-th order Runge-Kyutta-Nyström method as described in the DiffEq 2nd order ODE solver [docs](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dynamical_solve/). This is an adaptive solver, thus the only solver arguments are the absolute and relative error tolerances, which can be set using `abs_tol` and `rel_tol` respectively. By default they are both set to  `1e-10`.
- `dt`: time step. Only used if the solver is not adaptive, such as `VerletLeapfrog`, `McAte5`, `Yoshida6`, etc.. Can either be a `Real` or a `Quantity`, where a `Real` would specify a timestep equal to `dt` multiplied by the innermost orbital period.
- `showprogress`: whether to display a progress tracker, showing the current absolute system time and its percentage of `t_sim`. Will slow down the simulation.
- `verbose`: whether to output information about the setup and the final outcome, including total runtime and energy loss.
- `max_cpu_time`: Maximum amount of CPU time a simulation should run. Default is Inf.
- `callbacks`: see section on callbacks.
- `potential`: see section on potentials. -->

### Examples
```julia
sim = simulation(triple, t_sim=10, npoints=10_000) # simulate for 10 outer orbits, and save 10 000 snapshots uniformly distributed
sim = simulation(triple, t_sim=10, saveat=1u"d") # save 1 day
sim = simulation(triple, t_sim=10, save_every=10) # save every 10th timestep
sim = simulation(triple, t_sim=10u"kyr", save_everystep=false) # simulate for 10 000 years, and only save the initial and final state of the system
sim = simulation(triple, t_sim=1, alg=Syzygy.McAte5(), dt=1e-5) # use a timestep of 1e-5 * P_in
sim = simulation(triple, t_sim=1, alg=Syzygy.McAte5(), dt=1.0u"s") # use a timestep of 1 second
```

## Saving controls

Syzygy provides flexibility in how to save the output of a simulation. By default, the solver takes a snapshot at each time step. The output will also be dense, meaning you can interpolate into new time points. To set other saving options, you can send the following keyword arguments to `simulation`:

- `npoints`: number of snapshots to save. With this parameter, a snapshot will be taken every `(t_end - t_start)/npoints`, and will therefore be uniformly distributed in time. 
- `saveat`: Times at which to save a snapshot. If given as a number with a time unit, a snapshot will be taken every `saveat`. If given as an array, a snapshot will be taken at those times.
- `save_every`: Specify that you want to save a snapshot every `save_every`-th time step. A value of e.g. 2 will save a snapshot every 2nd time step. 
- `save_everystep`: If `false`, only save the initial and final steps of the solution, otherwise save every step. Default is `true`.

## Callbacks

The package contains several pre-defined callbacks that can be specified when setting up the simulation. These mostly define stopping conditions, but can also be e.g, flags or events. Callbacks can be specified using the `callbacks` keyword argument when calling `simulation`, and be given as a `Vector{<:AbstractSyzygyCallback}`. There are number of pre-defined callbacks, but you can also define your own. (see Advanced usage for more details). Some of the more important pre-defined callbacks are:

- `CollisionCB`: a stopping condition that terminates the integration when two objects gets close enough. The exact condition for specifying a collision depends on the stellar types of the objects. If the two objects are both stars, then the code simply checks if the radii overlap. If the objects are stellar remnants (black holes, neutron stars, white dwarfs), then the code checks if the distance is smaller than some multiple of their mutual gravitational radius. By default, this multiple value is set to 10000. If the two objects contain one star and one remnant, then the code checks if the star is within the tidal disruption radius of the remnant object.
- `EscapeCB`: a stopping condition for checking whether an object has been ejected from the system. The callback has three checks that need to be passed in order for the condition to be invoked. See Standish & Myles 1971 for more details. Currently this callback is only possible to use for a triple system.


Certain callbacks accept arguments when instantiating them. This can include how often the solver should check for the condition (for every n-th step). This is implemented as the field `check_every`. 


### Examples
```julia
callbacks = [CollisionCB(5, 10_000), EscapeCB()]
sim = simulation(system, t_sim=1, callbacks=callbacks)
```

## Potentials

`Syzygy.jl` has support for including/varying the potentials and corresponding acceleration functions in a simulation.  In this package, you can specify the acceleration functions to include by setting the `potential` keyword argument when calling `simulate` or `simulation`. This argument accepts a `Vector{MultiBodyPotential}` where `MultiBodyPotential` is a supertype of all the potentials defined in the package, and for each of which there exists a unique acceleration function. The total acceleration of a particle at each time step is thus a sum of all the acceleration functions defined by the `potential` vector.

Currently, there are a number of possible potentials:

- `PureGravitationalAcceleration`: Newtonian gravitational acceleration
- `DynamicalTidalAcceleration`: energy dissipation from dynamical tides, following the prescription of [Samsing et al. 2018](http://arxiv.org/abs/1803.08215).
- `EquilibriumTidalAcceleration`: tidal drag force from equilibrium tides, as described by [Hurley at al. 2002](http://arxiv.org/abs/astro-ph/0201220), using the model of [Hut 1981](https://ui.adsabs.harvard.edu/abs/1981A&A....99..126H).
- `TimeDependentEquilbriumTidalAcceleration`: same as above, but with the assumption that masses, radii, and other structural parameters of the particles can change during the simulation, which means that structural properties such as envelope mass and radius are calculated during the simulation.

There are also potentials for various general relatistic effects via post-Newtonian expansion. See the following section on post-Newtonian physics.

Each of these accepts a specific set of arguments related to the acceleration functions. To get an overview of these, enter the `help`-mode in the REPL with `?` and type any of the above names. 

!!! important
    The numerical implementation of potentials and acceleration functions in this package is heavily inspired by or completely taken from [NbodySimulator.jl](https://github.com/SciML/NBodySimulator.jl). Credit go the authors.

### Example

```julia
grav_pot = PureGravitationalAcceleration()
tidal_pot = DynamicalTidalAcceleration(;kwargs...)
res = simulate(system, potential=[grav_pot, tidal_pot]) 
```

## Post-Newtonian expansion

For systems where general relativistic effects can become important, `Syzygy.jl` supports post-Newtonian acceleration potentials up to order 2.5, following the formalism of [Blanchet 2015](http://arxiv.org/abs/1310.1528).

- `PN1Potential`: post-Newtonian acceleration of order 1.
- `PN2Potential`: post-Newtonian acceleration of order 2.
- `PN2p5Potential`: post-Newtonian acceleration of order 2.5.
- `PNPotential`: post-Newtonian acceleration of order 1 to 2.5. This is more efficient than including all of the above potentials individually.

!!! note
    The default ODE algorithm (`DPRKN8`) does not support velocity-dependent potentials, and will therefore not 
    produce correct trajectories when simulating system with post-Newtonian terms. Other algorithms that do support 
    velocity dependent potentials are `Vern7`, `Vern8`, `Vern9`, and `FineRKN5`.

```julia
sim = simulation(system, potential=[PureGravitationalPotential(), 
                                    PNPotential()], 
                        alg=Syzygy.FineRKN5())
```

## Running a simulation

Once the simulation has been instantiated, it can be run using `simulate(sim)`. You can also skip the previous step and directly call `simulate` with the exact same arguments as for `simulation`.

```julia
res = simulate(triple, t_sim=10, npoints=1000)
```


## Analysis and visualization

Once a simulation has finished, it can be analyzed and converted to a `MultiBodySolution` type, which allows for easy access to the state vectors and other parameters. This is done using the `to_solution` function, which only takes in the resulting output from `simulate` as its argument. In the `MultiBodySolution` type, the state vectors are stored in 3-dimensional `AxisArray` from [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl), with the axes being (dimension, particle, time). 

```julia
sol = to_solution(res)
N = length(sol.t)
r1 = sol.r[particle=1] # a (3, N) array, alias for sol.r[:,1,:]
v1 = sol.v[particle=1]
```

`Syzygy.jl` contains a few plot recipes for [Plots.jl](https://docs.juliaplots.org/latest/tutorial/), which can be used to create quick and simple visualizations of the orbits. To use this, simply import `Plots.jl`, and use the function `orbitplot`, which takes in a `MultiBodySolution` as its first argument, followed by keyword arguments:

- `dims`: which dimensions to plot. Should be a `Vector{Int}`, with any combination of `(1, 2, 3)` (x, y, z).
- `bodies`: which bodies (particles) to plot. Should be a `Vector{Int}`, with any combination of `1:N_particles`.
- `tspan`: the temporal selection to plot. Should be a `Tuple{Quantity, Quantity}`, with the start- and end-points.

```julia
using Plots
orbitplot(sol, bodies=[1, 2], dims=[1, 2]) # plot only particles 1 and 2 in the x-y plane
```

## Higher precision

By default `Syzygy.jl` uses `Float64` datatypes for the variables of integration, which means that we can not set the error tolerances for the ODE solver lower than $\sim 10^{-14}$ (see [DiffEq docs](https://docs.sciml.ai/DiffEqDocs/stable/basics/faq/#How-to-get-to-zero-error)). While this is generally pretty good, for longer simulation times or systems with very close passages, this might not be enough. Therefore, `Syzygy.jl` supports higher and arbitracy arithmetic, which means you can essentially get as small error as you want tT the cost of speed. This is done through the (DoubleFloats.jl)[] and (ArbNumerics.jl)[] packages. You get the best performance using the `Double64` type from the former package, which has twice the precision of `Float64`s. If you need even higher precision, you can use the arbitrary precision floating point numbers from (ArbNumerics.jl)[https://github.com/JeffreySarnoff/ArbNumerics.jl]. However, using these will dramatically slow down the code. To use this functionality, simply specify the `precision` keyword when calling either `simulation` or `simulate`. Fixed-bit types are specified with symbols, while arbitrary precision are specified with an integer, determining the number of bits to use. 

```julia
binary = multibodysystem([2.0, 1.0]u"Msun", a=1.0u"AU")
result = simulate(binary, t_sim=1, save_everystep=false, precision=:Double64) # will use Double64 precision (106 bits)
result = simulate(binary, t_sim=1, save_everystep=false, precision=254) # will use ArbFloat type with 254 bit precision
```
