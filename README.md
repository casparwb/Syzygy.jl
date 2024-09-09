# Syzygy.jl

<div align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" 
        srcset="logo-big-svg.svg" >
      <img alt="Syzygy.jl logo" 
        src="logo.svg" width="450">
    </picture>
</div>

> In astronomy, a syzygy (/ˈsɪzədʒi/ SIZ-ə-jee; from Ancient Greek συζυγία (suzugía) 'union, yoke') is a roughly straight-line configuration of three or more celestial bodies in a gravitational system.

`Syzygy.jl` is a high-performance few-body simulator written in [Julia](https://julialang.org/). This code is mainly aimed at simulating and visualizing the dynamics of hierarchical multistar systems, but it also supports non-hierarchical few-body systems. 

> [!CAUTION] 
> This package is still being developed, so certain features may not work as expected. If you want to use this package, please ensure that you perform tests to ensure that your results are to be expected.

The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations. This makes the package highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances, and the flexibility of callbacks for code injection.

# Usage

## Installation

`Syzygy.jl` is not a registered package, which means it can be installed using the link to this Github repo. To install it in the REPL, enter the Pkg mode with `]` and just type

```julia
add https://github.com/casparwb/Syzygy.jl
```

Alternatively, from you terminal run

```bash
julia -e 'import Pkg;  Pkg.add(url='https://github.com/casparwb/Syzygy.jl")'
```

## System initialization

A simulation in `Syzygy.jl` begins by setting up the system you want to simulate. This is done by the `multibodysystem` function, which takes in the structural arguments of the bodies in the system - masses, radii, stellar types, luminositites, etc... -, and the orbital parameters of the binaries - semi-major axes, eccentricities, etc.. . The masses are set as the first positional argument, while all other parameters are set using keyword arguments. `Syzygy.jl` uses units by utilizing [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro.jl/stable/), which of which are re-exported upon loading `Syzygy.jl`. Arguments that are not unitless need to be defined with a unit when initializing the system. The system is set up following [Hamers & Portegies Zwart 2016](https://doi.org/10.1093/mnras/stw784), with labeling being done as shown in this figure (Evans 1968):


<img src="https://github.com/casparwb/Syzygy.jl/assets/42641901/971bfd0d-d206-4912-963c-5edd2eeee186" width="500" />

> [!IMPORTANT]
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


> [!NOTE]
> `Syzygy.jl` uses units when initializing a system and when postprocessing a simulation results. The units are discarded before the actual simulation begins in order to ensure performance. The initial conditions and parameters that are sent to the ODE solver are always the preferred units of the respective dimensions. The default unit system is solar mass ($M_\odot$), solar radii ($R_\odot$), and year. If you want to change the default unit system, you have to call the `preferunits` function from `Unitful.jl` BEFORE importing `Syzygy.jl`.

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

`Syzygy.jl` has support for including/varying the potentials and corresponding acceleration functions in a simulation.  In this package, you can specify the acceleration functions to include by setting the `potentials` keyword argument when calling `simulate` or `simulation`. This argument accepts a `Vector{MultiBodyPotential}` where `MultiBodyPotential` is a supertype of all the potentials defined in the package, and for each of which there exists a unique acceleration function. The total acceleration of a particle at each time step is thus a sum of all the acceleration functions defined by the `potential` vector.

Currently, there are 4 possible potentials, with post-newtonian acceleration being being implemented. These are

- `PureGravitationalAcceleration`: Newtonian gravitational acceleration
- `DynamicalTidalAcceleration`: energy dissipation from dynamical tidal, following the prescription of [Samsing et al. 2018](http://arxiv.org/abs/1803.08215).
- `EquilbriumTidalAcceleration`: tidal drag force from equilibrium tides, as described by [Hurley at al. 2002](http://arxiv.org/abs/astro-ph/0201220)
- `StaticEquilibriumTidalAcceleration`: same as above, but with the assumption that masses, radii, and other structural parameters of the particles do not change during the simulation. This is more efficient as it only calculated certain quantities once, and re-uses them throughout the simulation.

Each of these accepts a specific set of arguments related to the acceleration functions. To get an overview of these, enter the `help`-mode in the REPL with `?` and type any of the above names. 

> [!IMPORTANT]
> The numerical implementation of potentials and acceleration functions in this package is heavily inspired by or completely taken from [NbodySimulator.jl](https://github.com/SciML/NBodySimulator.jl). All credit go the authors.

### Example

```julia
grav_pot = PureGravitationalAcceleration()
tidal_pot = DynamicalTidalAcceleration(;kwargs...)
res = simulate(system, potential=[grav_pot, tidal_pot]) 
```

## Running a simulation

Once the simulation has been instantiated, it can be run using `simulate(sim)`. You can also skip the previous step and directly call `simulate` with the exact same arguments as for `simulation`.

```julia
res = simulate(triple, t_sim=10, npoints=1000)
```


## Analysis and visualization

Once a simulation has finished, it can be analyzed and converted to a `MultiBodySolution` type, which allows for easy access to the state vectors and other parameters. This is done using the `analyse_simulation` function, which only takes in the resulting output from `simulate` as its argument. In the `MultiBodySolution` type, the state vectors are stored in 3-dimensional `AxisArray` from [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl), with the axes being (dimension, particle, time). 

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

## Arbitrary precision

By default `Syzygy.jl` uses `Float64` datatypes for the variables of integration, which means that we can not set the error tolerances for the ODE solver lower than $\sim 10^{-14}$ (see [DiffEq docs](https://docs.sciml.ai/DiffEqDocs/stable/basics/faq/#How-to-get-to-zero-error)). While this is generally pretty good, for longer simulation times, this might not be enough. Therefore, `Syzygy.jl` supports arbitrary precision arithmetic, which means you can essentially get as small error as you want, as the cost of speed of course. You can still get pretty decent speed if you use the `Double64` type from (DoubleFloats.jl)[https://github.com/JuliaMath/DoubleFloats.jl]. These have twice as many bits as a standard `Float64`, while still being somewhat performant. If you need even higher precision, you can use the arbitrary precision floating point numbers from (ArbNumerics.jl)[https://github.com/JeffreySarnoff/ArbNumerics.jl]. However, using these will dramatically slow down the code (using Double64 gives a slowdown of $\sim 70\times$, while using ArbFloat with the same number of bits makes the code run over $2000\times$ slower.) To use this functionality, simply specify the `precision` keyword when calling either `simulation` or `simulate`. Fixed-bit types are specified with symbols, while arbitrary precision are specified with an integer, determining the number of bits to use. 

```julia
binary = multibodysystem([2.0, 1.0]u"Msun", a=1.0u"AU")
result = simulate(binary, t_sim=1, save_everystep=false, precision=:Double64) # will use Double64 precision (106 bits)
result = simulate(binary, t_sim=1, save_everystep=false, precision=254) # will use ArbFloat type with 254 bit precision
```

# Advanced Usage

`Syzygy.jl` is designed to be highly composable and flexible, and not only includes additional acceleration functions for including other forces, but also allows the users to define their own potentials and callbacks. 

## Defining your own callbacks

To define you own callback, simply set it up following the the [DiffEq docs](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions), and include it in the `callbacks` vector when initializing a simulation. 


## Adding another potential

To add a new potential and acceleration function, you first need to define the acceleration function itself, which should be an in-place function that calculates the pairwise acceleration for two particles $(i, j)$ and adds it to an acceleration vectors (dvi and dvj). The function must return `nothing`. 
Example:

```julia
function my_acceleration_function!(dvi, dvj, rs, vs, pair, params, potential)

    i, j = pair
    mi, mj = params.M[i], params.M[j]
    ai = acceleration(i, mj)
    aj = acceleration(j, mi)

    dvi .+= ai
    dvj .+= aj
    nothing
end
```

You then need to define a potential type, which must be a subtype of the `Syzygy.MultiBodyPotential` supertype.  

```julia
struct MyPotential <: Syzygy.MultiBodyPotential end
```

You can optionally include parameters in your potential

```julia
struct MyPotential{T1, T2} <: Syzygy.MultiBodyPotential
    parameter1::T1
    parameter2::T2
end
```

Finally, you add a method to the `Syzygy.get_accelerating_function` such that is returns a wrapper around your acceleration function. The wrapper has to have the following signature:

```julia
function Syzygy.get_accelerating_function(parameters::MyPotential, n)
    (dvi, dvj, rs, vs, pair, time, params)-> my_acceleration_function!(dvi, dvj, rs, vs, pair, params, potential)
end
```

To use this potential in a simulation, create an instance of your new `MyPotential` type, and include it in the `potential` vector as a keyword argument when initializing the simulation. 





