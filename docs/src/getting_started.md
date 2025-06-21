# Getting started

## Installation

`Syzygy.jl` is not a registered package, which means it can be installed using the link to this Github repo. To install it in the REPL, enter the Pkg mode with `]` and just type

```julia
julia> ] # enter Pkg mode
pkg> add https://github.com/casparwb/Syzygy.jl
```

Alternatively, from your terminal run

```bash
> julia -e 'import Pkg;  Pkg.add(url="https://github.com/casparwb/Syzygy.jl")'
```

Or, if you want to make changes to the code, just clone the repo and do

```julia
julia> ]
pkg> dev path/to/Syzygy.jl
```

## Quick-start guide

A system can be set up using the [`multibodysystem`](@ref) function, which takes a vector of masses as the first argument, followed by keyword argument specifying orbital and stellar properties. Parameters have to be specified using units from [Unitful.jl](https://github.com/PainterQubits/Unitful.jl/), which are automatically imported when you load Syzygy.

```julia
binary = multibodysystem([2.0, 1.0]u"Msun", semi_major_axis=1.0u"AU", eccentricity=0.1)
triple = multibodysystem([2.0, 1.0, 3.0]u"Msun", a=[0.1, 1.0]u"AU", e=[0.1, 0.2])
```

!!! tip "Argument aliases"
    Many of the keyword arguments for `multibodysystem` have several alises. For example, for the semi-major axis, you can specify any of `sma, a, semi_major_axis, semi_major_axes, semimajor_axis`, or `semimajor_axes`. To see all the alises for all the parameters, call `Syzygy.multibodysystem_parameter_aliases`.

You can also specify a system using only masses, positions, and velocities as positional arguments, and stellar properties as keyword arguments.

```julia
three_body_system = multibodysystem([1.0, 1.0]u"kg", [rand(3), rand(3)]u"m", [rand(3), rand(3)]u"m/s", radii=[1.0, 1.0]u"m")
```

You can simulate the system by calling [`simulate`](@ref).

```julia
res = simulate(binary, t_sim=1) # simulate for 1 period
res = simulate(binary, t_sim=1u"yr") # simulate for 1 year
```

Finally, you can postprocess the result to make it easier for further analysis by calling `to_solution`. This allows you to access the positions and velocities of each body using the fields `r` and `v`, in addition to other structural arguments. 

```julia
sol = to_solution(res)
# r and v is stored as (N_dims, N_bodies, N_timesteps)
r1 = sol.r[:,1,:] # return an (N_dims, N_timesteps) matrix
r1 == sol.r[particle=1] 
```

`Syzygy.jl` has a few plotting recipes for `Plots.jl`, which allows you to quickly look at orbits, energy evolution, accelerations, and others. 

```julia
using Plots

orbitplot(sol) # plot all bodies at all timesteps
orbitplot(sol, dims=[1, 2]) # only plot the x-y coordinates
orbitplot(sol, tspan=(0.0, 1.0)u"d") # only plot the orbit in the given timespan
orbitplot(sol, bodies=[1, 3]) # only plot bodies 1 and 3

kozailidovplot(sol) # plot eccentricities and mutual inclination (only works for triples)

```


