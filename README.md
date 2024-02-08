# Syzygy.jl

> In astronomy, a syzygy (/ˈsɪzədʒi/ SIZ-ə-jee; from Ancient Greek συζυγία (suzugía) 'union, yoke') is a roughly straight-line configuration of three or more celestial bodies in a gravitational system.

`Syzygy.jl` is a [Julia](https://julialang.org/) code for simulating and visualizing the dynamics of hierarchical multistar systems. This package is being developed by someone with very little experience in software development, and is not meant for use by the general public. 

The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations. This makes the package highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances. 

# Usage

## System initialization

A simulation in `Syzygy.jl` begins by setting up the system you want to simulate. This is done by the `multibodysystem` function, which takes in the structural arguments of the bodies in the system - masses, radii, stellar types, luminositites, etc... -, and the orbital parameters of the binaries - semi-major axes, eccentricities, etc.. . The masses are set as the first positional argument, while all other parameters are set using keyword arguments. `Syzygy.jl` uses units by utilizing [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro.jl/stable/), which of which are re-exported upon loading `Syzygy.jl`. Arguments that are not unitless need to be defined with a unit when initializing the system. The system is set up following (Hamers & Portegies Zwart 2016)[https://doi.org/10.1093/mnras/stw784], with labeling being done as shown in this figure (Evans 1968)

![Screenshot 2024-02-08 143333](https://github.com/casparwb/Syzygy.jl/assets/42641901/971bfd0d-d206-4912-963c-5edd2eeee186)

### Examples


```julia

binary = multibodysystem([1.0, 1.0]u"Msun", a=1.0u"AU", e=0.4) # set up a binary system with two 1 solar-mass stars, in an orbit with 1 semi-major axis of 1 AU and an eccentricity of 0.4
triple = multibodysystem([2.0, 1.0, 3.0]u"Msun", a=[0.1, 0.5]u"Rsun", e=[0.1, 0.4], i=[90.0, 0.0]u"degree") # hierarchical triple
quadruple = multibodysystem([1.0, 1.0, 1.0, 1.0]u"Msun", a=[0.1, 0.5, 10.0]u"Rsun", e=[0.1, 0.4, 0.2], i=[90.0, 45.0, 0.0]u"degree", hierarchy=[4, 2, 1]) # 2+2 quadruple
```

The components of the system can be accessed via `system.binaries[index]`, with `index` being an integer from 1 to number of binaries, and `system.particles[index]`, for the individual particles.

## Simulation setup

## Running a simulation

