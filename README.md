# Syzygy.jl

> In astronomy, a syzygy (/ˈsɪzədʒi/ SIZ-ə-jee; from Ancient Greek συζυγία (suzugía) 'union, yoke') is a roughly straight-line configuration of three or more celestial bodies in a gravitational system.

`Syzygy.jl` is a [Julia](https://julialang.org/) code for simulating and visualizing the dynamics of hierarchical multistar systems. This package is being developed by someone with very little experience in software development, and is not meant for use by the general public. 

This package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations. This makes the package highly performant and flexible, allowing for both fixed and adaptive timestepping with adjustable error tolerances.  