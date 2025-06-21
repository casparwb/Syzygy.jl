```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Syzygy.jl"
  tagline: "A fast, flexible, direct n-body integrator written in pure Julia."
#   tagline: A Markdown backend designed to work with VitePress and Documenter.jl
  image:
    src: /logo.png
    alt: Syzygy
  actions:
    - theme: brand
      text: Get Started
      link: /getting_started
    - theme: alt
      text: View on Github
      link: https://github.com/casparwb/Syzygy.jl
    - theme: alt
      text: API
      link: /api

features:
  - 
    title: Fast and flexible
    details: Optimized for speed, and built on the DifferentialEquations.jl framework.
    link: /manual/simulating.md
  - 
    title: Composable
    details: Easily extend the code with your own potentials and conditions.
    link: /manual/advanced.md
  - 
    title: Feature-rich
    details: Supports high-precision numerics, post-Newtonian potentials, tidal prescriptions, stellar evolution parameters such as stellar types and stellar structure, and much more. 
    link: /manual/advanced.md
---
---
```

```@example
sleep(0.5) # hide
```

`Syzygy.jl` is a direct [N-body](https://en.wikipedia.org/wiki/N-body_simulation) simulator for astrophysical applications written in [Julia](https://julialang.org/). This code is mainly aimed at simulating and visualizing the dynamics of hierarchical multistar systems and planetary systems, but it also supports non-hierarchical configurations. The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations, making this code highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances, easy code injection with callbacks, and of the choice of a wide array of ODE solvers. 

```@example
using Plots, Syzygy
Plots.theme(:dracula) # hide

unstable_triple = multibodysystem([2.0, 1.0, 1.5]u"Msun", 
                                  semi_major_axes=[1.0, 4.0]u"AU", 
                                  eccentricities=[0.0, 0.4])
res = simulate(unstable_triple, t_sim=4)
sol = to_solution(res)
orbitplot(sol, dims=[1, 2])
```


!!! warning
    This package is still being developed, so certain features may not work as expected. If you want to use this package, please make sure that you perform tests to ensure that your results are as you expect.

