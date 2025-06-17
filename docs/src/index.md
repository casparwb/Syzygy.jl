```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Syzygy.jl"
  tagline: "A fast, flexible, direct n-nody integrator written in pure Julia."
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


---
```

`Syzygy.jl` is a fast and flexible few-body simulator for astrophysical applications written in [Julia](https://julialang.org/). This code is mainly aimed at simulating and visualizing the dynamics of hierarchical multistar systems and planetary systems, but it also supports non-hierarchical configurations. The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations, making this code highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances, the flexibility of callbacks for code injection, and of the choice of a wide array of ODE solvers. 



> [!CAUTION] 
> This package is still being developed, so certain features may not work as expected. If you want to use this package, please make sure that you perform tests to ensure that your results are to be expected.