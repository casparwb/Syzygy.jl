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
  - icon: <img width="64" height="64" src="https://docs.sciml.ai/DiffEqDocs/stable/assets/logo.png" alt="markdown"/>
    title: Fast and flexible
    details: Optimized for speed, and built on the DifferentialEquations.jl framework.
    link: /manual/simulating.md
  - icon: <img width="64" height="64" src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/refs/heads/master/images/juliadots.iconset/icon_1024x1024.png" />
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

`Syzygy.jl` is a direct [N-body](https://en.wikipedia.org/wiki/N-body_simulation) simulator for astrophysical applications written in [Julia](https://julialang.org/). This code is mainly aimed at simulating and visualizing the dynamics of hierarchical multistar systems and planetary systems, but it also supports non-hierarchical configurations. The package uses the [DifferentialEquations.jl](https://diffeq.sciml.ai/) ecosystem to solve the governing differential equations, making this code highly performant and flexible, allowing for, e.g., both fixed and adaptive timestepping with adjustable error tolerances, easy code injection with callbacks, and of the choice of a wide array of ODE solvers. 



!!! warning
    This package is still being developed, so certain features may not work as expected. If you want to use this package, please make sure that you perform tests to ensure that your results are to be expected.