
# Advanced Usage

`Syzygy.jl` is designed to be highly composable and flexible, and not only includes additional acceleration functions for including other forces, but also allows the users to define their own potentials and callbacks. 

## Defining your own callbacks

To define you own callback, first define a struct which is a subtype of `Syzygy.AbstractSyzygyCallback`. For example:

```julia
struct MyCallback <: Syzygy.AbstractSyzygyCallback end
```

You can then add a method to the `Syzygy.get_callback` function, which has the signature `get_callback(cb, system, retcodes, args)`, where `cb` is an object with a type of your callback struct, `system` is the system that is simulated, retcodes is a dictionary in which you can add return codes/flags, and args is a dictionary with additional arguments. Inside this function, you set up your `affect!` and `condition` functions, which determine when the callback is triggered, and what should happen when it does. It should then return a callback type from the SciMLBase package.

For example, if you want to add a callback that terminates the integration if it runs for a certain amount of time, you would do:

```julia

function Syzygy.get_callback(cb::MyCallback, system, retcodes, args)
    max_cpu_time = args[:max_cpu_time]
    t_start = time()
    
    condition_cpu_time(u, t, integrator) = true # check every step
    
    function max_cpu_time_callback!(integrator)
        if (time() - t_start) >= max_cpu_time
            retcode[:MaxCPUTime] = true
            terminate!(integrator)
        end
    end

    affect_cpu_time!(integrator) = max_cpu_time_callback!(integrator)
    
    callback_cpu_time = DiscreteCallback(condition_cpu_time, affect_cpu_time!, save_positions=(false, false))

    return callback_cpu_time
end
```


## Adding another potential

To add a new potential and acceleration function, you first need to define a potential type which will dispatch to the corresponding acceleration function. This which must be a subtype of the `Syzygy.MultiBodyPotential` supertype. If your potential has inherent constants (like the gravitational constant), the add a method for calculating the value of the given constant in the units of the system  For example:

```julia
struct MyPotential{T} <: Syzygy.MultiBodyPotential 
    G::T
    my_const::T
end

function MyPotential(system)
    G = Syzygy.get_G_in_system_units(system) # built-in function from Syzygy
    my_const = calculate_constant(system) # user-defined function

    return MyPotential(G, my_const)
end
```


Then you need to define the acceleration function itself, which should be an in-place function that calculates the pairwise acceleration for two particles `(i, j)`, whose indices are input as the tuple `pair`, and adds it to the acceleration matrix `dv`. The function must return `nothing`.

```julia
function my_acceleration_function!(dv, rs, vs, pair, params, pot::MyPotential)

    i, j = pair
    mi, mj = params.masses[i], params.masses[j]
    
    ai = some_acceleration(i, mj, pot.G, pot.my_const)
    aj = some_acceleration(j, mi, pot.G, pot.my_const)

    dv[1, i] += ai[1]
    dv[1, j] += aj[1]

    dv[2, i] += ai[2]
    dv[2, j] += aj[2]

    dv[3, i] += ai[3]
    dv[3, j] += aj[3]
    nothing
end
```


Finally, you add a method to the `Syzygy.get_accelerating_function` such that it returns a wrapper around your acceleration function. The wrapper has to have the following signature:

```julia
function Syzygy.get_accelerating_function(potential::MyPotential)
    (dv, rs, vs, pair, time, params)-> my_acceleration_function!(dv, rs, vs, pair, params, potential)
end
```

To use this potential in a simulation, create an instance of your new `MyPotential` type, and include it in the `potential` vector as a keyword argument when initializing the simulation. 

```julia-repl
mypotential = MyPotential(system)

res = simulate(system, potential=[mypotential])
```

### Parameters

You may have noticed that the acceleration function also accepts a `params` argument. This is an object which should contain properties of the particles, such as their masses and other parameters that are used by the acceleration function and/or callbacks. By default, the param object used is `Syzygy.DefaultSimulationParams`, which has the fields `radii`, `masses`, `stellar_types`, and `stellar_type_numbers`. You can define your own parameters, which can then be accessed in your own acceleration function. To do this, first define a new struct as a subtype of `Syzygy.SimulationParams`. For example:

```julia
struct MySimulationParams{T, T2} <: Syzygy.SimulationParams
    radii::T
    masses::T
    my_other_properties::T2
end
```

You then need to define a function which initializes this struct with the correct values for your given system. To do this, add a new method to the function `Syzygy.setup_params` which dispatches on the type of your previously defined struct:

```julia
function Syzygy.setup_params(::Type{<:MySimulationParams}, system, datatype=Float64; options)

    # preprocess the values that will go into the fields
    # return MySimulationParams(....)
end
```

Here, `options` should be a dictionary or named tuple, which can contain switches and other options for fine-tuning the parameter. Options can be specified using the `param_options` keyword in `simulate`.

To actually ensure that your simulation uses these parameter as opposed to the default ones, set the `params` keyword in `simulation` to the type of your parameter struct.

```julia-repl
res = simulate(system, potential=[mypotential], 
                       params=MySimulationParams, 
                       params_options=Dict(:some_switch => true))
```

