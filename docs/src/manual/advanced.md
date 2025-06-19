
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

To add a new potential and acceleration function, you first need to define the acceleration function itself, which should be an in-place function that calculates the pairwise acceleration for two particles $(i, j)$ and adds it to the acceleration vectors (dvi and dvj). The function must return `nothing`. 
Example:

```julia
function my_acceleration_function!(dvi, dvj, rs, vs, pair, params)

    i, j = pair
    mi, mj = params.M[i], params.M[j]
    ai = some_acceleration(i, mj)
    aj = some_acceleration(j, mi)

    dvi .+= ai
    dvj .+= aj
    nothing
end
```

You then need to define a potential type, which must be a subtype of the `Syzygy.MultiBodyPotential` supertype.

```julia
struct MyPotential{T1, T2} <: Syzygy.MultiBodyPotential end
```

Finally, you add a method to the `Syzygy.get_accelerating_function` such that it returns a wrapper around your acceleration function. The wrapper has to have the following signature:

```julia
function Syzygy.get_accelerating_function(potential::MyPotential)
    (dvi, dvj, rs, vs, pair, time, params)-> my_acceleration_function!(dvi, dvj, rs, vs, pair, params)
end
```

To use this potential in a simulation, create an instance of your new `MyPotential` type, and include it in the `potential` vector as a keyword argument when initializing the simulation. 

```julia-repl
mypotential = MyPotential()

res = simulate(system, potential=[PureGravitationalAcceleration(), MyPotential()])
```



