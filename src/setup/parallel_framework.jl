function DiffEqBase.SecondOrderODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        u0, v0, ::MultiThreading
    )

    accelerations = acc_funcs.fs
    println("Hello")
    n_particles = size(u0, 2)
    ais = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_particles]
    function soode_system!(dv, v, u, p, t)
        for acc in accelerations
            Threads.@threads for i = 1:n_particles
                ai = ais[i]
                fill!(ai, 0.0)
                for j = 1:n_particles
                    i == j && continue
                    ai += acc(u, i, j, p)
                end

                dv[1,i] = ai[1]
                dv[2,i] = ai[2]
                dv[3,i] = ai[3]
            end
        end
        return nothing
    end

    return SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end
