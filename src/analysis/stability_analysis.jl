using StaticArrays
using LinearAlgebra: norm

function lyapunov_exponent(system; t_sim, epsilon=1e-6, simulation_args...)

    masses = system.particles.mass
    initial_positions = system.particles.position
    initial_velocities = system.particles.velocity

    sim_1 = simulate(system, t_sim=t_sim; simulation_args...)

    t_final = sim_1.solution.t[end]

    final_positions = eachcol(sim_1.solution.u[end].x[2])
    final_velocities = eachcol(sim_1.solution.u[end].x[1])

    initial_positions_perturbed = initial_positions .* (1 - epsilon)
    initial_velocities_perturbed = initial_velocities .* (1 - epsilon)
    system_perturbed = multibodysystem(masses, initial_positions_perturbed, initial_velocities_perturbed)

    sim_2 = simulate(system_perturbed, t_sim=t_sim; simulation_args...)
    sol_2_at_final_t = sim_2.solution(t_final)
    final_positions_perturbed = eachcol(sol_2_at_final_t.x[2])
    final_velocities_perturbed = eachcol(sol_2_at_final_t.x[1])

    initial_dr_vector = map(x -> ustrip.(Syzygy.unit_length, x),                  initial_positions .- initial_positions_perturbed)
    initial_dv_vector = map(x -> ustrip.(Syzygy.unit_length/Syzygy.unit_time, x), initial_velocities .- initial_velocities_perturbed)
    initial_separation_vector = reduce(vcat, [initial_dr_vector..., initial_dv_vector...])

    final_dr_vector = final_positions .- final_positions_perturbed
    final_dv_vector = final_velocities .- final_velocities_perturbed
    final_separation_vector = reduce(vcat, [final_dr_vector..., final_dv_vector...])


    λ = 1/t_final*log(norm(final_separation_vector)/norm(initial_separation_vector))
end