@testset "IC to ODE" begin

    ain = 0.1au
    aout = 1.0au

    ein = 0.6
    eout = 0.2

    i = [0.0, π/4]

    νin = (π/4)
    νout = (2π/3)

    ωout = (π/3)
    ωin = (0.0)

    Ωout = ωout
    Ωin  = ωin
    
    masses = [2.0, 1.0, 3.0]Msun

    triple = multibodysystem(masses, a=[ain, aout], 
                                     e=[ein, eout],
                                     ν=[νin, νout],
                                     ω=[ωin, ωout],
                                     Ω=[Ωin, Ωout],
                                     i=i, nbody_units=true)


    sim = simulation(triple, t_sim=10)
    bodies = sim.bodies

    unit_mass = triple.units.u_mass
    unit_length = triple.units.u_length
    unit_speed = unit_length/triple.units.u_time

    @testset "Particle $i to ODE" for i = 1:3
        @test bodies[i].mass ≈ ustrip(unit_mass, triple.particles[i].mass)
        @test bodies[i].position ≈ ustrip.(unit_length, triple.particles[i].position)
        @test bodies[i].velocity ≈ ustrip.(unit_speed, triple.particles[i].velocity)
        # @test bodies[i].spin ≈ triple.particles[i].structure.S .|> ustrip

    end

end
