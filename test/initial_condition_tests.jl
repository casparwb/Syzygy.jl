@testset "IC to ODE" begin

    ain = 0.1u"AU"
    aout = 1.0u"AU"

    ein = 0.6
    eout = 0.2

    i = [0.0, π/4]u"rad"

    νin = (π/4)u"rad"
    νout = (2π/3)u"rad"

    ωout = (π/3)u"rad"
    ωin = (0.0)u"rad"

    Ωout = ωout
    Ωin  = ωin
    
    masses = [2.0, 1.0, 3.0]u"Msun"

    triple = multibodysystem(masses, a=[ain, aout], 
                                     e=[ein, eout],
                                     ν=[νin, νout],
                                     ω=[ωin, ωout],
                                     Ω=[Ωin, Ωout],
                                     i=i)


    sim = simulation(triple, t_sim=10)
    bodies = sim.bodies

    @testset "Particle $i to ODE" for i = 1:3
        @test bodies[i].mass ≈ triple.particles[i].mass |> ustrip
        @test bodies[i].position ≈ triple.particles[i].position .|> ustrip
        @test bodies[i].velocity ≈ triple.particles[i].velocity .|> ustrip
        @test bodies[i].spin ≈ triple.particles[i].structure.S .|> ustrip

    end

end
