@testset "Solver" begin

    @testset "IC to ODE" begin

        ain = 0.1u"AU"
        aout = 1.0u"AU"

        ein = 0.6
        eout = 0.2

        i = [0.0, π/4]

        νin = (π/4)
        νout = (2π/3)

        ωout = (π/3)
        ωin = (0.0)

        Ωout = ωout
        Ωin  = ωin
        
        masses = [2.0, 1.0, 3.0]u"Msun"

        triple = multibodysystem(masses, a=[ain, aout], 
                                        e=[ein, eout],
                                        ν=[νin, νout],
                                        ω=[ωin, ωout],
                                        Ω=[Ωin, Ωout],
                                        i=i)

        unit_length, unit_mass, unit_time = triple.units.u_length, triple.units.u_mass, triple.units.u_time
        sim = simulation(triple, t_sim=10)
        bodies = sim.bodies

        @testset "Particle $i to ODE" for i = 1:3
            @test bodies[i].mass ≈ ustrip(unit_mass, triple.particles[i].mass)
            @test bodies[i].position ≈ ustrip.(unit_length, triple.particles[i].position)
            @test bodies[i].velocity ≈ ustrip.(unit_length/unit_time, triple.particles[i].velocity) 
            # @test bodies[i].spin ≈ triple.particles[i].structure.S .|> ustrip
        end

    end

    
    @testset "Callbacks" begin

        callbacks = [EscapeCB(), CollisionCB(), RocheLobeOverflowCB(), DemocraticCheckCB(), CentreOfMassCB(), IonizationCB()]

        
        @testset "Collision" begin
            
            a = 5.0u"Rsun"
            e = 0.6
            ν = (π/2)
            rp = a*(1 - e)
            R = [rp/2, rp/2] .* 1.01
            masses = [2.0, 1.0]u"Msun"
            binary = multibodysystem(masses, a=a, e=e, R = R)

            res = simulate(binary, t_sim=1.0, callbacks=[CollisionCB()])

            @test :Collision in keys(res.retcode)
            @test res.retcode[:Collision][1] == [1,2]
        end

        # @testset "New collision" begin
            
        #     a = 5.0u"Rsun"
        #     e = 0.6
        #     ν = (π/2)
        #     rp = a*(1 - e)
        #     R = [rp/2, rp/2] .* 1.01
        #     masses = [2.0, 1.0]u"Msun"
        #     binary = multibodysystem(masses, a=a, e=e, R = R)

        #     res = simulate(binary, t_sim=1.0, callbacks=[Syzygy.NewCollisionCB()])

        #     @test :Collision in keys(res.retcode)
        #     @test res.retcode[:Collision][1] == [1,2]
        # end

        @testset "Running with callbacks" begin
            a = [1.0, 5.0]u"AU"
            e = [0.6, 0.1]
            ν = (π/2)

            masses = [2.0, 0.5, 1.0]u"Msun"
            triple = multibodysystem(masses, a=a, e=e)

            @testset "$cb" for cb in callbacks
                res = simulate(triple, t_sim=1, save_everystep=false, callbacks=[cb])
            end

        end

    end

        
    @testset "Arbitrary precision" begin

        binary = multibodysystem(ones(2)u"Msun")

        res = simulate(binary, t_sim=1)

    end

end