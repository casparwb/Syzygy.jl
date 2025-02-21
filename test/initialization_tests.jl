using LinearAlgebra: norm
@testset "Setup" begin

    # @testset "State vectors to elements" begin

    #     ain = 0.1u"AU"
    #     aout = 1.0u"AU"

    #     ein = 0.6
    #     eout = 0.2

    #     νin = (π/4)u"rad"
    #     νout = (2π/3)u"rad"

    #     ωout = (π/3)u"rad"
    #     ωin = (0.0)u"rad"

    #     Ωout = ωout
    #     Ωin  = ωin
        
    #     masses = [2.0, 1.0, 3.0]u"Msun"

    #     triple = multibodysystem(masses, a=[ain, aout], 
    #                                     e=[ein, eout],
    #                                     ν=[νin, νout],
    #                                     ω=[ωin, ωout],
    #                                     Ω=[Ωin, Ωout])

    #     Pin = triple.binaries[1].elements.P
    #     Pout = triple.binaries[2].elements.P


    #     rin = triple.particles.position[[1, 2]]
    #     vin = triple.particles.velocity[[1, 2]]
    #     min = masses[[1, 2]]

    #     r_com_in = Syzygy.centre_of_mass(rin, min)
    #     v_com_in = Syzygy.centre_of_mass(vin, min)

    #     rout = [r_com_in, triple.particles.position[3]]
    #     vout = [v_com_in, triple.particles.velocity[3]]
    #     mout = [masses[3], sum(min)]

    #     calculated_inner_binary_elements = Syzygy.binary_elements(rin, vin, min)
    #     calculated_outer_binary_elements = Syzygy.binary_elements(rout, vout, mout)

    #     @testset "$el" for el in propertynames(calculated_inner_binary_elements)
    #         @test getproperty(calculated_inner_binary_elements, el) ≈ getproperty(triple.binaries[1].elements, el)
    #         @test getproperty(calculated_outer_binary_elements, el) ≈ getproperty(triple.binaries[2].elements, el)
    #     end

    # end

    @testset "Triple initilization" begin

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

        Pin = triple.binaries[1].elements.P
        Pout = triple.binaries[2].elements.P


        rin = triple.particles.position[[1, 2]]
        vin = triple.particles.velocity[[1, 2]]
        min = masses[[1, 2]]

        r_com_in = Syzygy.centre_of_mass(rin, min)
        v_com_in = Syzygy.centre_of_mass(vin, min)

        rout = [r_com_in, triple.particles.position[3]]
        vout = [v_com_in, triple.particles.velocity[3]]
        mout = [masses[3], sum(min)]

        calculated_inner_binary_elements = Syzygy.binary_elements(rin, vin, min)
        calculated_outer_binary_elements = Syzygy.binary_elements(rout, vout, mout)

        @testset "$el" for el in propertynames(calculated_inner_binary_elements)
            @test getproperty(calculated_inner_binary_elements, el) ≈ getproperty(triple.binaries[1].elements, el)
            @test getproperty(calculated_outer_binary_elements, el) ≈ getproperty(triple.binaries[2].elements, el)
        end

        @testset "Mutual inclination" begin
            h1 = Syzygy.angular_momentum(rin[1] - rin[2], vin[1] - vin[2])
            h123 = Syzygy.angular_momentum(r_com_in - triple.particles.position[3], v_com_in - triple.particles.velocity[3])
            @test Syzygy.mutual_inclination(h1, h123) ≈ i[2]
        end

        νin = (0)u"rad"
        νout = (π)u"rad"
        
        masses = [2.0, 1.0, 3.0]u"Msun"

        triple = multibodysystem(masses, a=[ain, aout], 
                                        e=[ein, eout],
                                        ν=[νin, νout])

        rin = triple.particles.position[[1, 2]]
        min = masses[[1, 2]]
        r_com_in = Syzygy.centre_of_mass(rin, min)

        d12 = norm(triple.particles.position[2] - triple.particles.position[1])
        d312 = norm(r_com_in - triple.particles.position[3])

        @testset "Hierarchy" begin
            @test d12 < d312 
            @test triple.particles.mass ≈ masses
        end

        @testset "Separations" begin
            @test d12 ≈ ain*(1 - ein)
            @test d312 ≈ aout*(1 + eout)
        end

        @testset "Velocities" begin
            v_in = GRAVCONST*sum(min)*(2/d12 - 1/ain) |> sqrt
            v_out = GRAVCONST*sum(masses)*(2/d312 - 1/aout) |> sqrt


            @test v_in ≈ norm(vin[2] - vin[1])
            @test v_out ≈ norm(vout[2] - vout[1])
        end

        # @testset "Initialize with state vectors" begin

        # end

    end

    @testset "Arbitrary hierarchy" begin
    
        a = 1.0u"AU"
        # binary = multibodysystem(ones(2)u"Msun", a=a)
        # quint_2p2p1 = multibodysystem(ones(5)u"Msun", a=a .* [0.1, 1.0, 5.0, 10.0], hierarchy=[5, 2, 2, 1])
        # quint_3p2   = multibodysystem(ones(5)u"Msun", a=a .* [0.1, 1.0, 5.0, 10.0], hierarchy=[5, 1, 2, 1])
        
        
        @testset "Triple" begin
            triple = multibodysystem(ones(3)u"Msun", a=a .* [0.1, 1.0])
            p = triple.particles.position
            m = triple.particles.mass
    
            r12 = p[1] .- p[2]
            r13 = p[3] .- p[1]
            r23 = p[3] .- p[2]
    
            d12 = norm(r12) 
            d13 = norm(r13) 
            d23 = norm(r23) 
    
            com_12 = Syzygy.centre_of_mass(p[[1, 2]], m[[1, 2]])
    
            d123 = norm(p[3] .- com_12)
    
            @test d12 < d13
            @test d12 < d23
            @test d12 < d123
        end
    
        @testset "Quadruple" begin
            quad_3p1 = multibodysystem(ones(4)u"Msun", a=a .* [1.0, 1.0, 10.0], hierarchy=[4, 1, 1, 1])
            quad_2p2 = multibodysystem(ones(4)u"Msun", a=a .* [1.0, 1.0, 10.0], hierarchy=[4, 2, 1])
        
            let p = quad_3p1.particles.position, m = quad_3p1.particles.mass
                r12 = p[1] .- p[2]
                d12 = norm(r12) 
    
    
                com_12 = Syzygy.centre_of_mass(p[[1, 2]], m[[1, 2]])
                com_123 = Syzygy.centre_of_mass(p[[1, 2, 3]], m[[1, 2, 3]])
    
                d312 = norm(p[3] .- com_12)
                d4123 = norm(p[4] .- com_123)
    
                @test d12 < d312 < d4123
            end
    
            let p = quad_2p2.particles.position, m = quad_2p2.particles.mass
                r12 = p[1] .- p[2]
                r34 = p[3] .- p[4]
                d12 = norm(r12) 
                d34 = norm(r34)
    
                com_12 = Syzygy.centre_of_mass(p[[1, 2]], m[[1, 2]])
                com_34 = Syzygy.centre_of_mass(p[[3, 4]], m[[3, 4]])
    
                d = norm(com_12 .- com_34)
    
                @test d12 ≈ d34
                @test d12 < d
                @test d34 < d
            end
        end
    
        @testset "Quintuple" begin
            quint_4p1 = multibodysystem(ones(5)u"Msun", a=a .* [0.1, 1.0, 5.0, 10.0], hierarchy=[5, 1, 1, 1, 1])
            quint_2p2p1 = multibodysystem(ones(5)u"Msun", a=a .* [0.1, 1.0, 5.0, 10.0], hierarchy=[5, 2, 1, 1])
    
            let p = quint_4p1.particles.position, m = quint_4p1.particles.mass
                r12 = p[1] .- p[2]
                d12 = norm(r12) 
    
    
                com_12 = Syzygy.centre_of_mass(p[[1, 2]], m[[1, 2]])
                com_123 = Syzygy.centre_of_mass(p[[1, 2, 3]], m[[1, 2, 3]])
                com_1234 = Syzygy.centre_of_mass(p[[1, 2, 3, 4]], m[[1, 2, 3, 4]])
                d312 = norm(p[3] .- com_12)
                d4123 = norm(p[4] .- com_123)
                d51234 = norm(p[5] .- com_1234)
    
                @test d12 < d312 < d4123 < d51234
    
            end
    
            let p = quint_2p2p1.particles.position, m = quint_2p2p1.particles.mass
                r12 = p[1] .- p[2]
                r34 = p[3] .- p[4]
                d12 = norm(r12) 
                d34 = norm(r34)
    
                com_12 = Syzygy.centre_of_mass(p[[1, 2]], m[[1, 2]])
                com_34 = Syzygy.centre_of_mass(p[[3, 4]], m[[3, 4]])
                com_1234 = Syzygy.centre_of_mass(p[[1, 2, 3, 4]], m[[1, 2, 3, 4]])
                
                d_com_12_34 = norm(com_12 .- com_34)
                d5_com = norm(p[5] .- d_com_12_34)
    
    
                @test d12 ≈ d34
                @test d12 < d_com_12_34
                @test d34 < d_com_12_34
                @test d_com_12_34 < d5_com
            end
    
        end
    
    end

end