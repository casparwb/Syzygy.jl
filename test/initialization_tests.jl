
@testset "State vectors to elements" begin

    ain = 0.1u"AU"
    aout = 1.0u"AU"

    ein = 0.6
    eout = 0.2

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
                                     Ω=[Ωin, Ωout])

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

end

@testset "Initialization" begin

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
        @test triple.partices.mass ≈ masses
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

    @testset "Initialize with state vectors" begin

    end

end

