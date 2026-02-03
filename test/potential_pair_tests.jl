using LinearAlgebra: norm
using Test
# const G = Syzygy.UNITLESS_G
# const c = Syzygy.UNITLESS_c
@testset "Potential pairs" begin
    G = Syzygy.UNITLESS_G
    c = Syzygy.UNITLESS_c
    gravpot = Syzygy.PureGravitationalPotential(G, 0.0)
    PN1pot = Syzygy.PN1Potential(G, c)
    PN2pot = Syzygy.PN2Potential(G, c)
    PN2p5pot = Syzygy.PN2p5Potential(G, c)
    PNpot = Syzygy.PNPotential(G, c)

    unit_mass, unit_length, unit_time = Syzygy.default_unit_mass, Syzygy.default_unit_length, Syzygy.default_unit_time
    unit_speed = unit_length/unit_time
    Rsun_yr = Rsun/u"yr"
    @testset "Pure gravity" begin
        # assuming solar units

        masses =ustrip(unit_mass, 1.0*ones(3)Msun)
        radii = ustrip(unit_length, 1.0*ones(3)Rsun)
        stellar_types = [1, 1]

        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)
        y = 2*sin(deg2rad(60))
        r1 = ustrip(unit_length, [-1.0, 0.0, 0.0]Rsun)
        r2 = ustrip(unit_length, [1.0, 0.0, 0.0]Rsun)
        r3 = ustrip(unit_length, [0.0, y, 0.0]Rsun)
        rs = [r1 r2 r3]

        pair12 = (1, 2)
        pair13 = (1, 3)
        pair23 = (2, 3)
        pairs = (pair12, pair13, pair23)
        for p in pairs
            dv = zeros(3, 3)
            i, j = p
            Syzygy.pure_gravitational_acceleration!(dv, rs, p, params, gravpot)
            dvi = dv[:,i]
            dvj = dv[:,j]
            @test abs.(dvi) ≈ abs.(dvj)
        end
    end

    @testset "PN-1" begin

        masses =ustrip(unit_mass, 1.0*ones(3)Msun)
        radii = ustrip(unit_length, 1.0*ones(3)Rsun)
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = ustrip(unit_length, [-1.0, 0.0, 0.0]Rsun)
        r2 = ustrip(unit_length, [1.0, 0.0, 0.0]Rsun)
        rs = [r1 r2]

        v1 = ustrip(unit_speed, [0.0, -100.0, 0.0]Rsun_yr)
        v2 = ustrip(unit_speed, [0.0, 100.0, 0.0]Rsun_yr)
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN1_acceleration!(dv, rs, vs, pair, params, PN1pot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2" begin

        masses = ustrip(unit_mass, 1.0*ones(3)Msun)
        radii = ustrip(unit_length, 1.0*ones(3)Rsun)
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = ustrip(unit_length, [-1.0, 0.0, 0.0]Rsun)
        r2 = ustrip(unit_length, [1.0, 0.0, 0.0]Rsun)
        rs = [r1 r2]

        v1 = ustrip(unit_speed, [0.0, -100.0, 0.0]Rsun_yr)
        v2 = ustrip(unit_speed, [0.0, 100.0, 0.0]Rsun_yr)
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN2_acceleration!(dv, rs, vs, pair, params, PN2pot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2.5" begin

        masses =ustrip(unit_mass, 1.0*ones(3)Msun)
        radii = ustrip(unit_length, 1.0*ones(3)Rsun)
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = ustrip(unit_length, [-1.0, 0.0, 0.0]Rsun)
        r2 = ustrip(unit_length, [1.0, 0.0, 0.0]Rsun)
        rs = [r1 r2]

        v1 = ustrip(unit_speed, [0.0, -100.0, 0.0]Rsun_yr)
        v2 = ustrip(unit_speed, [0.0, 100.0, 0.0]Rsun_yr)
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN2p5_acceleration!(dv, rs, vs, pair, params, PN2p5pot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN1 to 2.5" begin

        masses = ustrip(unit_mass, 1.0*ones(3)Msun)
        radii = ustrip(unit_length, 1.0*ones(3)Rsun)
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = ustrip(unit_length, [-1.0, 0.0, 0.0]Rsun)
        r2 = ustrip(unit_length, [1.0, 0.0, 0.0]Rsun)
        rs = [r1 r2]

        v1 = ustrip(unit_speed, [0.0, -100.0, 0.0]Rsun_yr)
        v2 = ustrip(unit_speed, [0.0, 100.0, 0.0]Rsun_yr)
        vs = [v1 v2]

        pair = (1, 2)
        Syzygy.PN1_to_2p5_acceleration!(dv, rs, vs, pair, params, PNpot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)


        @testset "PN1_to_2p5_acceleration! equals PN1+PN2+PN2.5" begin
            dvv = zeros(3, 2)

            Syzygy.PN1_acceleration!(dvv, rs, vs, pair, params, PN1pot)
            Syzygy.PN2_acceleration!(dvv, rs, vs, pair, params, PN2pot)
            Syzygy.PN2p5_acceleration!(dvv, rs, vs, pair, params, PN2p5pot)
            dvv1 = dvv[:,1]
            dvv2 = dvv[:,2]

            @test dv1 ≈ dvv1
            @test dv2 ≈ dvv2
        end
    end


    # @testset "Static equilibrium tides" begin

    #     masses = 1.0*ones(2)Msun
    #     radii = 1.0*ones(2)Rsun#5.0*ones(2)Rsun
    #     luminosities = 1.0*ones(2)u"Lsun"
    #     stellar_types = [1, 1]

    #     M_cores = 0.1 .* masses
    #     R_cores = 0.01 .* radii

    #     ages = zeros(2)u"yr"

    #     binary = multibodysystem(masses, a=10.0Rsun, R=radii, m_core = M_cores, R_core = R_cores)

    #     params = Syzygy.TidalSimulationParams(Float64.(ustrip.(radii)), Float64.(ustrip.(masses)), 
    #                                           Float64.(ustrip.(luminosities)), [Syzygy.stellar_types[s] for s in stellar_types], stellar_types, 
    #                                           Float64.(ustrip.(M_cores)), Float64.(ustrip.(R_cores)), Float64.(ustrip.(ages)))

    #     dv = zeros(3, 2)

    #     rs = reduce(hcat, binary.particles.position) .|> ustrip
    #     vs = reduce(hcat, binary.particles.velocity) .|> ustrip

    #     pot = EquilibriumTidalPotential(binary, set_spin=true)
    #     pair = (1, 2)
    #     Syzygy.Syzygy.equilibrium_tidal_acceleration!(dv, rs, vs, pair, params, pot)
    #     dv1 = dv[:,1]
    #     dv2 = dv[:,2]
    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "Dynamical tides" begin

    #     masses = 1.0*ones(2)Msun
    #     radii = [1.0, 1.0]Rsun
    #     luminosities = 1.0*ones(2)u"Lsun"
    #     stellar_types = [1, 1]

    #     M_cores = [0.1 * masses[1], masses[2]]
    #     R_cores = [0.01*radii[1], radii[2]]

    #     ages = zeros(2)u"yr"

    #     binary = multibodysystem(masses, a=5.0Rsun, e=0.4, R=radii, m_core = M_cores, R_core = R_cores)

    #     params = Syzygy.TidalSimulationParams(Float64.(ustrip.(radii)), 
    #                                             Float64.(ustrip.(masses)), 
    #                                             Float64.(ustrip.(luminosities)), 
    #                                             [Syzygy.stellar_types[s] for s in stellar_types], 
    #                                             stellar_types, 
    #                                             Float64.(ustrip.(M_cores)), 
    #                                             Float64.(ustrip.(R_cores)), 
    #                                             Float64.(ustrip.(ages)))

    #     dv = zeros(3, 2)
        
    #     rs = reduce(hcat, binary.particles.position) .|> Rsun .|> ustrip
    #     vs = reduce(hcat, (binary.particles.velocity)) .|> Rsun_yr .|> ustrip

    #     pot = DynamicalTidalPotential(4, [1.5, 1.5])

    #     pair = (1, 2)
    #     Syzygy.Syzygy.dynamical_tidal_acceleration!(dv, rs, vs, pair, params, pot)
    #     dv1 = dv[:,1]
    #     dv2 = dv[:,2]
    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

end

