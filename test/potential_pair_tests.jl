using LinearAlgebra: norm
using Test

@testset "Potential pairs" begin

    @testset "Pure gravity" begin

        masses = 1.0*ones(3)#u"Msun"
        radii = 1.0*ones(3)#u"Rsun"
        stellar_types = [14, 14]

        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv1 = zeros(3)
        dv2 = zeros(3)
        dv3 = zeros(3)
        dvs = [dv1, dv2, dv3]

        y = 2*sin(deg2rad(60))
        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        r3 = [0.0, y, 0.0]#u"Rsun"
        rs = [r1 r2 r3]

        pair12 = (1, 2)
        pair13 = (1, 3)
        pair23 = (2, 3)
        pairs = (pair12, pair13, pair23)
        for p in pairs
            for dv in dvs
                fill!(dv, 0.0)
            end
            i, j = p
            Syzygy.pure_gravitational_acceleration!(dvs[i], dvs[j], rs, p, params)
            @test abs.(dvs[i]) ≈ abs.(dvs[j])
        end
    end

    @testset "PN-1" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv1 = zeros(3)
        dv2 = zeros(3)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN1_acceleration!(dv1, dv2, rs, vs, pair, params)
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv1 = zeros(3)
        dv2 = zeros(3)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN2_acceleration!(dv1, dv2, rs, vs, pair, params)
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2.5" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv1 = zeros(3)
        dv2 = zeros(3)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)
        Syzygy.PN2p5_acceleration!(dv1, dv2, rs, vs, pair, params)

        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN1 to 2.5" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv1 = zeros(3)
        dv2 = zeros(3)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)
        Syzygy.PN1_to_2p5_acceleration!(dv1, dv2, rs, vs, pair, params)

        @test abs.(dv1) ≈ abs.(dv2)


        @testset "PN1_to_2p5_acceleration! equals PN1+PN2+PN2.5" begin
            dv1_2 = zeros(3)
            dv2_2 = zeros(3)

            Syzygy.PN1_acceleration!(dv1_2, dv2_2, rs, vs, pair, params)
            Syzygy.PN2_acceleration!(dv1_2, dv2_2, rs, vs, pair, params)
            Syzygy.PN2p5_acceleration!(dv1_2, dv2_2, rs, vs, pair, params)


            @test dv1 ≈ dv1_2
            @test dv2 ≈ dv2_2
        end
    end


    @testset "Static equilibrium tides" begin

        masses = 1.0*ones(2)u"Msun"
        radii = 1.0*ones(2)u"Rsun"#5.0*ones(2)u"Rsun"
        luminosities = 1.0*ones(2)u"Lsun"
        stellar_types = [1, 1]

        M_cores = 0.1 .* masses
        R_cores = 0.01 .* radii

        ages = zeros(2)u"yr"

        binary = multibodysystem(masses, a=10.0u"Rsun", R=radii, m_core = M_cores, R_core = R_cores)

        params = Syzygy.TidalSimulationParams(Float64.(ustrip.(radii)), Float64.(ustrip.(masses)), 
                                              Float64.(ustrip.(luminosities)), [Syzygy.stellar_types[s] for s in stellar_types], stellar_types, 
                                              Float64.(ustrip.(M_cores)), Float64.(ustrip.(R_cores)), Float64.(ustrip.(ages)))

        dv1 = zeros(3)
        dv2 = zeros(3)

        rs = reduce(hcat, binary.particles.position) .|> ustrip
        vs = reduce(hcat, binary.particles.velocity) .|> ustrip

        pot = EquilibriumTidalPotential(binary, set_spin=true)
        pair = (1, 2)
        Syzygy.Syzygy.equilibrium_tidal_acceleration!(dv1, dv2, rs, vs, pair, params, pot)

        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "Dynamical tides" begin

        masses = 1.0*ones(2)u"Msun"
        radii = [1.0, 1.0]u"Rsun"
        luminosities = 1.0*ones(2)u"Lsun"
        stellar_types = [1, 1]

        M_cores = [0.1 * masses[1], masses[2]]
        R_cores = [0.01*radii[1], radii[2]]

        ages = zeros(2)u"yr"

        binary = multibodysystem(masses, a=5.0u"Rsun", e=0.4, R=radii, m_core = M_cores, R_core = R_cores)

        params = Syzygy.TidalSimulationParams(Float64.(ustrip.(radii)), 
                                                Float64.(ustrip.(masses)), 
                                                Float64.(ustrip.(luminosities)), 
                                                [Syzygy.stellar_types[s] for s in stellar_types], 
                                                stellar_types, 
                                                Float64.(ustrip.(M_cores)), 
                                                Float64.(ustrip.(R_cores)), 
                                                Float64.(ustrip.(ages)))

        dv1 = zeros(3)
        dv2 = zeros(3)
        
        rs = reduce(hcat, binary.particles.position) .|> u"Rsun" .|> ustrip
        vs = reduce(hcat, (binary.particles.velocity)) .|> u"Rsun/yr" .|> ustrip

        pot = DynamicalTidalPotential(4, [1.5, 1.5])

        pair = (1, 2)
        Syzygy.Syzygy.dynamical_tidal_acceleration!(dv1, dv2, rs, vs, pair, params, pot)
        @test abs.(dv1) ≈ abs.(dv2)

    end

end

