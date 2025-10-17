using LinearAlgebra: norm
using Test

@testset "Potential pairs" begin

    @testset "Pure gravity" begin

        masses = 1.0*ones(3)#u"Msun"
        radii = 1.0*ones(3)#u"Rsun"
        stellar_types = [14, 14]

        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

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
            dv = zeros(3, 3)
            i, j = p
            Syzygy.pure_gravitational_acceleration!(dv, rs, p, params)
            dvi = dv[:,i]
            dvj = dv[:,j]
            @test abs.(dvi) ≈ abs.(dvj)
        end
    end

    @testset "PN-1" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN1_acceleration!(dv, rs, vs, pair, params)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN2_acceleration!(dv, rs, vs, pair, params)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN-2.5" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)

        Syzygy.PN2p5_acceleration!(dv, rs, vs, pair, params)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

    @testset "PN1 to 2.5" begin

        masses = 10.0*ones(2)
        radii = zeros(2)#5.0*ones(2)u"Rsun"
        stellar_types = [14, 14]


        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dv = zeros(3, 2)

        r1 = [-1.0, 0.0, 0.0]#u"Rsun"
        r2 = [1.0, 0.0, 0.0]#u"Rsun"
        rs = [r1 r2]

        v1 = [0.0, -100.0, 0.0]#u"Rsun"
        v2 = [0.0, 100.0, 0.0]#u"Rsun"
        vs = [v1 v2]

        pair = (1, 2)
        Syzygy.PN1_to_2p5_acceleration!(dv, rs, vs, pair, params)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)


        @testset "PN1_to_2p5_acceleration! equals PN1+PN2+PN2.5" begin
            dvv = zeros(3, 2)

            Syzygy.PN1_acceleration!(dvv, rs, vs, pair, params)
            Syzygy.PN2_acceleration!(dvv, rs, vs, pair, params)
            Syzygy.PN2p5_acceleration!(dvv, rs, vs, pair, params)
            dvv1 = dvv[:,1]
            dvv2 = dvv[:,2]

            @test dv1 ≈ dvv1
            @test dv2 ≈ dvv2
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

        dv = zeros(3, 2)

        rs = reduce(hcat, binary.particles.position) .|> ustrip
        vs = reduce(hcat, binary.particles.velocity) .|> ustrip

        pot = EquilibriumTidalPotential(binary, set_spin=true)
        pair = (1, 2)
        Syzygy.Syzygy.equilibrium_tidal_acceleration!(dv, rs, vs, pair, params, pot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
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

        dv = zeros(3, 2)
        
        rs = reduce(hcat, binary.particles.position) .|> u"Rsun" .|> ustrip
        vs = reduce(hcat, (binary.particles.velocity)) .|> u"Rsun/yr" .|> ustrip

        pot = DynamicalTidalPotential(4, [1.5, 1.5])

        pair = (1, 2)
        Syzygy.Syzygy.dynamical_tidal_acceleration!(dv, rs, vs, pair, params, pot)
        dv1 = dv[:,1]
        dv2 = dv[:,2]
        @test abs.(dv1) ≈ abs.(dv2)

    end

end

