@testset "Units" begin

    struct UnitfulSimulationParams{LengthVecType, MassVecType, WattVecType, IntVecType, TimeVecType} <: Syzygy.SimulationParams
        R::LengthVecType # radii
        M::MassVecType # masses
        L::WattVecType # luminosities
        stellar_types::IntVecType 
        M_cores::MassVecType # core masses
        R_cores::LengthVecType # core radii
        ages::TimeVecType 
    end

    @testset "Pure gravity" begin

        masses = ones(3)u"Msun"
        radii = ones(3)u"Rsun"
        luminosities = ones(3)u"Lsun"
        stellar_types = [14, 14]
        M_cores = zeros(3)u"Msun"
        R_cores = zeros(3)u"Rsun"

        ages = zeros(3)u"Myr"

        params = UnitfulSimulationParams(radii, masses, luminosities, 
                                                stellar_types, M_cores, R_cores, ages)

        dv1 = zeros(3)u"km/s^2"
        dv2 = zeros(3)u"km/s^2"
        dv3 = zeros(3)u"km/s^2"
        dvs = [dv1, dv2, dv3]

        y = 2*sin(deg2rad(60))
        r1 = [-1.0, 0.0, 0.0]u"Rsun"
        r2 = [1.0, 0.0, 0.0]u"Rsun"
        r3 = [0.0, y, 0.0]u"Rsun"
        rs = [r1 r2 r3]

        pair12 = (1, 2)
        pair13 = (1, 3)
        pair23 = (2, 3)
        pairs = (pair12, pair13, pair23)
        for p in pairs
            for dv in dvs
                fill!(dv, zero(dv1[1]))
            end
            i, j = p
            Syzygy.pure_gravitational_acceleration!(dvs[i], dvs[j], rs, p, params)
            @test abs.(dvs[i]) ≈ abs.(dvs[j])
        end
    end

    # @testset "PN-1" begin

    #     masses = 10.0*ones(2)
    #     radii = zeros(2)#5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)

    #     ages = zeros(2)

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [1.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     pair = (1, 2)

    #     Syzygy.PN1_acceleration!(dv1, dv2, rs, vs, pair, params)
    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "PN-2" begin

    #     masses = 10.0*ones(2)
    #     radii = zeros(2)#5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)

    #     ages = zeros(2)

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [1.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     pair = (1, 2)

    #     Syzygy.PN2_acceleration!(dv1, dv2, rs, vs, pair, params)
    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "PN-2.5" begin

    #     masses = 10.0*ones(2)
    #     radii = zeros(2)#5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)

    #     ages = zeros(2)

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [1.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     pair = (1, 2)
    #     Syzygy.PN2_5_acceleration!(dv1, dv2, rs, vs, pair, params)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "PN-1.5 spin contribution" begin

    #     masses = 10.0*ones(2)
    #     radii = Syzygy.gravitational_radius.(masses) .|> upreferred |> ustrip #5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)
    #     ages = zeros(2)

    #     m1, m2 = masses

    #     Χ = 0.5 # dimensionless spin paramter (1 = maximum spinning)
    #     dir = rand(3)
    #     dir = dir/norm(dir)
    #     S1 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip
    
    #     # dir = rand(3)
    #     # dir = dir/norm(dir)
    #     S2 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [0.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     dS1 = Syzygy.spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
    #     dS2 = Syzygy.spin_precession_velocity(S2, r2, r1, v2, v1, m2, m1)
        
    #     @test abs.(dS1) ≈ abs.(dS2)

    #     rs = cat(rs, [S1 S2], dims=1)
    #     vs = cat(vs, [dS1 dS2], dims=1)

    #     pair = (1, 2)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     # Syzygy.pure_gravitational_acceleration!(dv1, dv2, rs, pair, params)
    #     # @show dv1 dv2

    #     dvs = [dv1 dv2]
    #     Syzygy.PN1_5_spin_acceleration!(dv1, dv2, rs, vs, pair, params)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "PN-2 spin contribution" begin

    #     masses = 10.0*ones(2)
    #     radii = Syzygy.gravitational_radius.(masses) .|> upreferred |> ustrip #5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)
    #     ages = zeros(2)

    #     m1, m2 = masses

    #     Χ = 0.5 # dimensionless spin paramter (1 = maximum spinning)
    #     dir = rand(3)
    #     dir = dir/norm(dir)
    #     S1 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip
    
    #     # dir = rand(3)
    #     # dir = dir/norm(dir)
    #     S2 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [0.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     dS1 = Syzygy.spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
    #     dS2 = Syzygy.spin_precession_velocity(S2, r2, r1, v2, v1, m2, m1)
        
    #     @test abs.(dS1) ≈ abs.(dS2)

    #     rs = cat(rs, [S1 S2], dims=1)
    #     vs = cat(vs, [dS1 dS2], dims=1)

    #     pair = (1, 2)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     # Syzygy.pure_gravitational_acceleration!(dv1, dv2, rs, pair, params)
    #     # @show dv1 dv2

    #     dvs = [dv1 dv2]
    #     Syzygy.PN2_spin_acceleration!(dv1, dv2, rs, vs, pair, params)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "PN-2.5 spin contribution" begin

    #     masses = 10.0*ones(2)
    #     radii = Syzygy.gravitational_radius.(masses) .|> upreferred |> ustrip #5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)
    #     ages = zeros(2)

    #     m1, m2 = masses

    #     Χ = 0.5 # dimensionless spin paramter (1 = maximum spinning)
    #     dir = rand(3)
    #     dir = dir/norm(dir)
    #     S1 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip
    
    #     # dir = rand(3)
    #     # dir = dir/norm(dir)
    #     S2 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [0.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     dS1 = Syzygy.spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
    #     dS2 = Syzygy.spin_precession_velocity(S2, r2, r1, v2, v1, m2, m1)
        
    #     @test abs.(dS1) ≈ abs.(dS2)

    #     rs = cat(rs, [S1 S2], dims=1)
    #     vs = cat(vs, [dS1 dS2], dims=1)

    #     pair = (1, 2)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     # Syzygy.pure_gravitational_acceleration!(dv1, dv2, rs, pair, params)
    #     # @show dv1 dv2

    #     dvs = [dv1 dv2]
    #     Syzygy.PN2_5_spin_acceleration!(dv1, dv2, rs, vs, pair, params)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "Spin precession" begin

    #     masses = 10.0*ones(2)
    #     radii = Syzygy.gravitational_radius.(masses) .|> upreferred |> ustrip #5.0*ones(2)u"Rsun"
    #     luminosities = zeros(2)
    #     stellar_types = [14, 14]
    #     M_cores = zeros(2)
    #     R_cores = zeros(2)
    #     ages = zeros(2)

    #     m1, m2 = masses

    #     Χ = 0.5 # dimensionless spin paramter (1 = maximum spinning)
    #     dir = rand(3)
    #     dir = dir/norm(dir)
    #     S1 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip
    
    #     # dir = rand(3)
    #     # dir = dir/norm(dir)
    #     S2 = GRAVCONST*Χ*dir/Syzygy.c .|> upreferred .|> ustrip

    #     params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
    #                                             stellar_types, M_cores, R_cores, ages)

    #     r1 = [-1.0, 0.0, 0.0]#u"Rsun"
    #     r2 = [1.0, 0.0, 0.0]#u"Rsun"
    #     rs = [r1 r2]

    #     v1 = [0.0, -100.0, 0.0]#u"Rsun"
    #     v2 = [0.0, 100.0, 0.0]#u"Rsun"
    #     vs = [v1 v2]

    #     dS1 = Syzygy.spin_precession_velocity(S1, r1, r2, v1, v2, m1, m2)
    #     dS2 = Syzygy.spin_precession_velocity(S2, r2, r1, v2, v1, m2, m1)
        
    #     @test abs.(dS1) ≈ abs.(dS2)

    #     rs = cat(rs, [S1 S2], dims=1)
    #     vs = cat(vs, [dS1 dS2], dims=1)

    #     pair = (1, 2)

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     # Syzygy.pure_gravitational_acceleration!(dv1, dv2, rs, pair, params)
    #     # @show dv1 dv2

    #     dvs = [dv1 dv2]
    #     Syzygy.spin_precession!(dv1, dv2, dvs, rs, vs, pair, params)

    #     @test dv1 ≈ dv2

    # end

    # @testset "Static equilibrium tides" begin

    #     masses = 1.0*ones(2)u"Msun"
    #     radii = 1.0*ones(2)u"Rsun"#5.0*ones(2)u"Rsun"
    #     luminosities = 1.0*ones(2)u"Lsun"
    #     stellar_types = [1, 1]

    #     M_cores = 0.1 .* masses
    #     R_cores = 0.01 .* radii

    #     ages = zeros(2)u"yr"

    #     binary = multibodysystem(masses, a=10.0u"Rsun", R=radii, m_core = M_cores, R_core = R_cores)

    #     params = Syzygy.DefaultSimulationParams(Float64.(ustrip(radii)), Float64.(ustrip(masses)), 
    #                                             Float64.(ustrip(luminosities)), stellar_types, 
    #                                             Float64.(ustrip(M_cores)), Float64.(ustrip(R_cores)), Float64.(ustrip(ages)))

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)

    #     rs = reduce(hcat, binary.particles.position) .|> ustrip
    #     vs = reduce(hcat, binary.particles.velocity) .|> ustrip

    #     Ss = Syzygy.stellar_spin.(masses, radii)
    #     Ss = reduce(hcat, [[upreferred(S).val, 0.0, 0.0] for S in Ss])
    #     vSs = reduce(hcat, [zeros(3) for i = 1:2])

    #     rs = [rs; Ss]
    #     vs = [vs; vSs]

    #     pot = Syzygy.TimeDependentEquilibriumTidalPotential()
    #     pair = (1, 2)
    #     Syzygy.Syzygy.equilibrium_tidal_acceleration!(dv1, dv2, rs, vs, pair, params, pot)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "Static equilibrium tides" begin

    #     masses = 1.0*ones(2)u"Msun"
    #     radii = 1.0*ones(2)u"Rsun"#5.0*ones(2)u"Rsun"
    #     luminosities = 1.0*ones(2)u"Lsun"
    #     stellar_types = [1, 1]

    #     M_cores = 0.1 .* masses
    #     R_cores = 0.01 .* radii

    #     ages = zeros(2)u"yr"

    #     binary = multibodysystem(masses, a=10.0u"Rsun", R=radii, m_core = M_cores, R_core = R_cores)

    #     params = Syzygy.DefaultSimulationParams(Float64.(ustrip(radii)), Float64.(ustrip(masses)), 
    #                                             Float64.(ustrip(luminosities)), stellar_types, 
    #                                             Float64.(ustrip(M_cores)), Float64.(ustrip(R_cores)), Float64.(ustrip(ages)))

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)


    #     rs = reduce(hcat, binary.particles.position) .|> ustrip
    #     vs = reduce(hcat, binary.particles.velocity) .|> ustrip

    #     Ss = Syzygy.stellar_spin.(masses, radii)
    #     Ss = reduce(hcat, [[upreferred(S).val, 0.0, 0.0] for S in Ss])
    #     vSs = reduce(hcat, [zeros(3) for i = 1:2])

    #     rs = [rs; Ss]
    #     vs = [vs; vSs]

    #     pot = Syzygy.EquilibriumTidalPotential(binary)
    #     pair = (1, 2)
    #     Syzygy.Syzygy.equilibrium_tidal_acceleration!(dv1, dv2, rs, vs, pair, params, pot)

    #     @test abs.(dv1) ≈ abs.(dv2)

    # end

    # @testset "Dynamical tides" begin

    #     masses = 1.0*ones(2)u"Msun"
    #     radii = [1.0, 1.0]u"Rsun"
    #     luminosities = 1.0*ones(2)u"Lsun"
    #     stellar_types = [1, 1]

    #     M_cores = [0.1 * masses[1], masses[2]]
    #     R_cores = [0.01*radii[1], radii[2]]

    #     ages = zeros(2)u"yr"

    #     binary = multibodysystem(masses, a=5.0u"Rsun", e=0.4, R=radii, m_core = M_cores, R_core = R_cores)

    #     params = Syzygy.DefaultSimulationParams(Float64.(ustrip(radii)), 
    #                                             Float64.(ustrip(masses)), 
    #                                             Float64.(ustrip(luminosities)), 
    #                                             stellar_types, 
    #                                             Float64.(ustrip(M_cores)), 
    #                                             Float64.(ustrip(R_cores)), 
    #                                             Float64.(ustrip(ages)))

    #     dv1 = zeros(3)
    #     dv2 = zeros(3)
        
    #     rs = reduce(hcat, binary.particles.position) .|> u"Rsun" .|> ustrip
    #     vs = reduce(hcat, (binary.particles.velocity)) .|> u"Rsun/yr" .|> ustrip

    #     pot = DynamicalTidalPotential(n=4, γ=[1.5, 1.5])

    #     pair = (1, 2)
    #     Syzygy.Syzygy.dynamical_tidal_acceleration!(dv1, dv2, rs, vs, pair, params, pot)
    #     @test abs.(dv1) ≈ abs.(dv2)

    # end
end