using LinearAlgebra: norm
using Test

@testset "Potentials" begin

    @testset "Pure gravity" begin

        masses = 1.0*ones(3)
        radii = zeros(3)#5.0*ones(2)u"Rsun"
        luminosities = zeros(3)
        stellar_types = [14, 14]
        M_cores = zeros(3)
        R_cores = zeros(3)

        ages = zeros(3)

        params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
                                                stellar_types, M_cores, R_cores, ages)

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

    
end