using LinearAlgebra: norm
using Test

@testset "Potential values" begin

    @testset "Pure gravity" begin
        """
        Binary sysem with masses [2, 1]kg at a distance of 2m.

        """

        n = 2
        masses = [2.0, 1.0]
        radii = zeros(2)
        stellar_types = [14, 14]

        ages = zeros(2)

        params = Syzygy.DefaultSimulationParams(radii, masses, [Syzygy.stellar_types[s] for s in stellar_types], stellar_types)

        dvs = zeros(3, 2)

        r1 = [-1.0, 0.0, 0.0]
        r2 = [1.0, 0.0, 0.0]
        rs = [r1 r2]

        pot = Syzygy.PureGravitationalPotential(Syzygy.UNITLESS_G)
        Syzygy.pure_gravitational_acceleration!(dvs, rs, (1, 2), params, pot)
       
        expected_force = ustrip(3.337e-11u"N")

        calculated_force_1 = norm(dvs[:,1])*masses[1]
        calculated_force_2 = norm(dvs[:,2])*masses[2]

        ΔF1 = abs(calculated_force_1 - expected_force)/expected_force
        ΔF2 = abs(calculated_force_2 - expected_force)/expected_force
        @test ΔF1 < 1e-4
        @test ΔF2 < 1e-4
    end

    
end