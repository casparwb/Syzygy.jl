using LinearAlgebra: norm
using Test

@testset "Potential values" begin

    @testset "Pure gravity" begin
        """
        Binary sysem with masses [2, 1]M⊙ at a distance of 2 R⊙.

        """

        n = 2
        masses = [2.0, 1.0]
        radii = zeros(2)
        luminosities = zeros(2)
        stellar_types = [14, 14]
        core_masses = zeros(2)
        core_radii = zeros(2)

        ages = zeros(2)

        params = Syzygy.DefaultSimulationParams(radii, masses, luminosities, 
                                                stellar_types, core_masses, core_radii, ages)

        dv1 = zeros(3)
        dv2 = zeros(3)

        r1 = [-1.0, 0.0, 0.0]
        r2 = [1.0, 0.0, 0.0]
        rs = [r1 r2]

        Syzygy.pure_gravitational_acceleration!(dv1, dv2, rs, (1, 2), params)
       
        expected_force = ustrip(u"Msun*Rsun/yr^2", 2.726e+32u"N")#2.72772941815376324831096063956728e32u"N")
        # expected_force_1 = expected_force#/(masses[1]u"Msun") |> u"Msun*Rsun/yr^2" |> ustrip
        # expected_force_1 = expected_force#/(masses[2]u"Msun") |> u"Msun*Rsun/yr^2" |> ustrip
        @show dv1 dv2
        calculated_force_1 = norm(dv1)*masses[1]
        calculated_force_2 = norm(dv2)*masses[2]

        @test expected_force ≈ calculated_force_1
        @test expected_force ≈ calculated_force_2
    end

    
end