@testset "Physics begin" begin

    @testset "Orbital elements" begin

    end

    @testset "Centre of mass" begin

    end

    @testset "Spin precession speed" begin

        
        masses = 10.0*ones(2)u"Msun"
        radii = Syzygy.schwarzschild_radius.(masses)
        luminosities = zeros(2)u"Lsun"
        stellar_types = [14, 14]
        ages = zeros(2)u"yr"
        
        sma = 1000*Syzygy.gravitational_radius(masses[1])
        e = 0.4

        m1, m2 = masses

        Χ = 0.5 # dimensionless spin parameter (1 = maximum spinning)
        dir = rand(3)
        dir = dir/norm(dir)
        S1 = GRAVCONST*Χ*dir*m1^2
        S2 = GRAVCONST*Χ*dir*m2^2

        binary = multibodysystem(masses, a=sma, e=e)

        p1, p2 = binary.particles[1], binary.particles[2]
        @test Syzygy.PN1_spin_precession_velocity(p1, p2) ≈ Syzygy.PN1_spin_precession_velocity(p2, p1)
        @test Syzygy.PN1p5_spin_precession_velocity(p1, p2) ≈ Syzygy.PN1p5_spin_precession_velocity(p2, p1)
        @test Syzygy.PN2_spin_precession_velocity(p1, p2) ≈ Syzygy.PN2_spin_precession_velocity(p2, p1)
    end

end