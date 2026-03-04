##########################################################

struct PN1p5SpinPotential <: SpinPotential end

struct PN2SpinPotential <: SpinPotential end

struct PN2p5SpinPotential <: SpinPotential end

struct PN1SpinPrecessionPotential <: SpinPotential end

struct PN1p5SpinPrecessionPotential <: SpinPotential end

struct PN2SpinPrecessionPotential <: SpinPotential end

struct SpinPrecessionPotential <: SpinPotential end

function PN1_spin_precession!(
        dvi,
        dvj,
        dvs,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )
    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    ā₁ = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]


    m₁ = params.M[i]
    m₂ = params.M[j]

    # add @fastmath?

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r² = r * r
    r⁻¹ = 1 / r
    r⁻² = 1 / r²

    n = r̄ / r
    nv = dot(n, v̄)

    dr_dt = nv
    dn_dt = (r * v̄ - r̄ * dr_dt) * r⁻²

    Gm₁ = UNITLESS_G * m₁
    Gm₂ = UNITLESS_G * m₂

    dT1PN_dt₁ = let

        nS₁ = dot(n, S̄₁)
        vS₁ = dot(v̄, S̄₁)

        dnS₁_dt = dot(n, dS̄₁) + dot(S̄₁, dn_dt)
        dvS₁_dt = dot(v̄, dS̄₁) + dot(S̄₁, ā)
        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)

        Gm₂r⁻² = Gm₂ * r⁻²
        dGm₂r⁻²_dt = -2 * Gm₂ * dr_dt * r⁻¹ * r⁻²
        fac = (v̄₁ - 2 * v̄₂) * nS₁ + S̄₁ * nv - 2 * n * vS₁
        dfac_dt = (v̄₁ - 2 * v̄₂) * dnS₁_dt + (ā₁ - 2 * ā₂) * nS₁ +
            S̄₁ * dnv_dt + nv * dS̄₁ -
            2 * n * dvS₁_dt - 2 * vS₁ * dn_dt
        fac * dGm₂r⁻²_dt + Gm₂r⁻² * dfac_dt
    end

    dT1PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā,

            dn_dt = -dn_dt
        nS₂ = dot(n, S̄₂)
        vS₂ = dot(v̄, S̄₂)

        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)
        dnS₂_dt = dot(n, dS̄₂) + dot(S̄₂, dn_dt)
        dvS₂_dt = dot(v̄, dS̄₂) + dot(S̄₂, ā)


        Gm₁r⁻² = Gm₁ * r⁻²
        dGm₁r⁻²_dt = -2 * Gm₁ * dr_dt * r⁻¹ * r⁻²
        fac = (v̄₂ - 2 * v̄₁) * nS₂ + S̄₂ * nv - 2 * n * vS₂
        dfac_dt = (v̄₂ - 2 * v̄₁) * dnS₂_dt + (ā₂ - 2 * ā₁) * nS₂ +
            S̄₂ * dnv_dt + nv * dS̄₂ -
            2 * n * dvS₂_dt - 2 * vS₂ * dn_dt
        fac * dGm₁r⁻²_dt + Gm₁r⁻² * dfac_dt
    end

    dvi .+= dT1PN_dt₁ * c⁻²
    dvj .+= dT1PN_dt₂ * c⁻²
    return nothing
end

function PN1p5_spin_precession!(
        dvi,
        dvj,
        dvs,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )
    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1 / r
    r⁻² = r⁻¹ * r⁻¹
    r⁻³ = r⁻² * r⁻¹

    n = r̄ / r
    nv = dot(n, v̄)

    dr_dt = nv
    dn_dt = (r * v̄ - r̄ * dr_dt) * r⁻²

    Gr⁻³ = -UNITLESS_G * r⁻³
    dGr⁻³_dt = 3 * UNITLESS_G * dr_dt * r⁻² * r⁻²

    dT1p5PN_dt₁ = let
        nS₂ = dot(n, S̄₂)
        dnS₂_dt = dot(n, dS̄₂) + dot(S̄₂, dn_dt)

        # F1p5PN = S̄₂ - 3nS₂*n
        # dF1p5PN_dt = dS̄₂ - 3*n*dnS₂_dt - 3nS₂*dn_dt
        # dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₁) + (F1p5PN × dS̄₁))

        F = (3nS₂ * n × S̄₁) - (S̄₂ × S̄₁)
        dF_dt = 3 * ((nS₂ * dn_dt + n * dnS₂_dt) × S̄₁ + nS₂ * n × dS̄₁) - (S̄₁ × dS̄₂) + (S̄₂ × dS̄₁)

        Gr⁻³ * dF_dt + F * dGr⁻³_dt
    end

    dT1p5PN_dt₂ = let n = -n, dn_dt = -dn_dt
        nS₁ = dot(n, S̄₁)
        dnS₁_dt = dot(n, dS̄₁) + dot(S̄₁, dn_dt)

        # F1p5PN = S̄₁ - 3nS₁*n
        # dF1p5PN_dt = dS̄₁ - 3*n*dnS₁_dt - 3nS₁*dn_dt
        # dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₂) + (F1p5PN × dS̄₂))
        # ₁ ₂ ₃

        F = (3nS₁ * n × S̄₂) - (S̄₁ × S̄₂)
        dF_dt = 3 * ((nS₁ * dn_dt + n * dnS₁_dt) × S̄₂ + nS₁ * n × dS̄₂) - (S̄₂ × dS̄₁) + (S̄₁ × dS̄₂)

        Gr⁻³ * dF_dt + F * dGr⁻³_dt
    end

    dvi .+= dT1p5PN_dt₁ * c⁻³
    dvj .+= dT1p5PN_dt₂ * c⁻³
    return nothing
end

function PN2_spin_precession!(
        dvi,
        dvj,
        dvs,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )

    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    ā₁ = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]

    m₁ = params.M[i]
    m₂ = params.M[j]
    δm = m₁ - m₂

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂
    r = norm(r̄) # r₁₂

    r² = r * r
    r⁻¹ = 1 / r
    r⁻² = 1 / r²

    n = r̄ / r
    nv = dot(n, v̄)
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    nv₁² = nv₁^2
    nv₂² = nv₂^2

    dr_dt = nv
    dn_dt = (r * v̄ - r̄ * dr_dt) * r⁻²

    Gm₁ = UNITLESS_G * m₁
    Gm₂ = UNITLESS_G * m₂

    Gm₁r⁻² = Gm₁ * r⁻²
    dGm₁r⁻²_dt = -2 * Gm₁ * dr_dt * r⁻¹ * r⁻²

    Gm₂r⁻² = Gm₂ * r⁻²
    dGm₂r⁻²_dt = -2 * Gm₂ * dr_dt * r⁻¹ * r⁻²

    dT2PN_dt₁ = let

        vv₂ = dot(v̄, v̄₂)
        nS₁ = dot(n, S̄₁)
        vS₁ = dot(v̄, S̄₁)
        v₁S₁ = dot(v̄₁, S̄₁)
        v₂S₁ = dot(v̄₂, S̄₁)


        dnS₁_dt = dot(n, dS̄₁) + dot(S̄₁, dn_dt)
        dvS₁_dt = dot(v̄, dS̄₁) + dot(S̄₁, ā)
        dv₁S₁_dt = dot(v̄₁, dS̄₁) + dot(S̄₁, ā₁)
        dv₂S₁_dt = dot(v̄₂, dS̄₁) + dot(S̄₁, ā₂)
        dvv₂_dt = dot(v̄, ā₂) + dot(v̄₂, ā)
        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)
        dnv₁_dt = dot(n, ā₁) + dot(v̄₁, dn_dt)
        dnv₂_dt = dot(n, ā₂) + dot(v̄₂, dn_dt)


        # T2PN = m₂*r⁻²*(S̄₁*(nv₂*vv₂ - 3/2*nv₂²*nv + Gm₁*r⁻¹*nv₁ - Gm₂*r⁻¹*nv) +
        #                 n*(vS₁*(3*nv₂² + 2*vv₂) + Gm₁*r⁻¹*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁) +
        #                     2nS₁*Gm₂*r⁻¹*nv) - v̄₁*(3/2*nS₁*nv₂² + vS₁*nv₂ -
        #                     nS₁*UNITLESS_G*r⁻¹*(6m₁ - m₂)) + v̄₂*(nS₁*(2vv₂ + 3nv₂²) +
        #                     2nv*(v₁S₁ + v₂S₁) - 5nS₁*UNITLESS_G*r⁻¹*δm)
        #                 )

        # dF2PN_dt = dn_dt*(Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*r⁻¹ + 2*Gm₂*nS₁*nv*r⁻¹ + (3*nv₂² + 2*vv₂)*vS₁) +
        #            ā₂*(-5*UNITLESS_G*δm*nS₁*r⁻¹ + (3*nv₂² + 2*vv₂)*nS₁ + 2*(v₁S₁ + v₂S₁)*nv) -
        #            (-UNITLESS_G*(6*m₁ - m₂)*nS₁*r⁻¹ + 3*nS₁*nv₂²/2 + v̄₂*vS₁)*ā₁ +
        #            (Gm₁*nv₁*r⁻¹ - Gm₂*nv*r⁻¹ - 3*nv*nv₂²/2 + nv₂*vv₂)*dS₁_dt +
        #            (5*UNITLESS_G*δm*nS₁*dr_dt*r⁻² - 5*UNITLESS_G*δm*dnS₁_dt*r⁻¹ +
        #             (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*nS₁ + (3*nv₂² + 2*vv₂)*dnS₁_dt +
        #             2*(v₁S₁ + v₂S₁)*dnv_dt + 2*(dv₁S₁_dt + dv₂S₁_dt)*nv)*v̄₂ -
        #            (UNITLESS_G*(6*m₁ - m₂)*nS₁*dr_dt*r⁻² - UNITLESS_G*(6*m₁ - m₂)*dnS₁_dt*r⁻¹ +
        #             3*nS₁*nv₂*dnv₂_dt + 3*nv₂²*dnS₁_dt/2 + nv₂*dvS₁_dt + vS₁*dnv₂_dt)*v̄₁ +
        #             (-Gm₁*(-16*nS₁*nv + 3*v₁S₁ - 7*v₂S₁)*dr_dt*r⁻² + Gm₁*(-16*nS₁*dnv_dt -
        #              16*nv*dnS₁_dt + 3*dv₁S₁_dt - 7*dv₂S₁_dt)*r⁻¹ - 2*Gm₂*nS₁*nv*dr_dt*r⁻² +
        #              2*Gm₂*nS₁*dnv_dt*r⁻¹ + 2*Gm₂*nv*dnS₁_dt*r⁻¹ + (6*nv₂*dnv₂_dt + 2*dvv₂_dt)*vS₁ +
        #             (3*nv₂² + 2*vv₂)*dvS₁_dt)*n +
        #             (-Gm₁*nv₁*dr_dt*r⁻² + Gm₁*dnv₁_dt*r⁻¹ + Gm₂*nv*dr_dt*r⁻² -
        #              Gm₂*dnv_dt*r⁻¹ - 3*nv*nv₂*dnv₂_dt - 3*nv₂²*dnv_dt/2 + nv₂*dvv₂_dt + vv₂*dnv₂_dt)*S₁

        num = (Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * r⁻¹ + 2 * Gm₂ * nS₁ * nv * r⁻¹ + (3 * nv₂² + 2 * vv₂) * vS₁) * n
        num += (-5 * UNITLESS_G * δm * nS₁ * r⁻¹ + (3 * nv₂² + 2 * vv₂) * nS₁ + 2 * (v₁S₁ + v₂S₁) * nv) * v̄₂
        num += -(-UNITLESS_G * (6 * δm) * nS₁ * r⁻¹ + 3 * nS₁ * nv₂² / 2 + nv₂ * vS₁) * v̄₁
        num += (Gm₁ * nv₁ * r⁻¹ - Gm₂ * nv * r⁻¹ - 3 * nv * nv₂² / 2 + nv₂ * vv₂) * S̄₁
        num *= -2 * Gm₂ * dr_dt * r⁻² * r⁻¹
        num += (Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * r⁻¹ + 2 * Gm₂ * nS₁ * nv * r⁻¹ + (3 * nv₂² + 2 * vv₂) * vS₁) * dn_dt
        num += (-5 * UNITLESS_G * δm * nS₁ * r⁻¹ + (3 * nv₂² + 2 * vv₂) * nS₁ + 2 * (v₁S₁ + v₂S₁) * nv) * ā₂
        num += -(-UNITLESS_G * (6 * δm) * nS₁ * r⁻¹ + 3 * nS₁ * nv₂² / 2 + nv₂ * vS₁) * ā₁
        num += (Gm₁ * nv₁ * r⁻¹ - Gm₂ * nv * r⁻¹ - 3 * nv * nv₂² / 2 + nv₂ * vv₂) * dS̄₁
        num += (
            5 * UNITLESS_G * δm * nS₁ * dr_dt * r⁻² - 5 * UNITLESS_G * δm * dnS₁_dt * r⁻¹ + (6 * nv₂ * dnv₂_dt + 2 * dvv₂_dt) * nS₁ +
                (3 * nv₂² + 2 * vv₂) * dnS₁_dt + 2 * (v₁S₁ + v₂S₁) * dnv_dt + 2 * (dv₁S₁_dt + dv₂S₁_dt) * nv
        ) * v̄₂
        num += -(
            UNITLESS_G * (6 * δm) * nS₁ * dr_dt * r⁻² - UNITLESS_G * (6 * δm) * dnS₁_dt * r⁻¹ +
                3 * nS₁ * nv₂ * dnv₂_dt + 3 * nv₂² * dnS₁_dt / 2 + nv₂ * dvS₁_dt + vS₁ * dnv₂_dt
        ) * v̄₁
        num += (
            -Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * dr_dt * r⁻² +
                Gm₁ * (
                -16 * nS₁ * dnv_dt - 16 * nv * dnS₁_dt +
                    3 * dv₁S₁_dt - 7 * dv₂S₁_dt
            ) * r⁻¹ -
                2 * Gm₂ * nS₁ * nv * dr_dt * r⁻² + 2 * Gm₂ * nS₁ * dnv_dt * r⁻¹ +
                2 * Gm₂ * nv * dnS₁_dt * r⁻¹ +
                (6 * nv₂ * dnv₂_dt + 2 * dvv₂_dt) * vS₁ +
                (3 * nv₂² + 2 * vv₂) * dvS₁_dt
        ) * n
        num += (
            -Gm₁ * nv₁ * dr_dt * r⁻² + Gm₁ * dnv₁_dt * r⁻¹ +
                Gm₂ * nv * dr_dt * r⁻² - Gm₂ * dnv_dt * r⁻¹ -
                3 * nv * nv₂ * dnv₂_dt - 3 * nv₂² * dnv_dt / 2 +
                nv₂ * dvv₂_dt + vv₂ * dnv₂_dt
        ) * S̄₁
        num *= Gm₂ * r⁻²

        # Gm₁r⁻²*dF2PN_dt + dGm₁r⁻²_dt*T2PN
        num * Gm₂ * r⁻²
    end

    dT2PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā, δm = -δm

        dn_dt = -dn_dt
        nv₁ = -nv₁
        nv₂ = -nv₂
        vv₁ = dot(v̄, v̄₁)
        nS₂ = dot(n, S̄₂)
        vS₂ = dot(v̄, S̄₂)
        v₂S₂ = dot(v̄₂, S̄₂)
        v₁S₂ = dot(v̄₁, S̄₂)

        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)
        dnv₁_dt = dot(n, ā₁) + dot(v̄₁, dn_dt)
        dnv₂_dt = dot(n, ā₂) + dot(v̄₂, dn_dt)
        dnS₂_dt = dot(n, dS̄₂) + dot(S̄₂, dn_dt)
        dvS₂_dt = dot(v̄, dS̄₂) + dot(S̄₂, ā)
        dv₂S₂_dt = dot(v̄₂, dS̄₂) + dot(S̄₂, ā₂)
        dv₁S₂_dt = dot(v̄₁, dS̄₂) + dot(S̄₂, ā₁)
        dvv₁_dt = dot(v̄, ā₁) + dot(v̄₁, ā)

        # T2PN = m₁*r⁻²*(S̄₂*(nv₁*vv₁ - 3/2*nv₁²*nv + Gm₂*r⁻¹*nv₂ - Gm₁*r⁻¹*nv) +
        #                 n*(vS₂*(3*nv₁² + 2*vv₁) + Gm₂*r⁻¹*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂) +
        #                     2nS₂*Gm₁*r⁻¹*nv) - v̄₂*(3/2*nS₂*nv₁² + vS₂*nv₁ -
        #                     nS₂*UNITLESS_G*r⁻¹*(6m₂ - m₁)) + v̄₁*(nS₂*(2vv₁ + 3nv₁²) +
        #                     2nv*(v₂S₂ + v₁S₂) - 5nS₂*UNITLESS_G*r⁻¹*δm)
        #                 )

        # dF2PN_dt = dn_dt*(Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*r⁻¹ + 2*Gm₁*nS₂*nv*r⁻¹ + (3*nv₁² + 2*vv₁)*vS₂) +
        #            ā₁*(-5*UNITLESS_G*δm*nS₂*r⁻¹ + (3*nv₁² + 2*vv₁)*nS₂ + 2*(v₂S₂ + v₁S₂)*nv) -
        #            (-UNITLESS_G*(6*m₂ - m₁)*nS₂*r⁻¹ + 3*nS₂*nv₁²/2 + v̄₁*vS₂)*ā₂ +
        #            (Gm₂*nv₂*r⁻¹ - Gm₁*nv*r⁻¹ - 3*nv*nv₁²/2 + nv₁*vv₁)*dS₂_dt +
        #            (5*UNITLESS_G*δm*nS₂*dr_dt*r⁻² - 5*UNITLESS_G*δm*dnS₂_dt*r⁻¹ +
        #             (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*nS₂ + (3*nv₁² + 2*vv₁)*dnS₂_dt +
        #             2*(v₂S₂ + v₁S₂)*dnv_dt + 2*(dv₂S₂_dt + dv₁S₂_dt)*nv)*v̄₁ -
        #            (UNITLESS_G*(6*m₂ - m₁)*nS₂*dr_dt*r⁻² - UNITLESS_G*(6*m₂ - m₁)*dnS₂_dt*r⁻¹ +
        #             3*nS₂*nv₁*dnv₁_dt + 3*nv₁²*dnS₂_dt/2 + nv₁*dvS₂_dt + vS₂*dnv₁_dt)*v̄₂ +
        #             (-Gm₂*(-16*nS₂*nv + 3*v₂S₂ - 7*v₁S₂)*dr_dt*r⁻² + Gm₂*(-16*nS₂*dnv_dt -
        #              16*nv*dnS₂_dt + 3*dv₂S₂_dt - 7*dv₁S₂_dt)*r⁻¹ - 2*Gm₁*nS₂*nv*dr_dt*r⁻² +
        #              2*Gm₁*nS₂*dnv_dt*r⁻¹ + 2*Gm₁*nv*dnS₂_dt*r⁻¹ + (6*nv₁*dnv₁_dt + 2*dvv₁_dt)*vS₂ +
        #             (3*nv₁² + 2*vv₁)*dvS₂_dt)*n +
        #             (-Gm₂*nv₂*dr_dt*r⁻² + Gm₂*dnv₂_dt*r⁻¹ + Gm₁*nv*dr_dt*r⁻² -
        #              Gm₁*dnv_dt*r⁻¹ - 3*nv*nv₁*dnv₁_dt - 3*nv₁²*dnv_dt/2 + nv₁*dvv₁_dt + vv₁*dnv₁_dt)*S₂
        # Gm₂r⁻²*dF2PN_dt + dGm₂r⁻²_dt*T2PN

        num = (Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * r⁻¹ + 2 * Gm₁ * nS₂ * nv * r⁻¹ + (3 * nv₁² + 2 * vv₁) * vS₂) * n
        num += (-5 * UNITLESS_G * δm * nS₂ * r⁻¹ + (3 * nv₁² + 2 * vv₁) * nS₂ + 2 * (v₂S₂ + v₁S₂) * nv) * v̄₁
        num += -(-UNITLESS_G * (6 * δm) * nS₂ * r⁻¹ + 3 * nS₂ * nv₁² / 2 + nv₁ * vS₂) * v̄₂
        num += (Gm₂ * nv₂ * r⁻¹ - Gm₁ * nv * r⁻¹ - 3 * nv * nv₁² / 2 + nv₁ * vv₁) * S̄₂
        num *= -2 * Gm₁ * dr_dt * r⁻² * r⁻¹
        num += (Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * r⁻¹ + 2 * Gm₁ * nS₂ * nv * r⁻¹ + (3 * nv₁² + 2 * vv₁) * vS₂) * dn_dt
        num += (-5 * UNITLESS_G * δm * nS₂ * r⁻¹ + (3 * nv₁² + 2 * vv₁) * nS₂ + 2 * (v₂S₂ + v₁S₂) * nv) * ā₁
        num += -(-UNITLESS_G * (6 * δm) * nS₂ * r⁻¹ + 3 * nS₂ * nv₁² / 2 + nv₁ * vS₂) * ā₂
        num += (Gm₂ * nv₂ * r⁻¹ - Gm₁ * nv * r⁻¹ - 3 * nv * nv₁² / 2 + nv₁ * vv₁) * dS̄₂
        num += (
            5 * UNITLESS_G * δm * nS₂ * dr_dt * r⁻² - 5 * UNITLESS_G * δm * dnS₂_dt * r⁻¹ + (6 * nv₁ * dnv₁_dt + 2 * dvv₁_dt) * nS₂ +
                (3 * nv₁² + 2 * vv₁) * dnS₂_dt + 2 * (v₂S₂ + v₁S₂) * dnv_dt + 2 * (dv₂S₂_dt + dv₁S₂_dt) * nv
        ) * v̄₁
        num += -(
            UNITLESS_G * (6 * δm) * nS₂ * dr_dt * r⁻² - UNITLESS_G * (6 * δm) * dnS₂_dt * r⁻¹ +
                3 * nS₂ * nv₁ * dnv₁_dt + 3 * nv₁² * dnS₂_dt / 2 + nv₁ * dvS₂_dt + vS₂ * dnv₁_dt
        ) * v̄₂
        num += (
            -Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * dr_dt * r⁻² +
                Gm₂ * (
                -16 * nS₂ * dnv_dt - 16 * nv * dnS₂_dt +
                    3 * dv₂S₂_dt - 7 * dv₁S₂_dt
            ) * r⁻¹ -
                2 * Gm₁ * nS₂ * nv * dr_dt * r⁻² + 2 * Gm₁ * nS₂ * dnv_dt * r⁻¹ +
                2 * Gm₁ * nv * dnS₂_dt * r⁻¹ +
                (6 * nv₁ * dnv₁_dt + 2 * dvv₁_dt) * vS₂ +
                (3 * nv₁² + 2 * vv₁) * dvS₂_dt
        ) * n
        num += (
            -Gm₂ * nv₂ * dr_dt * r⁻² + Gm₂ * dnv₂_dt * r⁻¹ +
                Gm₁ * nv * dr_dt * r⁻² - Gm₁ * dnv_dt * r⁻¹ -
                3 * nv * nv₁ * dnv₁_dt - 3 * nv₁² * dnv_dt / 2 +
                nv₁ * dvv₁_dt + vv₁ * dnv₁_dt
        ) * S̄₂
        num *= Gm₁ * r⁻²

        num * Gm₁ * r⁻²
    end

    dvi .+= dT2PN_dt₁ * c⁻⁴
    dvj .+= dT2PN_dt₂ * c⁻⁴
    return nothing
end

function spin_precession!(
        dvi,
        dvj,
        dvs,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )

    # i = 1, j = 2
    i, j = pair

    if any(x -> params.stellar_types[x] < 13, pair)
        return nothing
    end

    ā₁ = @SVector [dvs[1, i], dvs[2, i], dvs[3, i]]
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]

    ā₂ = @SVector [dvs[1, j], dvs[2, j], dvs[3, j]]
    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]

    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]
    dS̄₁ = @SVector [vs[4, i], vs[5, i], vs[6, i]]

    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]
    dS̄₂ = @SVector [vs[4, j], vs[5, j], vs[6, j]]


    m₁ = params.M[i]
    m₂ = params.M[j]
    δm = m₁ - m₂

    ā = ā₁ - ā₂
    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r² = r * r
    r⁻¹ = 1 / r
    r⁻² = 1 / r²

    n = r̄ / r
    nv = dot(n, v̄)
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    nS₁ = dot(n, S̄₁)
    nS₂ = dot(n, S̄₂)

    dr_dt = nv
    dn_dt = (r * v̄ - r̄ * dr_dt) * r⁻²

    Gm₁ = UNITLESS_G * m₁
    Gm₂ = UNITLESS_G * m₂


    Gr⁻³ = UNITLESS_G * r⁻² * r⁻¹
    dGr⁻³_dt = 3 * dr_dt * Gr⁻³ * r⁻¹

    dT1PN_dt₁, dT1p5PN₁, dT2PN_dt₁ = let

        vv₂ = dot(v̄, v̄₂)
        vS₁ = dot(v̄, S̄₁)
        v₁S₁ = dot(v̄₁, S̄₁)
        v₂S₁ = dot(v̄₂, S̄₁)

        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)
        dnv₁_dt = dot(n, ā₁) + dot(v̄₁, dn_dt)
        dnv₂_dt = dot(n, ā₂) + dot(v̄₂, dn_dt)
        dnS₁_dt = dot(n, dS̄₁) + dot(S̄₁, dn_dt)
        dnS₂_dt = dot(n, dS̄₂) + dot(S̄₂, dn_dt)
        dvS₁_dt = dot(v̄, dS̄₁) + dot(S̄₁, ā)
        dv₁S₁_dt = dot(v̄₁, dS̄₁) + dot(S̄₁, ā₁)
        dv₂S₁_dt = dot(v̄₂, dS̄₁) + dot(S̄₁, ā₂)
        dvv₂_dt = dot(v̄, ā₂) + dot(v̄₂, ā)

        nv₂² = nv₂^2

        ######################## PN-1 ########################
        Gm₂r⁻² = Gm₂ * r⁻²
        dGm₂r⁻²_dt = -2 * Gm₂ * dr_dt * r⁻¹ * r⁻²
        fac = (v̄₁ - 2 * v̄₂) * nS₁ + S̄₁ * nv - 2 * n * vS₁
        dfac_dt = (v̄₁ - 2 * v̄₂) * dnS₁_dt + (ā₁ - 2 * ā₂) * nS₁ +
            S̄₁ * dnv_dt + nv * dS̄₁ -
            2 * n * dvS₁_dt - 2 * vS₁ * dn_dt
        dT1PN = fac * dGm₂r⁻²_dt + Gm₂r⁻² * dfac_dt
        #######################################################

        ####################### PN-1.5 ########################
        # F1p5PN = S̄₂ .- 3nS₂*n
        # dF1p5PN_dt = dS̄₂ .- 3*n*dnS₂_dt .- 3nS₂*dn_dt
        # dT1p5PN = dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₁) + (F1p5PN × dS̄₁))
        F = (3nS₂ * n × S̄₁) - (S̄₂ × S̄₁)
        dF_dt = 3 * ((nS₂ * dn_dt + n * dnS₂_dt) × S̄₁ + nS₂ * n × dS̄₁) - (S̄₁ × dS̄₂) + (S̄₂ × dS̄₁)

        dT1p5PN = Gr⁻³ * dF_dt + F * dGr⁻³_dt

        #######################################################

        ######################## PN-2 #########################
        num = (Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * r⁻¹ + 2 * Gm₂ * nS₁ * nv * r⁻¹ + (3 * nv₂² + 2 * vv₂) * vS₁) * n
        num += (-5 * UNITLESS_G * δm * nS₁ * r⁻¹ + (3 * nv₂² + 2 * vv₂) * nS₁ + 2 * (v₁S₁ + v₂S₁) * nv) * v̄₂
        num += -(-UNITLESS_G * (6 * δm) * nS₁ * r⁻¹ + 3 * nS₁ * nv₂² / 2 + nv₂ * vS₁) * v̄₁
        num += (Gm₁ * nv₁ * r⁻¹ - Gm₂ * nv * r⁻¹ - 3 * nv * nv₂² / 2 + nv₂ * vv₂) * S̄₁
        num *= -2 * Gm₂ * dr_dt * r⁻² * r⁻¹
        num += (Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * r⁻¹ + 2 * Gm₂ * nS₁ * nv * r⁻¹ + (3 * nv₂² + 2 * vv₂) * vS₁) * dn_dt
        num += (-5 * UNITLESS_G * δm * nS₁ * r⁻¹ + (3 * nv₂² + 2 * vv₂) * nS₁ + 2 * (v₁S₁ + v₂S₁) * nv) * ā₂
        num += -(-UNITLESS_G * (6 * δm) * nS₁ * r⁻¹ + 3 * nS₁ * nv₂² / 2 + nv₂ * vS₁) * ā₁
        num += (Gm₁ * nv₁ * r⁻¹ - Gm₂ * nv * r⁻¹ - 3 * nv * nv₂² / 2 + nv₂ * vv₂) * dS̄₁
        num += (
            5 * UNITLESS_G * δm * nS₁ * dr_dt * r⁻² - 5 * UNITLESS_G * δm * dnS₁_dt * r⁻¹ + (6 * nv₂ * dnv₂_dt + 2 * dvv₂_dt) * nS₁ +
                (3 * nv₂² + 2 * vv₂) * dnS₁_dt + 2 * (v₁S₁ + v₂S₁) * dnv_dt + 2 * (dv₁S₁_dt + dv₂S₁_dt) * nv
        ) * v̄₂
        num += -(
            UNITLESS_G * (6 * δm) * nS₁ * dr_dt * r⁻² - UNITLESS_G * (6 * δm) * dnS₁_dt * r⁻¹ +
                3 * nS₁ * nv₂ * dnv₂_dt + 3 * nv₂² * dnS₁_dt / 2 + nv₂ * dvS₁_dt + vS₁ * dnv₂_dt
        ) * v̄₁
        num += (
            -Gm₁ * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) * dr_dt * r⁻² +
                Gm₁ * (
                -16 * nS₁ * dnv_dt - 16 * nv * dnS₁_dt +
                    3 * dv₁S₁_dt - 7 * dv₂S₁_dt
            ) * r⁻¹ -
                2 * Gm₂ * nS₁ * nv * dr_dt * r⁻² + 2 * Gm₂ * nS₁ * dnv_dt * r⁻¹ +
                2 * Gm₂ * nv * dnS₁_dt * r⁻¹ +
                (6 * nv₂ * dnv₂_dt + 2 * dvv₂_dt) * vS₁ +
                (3 * nv₂² + 2 * vv₂) * dvS₁_dt
        ) * n
        num += (
            -Gm₁ * nv₁ * dr_dt * r⁻² + Gm₁ * dnv₁_dt * r⁻¹ +
                Gm₂ * nv * dr_dt * r⁻² - Gm₂ * dnv_dt * r⁻¹ -
                3 * nv * nv₂ * dnv₂_dt - 3 * nv₂² * dnv_dt / 2 +
                nv₂ * dvv₂_dt + vv₂ * dnv₂_dt
        ) * S̄₁
        num *= Gm₂ * r⁻²
        dT2PN = num * Gm₂ * r⁻²
        #######################################################

        dT1PN, dT1p5PN, dT2PN
    end

    dT1PN_dt₂, dT1p5PN₂, dT2PN_dt₂ = let n = -n, v̄ = -v̄, ā = -ā, δm = -δm

        dn_dt = -dn_dt
        nv₁ = -nv₁
        nv₂ = -nv₂
        nS₁ = -nS₁
        nS₂ = -nS₂

        vv₁ = dot(v̄, v̄₁)
        vS₂ = dot(v̄, S̄₂)
        v₂S₂ = dot(v̄₂, S̄₂)
        v₁S₂ = dot(v̄₁, S̄₂)

        dnv_dt = dot(n, ā) + dot(v̄, dn_dt)
        dnv₁_dt = dot(n, ā₁) + dot(v̄₁, dn_dt)
        dnv₂_dt = dot(n, ā₂) + dot(v̄₂, dn_dt)
        dnS₁_dt = dot(n, dS̄₁) + dot(S̄₁, dn_dt)
        dnS₂_dt = dot(n, dS̄₂) + dot(S̄₂, dn_dt)
        dvS₂_dt = dot(v̄, dS̄₂) + dot(S̄₂, ā)
        dv₂S₂_dt = dot(v̄₂, dS̄₂) + dot(S̄₂, ā₂)
        dv₁S₂_dt = dot(v̄₁, dS̄₂) + dot(S̄₂, ā₁)
        dvv₁_dt = dot(v̄, ā₁) + dot(v̄₁, ā)


        nv₁² = nv₁^2

        ######################## PN-1 ########################
        Gm₁r⁻² = Gm₁ * r⁻²
        dGm₁r⁻²_dt = -2 * Gm₁ * dr_dt * r⁻¹ * r⁻²
        fac = (v̄₂ - 2 * v̄₁) * nS₂ + S̄₂ * nv - 2 * n * vS₂
        dfac_dt = (v̄₂ - 2 * v̄₁) * dnS₂_dt + (ā₂ - 2 * ā₁) * nS₂ +
            S̄₂ * dnv_dt + nv * dS̄₂ -
            2 * n * dvS₂_dt - 2 * vS₂ * dn_dt
        dT1PN = fac * dGm₁r⁻²_dt + Gm₁r⁻² * dfac_dt
        #######################################################

        ###################### PN-1.5 #########################
        # F1p5PN = S̄₁ .- 3nS₁*n
        # dF1p5PN_dt = dS̄₁ - 3*n*dnS₁_dt - 3nS₁*dn_dt
        # dT1p5PN = dGr⁻³_dt*F1p5PN - Gr⁻³*((dF1p5PN_dt × S̄₂) + (F1p5PN × dS̄₂))

        F = (3nS₁ * n × S̄₂) - (S̄₁ × S̄₂)
        dF_dt = 3 * ((nS₁ * dn_dt + n * dnS₁_dt) × S̄₂ + nS₁ * n × dS̄₂) - (S̄₂ × dS̄₁) + (S̄₁ × dS̄₂)

        dT1p5PN = Gr⁻³ * dF_dt + F * dGr⁻³_dt
        #######################################################


        ###################### PN-2 #########################
        num = (Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * r⁻¹ + 2 * Gm₁ * nS₂ * nv * r⁻¹ + (3 * nv₁² + 2 * vv₁) * vS₂) * n
        num += (-5 * UNITLESS_G * δm * nS₂ * r⁻¹ + (3 * nv₁² + 2 * vv₁) * nS₂ + 2 * (v₂S₂ + v₁S₂) * nv) * v̄₁
        num += -(-UNITLESS_G * (6 * δm) * nS₂ * r⁻¹ + 3 * nS₂ * nv₁² / 2 + nv₁ * vS₂) * v̄₂
        num += (Gm₂ * nv₂ * r⁻¹ - Gm₁ * nv * r⁻¹ - 3 * nv * nv₁² / 2 + nv₁ * vv₁) * S̄₂
        num *= -2 * Gm₁ * dr_dt * r⁻² * r⁻¹
        num += (Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * r⁻¹ + 2 * Gm₁ * nS₂ * nv * r⁻¹ + (3 * nv₁² + 2 * vv₁) * vS₂) * dn_dt
        num += (-5 * UNITLESS_G * δm * nS₂ * r⁻¹ + (3 * nv₁² + 2 * vv₁) * nS₂ + 2 * (v₂S₂ + v₁S₂) * nv) * ā₁
        num += -(-UNITLESS_G * (6 * δm) * nS₂ * r⁻¹ + 3 * nS₂ * nv₁² / 2 + nv₁ * vS₂) * ā₂
        num += (Gm₂ * nv₂ * r⁻¹ - Gm₁ * nv * r⁻¹ - 3 * nv * nv₁² / 2 + nv₁ * vv₁) * dS̄₂
        num += (
            5 * UNITLESS_G * δm * nS₂ * dr_dt * r⁻² - 5 * UNITLESS_G * δm * dnS₂_dt * r⁻¹ + (6 * nv₁ * dnv₁_dt + 2 * dvv₁_dt) * nS₂ +
                (3 * nv₁² + 2 * vv₁) * dnS₂_dt + 2 * (v₂S₂ + v₁S₂) * dnv_dt + 2 * (dv₂S₂_dt + dv₁S₂_dt) * nv
        ) * v̄₁
        num += -(
            UNITLESS_G * (6 * δm) * nS₂ * dr_dt * r⁻² - UNITLESS_G * (6 * δm) * dnS₂_dt * r⁻¹ +
                3 * nS₂ * nv₁ * dnv₁_dt + 3 * nv₁² * dnS₂_dt / 2 + nv₁ * dvS₂_dt + vS₂ * dnv₁_dt
        ) * v̄₂
        num += (
            -Gm₂ * (-16 * nS₂ * nv + 3 * v₂S₂ - 7 * v₁S₂) * dr_dt * r⁻² +
                Gm₂ * (
                -16 * nS₂ * dnv_dt - 16 * nv * dnS₂_dt +
                    3 * dv₂S₂_dt - 7 * dv₁S₂_dt
            ) * r⁻¹ -
                2 * Gm₁ * nS₂ * nv * dr_dt * r⁻² + 2 * Gm₁ * nS₂ * dnv_dt * r⁻¹ +
                2 * Gm₁ * nv * dnS₂_dt * r⁻¹ +
                (6 * nv₁ * dnv₁_dt + 2 * dvv₁_dt) * vS₂ +
                (3 * nv₁² + 2 * vv₁) * dvS₂_dt
        ) * n
        num += (
            -Gm₂ * nv₂ * dr_dt * r⁻² + Gm₂ * dnv₂_dt * r⁻¹ +
                Gm₁ * nv * dr_dt * r⁻² - Gm₁ * dnv_dt * r⁻¹ -
                3 * nv * nv₁ * dnv₁_dt - 3 * nv₁² * dnv_dt / 2 +
                nv₁ * dvv₁_dt + vv₁ * dnv₁_dt
        ) * S̄₂
        num *= Gm₁ * r⁻²
        dT2PN = num * Gm₁ * r⁻²
        #######################################################

        dT1PN, dT1p5PN, dT2PN
    end

    dvi .+= dT1PN_dt₁ * c⁻² + dT1p5PN₁ * c⁻³ + dT2PN_dt₁ * c⁻⁴
    dvj .+= dT1PN_dt₂ * c⁻² + dT1p5PN₂ * c⁻³ + dT2PN_dt₂ * c⁻⁴
    return nothing
end

function PN2p5_spin_acceleration!(
        dvi,
        dvj,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )

    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]

    m₁ = params.M[i]
    m₂ = params.M[j]

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    v₁xv₂ = v̄₁ × v̄₂
    v₂xv₁ = v̄₂ × v̄₁

    r = norm(r̄) # r₁₂

    r⁻¹ = 1 / r

    n = r̄ * r⁻¹
    nv = dot(n, v̄)
    nv₁ = dot(n, v̄₁)
    nv₂ = dot(n, v̄₂)

    vv₁ = dot(v̄, v̄₁)
    vv₂ = dot(v̄, v̄₂)

    nS₁ = dot(n, S̄₁)
    nS₂ = dot(n, S̄₂)

    nxv = n × v̄
    nxv₁ = n × v̄₁
    nxv₂ = n × v̄₂

    nxS₁ = n × S̄₁
    nxS₂ = n × S̄₂

    vxS₁ = v̄ × S̄₁
    vxS₂ = v̄ × S̄₂

    v₁S₁ = dot(v̄₁, S̄₁)
    v₂S₂ = dot(v̄₂, S̄₂)

    nv₁v₂ = dot(n, v₁xv₂)
    nv₂v₁ = dot(n, v₂xv₁)
    S₁nv = dot(S̄₁, nxv)
    S₂nv = dot(S̄₂, nxv)
    S₁nv₁ = dot(S̄₁, nxv₁)
    S₁nv₂ = dot(S̄₁, nxv₂)
    S₂nv₁ = dot(S̄₂, nxv₁)
    S₂nv₂ = dot(S̄₂, nxv₂)

    G_r = UNITLESS_G * r⁻¹
    G_r³ = G_r * r⁻¹ * r⁻¹

    Gm₁ = UNITLESS_G * m₁
    Gm₂ = UNITLESS_G * m₂

    ai = n * (
        -6 * nv₁v₂ * (v₁S₁ / m₁ + v₂S₂ / m₂) -
            S₁nv / m₁ * (15nv₂^2 + 6vv₂ + 26Gm₁ / r + 18Gm₂ / r) -
            S₂nv / m₂ * (15nv₂^2 + 6vv₂ + 49 / 2 * Gm₁ / r + 20Gm₂ / r)
    )
    ai += v̄₁ * (
        -3 * S₁nv₁ / m₁ * (nv₁ + nv₂) + 6nv₁ * S₁nv₂ / m₁ - 3 * dot(S̄₁, v₁xv₂) / m₁ -
            6nv₁ * S₂nv₁ / m₂ + S₂nv₂ / m₂ * (12nv₁ - 6nv₂) - 4 * dot(S̄₂, v₁xv₂) / m₂
    )
    ai += v̄₂ * (6nv₁ * S₁nv / m₁ + 6nv₁ * S₂nv / m₂)
    ai -= nxv₁ * (3nv * v₁S₁ / m₁ + 4Gm₁ / r * nS₂ / m₂)
    ai -= nxv₂ * (6nv * v₂S₂ / m₂ - 4Gm₁ / r * nS₂ / m₂)
    ai += v₁xv₂ * (3v₁S₁ / m₁ + 4v₂S₂ / m₂)
    ai += nxS₁ / m₁ * (
        -15 / 2 * nv * nv₂^2 + 3nv₂ * vv₂ -
            14Gm₁ / r * nv - 9Gm₂ / r * nv
    )
    ai += nxS₂ / m₂ * (-15nv * nv₂^2 - 6nv₁ * vv₂ + 12nv₂ * vv₂ + Gm₁ / r * (-35 / 2 * nv₁ + 39 / 2 * nv₂) - 16Gm₂ / r * nv)
    ai += vxS₁ / m₁ * (-3nv₁ * nv₂ + 15 / 2 * nv₂^2 + UNITLESS_G / r * (14m₁ + 9m₂) + 3vv₂)
    ai += vxS₂ / m₂ * (6nv₂^2 + 4vv₂ + 23 / 2 * Gm₁ / r + 12Gm₂ / r)
    ai *= G_r³ * m₂

    aj = let n = -n

        nv₁ = -nv₁
        nv₂ = -nv₂
        vv₁ = -vv₁
        nS₁ = -nS₁

        nxv₁ = -nxv₁
        nxv₂ = -nxv₂

        nxS₁ = -nxS₁
        nxS₂ = -nxS₂

        vxS₁ = -vxS₁
        vxS₂ = -vxS₂

        nv₂v₁ = -nv₂v₁
        S₁nv₁ = -S₁nv₁
        S₁nv₂ = -S₁nv₂
        S₂nv₁ = -S₂nv₁
        S₂nv₂ = -S₂nv₂

        a = n * (
            -6 * nv₂v₁ * (v₂S₂ / m₂ + v₁S₁ / m₁) -
                S₂nv / m₂ * (15nv₁^2 + 6vv₁ + 26Gm₂ / r + 18Gm₁ / r) -
                S₁nv / m₁ * (15nv₁^2 + 6vv₁ + 49 / 2 * Gm₂ / r + 20Gm₁ / r)
        )
        a += v̄₂ * (
            -3 * S₂nv₂ / m₂ * (nv₂ + nv₁) + 6nv₂ * S₂nv₁ / m₂ - 3 * dot(S̄₂, v₂xv₁) / m₂ -
                6nv₂ * S₁nv₂ / m₁ + S₁nv₁ / m₁ * (12nv₂ - 6nv₁) - 4 * dot(S̄₁, v₂xv₁) / m₁
        )
        a += v̄₁ * (6nv₂ * S₂nv / m₂ + 6nv₂ * S₁nv / m₁)
        a -= nxv₂ * (3nv * v₂S₂ / m₂ + 4Gm₂ / r * nS₁ / m₁)
        a -= nxv₁ * (6nv * v₁S₁ / m₁ - 4Gm₂ / r * nS₁ / m₁)
        a += v₂xv₁ * (3v₂S₂ / m₂ + 4v₁S₁ / m₁)
        a += nxS₂ / m₂ * (
            -15 / 2 * nv * nv₁^2 + 3nv₁ * vv₁ -
                14Gm₂ / r * nv - 9Gm₁ / r * nv
        )
        a += nxS₁ / m₁ * (-15nv * nv₁^2 - 6nv₂ * vv₁ + 12nv₁ * vv₁ + Gm₂ / r * (-35 / 2 * nv₂ + 39 / 2 * nv₁) - 16Gm₁ / r * nv)
        a += vxS₂ / m₂ * (-3nv₂ * nv₁ + 15 / 2 * nv₁^2 + UNITLESS_G / r * (14m₂ + 9m₁) + 3vv₁)
        a += vxS₁ / m₁ * (6nv₁^2 + 4vv₁ + 23 / 2 * Gm₂ / r + 12Gm₁ / r)
        a *= G_r³ * m₁

        a
    end

    dvi .+= ai * c⁻⁵
    dvj .+= aj * c⁻⁵
    return nothing
end

function PN2_spin_acceleration!(
        dvi,
        dvj,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )

    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    # v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    # v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]

    m₁ = params.M[i]
    m₂ = params.M[j]

    r̄ = r̄₁ - r̄₂
    # v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1 / r

    n = r̄ * r⁻¹
    S₁S₂ = dot(S̄₁, S̄₂)
    S₂S₁ = S₁S₂
    nS₁ = dot(n, S̄₁)
    nS₂ = dot(n, S̄₂)

    G_r = UNITLESS_G * r⁻¹
    G_r⁴ = G_r * r⁻¹ * r⁻¹ * r⁻¹

    # PN-2 spin-spin interaction
    ai = -G_r⁴ * 3 / m₁ * (n * S₁S₂ + S̄₁ * nS₁ + S̄₂ * nS₁ - 5n * nS₁ * nS₂)

    aj = -G_r⁴ * 3 / m₂ * (n * S₂S₁ + S̄₂ * nS₂ + S̄₁ * nS₂ - 5n * nS₂ * nS₁)

    dvi .+= ai * c⁻⁴
    dvj .+= aj * c⁻⁴
    return nothing
end

function PN1p5_spin_acceleration!(
        dvi,
        dvj,
        rs,
        vs,
        pair::Tuple{Int, Int},
        params::SimulationParams
    )

    i, j = pair # i = 1, j = 2
    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    S̄₁ = @SVector [rs[4, i], rs[5, i], rs[6, i]]

    r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
    v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
    S̄₂ = @SVector [rs[4, j], rs[5, j], rs[6, j]]

    m₁ = params.M[i]
    m₂ = params.M[j]

    r̄ = r̄₁ - r̄₂
    v̄ = v̄₁ - v̄₂

    r = norm(r̄) # r₁₂

    r⁻¹ = 1 / r
    r⁻³ = r⁻¹ / r^2

    n = r̄ * r⁻¹
    nv = dot(n, v̄)
    nxv = n × v̄

    # @show (n × S̄₂) (v̄ × S̄₁)

    G_r = UNITLESS_G * r⁻¹
    G_r³ = G_r * r⁻³

    # PN-1.5 acceleration from spin contribution
    ai = G_r³ * m₂ * (
        (6 * dot(S̄₁, nxv) / m₁ + 6 * dot(S̄₂, nxv) / m₂) * n + 3nv * (n × S̄₁) / m₁ +
            6nv * (n × S̄₂) / m₂ - 3 * (v̄ × S̄₁) / m₁ - 4 * (v̄ × S̄₂) / m₂
    )

    aj = G_r³ * m₁ * (
        (6 * dot(S̄₂, nxv) / m₂ + 6 * dot(S̄₁, nxv) / m₁) * n + 3nv * (n × S̄₂) / m₂ +
            6nv * (n × S̄₁) / m₁ - 3 * (v̄ × S̄₂) / m₂ - 4 * (v̄ × S̄₁) / m₁
    )

    dvi .+= ai * c⁻³
    dvj .+= aj * c⁻³
    return nothing
end


########################### All PN spin velocity terms ###########################
get_spin_precession_velocity(object1, object2, potential::PN1SpinPrecessionPotential) = PN1_spin_precession_velocity(object1, object2)
get_spin_precession_velocity(object1, object2, potential::PN1p5SpinPrecessionPotential) = PN1p5_spin_precession_velocity(object1, object2)
get_spin_precession_velocity(object1, object2, potential::PN2SpinPrecessionPotential) = PN2_spin_precession_velocity(object1, object2)
get_spin_precession_velocity(object1, object2, potential::SpinPrecessionPotential) = spin_precession_velocity(object1, object2)

"""
    spin_precession_velocity(particle1::Particle, particle2::Particle)


"""
function spin_precession_velocity(particle1::Particle, particle2::Particle)
    S1 = particle1.structure.S
    S2 = particle2.structure.S
    r1 = particle1.position
    r2 = particle2.position
    v1 = particle1.velocity
    v2 = particle2.velocity

    m1 = particle1.mass
    m2 = particle2.mass
    return spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function spin_precession_velocity(body1::MassBody, body2::MassBody)
    S1 = body1.spin
    S2 = body2.spin
    r1 = body1.position
    r2 = body2.position
    v1 = body1.velocity
    v2 = body2.velocity

    m1 = body1.mass
    m2 = body2.mass
    return spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::Quantity, m2::Quantity)
    T1PN = PN1_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    T1p5PN = PN1p5_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    T2PN = PN2_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2, GRAVCONST)
    return GRAVCONST * (T1PN / c^2 + T1p5PN / c^3 + T2PN / c^4)
end

function spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::AbstractFloat, m2::AbstractFloat)
    T1PN = PN1_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    T1p5PN = PN1p5_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    T2PN = PN2_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2, UNITLESS_G)
    return UNITLESS_G * (T1PN / c² + T1p5PN / c³ + T2PN / c⁴)
end
####################################################################################

############################# PN-1 spin velocity terms #############################
function PN1_spin_precession_velocity(particle1::Particle, particle2::Particle)
    S1 = particle1.structure.S
    S2 = particle2.structure.S
    r1 = particle1.position
    r2 = particle2.position
    v1 = particle1.velocity
    v2 = particle2.velocity

    m1 = particle1.mass
    m2 = particle2.mass
    return PN1_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN1_spin_precession_velocity(body1::MassBody, body2::MassBody)
    S1 = body1.spin
    S2 = body2.spin
    r1 = body1.position
    r2 = body2.position
    v1 = body1.velocity
    v2 = body2.velocity

    m1 = body1.mass
    m2 = body2.mass
    return PN1_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN1_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::Quantity, m2::Quantity)
    T1PN = PN1_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    return GRAVCONST * T1PN / c^2
end

function PN1_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::AbstractFloat, m2::AbstractFloat)
    T1PN = PN1_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    return UNITLESS_G * T1PN / c²
end
####################################################################################

############################# PN-1.5 spin velocity terms ###########################
function PN1p5_spin_precession_velocity(particle1::Particle, particle2::Particle)
    S1 = particle1.structure.S
    S2 = particle2.structure.S
    r1 = particle1.position
    r2 = particle2.position
    v1 = particle1.velocity
    v2 = particle2.velocity

    m1 = particle1.mass
    m2 = particle2.mass
    return PN1p5_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN1p5_spin_precession_velocity(body1::MassBody, body2::MassBody)
    S1 = body1.spin
    S2 = body2.spin
    r1 = body1.position
    r2 = body2.position
    v1 = body1.velocity
    v2 = body2.velocity

    m1 = body1.mass
    m2 = body2.mass
    return PN1p5_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN1p5_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::Quantity, m2::Quantity)
    T1p5PN = PN1p5_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    return GRAVCONST * T1p5PN / c^3
end

function PN1p5_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::AbstractFloat, m2::AbstractFloat)
    T1p5PN = PN1p5_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2)
    return UNITLESS_G * T1p5PN / c³
end
####################################################################################

############################# PN-2 spin velocity terms #############################
function PN2_spin_precession_velocity(particle1::Particle, particle2::Particle)
    S1 = particle1.structure.S
    S2 = particle2.structure.S
    r1 = particle1.position
    r2 = particle2.position
    v1 = particle1.velocity
    v2 = particle2.velocity

    m1 = particle1.mass
    m2 = particle2.mass
    return PN2_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN2_spin_precession_velocity(body1::MassBody, body2::MassBody)
    S1 = body1.spin
    S2 = body2.spin
    r1 = body1.position
    r2 = body2.position
    v1 = body1.velocity
    v2 = body2.velocity

    m1 = body1.mass
    m2 = body2.mass
    return PN2_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1, m2)
end

function PN2_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::Quantity, m2::Quantity)
    T2PN = PN2_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2, GRAVCONST)
    return GRAVCONST * T2PN / c^4
end

function PN2_spin_precession_velocity(S1, S2, r1, r2, v1, v2, m1::AbstractFloat, m2::AbstractFloat)
    T2PN = PN2_spin_precession_velocity_factor(S1, S2, r1, r2, v1, v2, m1, m2, UNITLESS_G)
    return UNITLESS_G * T2PN / c⁴
end
####################################################################################


function PN1_spin_precession_velocity_factor(S̄₁, S̄₂, r₁, r₂, v₁, v₂, m₁, m₂)
    r̄ = r₁ - r₂
    v̄ = v₁ - v₂
    r = norm(r̄)

    n̄ = r̄ / r
    return m₂ / r^2 * (S̄₁ * (n̄ ⋅ v̄) - 2n̄ * (v̄ ⋅ S̄₁) + (v₁ - 2v₂) * (n̄ ⋅ S̄₁))
end

function PN1p5_spin_precession_velocity_factor(S̄₁, S̄₂, r₁, r₂, v₁, v₂, m₁, m₂)
    r̄ = r₁ - r₂
    r = norm(r̄)

    n̄ = r̄ / r

    # with a minus or not???
    return 1 / r^3 * (S̄₂ - 3 * (n̄ ⋅ S̄₂) * n̄) × S̄₁
end

function PN2_spin_precession_velocity_factor(S̄₁, S̄₂, r₁, r₂, v₁, v₂, m₁, m₂, G)

    r̄ = r₁ - r₂
    v̄ = v₁ - v₂

    r = norm(r̄)

    n̄ = r̄ / r
    nS₁ = dot(n̄, S̄₁)
    nv = dot(n̄, v̄)
    nv₁ = dot(n̄, v₁)
    nv₂ = dot(n̄, v₂)
    vv₂ = dot(v̄, v₂)
    vS₁ = dot(v̄, S̄₁)
    v₁S₁ = dot(v₁, S̄₁)
    v₂S₁ = dot(v₂, S̄₁)

    return m₂ / r^2 * (
        S̄₁ * (nv₂ * vv₂ - 3 / 2 * nv₂^2 * nv + G * m₁ / r * nv₁ - G * m₂ / r * nv) +
            n̄ * (
            vS₁ * (3 * nv₂^2 + 2 * vv₂) + G * m₁ / r * (-16 * nS₁ * nv + 3 * v₁S₁ - 7 * v₂S₁) +
                2nS₁ * G * m₂ / r * nv
        ) - v₁ * (
            3 / 2 * nS₁ * nv₂^2 + vS₁ * nv₂ -
                nS₁ * G / r * (6m₁ - m₂)
        ) + v₂ * (
            nS₁ * (2vv₂ + 3nv₂^2) +
                2nv * (v₁S₁ + v₂S₁) - 5nS₁ * G / r * (m₁ - m₂)
        )
    )
end

function get_initial_conditions(simulation::MultiBodySimulation, dtype, ::Type{SpinPotential})
    bodies = simulation.bodies

    spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    if any(x -> x isa SpinPotential, values(simulation.potential))
        for pot in values(simulation.potential)
            !(pot isa SpinPotential) && continue

            system = simulation.ic
            for pair in system.pairs
                i, j = pair
                b1 = bodies[i]
                b2 = bodies[j]
                spinvelocity[i] += get_spin_precession_velocity(b1, b2, pot)
                spinvelocity[j] += get_spin_precession_velocity(b2, b1, pot)
            end
        end
    end

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]
    ss = [b.spin for b in bodies]

    u0, v0 = make_initial_conditions(us, vs, ss, spinvelocity, dtype)

    return u0, v0
end

function DiffEqBase.SecondOrderODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        spin_acc_funcs::AccelerationFunctions,
        dtype::Type{ArbFloat}
    )

    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    return SecondOrderODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0)
end

function DiffEqBase.SecondOrderODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        spin_acc_funcs::AccelerationFunctions,
        dtype::Type{<:AbstractFloat}
    )

    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    return SecondOrderODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0)
end


function DiffEqBase.SecondOrderODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        spin_acc_funcs::AccelerationFunctions,
        u0, v0, ai, aj
    )
    pairs = simulation.ic.pairs

    accelerations = FunctionWrangler(acc_funcs.fs)
    output = Vector{Nothing}(undef, acc_funcs.N)

    spin_accelerations = FunctionWrangler(spin_acc_funcs.fs)
    spin_out = Vector{Nothing}(undef, spin_acc_funcs.N)


    dtype = eltype(u0)
    dtype_0 = zero(dtype)
    function soode_system!(dv, v, u, p, t)
        fill!(dv, dtype_0)
        return @inbounds for pair in pairs
            i, j = pair
            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(output, accelerations, ai, aj, u, v, pair, t, p)

            # @inbounds for k = 1:3
            #     dv[k, i] += ai[k]
            #     dv[k, j] += aj[k]
            # end

            dv[1, i] += ai[1]
            dv[1, j] += aj[1]

            dv[2, i] += ai[2]
            dv[2, j] += aj[2]

            dv[3, i] += ai[3]
            dv[3, j] += aj[3]

            fill!(ai, dtype_0)
            fill!(aj, dtype_0)

            smap!(spin_out, spin_accelerations, ai, aj, dv, u, v, pair, t, p)

            dv[4, i] += ai[4]
            dv[4, j] += aj[4]

            dv[5, i] += ai[5]
            dv[5, j] += aj[5]

            dv[6, i] += ai[6]
            dv[6, j] += aj[6]
        end

    end

    return SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan, simulation.params)
end

function get_initial_conditions_static(simulation::MultiBodySimulation)
    bodies = simulation.bodies
    n = length(bodies)

    spinvelocity = [zeros(eltype(b.spin), 3) for b in bodies]
    if any(x -> x isa SpinPotential, values(simulation.potential))
        system = simulation.ic
        for pair in system.pairs
            i, j = pair
            b1 = bodies[i]
            b2 = bodies[j]
            dS = spin_precession_velocity(b1, b2)
            spinvelocity[i] += dS
        end
    end

    us = [b.position for b in bodies]
    vs = [b.velocity for b in bodies]
    ss = [b.spin for b in bodies]

    u0 = SMatrix{6, n}([reduce(hcat, us); reduce(hcat, ss)])
    v0 = SMatrix{6, n}([reduce(hcat, vs); reduce(hcat, spinvelocity)])

    return (u0, v0, n)
end

function DiffEqBase.ODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        spin_acc_funcs::AccelerationFunctions,
        dtype::Type{ArbFloat}
    )

    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    ai = SizedVector{3, dtype}(zeros(dtype, 3)...)
    aj = SizedVector{3, dtype}(zeros(dtype, 3)...)
    return ODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0, ai, aj)
end

function DiffEqBase.ODEProblem(
        simulation::MultiBodySimulation,
        acc_funcs::AccelerationFunctions,
        spin_acc_funcs::AccelerationFunctions,
        dtype::Type{<:AbstractFloat}
    )

    u0, v0 = get_initial_conditions(simulation, dtype, SpinPotential)
    ai = MVector{3, dtype}(zeros(dtype, 3)...)
    aj = MVector{3, dtype}(zeros(dtype, 3)...)
    return ODEProblem(simulation, acc_funcs, spin_acc_funcs, u0, v0, ai, aj)
end


# function get_accelerating_function(potential::PN1p5SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN1p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN2_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2p5SpinPotential)
#     (dvi, dvj, rs, vs, pair, time, params) -> PN2p5_spin_acceleration!(dvi, dvj, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN1SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN1p5SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN1p5_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end

# function get_accelerating_function(potential::PN2SpinPrecessionPotential)
#     (dvi, dvj, dvs, rs, vs, pair, time, params) -> PN2_spin_precession!(dvi, dvj, dvs, rs, vs, pair, params)
# end
