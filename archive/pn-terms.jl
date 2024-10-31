function PN3_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN3Potential)

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

     
    
    m₁ = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    a₂ =  @MVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @inbounds for j = 1:n
        if j != i   
            r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂ = norm(v̄₂)

            v₂² = v₂^2

            r̄ = r̄₁ - r̄₂
            v̄ = v̄₁ - v̄₂

            r = norm(r̄) # r₁₂
            v = norm(v̄) # v₁₂

            r⁻¹ = 1/r

            r² = r^2
            r³ = r²*r
            r⁴ = r³*r
            r⁵ = r⁴*r
            
            v² = v^2
            n = r̄*r⁻¹

            v₁v₂ = dot(v̄₁, v̄₂) 
            nv₁ = dot(n, v̄₁)
            nv₂ = dot(n, v̄₂)
            nv = dot(n, v̄)
            # nv² = nv^2

            v₁v₂² = v₁v₂^2 

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            # nv₁³ = nv₁^3
            nv₂³ = nv₂^3

            # nv₁⁴ = nv₁^4
            nv₂⁴ = nv₂^4

            m₂ = params.M[j]
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2
            m₁²m₂² = m₁^2*m₂^2

            G_r = G*r⁻¹
            G_r² = G_r*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            r₁′ = r₂′ = 1.0
      
            # TO-DO: split up to avoid allocations; calculate the gauge constants r′
            a = @. n*(G_r²*m₂*(35/16*nv₂^6 - 15/8*nv₂⁴*v₁² + 15/2*nv₂⁴*v₁v₂ + 3*nv₂²*v₁v₂² -
                                     15/2*nv₂⁴*v₂² + 3/2*nv₂²*v₁²*v₂² - 12*nv₂²*v₁v₂*v₂² - 2*v₁v₂²*v₂² + 
                                     15/2*nv₂²*v₂^4 + 4*v₁v₂*v₂^4 - 2v₂^6
                                    ) +
                      G²_r³*m₁m₂*(-171/8*nv₁^4 + 171/2*nv₁³*nv₂ - 723/4*nv₁²*nv₂² +
                                          383/2*nv₁*nv₂³ - 455/8*nv₂⁴ + 229/4*nv₁²*v₁² - 
                                          205/2*nv₁*nv₁*v₁² + 191/4*nv₂²*v₁² - 91/8*v₁^4 - 229/2*nv₁²*v₁v₂ +
                                          244*nv₁*nv₂*v₁v₂ - 225/2*nv₂²*v₁v₂ + 91/2*v₁²*v₁v₂ -
                                          177/4*v₁v₂² + 229/4*nv₁²*v₂² - 283/2*nv₁*nv₂*v₂² +
                                          259/4*nv₂²*v₂² - 91/4*v₁²*v₂² + 43*v₁v₂*v₂² - 81/8*v₂^4
                                         ) +
                     G²_r³*m₂^2*(-6*nv₁²*nv₂² + 12*nv₁*nv₂³ + 6*nv₂⁴ + 
                                         4*nv₁*nv₂*v₁v₂ + 12*nv₂³*v₁v₂ + 4*v₁v₂ -
                                         4*nv₁*nv₂*v₂² - 12*nv₂²*v₂² - 8*v₁v₂*v₂² + 4v₂^4
                                       ) +
                     G³_r⁴*m₂^3*(-nv₁² + 2*nv₁*nv₂ + 43/2*nv₂² + 18*v₁v₂ - 9v₂²) +
                     G³_r⁴*m₁m₂²*(415/8*nv₁² - 375/4*nv₁*nv₂ + 1113/8*nv₂² - 615/64*nv²*π² +
                                           18v₁² + 123/64*π²*v² + 33*v₁v₂ - 33/2*v₂²) + 
                     G³_r⁴*m₁²m₂*(-45887/168*nv₁² + 24025/42*nv₁*nv₂ - 10469/42*nv₂² + 48197/840*v₁² -
                                           36227/420*v₁v₂ + 36227*v₂² + 110*nv²*log(r̄/r₁′) - 22*v²*log(r̄/r₁′)) + 
                     16G⁴*m₂^4/r⁵ + G⁴*m₁²m₂²/r⁵*(175 - 41/16*π² - 44/3*log(r̄/r₂′))) +
                     (G_r²*m₂*(15/2*nv₁*nv₂⁴ - 45/8*nv₂^5 - 3/2*nv₂³*v₁² + 6*nv₁*nv₂²*v₁v₂ -
                                     6*nv₂³*v₁v₂ - 2*nv₂*v₁v₂² - 12*nv₁*nv₂²*v₂² + 12*nv₂³*v₂² +
                                     nv₂*v₁²*v₂² - 4*nv₁*v₁v₂*v₂² + 8*nv₂*v₁v₂*v₂² + 4*nv₁*v₂^4 -
                                     7*nv₂*v₂^4) +
                      G²_r³*m₂^2*(-2*nv₁²*nv₂ + 8*nv₁*nv₂² + 2*nv₂³ + 2*nv₁*v₁v₂ + 
                                         4*nv₂*v₁v₂ - 2*nv₁*v₂² - 4*nv₂*v₂²) +
                      G²_r³*m₁m₂*(-243/4*nv₁³ + 565/4*nv₁²*nv₂ - 269/4*nv₁*nv₂² -
                                          95/12*nv₂³ + 207/8*nv₁*v₁² - 137/8*nv₂*v₁² - 36*nv₁*v₁v₂ + 
                                          27/4*nv₂*v₁v₂ + 81/8*nv₁*v₂² + 83/8*nv₂*v₂²) + 
                      G³_r⁴*m₂^3*(4*nv₁ + 5*nv₂) + 
                      G³_r⁴*m₁m₂²*(-307/8*nv₁ - 479/8*nv₂ + 123/32*nv*π²) + 
                      G³_r⁴*m₁²m₂*(31397/420*nv₁ - 36227/427*nv₂ - 44*nv*log(r/r₁′)))*v̄

            accel += a
        end
    end

    @. dv += accel * c⁻⁶
end

function PN3_5_acceleration!(dv,
                            rs,
                            vs,
                            params::SimulationParams,
                            i::Int,
                            n::Int,
                            potential::PN3_5Potential)

    r̄₁ = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    v̄₁ = @SVector [vs[1, i], vs[2, i], vs[3, i]]
    v₁ = norm(v̄₁)

    v₁² = v₁^2

     
    
    m₁ = params.M[i]
    accel = @SVector [0.0, 0.0, 0.0]
    
    # i = 1, j = 2
    # add @fastmath?
    @inbounds for j = 1:n
        if j != i   
            r̄₂ = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            v̄₂ = @SVector [vs[1, j], vs[2, j], vs[3, j]]
            v₂ = norm(v̄₂)

            v₂² = v₂^2

            r̄ = r̄₁ - r̄₂
            v̄ = v̄₁ - v̄₂

            r = norm(r̄) # r₁₂
            v = norm(v̄) # v₁₂

            r⁻¹ = 1/r

            r² = r^2
            r³ = r²*r
            r⁴ = r³*r
            r⁵ = r⁴*r
         
            v² = v^2


            n = r̄*r⁻¹

            v₁v₂ = dot(v̄₁, v̄₂) 
            nv₁ = dot(n, v̄₁)
            nv₂ = dot(n, v̄₂)
            nv = dot(n, v̄)
            # nv² = nv^2

            v₁v₂² = v₁v₂^2 

            nv₁² = nv₁^2
            nv₂² = nv₂^2

            # nv₁³ = nv₁^3
            nv₂³ = nv₂^3

            # nv₁⁴ = nv₁^4

            m₂ = params.M[j]
            m₁m₂ = m₁*m₂
            m₁²m₂ = m₁^2*m₂
            m₁m₂² = m₁*m₂^2
            m₁²m₂² = m₁^2*m₂^2

            G_r = G*r⁻¹
            G²_r³ = G²*r⁻¹^3
            G³_r⁴ = G³*r⁻¹^4

            # r₁′ = r₂′ = 1.0
      

            a = @. n*(G⁴*m₁^3*m₂/r⁵*(3992/105*nv₁ - 4328/105*nv₂) + 
                      G⁴*m₁²m₂²/r⁶*(-13576/105*nv₁ + 2872/21*nv₂) - 3172/21*G⁴*m₁*m₂^3/r⁶*nv +
                      G³_r⁴*m₁²m₂*(48*nv₁³ - 696/5*nv₁²*nv₂ + 744/5*nv₁*nv₂² - 288/5*nv₂³ -
                                            4888/105*nv₁*v₁² + 5056*nv₂*v₁² + 2056/21*nv₁*v₁v₂ -
                                            2224/21*nv₂*v₁v₂ - 1028/21*nv₁*v₂² + 5812/105*nv₂*v₂²) + 
                      G³_r⁴*m₁m₂²*(-582/5*nv₂³ + 1746/5*nv₁²*nv₂ - 1954/5*nv₁*nv₂² +
                                            158*nv₂³ + 3568/105*nv*v₁² - 2864/35*nv₁*v₁v₂ +
                                            10048/105*nv₂*v₁v₂ + 1432/35*nv₁*v₂² - 5752/105*nv₂*v₂²) +
                      G²_r³*m₁m₂*(-56*nv^5 + 60*nv₁³*v² - 180*nv₁²*nv₂*v² + 
                                          174*nv₁*nv₂²*v² - 54*nv₂³*v² - 246/35*nv*v₁^4 +
                                          1068/35*nv₁*v₁²*v₁v₂ - 984/35*nv₂*v₂²*v₁v₂ - 1068/35*nv₁*v₁v₂² +
                                          180/7*nv₂*v₁v₂² - 534/35*nv₁*v₁²*v₂² + 90/7*nv₂*v₁²*v₂² +
                                          984/35*nv₁*v₁v₂*v₂² - 732/35*nv₂*v₁v₂*v₂² - 204/35*nv₁*v₂^4 + 
                                          24/7*nv₂*v₂^4)) + 
                   v*(-184/21*G⁴*m₁^3*m₂/r⁵ + 6224/105*G⁴*m₁²m₂²/r⁶ + 6388/105*G⁴*m₁*m₂^3/r⁶ +
                      G³*m₁²m₂*(52/15*nv₁² - 56/15*nv₁*nv₂ - 44/15*nv₂² - 132/35*v₁² + 152/35*v₁v₂ - 48/35*v₂²) +
                      G³_r⁴*m₁m₂²*(454/15*nv₁² - 372/5*nv₁*nv₂ + 854/15*nv₂² - 152/21*v₁² + 
                                            2864/105*v₁v₂ - 1768/105*v₂²) +
                      G²_r³*m₁m₂*(60*nv^4 - 348/5*nv₁²*v² + 684/5*nv₁*nv₁*v² -
                                          66*nv₂²*v² + 334/35*v₁^4 - 1336/35*v₁²*v₁v₂ + 1308/35*v₁v₂² + 654/35*v₁²*v₂² -
                                          1252/35*v₁v₂*v₂² + 292/35*v₂^4))

            accel += a
        end
    end

    @. dv += accel * c⁻⁷
end
