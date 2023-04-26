


"""

hierarchy_vector should be an Integer vector with number of bodies as first
element, and number of binaries on each level as the remaining elements.
"""
function hierarchy_matrix(hierarchy_vector)

    n_bodies = hierarchy_vector[1]
    binaries = hierarchy_vector[2:end]

    level = ifelse(isone(binaries[end]), length(binaries), length(binaries) - 1)

    H = zeros(Int, n_bodies, n_bodies)

    bodies, binariesp, row = first_level!(H, n_bodies, binaries)
    if binaries[end] > 1
        binaries = binaries[1:end-1]
        nested_symmetric!(H, n_bodies, binaries, bodies, binariesp, row, 2)
    elseif 2*binariesp == n_bodies
        symmetric!(H, n_bodies, binaries, bodies, binariesp, row, 2)
    else
        n_level!(H, n_bodies, binaries, bodies, binariesp, row, 2)
    end

    H[end, :] .= -1
    return H
end

function first_level!(H, n_bodies, binaries)

    H[1,1] = -1.0
    H[1,2] = 1.0
    if binaries[1] > 1
        j = ifelse(isone(binaries[end]), 1, n_bodies/2 - 1)
        for i = 1:binaries[1] - 1
            if (j+3) <= n_bodies
                H[i+1,j+2] = -1
                H[i+1,j+3] = 1
                j += 2
            end
        end
    end
    binariesp = binaries[1]
    row = binariesp[1] + 1
    bodies = binaries[1] * 2

    return bodies, binariesp, row
end

function n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration)

    if iteration <= length(binaries)
        if n_bodies != bodies
            if binaries[iteration] == binariesp
                bodies = bodies + binaries[iteration]
            elseif binaries[iteration] > binariesp
                bodies = bodies + 2*(binaries[iteration]) - 1
            end
        elseif n_bodies == bodies
            if row == n_bodies-1
                if (bodies%4) > 2
                    H[row,1:row-2] .= -1
                    H[row,row-1:bodies] .= 1
                    return H
                elseif (bodies%4) <= 2
                    H[row,1:row-1] .= -1
                    H[row,row:bodies] .= 1
                    return H
                end
            end
        end
    
        # Looks at the previous and current binaries and calls nlevel again
        if n_bodies >= bodies
            if binariesp == 1
                if binaries[iteration] == 1
                    H[row,1:bodies-binariesp] .= -1
                    H[row,bodies-binariesp + 1:bodies] .= 1
                    row = row + binariesp
                    binariesp = binaries[iteration]
                    n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration+1)
                elseif binaries[iteration] == 2
                    H[row,1:bodies-(2*binariesp)-1] .= -1
                    H[row,bodies-(2*binariesp)] = 1
                    H[row+1,bodies-(2*binariesp)+1] = -1
                    H[row+1,bodies] = 1
                    row = row + 2
                    binariesp = binaries[iteration]
                    n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration+1)
                end
            elseif binariesp == 2
                if binaries[iteration] == 1
                    H[row,1:row-1] .= -1
                    H[row,row:bodies] .= 1
                    row = row + 1
                    binariesp = binaries[iteration]
                    n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration+1)
                elseif binaries[iteration] == 2
                    H[row,1:row-1] .= -1
                    H[row,row:bodies-2*(binariesp-1)] .= 1
                    H[row+1,bodies-2*(binariesp-1)+1] = -1
                    H[row+1,bodies] = 1
                    row = row + 2
                    binariesp = binaries[iteration]
                    n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration+1)
                end
            elseif binariesp == 3
                if binaries[iteration] == 2
                    H[row,1:Int64(bodies/2)-1] .= -1
                    H[row,Int64(bodies/2):2*Int64(bodies/3)] .= 1
                    H[row+1,2*Int64(bodies/3)+1:bodies] .= -1
                    H[row+1,bodies+1] = 1
                    row = row + 2
                    binariesp = binaries[iteration]
                    bodies = bodies + 1
                    n_level!(H, n_bodies, binaries, bodies, binariesp, row, iteration+1)
                end
            end
        end
    end
end

function symmetric!(H, n_bodies, binaries, bodies, binariesp, row, iteration)

    j = 0
    while row < n_bodies - 1
        H[row, 1 + j:2 + j] .= -1
        H[row, 3 + j:4 + j] .= 1
        j += 4
        row += 1
    end

    H[row, 1:binariesp] .= -1
    H[row, binariesp + 1:n_bodies] .= 1
end

function nested_symmetric!(H, n_bodies, binaries, bodies, binariesp, row, iteration)
    half = Int(n_bodies/2)
    j = 1
    while row < n_bodies - 1
        H[row    , 1:2 + j] .= -1
        H[row    , 3 + j] = 1
        H[row + 1, half + 1:half + 2 + j] .= -1
        H[row + 1, half + 3 + j] = 1
        j += 1
        row += 2
    end

    H[row, 1:half] .= -1
    H[row, half + 1:n_bodies] .= 1

end