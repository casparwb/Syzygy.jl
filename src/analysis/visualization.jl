#using Plots: plot, @animate, plot!, scatter!, scatter
# using Plots
using ProgressMeter
using RecipesBase

@userplot EnergyPlot
@userplot OrbitPlot
@userplot ElementPlot
@userplot StructurePlot
@userplot DistancePlot
@userplot AngularMomentumPlot
@userplot ICPlot

# function animate_threebody(sol::TripleSolution, name; 
#                            max_diff=100, nframes=200, fps=20,
#                            save_params=true)

#     r₁, r₂, r₃ = u"AU".(sol.r₁), u"AU".(sol.r₂), u"AU".(sol.r₃)

#     max_1, max_2 = maximum.((abs.(r₁[:, 1:end÷10]), abs.(r₂[:, 1:end÷10])), dims=2)
#     max_3 = maximum(abs.(r₃[:, 1:end÷10]), dims=2)
#     max_x12 = max(max_1[1], max_2[1]).val * 2.0
#     max_y12 = max(max_1[2], max_2[2]).val * 2.0
#     max_z12 = max(max_1[3], max_2[3]).val * 2.0

#     max_x3 = max_3[1].val * 1.5#1.1
#     max_y3 = max_3[2].val * 1.5#1.1
#     max_z3 = max_3[3].val * 1.5#1.5

#     if all(iszero.((max_z12, max_z3)))
#         max12 = (max_x12, max_y12)
#         max3 = (max_x3, max_y3)
#         anim =  animate_threebody_2d(r₁, r₂, r₃, sol.t, max12, max3, max_diff, nframes)
#     else
#         max12 = (max_x12, max_y12, max_z12)
#         max3 = (max_x3, max_y3, max_z3)
#         anim =  animate_threebody_3d(r₁, r₂, r₃, sol.t, max12, max3, max_diff, nframes)
#     end

#     if occursin(".gif", name)
#         gif(anim, joinpath(@__DIR__, "..", "..", "figures", name), fps=fps)
#     else
#         name = endswith(name, ".mp4") ? name : name*".mp4"
#         mp4(anim, joinpath(@__DIR__, "..", "..", "figures", name), fps=fps)
#     end

#     if save_params
#         outfile = replace(replace(name, ".gif" => ""), ".mp4"  => "")
#         outfile = joinpath(@__DIR__, "..", "..", "data", outfile*".jld2")
#         if !isfile(outfile)
#             @info "Saving parameters to" relpath(outfile)
#             JLD2.save(outfile, "params", sol.initial_conditions)
#         end
#     end

#     nothing
# end

# function animate_threebody_3d(r₁, r₂, r₃, t, max12, max3, max_diff=100, nframes=200)
#     max_x12, max_y12, max_z12 = max12
#     max_x3, max_y3, max_z3 = max3

#     max_x = max(max_x12, max_x3)
#     max_y = max(max_y12, max_y3)
#     max_z = max(max_z12, max_z3)

#     N_step = max_diff/10 |> x -> round(Int, x)
#     r₁_coords = [eltype(r₁)[r₁[1,1]], eltype(r₁)[r₁[1,1]], eltype(r₁)[r₁[1,1]]]
#     r₂_coords = [eltype(r₂)[r₂[1,1]], eltype(r₂)[r₂[1,1]], eltype(r₂)[r₂[1,1]]]
#     r₃_coords = [eltype(r₃)[r₃[1,1]], eltype(r₃)[r₃[1,1]], eltype(r₃)[r₃[1,1]]]
#     nframes = min(length(t), nframes) |> Int
#     progress = Progress(nframes)

#     c = [1, 2, 3]

#     times = round.(u"yr", t, digits=0)
#     anim = @animate for i in eachindex(t)[1:end÷nframes:end]
#         startidx = i > max_diff ? i-max_diff : 1

#         # p = plot(xlims=(-max_y, max_y), 
#         #          ylims=(-max_y, max_y), 
#         #          zlims=(-max_z, max_z), xlabel="x", ylabel="y", zlabel="z")

#         p123 = plot(xlims=(-max_x, max_x), 
#                     ylims = (-max_y, max_y), 
#                     zlims = (-max_z, max_z), xlabel="x", ylabel="y", zlabel="z", 
#                     title="t = $(times[i])")
#         p12 = plot(xlims=(-max_x12, max_x12), 
#                     ylims = (-max_y12, max_y12), 
#                     zlims = (-max_z12, max_z12), xlabel="x", ylabel="y", zlabel="z")

#         # p3 = plot(xlims=(-max_x3, max_x3), 
#         #          ylims=(-max_y3, max_y3), 
#         #          zlims=(-max_z3, max_z3), xlabel="x", ylabel="y", zlabel="z")

#         r₁_coords[1][1] = r₁[1,i]
#         r₁_coords[2][1] = r₁[2,i]
#         r₁_coords[3][1] = r₁[3,i]

#         r₂_coords[1][1] = r₂[1,i]
#         r₂_coords[2][1] = r₂[2,i]
#         r₂_coords[3][1] = r₂[3,i]

#         r₃_coords[1][1] = r₃[1,i]
#         r₃_coords[2][1] = r₃[2,i]
#         r₃_coords[3][1] = r₃[3,i]

#         scatter!(p12, r₁_coords..., c=c[1], markersize=2, label="Primary")
#         scatter!(p12, r₂_coords..., c=c[2], markersize=2, label="Secondary")
#         # scatter!(p3,  r₃_coords..., c=c[3], markersize=2, label="Tertiary")
#         # @show size(xs_1)
#         plot!(p12, r₁[1,startidx:i], r₁[2,startidx:i], r₁[3,startidx:i], 
#                    c=c[1], label=false, alpha=0.8)
#         plot!(p12, r₂[1,startidx:i], r₂[2,startidx:i], r₂[3,startidx:i], 
#                    c=c[2], label=false, alpha=0.8)
#         # plot!(p3,  r₃[1,startidx:i], r₃[2,startidx:i], r₃[3,startidx:i], 
#         #    c=c[3], label=false, alpha=0.8)


#         scatter!(p123, r₁_coords..., c=c[1], markersize=2, label="Primary")
#         scatter!(p123, r₂_coords..., c=c[2], markersize=2, label="Secondary")
#         scatter!(p123, r₃_coords..., c=c[3], markersize=2, label="Tertiary")
#         # @show size(xs_1)
#         plot!(p123, r₁[1,startidx:i], r₁[2,startidx:i], r₁[3,startidx:i], 
#                     c=c[1], label=false, alpha=0.8)
#         plot!(p123, r₂[1,startidx:i], r₂[2,startidx:i], r₂[3,startidx:i], 
#                     c=c[2], label=false, alpha=0.8)
#         plot!(p123, r₃[1,startidx:i], r₃[2,startidx:i], r₃[3,startidx:i], 
#                     c=c[3], label=false, alpha=0.8)


#         next!(progress)
#         plot(p123, p12, layout=(1, 2), size=(800, 400))
#     end
#     return anim
# end

# function animate_threebody_2d(r₁, r₂, r₃, t, max12, max3, max_diff=100, nframes=200)
#     max_x12, max_y12 = max12
#     max_x3, max_y3 = max3

#     max_x = max(max_x12, max_x3)
#     max_y = max(max_y12, max_y3)

#     c = [1, 2, 3]

#     r₁_coords = [eltype(r₁)[r₁[1,1]], eltype(r₁)[r₁[1,1]]]
#     r₂_coords = [eltype(r₂)[r₂[1,1]], eltype(r₂)[r₂[1,1]]]
#     r₃_coords = [eltype(r₃)[r₃[1,1]], eltype(r₃)[r₃[1,1]]]
#     nframes   = min(length(t), nframes) |> Int
#     progress  = Progress(nframes)

#     times = round.(u"yr", t, digits=0)
#     anim = @animate for i in eachindex(t)[1:end÷nframes:end]
#         startidx = i > max_diff ? i-max_diff : 1


#         p123 = plot(xlims=(-max_x, max_x), 
#                     ylims=(-max_y, max_y), 
#                     xlabel="x", ylabel="y", zlabel="z",
#                     title="t = $(times[i])")

#         p12 = plot(xlims=(-max_x12, max_x12), 
#                    ylims=(-max_y12, max_y12),
#                    xlabel="x", ylabel="y", zlabel="z")

#         r₁_coords[1][1] = r₁[1,i]
#         r₁_coords[2][1] = r₁[2,i]

#         r₂_coords[1][1] = r₂[1,i]
#         r₂_coords[2][1] = r₂[2,i]

#         r₃_coords[1][1] = r₃[1,i]
#         r₃_coords[2][1] = r₃[2,i]

#         scatter!(p123, r₁_coords..., c=c[1], markersize=2, label="Primary")
#         scatter!(p123, r₂_coords..., c=c[2], markersize=2, label="Secondary")
#         scatter!(p123, r₃_coords..., c=c[3], markersize=2, label="Tertiary")

#         scatter!(p12, r₁_coords..., c=c[1], markersize=2, label="Primary")
#         scatter!(p12, r₂_coords..., c=c[2], markersize=2, label="Secondary")
#         # @show size(xs_1)
#         plot!(p123, r₁[1,startidx:i], r₁[2,startidx:i], 
#                     c=c[1], label=false, alpha=0.8)
#         plot!(p123, r₂[1,startidx:i], r₂[2,startidx:i], 
#                     c=c[2], label=false, alpha=0.8)
#         plot!(p123, r₃[1,startidx:i], r₃[2,startidx:i], 
#                     c=c[3], label=false, alpha=0.8)

#         plot!(p12, r₁[1,startidx:i], r₁[2,startidx:i], 
#                     c=c[1], label=false, alpha=0.8)
#         plot!(p12, r₂[1,startidx:i], r₂[2,startidx:i], 
#                     c=c[2], label=false, alpha=0.8)

#         next!(progress)
#         plot(p123, p12, layout=(1, 2), size=(800, 400))
#     end

#     return anim
# end


@recipe function plot(ic::MultiBodySystem; 
                      bodies=collect(keys(ic.particles)), dims=[1,2])

    particles = [ic.particles[body] for body in sort(bodies)]

    aspect_ratio --> 1

    masses = [p.mass for p in particles]
    m_0_1 = (masses .- minimum(masses))/(maximum(masses) - minimum(masses))
    # mratio = masses ./ maximum(masses) .* 5
    markersizes = all(isnan, m_0_1) ? repeat([5], length(masses)) : (m_0_1 .* 5) .+ 8

    ylims --> :auto
    xlims --> :auto

    for particle in particles
        # @show particle
        @series begin
            seriestype --> :scatter
            markersize := markersizes[particle.key.i]
            label --> "Particle $(particle.key.i) ($(particle.mass))"
            data = Tuple([[u"AU"(r)] for r in particle.position[dims]])
            # rs = r.(particle.)
            data
        end
    end
end 

@recipe function plot(res::SimulationResult; bodies="all", dims=[1, 2, 3],
                     tspan=nothing, step=1, ref_frame="com")

    # labels = ["Primary", "Secondary"]
    if bodies isa String 
        bodies = keys(res.simulation.ic.particles) |> collect |> sort
    end

    @assert length(bodies) <= res.simulation.ic.n "Number of bodies to plot is greater than bodies in system."
    
    if all(isone, res.simulation.ic.hierarchy[2:end])
        labels = hierarchy_labels#[bodies]
    else
        labels = ["Partcle $i" for i in bodies]
    end

    time = res.solution.t  .* u"s"
    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]
    
    dimlabels = ["x", "y", "z"]
    xlabel := dimlabels[dims[1]]
    ylabel := dimlabels[dims[2]]
    zlabel := length(dims) == 3 ? dimlabels[dims[3]] : nothing
    # aspect_ratio --> 1

    r = [res.solution.u[i].x[2][:,bodies] .* u"m" for i ∈ indices]
    masses = get_masses(res.simulation.system)[bodies] .* u"kg"

    if ref_frame == "com"
        com = zeros(eltype(r[1]), 3, length(indices))
        for i in eachindex(indices)
            com[:,i] .= centre_of_mass(r[i], masses)
        end
        shift_coord = com
    elseif !isnothing(tryparse(Int, ref_frame))
        ref_frame_idx = parse(Int, ref_frame)
        shift_coord = reduce(hcat, [rr[:,ref_frame_idx] for rr in r])
        
    else
        shift_coord = fill(0.0u"m", (3, length(indices)))
    end


    for idx in eachindex(bodies)
        data = [r[i][dims,idx] .- shift_coord[dims,i] .|> u"AU" for i ∈ indices]
        data = Tuple([[d[dim] for d in data] for dim in dims])
        @series begin 
            label --> labels[idx]
            data
        end
    end

end

function getsibling(system::MultiBodySystem, particle::Particle)
    system[particle.sibling]
end

function getsibling(system::MultiBodySystem, binary::Binary)
    @assert binary.parent.i != -1 "Root binary has no sibling."
    parent = system[binary.parent]
    return [child for child in parent.children if child.key != binary.key][1]
end

@recipe function f(eO::OrbitPlot; bodies="all", dims=[1, 2, 3],
                     tspan=nothing, step=1, ref_frame="com", axis_units=u"AU")

    sol = eO.args[1]
    @assert sol isa MultiBodySolution "First argument must be a `MultiBodySolution`. Got $(typeof(sol))"
    # labels = ["Primary", "Secondary"]
    if bodies isa String 
        bodies = keys(sol.initial_conditions.particles) |> collect |> sort
    end

    @assert length(bodies) <= sol.initial_conditions.n "Number of bodies to plot is greater than bodies in system."
    
    if all(isone, sol.initial_conditions.hierarchy[2:end])
        labels = hierarchy_labels#[bodies]
    else
        labels = ["Partcle $i" for i in bodies]
    end

    time = sol.t#  .|> u"kyr"
    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]
    dimlabels = ["x", "y", "z"]
    xlabel := dimlabels[dims[1]]
    ylabel := dimlabels[dims[2]]
    zlabel := length(dims) == 3 ? dimlabels[dims[3]] : nothing
    aspect_ratio --> 1

    # data = 
    if ref_frame == "com"
        shift_coord = centre_of_mass(sol, bodies, tspan=tspan)#[:,indices]
    elseif !isnothing(tryparse(Int, ref_frame))
        ref_frame_idx = parse(Int, ref_frame)
        shift_coord = sol.r[:,ref_frame_idx, indices]
    else
        shift_coord = fill(0.0u"m", (3, length(indices)))
    end

    for idx in bodies
        # data = Tuple([sol.r[dim,idx,indices] .- shift_coord[dim,:]  .|> u"AU" for dim in dims])
        # @show size(sol.r[dim=1,particle=idx,time=tspan[1]..tspan[2]]), size(shift_coord[1,:])
        # data = Tuple([sol.r[dim=dim,particle=idx,time=tspan[1]..tspan[2]] .- shift_coord[dim,:]  .|> axis_units for dim in dims])
        data = Tuple([sol.r[dim=dim,particle=idx,time=indices] .- shift_coord[dim,1:step:end]  .|> axis_units for dim in dims])
        @series begin 
            # c --> idx
            label --> labels[idx]
            # sol.r[1,idx,indices], .- shift_coord[1,:] .|> u"AU", 
            # sol.r[2,idx,indices], .- shift_coord[2,:] .|> u"AU" 
            # sol.r[3,idx,indices] .- shift_coord[3,:] .|> u"AU"
            data
        end
    end

end

@recipe function f(eP::EnergyPlot; tspan=nothing, step=1, t_unit=u"kyr")

    sol = eP.args[1]
    time = sol.t  .|> t_unit

    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]

    t = time[indices]
    Etot = sol.quantities.E[indices] ./ first(sol.quantities.E)
    Ekin = sol.quantities.K[indices] ./ first(sol.quantities.K)
    Epot = sol.quantities.U[indices] ./ first(sol.quantities.U)

    if get(plotattributes, :yscale, :log10) == :log10
        Etot = 1 .- Etot .|> abs
        # Ekin = 1 .- Ekin .|> abs
        # Epot = 1 .- Epot .|> abs

        Etot[1] += Etot[2]
        # Ekin[1] += Ekin[2]
        # Epot[1] += Epot[2]

        ylabel = "1 - E/E₀"
    else
        ylabel = "E/E₀"
    end

    # set up subplots
    legend := :outerright
    tickfontsize --> 10
    labelfontsize --> 14
    # link   := :both
    layout := (3, 1)
    # label := false
    size --> (800, 600)

    lowerlim = minimum(minimum.([Ekin, Epot]))
    upperlim = maximum(maximum.([Ekin, Epot]))
    Elims    = extrema(Etot)
    
    # Subplot 1 (potential energy)
    @series begin
        seriestype --> :line
        yscale := :identity
        subplot := 1
        #ylims --> (0.9*lowerlim, 1.1*upperlim)
        alpha --> 0.5
        # color --> :green
        ylabel := "E/E₀"
        title --> "Kinetic"
        xticks := nothing
        label := false
        ustrip(t), Ekin
    end

    # Subplot 2 (total energy)
    @series begin
        seriestype --> :line
        yscale --> get(plotattributes, :yscale, :log10)
        subplot := 2
        linewidth --> 3
        ylabel --> ylabel
        label --> false
        #ylims --> (0.99Elims[1], 1.01Elims[2])
        title --> "Total"
        xticks := nothing
        ustrip(t), Etot
    end

    # Subplot 1 (potential energy)
    @series begin
        seriestype --> :line
        yscale := :identity
        subplot := 3
        #ylims --> (-1.1*upperlim, 0.9*lowerlim)
        alpha --> 0.5
        # color --> :red
        ylabel := "E/E₀"
        xlabel --> "Time"
        title --> "Potential"
        label := false
        t, -Epot
    end
end


@recipe function f(eP::ElementPlot; loge=false, tspan=nothing, step=1, t_unit=u"kyr")

    sol = eP.args[1]
    # @assert sol isa MultiBodySolution
    
    time = sol.t
    tspan = isnothing(tspan) ? extrema(time) : tspan
    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]
    t = t_unit.(time[indices])

    system = sol.initial_conditions
    n_binaries = length(system.binaries)
    n_particles = length(system.particles)

    # label := false
    # layout := @layout [a b; c d; e]
    layout --> (3, 1)
    size --> (900, 700)

    legend := :outerright
    legendfontsize --> 12
    titlefontsize --> 20
    tickfontsize --> 10


    # Semi-major axis plot
    for bin in 1:n_binaries

        @series begin
            seriestype --> :line
            label := "Binary $bin"
            subplot := 1
            ylabel --> "a / a₀"
            xticks --> nothing
            title --> "Semi-major axis"

            ustrip(t), sol.elements[bin].a ./ first(sol.elements[bin].a)

        end
    end

    # Eccentricity plot
    for bin in 1:n_binaries

        @series begin
            seriestype --> :line
            label := false#"Binary $bin"
            subplot := 2

            yscale --> ifelse(loge, :log10, :identity)
            # @show plotattributes[:yscale]
            edata = ifelse(loge, 1 .- sol.elements[bin].e, sol.elements[bin].e ./ first(sol.elements[bin].e))
            ylabel := ifelse(loge,"log (1 - e)", "e / e₀")
            xticks --> nothing
            title --> "Eccentricity"

            ustrip(t), edata

        end
    end


    # # Mutual inclination plot
    for bin in 1:n_binaries

        @series begin
            seriestype --> :line
            label --> false#"Binary $bin"
            subplot := 3

            ylabel := "i"
            # @show plotattributes[:yscale]
            # xticks --> nothing
            title --> "Inclination"

            t, sol.elements[bin].i

        end
    end

    # link := :x


end


@recipe function f(eD::DistancePlot; tspan=nothing, step=1)

    sol = eD.args[1]
    time = sol.t  .|> u"Myr"

    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]

    t = time[indices]

    r1 = sol.r[:,1,indices] .|> u"AU"
    r2 = sol.r[:,2,indices] .|> u"AU"
    r3 = sol.r[:,3,indices] .|> u"AU"
    
    d12 = distances(r2 .- r1)
    d312 = distances(r3 .- centre_of_mass(sol, (1, 2))[:,indices] .|> u"AU")

    layout --> (2, 1)

    @series begin
        subplot := 1
        title --> "Inner binary"
        xticks --> nothing
        label --> false
        ylabel --> "Distanc"
        ustrip(t), d12 #./ first(d12)
    end

    @series begin
        subplot := 2
        title --> "Outer binary"
        xlabel --> "Time"
        label := false
        ylabel --> "Distance"
        t, d312 #./ first(d312)
    end

    # p1 = plot(ustrip(t), d12 ./ first(d12), 
    #           title="Inner binary distance", xticks=nothing, 
    #           label=false, ylabel="d/d₀")
    # p2 = plot(t, d312 ./ first(d312), 
    #           title="Outer binary distance", label=false, 
    #           xlabel="Time", ylabel="d/d₀")
end


@recipe function f(eS::StructurePlot; tspan=nothing, step=1)

    sol = eS.args[1]
    time = sol.t  .|> u"Myr"

    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]

    t = time[indices]
    
    masses = sol.structure.m[:, indices] ./ sol.structure.m[:, 1] 
    radii = sol.structure.R[:, indices] ./ sol.structure.R[:, 1] 

    layout --> (2, 1)
    legend --> :bottomleft

    @series begin
        subplot := 1
        xticks --> nothing
        ylabel --> "m/m₀"
        label --> ["Primary" "Secondary" "Tertiary"]
        ylims --> (0.999*minimum(masses), maximum(masses)*1.001)
        linewidth --> 4
        ustrip(t), masses'
    end

    @series begin
        subplot := 2
        label --> false
        ylabel --> "R/R₀"
        linewidth --> 4
        ylims --> (0.9*minimum(radii), maximum(radii)*1.1)
        t, radii'
    end
end


@recipe function f(eS::AngularMomentumPlot; total=false, tspan=nothing, step=1)

    sol = eS.args[1]
    time = sol.t  .|> u"Myr"

    tspan = isnothing(tspan) ? extrema(time) : tspan

    indices = (argmin(abs.(time .- tspan[1])), argmin(abs.(time .- tspan[2])))
    indices = indices[1]:step:indices[2]

    t = time[indices]

    h12 = sol.quantities.h[:,1,indices]'
    h3  = sol.quantities.h[:,2,indices]'

    # if total
    #     h12 = norm.(eachrow(h12)) ./ norm(h12[1,:])
    #     h3 = norm.(eachrow(h3)) ./ norm(h3[1,:])
    # end

    nbodies = sol.initial_conditions.n
    total_angular_momentum = zeros(typeof(1.0u"m^2/s"), 3, length(t))
    for i in eachindex(t)
        for body in 1:nbodies
            total_angular_momentum[:, i] = sol.r[:,body,i] × (sol.v[:,body,i])
        end
    end

    htot = sum(total_angular_momentum, dims=1)[1,:] #./ 3
    @series begin
        t, htot ./ htot[1] .- 1
    end
    # layout --> (2, 1)
    # legend --> :bottomleft

    # @series begin
    #     subplot := 1
    #     xticks --> nothing
    #     ylabel --> "h₁₂"
    #     label --> total ? "Total" : ["Primary" "Secondary" "Tertiary"]
    #     # ylims --> (0.999*minimum(masses), maximum(masses)*1.001)
    #     # linewidth --> 4
    #     ustrip(t), h12
    # end

    # @series begin
    #     subplot := 2
    #     label --> false
    #     ylabel --> "h₃"
    #     # linewidth --> 4
    #     # ylims --> (0.9*minimum(radii), maximum(radii)*1.1)
    #     t, h3
    # end
end

# anim = @animate for i in eachindex(sol.t)[1:100:end]
#     sidx = i > 1000 ? i-1000 : 1
#     # p = plot(xlims=(-1.2e4, 1.2e4), ylims=(-1.2e4, 2.2e4), zlims=(-3.0e4, 1.3e4), legend=:outerright, aspect_ratio=1)
#     p = plot(xlims=(-0.2, 0.08), ylims=(-0.05, 0.1), zlims=(-0.15, 0.7), 
#              legend=:outerright, aspect_ratio=1)
#     scatter!(p, [sol.r[1,1,i] |> u"pc"], 
#                 [sol.r[2,1,i] |> u"pc"], 
#                 [sol.r[3,1,i] |> u"pc"], 
#                 label="Primary", c=:blue)
#     scatter!(p, [sol.r[1,2,i] |> u"pc"], 
#                 [sol.r[2,2,i] |> u"pc"], 
#                 [sol.r[3,2,i] |> u"pc"], 
#                 label="Secondary", c=:yellow)
#     scatter!(p, [sol.r[1,3,i] |> u"pc"], 
#                 [sol.r[2,3,i] |> u"pc"], 
#                 [sol.r[3,3,i] |> u"pc"], 
#                 label="Tertiary", c=:red)

#     plot!(p, sol.r[1,1,sidx:i] .|> u"pc", 
#              sol.r[2,1,sidx:i] .|> u"pc", 
#              sol.r[3,1,sidx:i] .|> u"pc", 
#              label=false, c=:blue, alpha=0.2)
#     plot!(p, sol.r[1,2,sidx:i] .|> u"pc", 
#              sol.r[2,2,sidx:i] .|> u"pc", 
#              sol.r[3,2,sidx:i] .|> u"pc", 
#              label=false, c=:yellow, alpha=0.2)
#     plot!(p, sol.r[1,3,sidx:i] .|> u"pc", 
#              sol.r[2,3,sidx:i] .|> u"pc", 
#              sol.r[3,3,sidx:i] .|> u"pc", 
#              label=false, c=:red, alpha=0.2)
#     println(i)
# end

# let sol = sol2
#     anim = @animate for i in eachindex(sol.t)[end-3000:2:end]
#         sidx = i > 1000 ? i-1000 : 1
#         # p = plot(xlims=(-1.2e4, 1.2e4), ylims=(-1.2e4, 2.2e4), zlims=(-3.0e4, 1.3e4), legend=:outerright, aspect_ratio=1)
#         p = plot(xlims=(-0.05, 0.08), ylims=(-0.06, 0.12), zlims=(-0.3, 0.2), 
#                  legend=:outerright, aspect_ratio=1)
#         scatter!(p, [sol.r[1,1,i] |> u"pc"], 
#                     [sol.r[2,1,i] |> u"pc"], 
#                     label="Primary", c=:blue)
#         scatter!(p, [sol.r[1,2,i] |> u"pc"], 
#                     [sol.r[2,2,i] |> u"pc"], 
#                     label="Secondary", c=:yellow)
#         scatter!(p, [sol.r[1,3,i] |> u"pc"], 
#                     [sol.r[2,3,i] |> u"pc"], 
#                     label="Tertiary", c=:red)

#         plot!(p, sol.r[1,1,sidx:i] .|> u"pc", 
#                  sol.r[2,1,sidx:i] .|> u"pc", 
#                  label=false, c=:blue, alpha=0.2)
#         plot!(p, sol.r[1,2,sidx:i] .|> u"pc", 
#                  sol.r[2,2,sidx:i] .|> u"pc", 
#                  label=false, c=:yellow, alpha=0.2)
#         plot!(p, sol.r[1,3,sidx:i] .|> u"pc", 
#                  sol.r[2,3,sidx:i] .|> u"pc", 
#                  label=false, c=:red, alpha=0.2)
#         println(i)
#     end
#     gif(anim, "blackhole_final.mp4", fps=60)
# end

# let sol = sol2
#     com = MultiBodySimulator.centre_of_mass(sol, [1,2], tspan=extrema(sol.t[end-3000:end]))
#     k = 1
#     step = 1
#     anim = @animate for i in eachindex(sol.t)[end-3000:step:end]
#         sidx = length(sol.t)-3000#i > end-3000 ? i-1000 : 1
#         # p = plot(xlims=(-1.2e4, 1.2e4), ylims=(-1.2e4, 2.2e4), zlims=(-3.0e4, 1.3e4), legend=:outerright, aspect_ratio=1)
#         p = plot(xlims=(-0.002, 0.002), ylims=(-0.002, 0.0015), 
#                  legend=:outerright, aspect_ratio=1)
#         scatter!(p, [sol.r[1,1,i] - com[1,k] |> u"pc"], 
#                     [sol.r[2,1,i] - com[2,k] |> u"pc"], 
#                     label="Primary", c=:blue)
#         scatter!(p, [sol.r[1,2,i] - com[1,k] |> u"pc"], 
#                     [sol.r[2,2,i] - com[2,k] |> u"pc"], 
#                     label="Secondary", c=:yellow)


#         plot!(p, sol.r[1,1,sidx:i] .- com[1,1:k] .|> u"pc", 
#                  sol.r[2,1,sidx:i] .- com[2,1:k] .|> u"pc", 
#                  label=false, c=:blue, alpha=0.2)
#         plot!(p, sol.r[1,2,sidx:i] .- com[1,1:k] .|> u"pc", 
#                  sol.r[2,2,sidx:i] .- com[2,1:k] .|> u"pc", 
#                  label=false, c=:yellow, alpha=0.2)
#         k += step
#         println(i)
#     end
#     gif(anim, "blackhole_final_inner.mp4", fps=60)
# end


# let
#     com = MultiBodySimulator.centre_of_mass(sol, [1,2], tspan=extrema(sol.t[end-1000:10:end]))

#     k = 1
#     anim = @animate for i in eachindex(sol.t)[end-1000:end]
#         sidx = length(sol.t)-1000#i > 1000 ? i-1000 : 1
#         # p = plot(xlims=(-1.2e4, 1.2e4), ylims=(-1.2e4, 2.2e4), zlims=(-3.0e4, 1.3e4), legend=:outerright, aspect_ratio=1)
#         p = plot(xlims=(-300, 100), ylims=(-100, 200), zlims=(-300, 200), legend=:outerright, aspect_ratio=1)
#         scatter!(p, [sol.r[1,1,i] - com[1,k] |> u"AU"], 
#                     [sol.r[2,1,i] - com[2,k] |> u"AU"], 
#                     [sol.r[3,1,i] - com[3,k] |> u"AU"], 
#                     label="Primary", c=:blue)
#         scatter!(p, [sol.r[1,2,i] - com[1,k] |> u"AU"], 
#                     [sol.r[2,2,i] - com[2,k] |> u"AU"], 
#                     [sol.r[3,2,i] - com[3,k] |> u"AU"], 
#                     label="Secondary", c=:yellow)

#         plot!(p, sol.r[1,1,sidx:i] .- com[1,1:k] .|> u"AU", 
#                  sol.r[2,1,sidx:i] .- com[2,1:k] .|> u"AU", 
#                  sol.r[3,1,sidx:i] .- com[3,1:k] .|> u"AU", 
#                  label=false, c=:blue)
#         plot!(p, sol.r[1,2,sidx:i] .- com[1,1:k] .|> u"AU", 
#                  sol.r[2,2,sidx:i] .- com[2,1:k] .|> u"AU", 
#                  sol.r[3,2,sidx:i] .- com[3,1:k] .|> u"AU", 
#                  label=false, c=:yellow)
#         k += 1
#         println(i)
#     end
#     gif(anim, "inner_blackhole.mp4", fps=60)
# end


# anim = @animate for i in eachindex(sol.t)[1:10:end]
#     sidx = i > 1000 ? i-1000 : 1
#     p = plot(xlims=(-4.2, 4.2), ylims=(-4.0, 4.0), legend=false, aspect_ratio=1)

#     dims = [1,2]
#     cs = [:blue, :red, :green, :purple, :purple]
#     for pidx in eachindex(sol.initial_conditions.particles)
#             scatter!(p, [[sol.r[dim,pidx,i] |> u"AU"] for dim in dims]..., label=hierarchy_labels[pidx], c=cs[pidx])
#             plot!(p, [sol.r[dim,pidx,sidx:i] .|> u"AU" for dim in dims]..., label=false, c=cs[pidx])

#     end
#     println(i)
#     p
# end

# function animate(sol; dims=[1,2,3], ids=eachindex(sol.t), step=1, fps=60, filename="animation.mp4", args...)
           
#     com = Syzygy.centre_of_mass(sol, [1,2], tspan=extrema(sol.t[ids]))

#     k = 1
#     anim = @animate for i in ids#eachindex(sol.t)[ids]
#         sidx = 1
#         p = plot(legend=:outerright, aspect_ratio=1; args...)
#         scatter!(p, [[sol.r[dim,1,i] - com[dim,k] |> u"AU"] for dim in dims]..., 
#                     label="Primary", c=:blue)

#         scatter!(p, [[sol.r[dim,2,i] - com[dim,k] |> u"AU"] for dim in dims]...,
#                     label="Secondary", c=:green)
#         # scatter!(p, [sol.r[1,1,i] - com[1,k] |> u"AU"], 
#         #             [sol.r[2,1,i] - com[2,k] |> u"AU"], 
#         #             3 in dims ? [sol.r[3,1,i] - com[3,k] |> u"AU"] : nothing, 
#         #             label="Primary", c=:blue)
#         # scatter!(p, [sol.r[1,2,i] - com[1,k] |> u"AU"], 
#         #             [sol.r[2,2,i] - com[2,k] |> u"AU"], 
#         #             3 in dims ? [sol.r[3,2,i] - com[3,k] |> u"AU"] : nothing, 
#         #             label="Secondary", c=:green)
#         # scatter!(p, [sol.r[1,3,i] - com[1,k] |> u"AU"], 
#         #             [sol.r[2,3,i] - com[2,k] |> u"AU"], 
#         #             3 in dims ? [sol.r[3,2,i] - com[3,k] |> u"AU"] : nothing, 
#         #             label="Tertiary", c=:red)

#         # plot!(p, sol.r[1,1,sidx:10:i] .|> u"AU", 
#         #          sol.r[2,1,sidx:10:i] .|> u"AU", 
#         #          3 in dims ? sol.r[3,1,sidx:10:i] .|> u"AU" : nothing, 
#         #          label=false, c=:blue, alpha=0.5)
#         # plot!(p, sol.r[1,2,sidx:10:i] .|> u"AU", 
#         #          sol.r[2,2,sidx:10:i] .|> u"AU", 
#         #          3 in dims ? sol.r[3,2,sidx:10:i] .|> u"AU" : nothing, 
#         #          label=false, c=:green, alpha=0.5)
#         # plot!(p, sol.r[1,3,sidx:10:i] .|> u"AU", 
#         #          sol.r[2,3,sidx:10:i] .|> u"AU", 
#         #          3 in dims ? sol.r[3,3,sidx:10:i] .|> u"AU" : nothing, 
#         #          label=false, c=:red, alpha=0.5)
#         k += 1
#     end every step

#     outpath = filename
#     gif(anim, outpath, fps=fps)
# end

# function animate(sol; dims=[1,2,3], ids=eachindex(sol.t), step=1, fps=60, filename="animation.mp4", args...)
                  
#     com = Syzygy.centre_of_mass(sol, [1,2, 3], tspan=extrema(sol.t[ids]))

#     k = 1
#     anim = @animate for i in ids#eachindex(sol.t)[ids]
#         sidx = 1
#         p = plot(legend=:outerright, aspect_ratio=1; args...)
#         scatter!(p, [[sol.r[dim,1,i] - com[dim,k] |> u"AU"] for dim in dims]..., 
#                     label="Primary", c=:blue)

#         scatter!(p, [[sol.r[dim,2,i] - com[dim,k] |> u"AU"] for dim in dims]...,
#                     label="Secondary", c=:green)

#         scatter!(p, [[sol.r[dim,3,i] - com[dim,k] |> u"AU"] for dim in dims]...,
#                     label="Tertiary", c=:red)

#         k += 1
#     end every step

#     outpath = filename
#     gif(anim, outpath, fps=fps)
# end