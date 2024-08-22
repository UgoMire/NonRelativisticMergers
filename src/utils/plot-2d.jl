# function plot_euler(
#         prob::FDProblem{Grid2D, Euler, <:Any, <:Any},
#         sol;
#         type = :heatmap
# )
#     (; xl, yl) = prob.grid

#     tspan = (sol.t[1], sol.t[end])

#     fig = Figure(; size = (900, 900))

#     if type == :heatmap
#         sg = SliderGrid(fig[2, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
#         tlift = sg.sliders[1].value

#         Na = 5 # Plot an arrow one every Na point. 
#         ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
#         heatmap!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))
#         arrows!(
#             ax1,
#             xl[1:Na:end],
#             yl[1:Na:end],
#             (@lift sol($tlift)[2, 1:Na:end, 1:Na:end] ./
#                    sol($tlift)[1, 1:Na:end, 1:Na:end]),
#             (@lift sol($tlift)[3, 1:Na:end, 1:Na:end] ./
#                    sol($tlift)[1, 1:Na:end, 1:Na:end]);
#             arrowsize = 10,
#             lengthscale = 0.5
#         )

#         ax4 = Axis(fig[1, 2]; title = "pressure", xlabel = "x")
#         heatmap!(
#             ax4,
#             xl,
#             yl,
#             (@lift (5 / 3 - 1) * (
#                 sol($tlift)[4, :, :] .-
#                 1 / 2 * sol($tlift)[1, :, :] .*
#                 (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
#             ));
#         )
#     elseif type == :surface
#         sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
#         tlift = sg.sliders[1].value

#         ax1 = Axis3(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
#         surface!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))

#         ax3 = Axis3(fig[2, 1]; title = "x velocity", xlabel = "x")
#         surface!(ax3, xl, yl, (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]))

#         ax3 = Axis3(fig[2, 2]; title = "y velocity", xlabel = "x")
#         surface!(ax3, xl, yl, (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]))

#         ax4 = Axis3(fig[1, 2]; title = "pressure", xlabel = "x")
#         surface!(
#             ax4,
#             xl,
#             yl,
#             (@lift (5 / 3 - 1) * (
#                 sol($tlift)[4, :, :] .-
#                 1 / 2 * sol($tlift)[1, :, :] .*
#                 (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
#             ));
#         )
#     end

#     return fig
# end

# function plot_euler(
#         prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
#         sol;
#         type = :heatmap
# )
#     (; xl, yl) = prob.grid

#     tspan = (sol.t[1], sol.t[end])

#     fig = Figure(; size = (900, 900))

#     if type == :heatmap
#         sg = SliderGrid(fig[2, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
#         tlift = sg.sliders[1].value

#         Na = 10 # Plot an arrow one every Na point. 
#         ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
#         contourf!(
#             ax1, xl, yl, (@lift sol($tlift)[1, :, :]); colormap = :inferno, levels = 20)
#         arrows!(
#             ax1,
#             xl[1:Na:end],
#             yl[1:Na:end],
#             (@lift sol($tlift)[2, 1:Na:end, 1:Na:end] ./
#                    sol($tlift)[1, 1:Na:end, 1:Na:end]),
#             (@lift sol($tlift)[3, 1:Na:end, 1:Na:end] ./
#                    sol($tlift)[1, 1:Na:end, 1:Na:end]);
#             arrowsize = 10,
#             lengthscale = 0.2,
#             color = (:grey, 0.5)
#         )

#         ax4 = Axis(fig[1, 2]; title = "pressure", xlabel = "x")
#         heatmap!(
#             ax4,
#             xl,
#             yl,
#             (@lift (5 / 3 - 1) * (
#                 sol($tlift)[4, :, :] .-
#                 1 / 2 * sol($tlift)[1, :, :] .*
#                 (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
#             ));
#             colormap = :inferno
#         )
#     elseif type == :surface
#         sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
#         tlift = sg.sliders[1].value

#         ax1 = Axis3(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
#         surface!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))

#         ax3 = Axis3(fig[2, 1]; title = "x velocity", xlabel = "x")
#         surface!(ax3, xl, yl, (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]))

#         ax3 = Axis3(fig[2, 2]; title = "y velocity", xlabel = "x")
#         surface!(ax3, xl, yl, (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]))

#         ax4 = Axis3(fig[1, 2]; title = "pressure", xlabel = "x")
#         surface!(
#             ax4,
#             xl,
#             yl,
#             (@lift (5 / 3 - 1) * (
#                 sol($tlift)[4, :, :] .-
#                 1 / 2 * sol($tlift)[1, :, :] .*
#                 (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
#             ));
#         )

#         ax5 = Axis3(fig[1, 3]; title = "potential", xlabel = "x")
#         surface!(ax5, xl, yl, (@lift solve_poisson(prob, sol($tlift)[1, :, :])))
#     end

#     return fig
# end

function add_main_layout!(fig, prob, sol; cmap = :lipari, Na = 10)
    (; xmin, xmax, Nx, xl, ymin, ymax, Ny, yl) = prob.grid

    tickformat = labels -> [rpad((Printf.@sprintf "%.0f" label), 4) for label in labels]

    tlift = Observable{Float64}(0)

    ρlift = Observable(zeros(Nx, Ny))
    Plift = Observable(zeros(Nx, Ny))
    ϕlift = Observable(zeros(Nx, Ny))
    elift = Observable(zeros(Nx, Ny))
    vxlift = Observable(zeros(length(1:Na:Nx), length(1:Na:Ny)))
    vylift = Observable(zeros(length(1:Na:Nx), length(1:Na:Ny)))

    for (i, j, quantity, title) in [(1, 1, ρlift, "density"), (1, 2, Plift, "pressure"),
        (2, 1, ϕlift, "potential"), (2, 2, elift, "internal energy")]
        ax = Axis(fig[i, j][1, 1]; title, xlabel = "x", ylabel = "y")
        limits!(ax, (xmin, xmax), (ymin, ymax))

        hm = heatmap!(ax, xl, yl, quantity; colormap = cmap, interpolate = true)
        arrows!(ax, xl[1:Na:end], yl[1:Na:end], vxlift, vylift;
            arrowsize = 6, lengthscale = 0.08, color = (:grey, 0.1), normalize = true)

        Colorbar(fig[i, j][1, 2], hm; tickformat)
    end

    on(tlift) do t
        ρ, vx, vy, P = get_primitive_variables(prob, sol(t))
        ϕ = solve_poisson(prob, ρ)

        ρlift[] = ρ
        Plift[] = P
        ϕlift[] = ϕ
        elift[] = P ./ ρ ./ (prob.model.γ - 1)
        vxlift[] = vx[1:Na:end, 1:Na:end]
        vylift[] = vy[1:Na:end, 1:Na:end]
    end

    return tlift
end

function plot_euler(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sol;
        cmap = :lipari,
        Na = 10,
        trange = nothing
)
    if isnothing(trange)
        trange = range(sol.t[1] + 0.001, sol.t[end], 100)
    end

    fig = Figure(; size = (950, 900))

    tlift = add_main_layout!(fig, prob, sol; cmap, Na)

    sg = SliderGrid(fig[3, 1:2], (label = "t", range = trange))

    on(sg.sliders[1].value) do t
        tlift[] = t
    end

    return fig
end

function record_euler(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sol,
        filename;
        cmap = :lipari,
        Na = 10,
        trange = nothing
)
    if isnothing(trange)
        trange = range(sol.t[1] + 0.001, sol.t[end], 100)
    end

    fig = Figure(; size = (950, 900), pt_per_unit = 0.1)

    tlift = add_main_layout!(fig, prob, sol; cmap, Na)

    record(fig, filename, trange; framerate = 30) do t
        tlift[] = t
    end
end
