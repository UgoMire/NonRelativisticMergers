function plot_euler(
        prob::FDProblem{Grid2D, Euler, <:Any, <:Any},
        sol;
        type = :heatmap
)
    (; xl, yl) = prob.grid

    tspan = (sol.t[1], sol.t[end])

    fig = Figure(; size = (900, 900))

    if type == :heatmap
        sg = SliderGrid(fig[2, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
        tlift = sg.sliders[1].value

        Na = 5 # Plot an arrow one every Na point. 
        ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
        heatmap!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))
        arrows!(
            ax1,
            xl[1:Na:end],
            yl[1:Na:end],
            (@lift sol($tlift)[2, 1:Na:end, 1:Na:end] ./
                   sol($tlift)[1, 1:Na:end, 1:Na:end]),
            (@lift sol($tlift)[3, 1:Na:end, 1:Na:end] ./
                   sol($tlift)[1, 1:Na:end, 1:Na:end]);
            arrowsize = 10,
            lengthscale = 0.5
        )

        ax4 = Axis(fig[1, 2]; title = "pressure", xlabel = "x")
        heatmap!(
            ax4,
            xl,
            yl,
            (@lift (5 / 3 - 1) * (
                sol($tlift)[4, :, :] .-
                1 / 2 * sol($tlift)[1, :, :] .*
                (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
            ));
        )
    elseif type == :surface
        sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
        tlift = sg.sliders[1].value

        ax1 = Axis3(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
        surface!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))

        ax3 = Axis3(fig[2, 1]; title = "x velocity", xlabel = "x")
        surface!(ax3, xl, yl, (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]))

        ax3 = Axis3(fig[2, 2]; title = "y velocity", xlabel = "x")
        surface!(ax3, xl, yl, (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]))

        ax4 = Axis3(fig[1, 2]; title = "pressure", xlabel = "x")
        surface!(
            ax4,
            xl,
            yl,
            (@lift (5 / 3 - 1) * (
                sol($tlift)[4, :, :] .-
                1 / 2 * sol($tlift)[1, :, :] .*
                (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
            ));
        )
    end

    return fig
end

function plot_euler(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sol;
        type = :heatmap
)
    (; xl, yl) = prob.grid

    tspan = (sol.t[1], sol.t[end])

    fig = Figure(; size = (900, 900))

    if type == :heatmap
        sg = SliderGrid(fig[2, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
        tlift = sg.sliders[1].value

        Na = 10 # Plot an arrow one every Na point. 
        ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
        contourf!(
            ax1, xl, yl, (@lift sol($tlift)[1, :, :]); colormap = :inferno, levels = 20)
        arrows!(
            ax1,
            xl[1:Na:end],
            yl[1:Na:end],
            (@lift sol($tlift)[2, 1:Na:end, 1:Na:end] ./
                   sol($tlift)[1, 1:Na:end, 1:Na:end]),
            (@lift sol($tlift)[3, 1:Na:end, 1:Na:end] ./
                   sol($tlift)[1, 1:Na:end, 1:Na:end]);
            arrowsize = 10,
            lengthscale = 0.2,
            color = (:grey, 0.5)
        )

        ax4 = Axis(fig[1, 2]; title = "pressure", xlabel = "x")
        heatmap!(
            ax4,
            xl,
            yl,
            (@lift (5 / 3 - 1) * (
                sol($tlift)[4, :, :] .-
                1 / 2 * sol($tlift)[1, :, :] .*
                (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
            ));
            colormap = :inferno
        )
    elseif type == :surface
        sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
        tlift = sg.sliders[1].value

        ax1 = Axis3(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
        surface!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]))

        ax3 = Axis3(fig[2, 1]; title = "x velocity", xlabel = "x")
        surface!(ax3, xl, yl, (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]))

        ax3 = Axis3(fig[2, 2]; title = "y velocity", xlabel = "x")
        surface!(ax3, xl, yl, (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]))

        ax4 = Axis3(fig[1, 2]; title = "pressure", xlabel = "x")
        surface!(
            ax4,
            xl,
            yl,
            (@lift (5 / 3 - 1) * (
                sol($tlift)[4, :, :] .-
                1 / 2 * sol($tlift)[1, :, :] .*
                (sol($tlift)[2, :, :] .^ 2 + sol($tlift)[3, :, :] .^ 2)
            ));
        )

        ax5 = Axis3(fig[1, 3]; title = "potential", xlabel = "x")
        surface!(ax5, xl, yl, (@lift solve_poisson(prob, sol($tlift)[1, :, :])))
    end

    return fig
end

function plot_euler2(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sol;
        cmap = :lipari,
        crange = Makie.MakieCore.Automatic(),
        Na = 10)
    (; xl, Nx, yl, Ny) = prob.grid

    fig = Figure(; size = (900, 900))

    ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
    ax2 = Axis(fig[1, 2]; title = "pressure", xlabel = "x", ylabel = "y")
    ax3 = Axis(fig[2, 1]; title = "potential", xlabel = "x", ylabel = "y")
    ax4 = Axis(fig[2, 2]; title = "internal energy", xlabel = "x", ylabel = "y")
    sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(sol.t[1], sol.t[end], 100)))

    ρlift = Observable(zeros(Nx, Ny))
    Plift = Observable(zeros(Nx, Ny))
    ϕlift = Observable(zeros(Nx, Ny))
    elift = Observable(zeros(Nx, Ny))
    vxlift = Observable(zeros(length(1:Na:Nx), length(1:Na:Ny)))
    vylift = Observable(zeros(length(1:Na:Nx), length(1:Na:Ny)))

    heatmap!(ax1, xl, yl, ρlift; colormap = cmap, colorrange = crange, interpolate = true)
    heatmap!(ax2, xl, yl, Plift; colormap = cmap, interpolate = true)
    heatmap!(ax3, xl, yl, ϕlift; colormap = cmap, interpolate = true)
    heatmap!(ax4, xl, yl, elift; colormap = cmap, interpolate = true)

    for ax in [ax1, ax2, ax3, ax4]
        arrows!(
            ax, xl[1:Na:end], yl[1:Na:end], vxlift, vylift;
            arrowsize = 6, lengthscale = 0.08, color = (:grey, 0.1), normalize = true)
    end

    tlift = sg.sliders[1].value

    on(tlift) do t
        E = sol(t)[4, :, :]
        ρ, vx, vy, P = get_primitive_variables(prob, sol(t))
        ϕ = solve_poisson(prob, ρ)

        ρlift[] = ρ
        Plift[] = P
        ϕlift[] = ϕ
        elift[] = P ./ ρ ./ (prob.model.γ - 1)
        vxlift[] = vx[1:Na:end, 1:Na:end]
        vylift[] = vy[1:Na:end, 1:Na:end]
    end

    tlift[] = sol.t[1]

    return fig
end
