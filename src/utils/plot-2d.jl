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

        ax5 = Axis3(fig[1, 3]; title = "potential", xlabel = "x")
        surface!(ax5, xl, yl, (@lift solve_poisson(prob, sol($tlift)[1, :, :])))
    end

    return fig
end
