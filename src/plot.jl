function plot_euler(sol, gd)
    u0 = sol.u[1]
    tspan = (sol.t[1], sol.t[end])

    fig = Figure(; size = (1200, 400))

    sg = SliderGrid(fig[2, 1], (label = "t", range = range(tspan[1], tspan[2], 100)))
    tlift = sg.sliders[1].value

    ax1 = Axis(fig[1, 1]; title = "density", xlabel = "x")
    lines!(ax1, gd.xl, u0[1, :])
    lines!(ax1, gd.xl, (@lift sol($tlift)[1, :]))

    ax2 = Axis(fig[1, 2]; title = "velocity", xlabel = "x")
    lines!(ax2, gd.xl, u0[2, :] ./ u0[1, :])
    lines!(ax2, gd.xl, (@lift sol($tlift)[2, :] ./ sol($tlift)[1, :]))

    ax3 = Axis(fig[1, 3]; title = "pressure", xlabel = "x")
    lines!(ax3, gd.xl, (5 / 3 - 1) * (u0[3, :] .- 1 / 2 * u0[1, :] .* u0[2, :] .^ 2))
    lines!(
        ax3,
        gd.xl,
        (@lift (5 / 3 - 1) *
               (sol($tlift)[3, :] .- 1 / 2 * sol($tlift)[1, :] .* sol($tlift)[2, :] .^ 2)),
    )

    on(tlift) do val
        autolimits!(ax1)
        autolimits!(ax2)
        autolimits!(ax3)
    end

    return fig
end

function plot_euler2d(prob, sol; interpolate = true)
    (; xl, yl) = prob.grid

    tspan = (sol.t[1], sol.t[end])

    fig = Figure(; size = (900, 900))

    sg = SliderGrid(fig[3, 1:2], (label = "t", range = range(tspan[1], tspan[2], 100)))
    tlift = sg.sliders[1].value

    ax1 = Axis3(fig[1, 1]; title = "density", xlabel = "x", ylabel = "y")
    surface!(ax1, xl, yl, (@lift sol($tlift)[1, :, :]); interpolate)
    # arrows!(
    #     ax1,
    #     xl,
    #     yl,
    #     (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]),
    #     (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]);
    #     arrowsize = 3,
    # )

    ax3 = Axis3(fig[2, 1]; title = "x velocity", xlabel = "x")
    surface!(ax3, xl, yl, (@lift sol($tlift)[2, :, :] ./ sol($tlift)[1, :, :]); interpolate)

    ax3 = Axis3(fig[2, 2]; title = "y velocity", xlabel = "x")
    surface!(ax3, xl, yl, (@lift sol($tlift)[3, :, :] ./ sol($tlift)[1, :, :]); interpolate)

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
        interpolate,
    )

    return fig
end

function plot_reconstruction(gd, model, u0)
    w_reconstruct = reconstruct(model.reconst, gd, u0)

    fig = Figure(; size = (1200, 400))

    ax1 = Axis(fig[1, 1]; title = "density")
    ax2 = Axis(fig[1, 2]; title = "velocity")
    ax3 = Axis(fig[1, 3]; title = "pressure")

    for i in 1:gd.Nx, (ax, j) in zip([ax1, ax2, ax3], 1:3)
        lines!(
            ax,
            [gd.xl[i] - gd.Δx / 2, gd.xl[i] + gd.Δx / 2],
            [u0[j, i], u0[j, i]];
            color = :black,
        )
        scatter!(ax, [gd.xl[i] - gd.Δx / 2], [w_reconstruct[j, i, 1]]; color = :red)
        scatter!(ax, [gd.xl[i] + gd.Δx / 2], [w_reconstruct[j, i, 2]]; color = :blue)
        lines!(
            ax,
            [gd.xl[i] - gd.Δx / 2, gd.xl[i] + gd.Δx / 2],
            [w_reconstruct[j, i, 1], w_reconstruct[j, i, 2]];
            color = :grey,
            linestyle = :dash,
        )
    end

    ax4 = Axis(fig[2, 1])
    ax5 = Axis(fig[2, 2])
    ax6 = Axis(fig[2, 3])

    for (ax, j) in zip([ax4, ax5, ax6], 1:3)
        F = solve_riemann(model.riemannsolver, gd, w_reconstruct)
        lines!(ax, gd.xl, F[j, :])
    end

    return fig
end