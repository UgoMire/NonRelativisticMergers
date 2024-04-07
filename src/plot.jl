function plot_euler(sol, gd)
    u0 = sol.u[1]
    tspan = (sol.t[1], sol.t[end])

    fig = Figure(size = (1200, 400))

    sg = SliderGrid(fig[2, 1], (label = "t", range = range(tspan[1], tspan[2], 100)))
    tlift = sg.sliders[1].value

    ax1 = Axis(fig[1, 1])
    lines!(ax1, gd.xl, u0[1, :])
    lines!(ax1, gd.xl, (@lift sol($tlift)[1, :]))

    ax2 = Axis(fig[1, 2])
    lines!(ax2, gd.xl, u0[2, :])
    lines!(ax2, gd.xl, (@lift sol($tlift)[2, :]))

    ax3 = Axis(fig[1, 3])
    lines!(ax3, gd.xl, u0[3, :])
    lines!(ax3, gd.xl, (@lift sol($tlift)[3, :]))

    on(tlift) do val
        autolimits!(ax1)
        autolimits!(ax2)
        autolimits!(ax3)
    end

    return fig
end