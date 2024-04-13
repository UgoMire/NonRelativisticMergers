function flux(model, ρ, v, p)
    γ = model.γ

    return (ρ * v, ρ * v^2 + p, (p / (γ - 1) + 1 / 2 * ρ * v^2 + p) * v)
end

function euler1d!(du, u, p, t)
    (; prob, wstore, fluxstore) = p
    gd = prob.gd

    wstore = get_tmp(wstore, u)
    fluxstore = get_tmp(fluxstore, u)

    reconstruct!(prob, wstore, u)

    solve_riemann_problem!(prob, fluxstore, wstore)

    for j in 1:3, i in 1:gd.Nx
        im = i == 1 ? gd.Nx : i - 1

        du[j, i] = -(fluxstore[j, i] - fluxstore[j, im]) / gd.Δx
    end
end

function solve(prob::FDProblem{<:Any,Euler,<:Any,<:Any}, ρ0l, v0l, p0l, tspan)
    Nx = prob.gd.Nx
    γ = prob.model.γ

    u0 = zeros(3, Nx)

    @. u0[1, :] = ρ0l
    @. u0[2, :] = v0l * ρ0l
    @. u0[3, :] = p0l / (γ - 1) + 1 / 2 * ρ0l * v0l^2

    wstore = DiffCache(zeros(3, Nx, 2))
    fluxstore = DiffCache(zeros(3, Nx))

    prob = ODEProblem(euler1d!, u0, tspan, (; prob, wstore, fluxstore))

    sol = DifferentialEquations.solve(
        prob,
        Tsit5();
        # TRBDF2();
        # AutoTsit5(Rosenbrock23());
        # QNDF();
        saveat = range(tspan[1], tspan[2]; length = 100),
        abstol = 1e-8,
        reltol = 1e-8,
        progress = true,
        progress_steps = 100,
    )

    return sol
end