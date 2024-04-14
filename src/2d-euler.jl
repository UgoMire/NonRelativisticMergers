function euler2d!(du, u, p, t)
    (; prob, wstore, fluxstore) = p

    return nothing
end

function setup_initial_state(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ0l, v0l, p0l)
    (; Nx, Ny) = prob.gd
    (; γ) = prob.model

    u0 = zeros(3, Nx, Ny)

    @. u0[1, :, :] = ρ0l
    @. u0[2, :, :] = v0l * ρ0l
    @. u0[3, :, :] = p0l / (γ - 1) + 1 / 2 * ρ0l * v0l^2

    return u0
end

function solve(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ0l, v0l, p0l, tspan)
    u0 = setup_initial_state(prob, ρ0l, v0l, p0l)

    wstore = DiffCache(zeros(3, Nx, Ny, 2))
    fluxstore = DiffCache(zeros(3, Nx, Ny))

    prob = ODEProblem(euler2d!, u0, tspan, (; prob, wstore, fluxstore))

    sol = OrdinaryDiffEq.solve(
        prob,
        # Tsit5();
        # TRBDF2();
        AutoTsit5(Rosenbrock23());
        # QNDF();
        saveat = range(tspan[1], tspan[2]; length = 100),
        abstol = 1e-8,
        reltol = 1e-8,
        progress = true,
        progress_steps = 100,
    )

    return sol
end