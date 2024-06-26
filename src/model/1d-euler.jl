@inline @fastmath function get_primitive_variables(
    prob::FDProblem{Grid1D,<:FDModel,<:Any,<:Any},
    u,
    i,
)
    (; γ) = prob.model

    ρ = u[1, i]
    v = u[2, i] / u[1, i]
    p = (γ - 1) * (u[3, i] - 1 / 2 * u[1, i] * u[2, i]^2)

    return ρ, v, p
end

function flux(prob::FDProblem{Grid1D,<:FDModel,<:Any,<:Any}, ρ, v, p)
    (; γ) = prob.model

    return (ρ * v, ρ * v^2 + p, (p / (γ - 1) + 1 / 2 * ρ * v^2 + p) * v)
end

function euler1d!(du, u, p, t)
    (; prob, wstore, fluxstore) = p
    (; Nx, Δx) = prob.grid

    wstore = get_tmp(wstore, u)
    fluxstore = get_tmp(fluxstore, u)

    reconstruct!(prob, wstore, u)
    solve_riemann_problem!(prob, fluxstore, wstore)

    for j in 1:3, i in 1:Nx
        im = i == 1 ? Nx : i - 1

        du[j, i] = -(fluxstore[j, i] - fluxstore[j, im]) / Δx
    end
end

function setup_initial_state(prob::FDProblem{Grid1D,<:FDModel,<:Any,<:Any}, ρ0l, v0l, p0l)
    (; Nx) = prob.grid
    (; γ) = prob.model

    u0 = zeros(3, Nx)

    @. u0[1, :] = ρ0l
    @. u0[2, :] = v0l * ρ0l
    @. u0[3, :] = p0l / (γ - 1) + 1 / 2 * ρ0l * v0l^2

    return u0
end

function solve(prob::FDProblem{Grid1D,Euler,<:Any,<:Any}, ρ0l, v0l, p0l, tspan)
    (; Nx) = prob.grid

    u0 = setup_initial_state(prob, ρ0l, v0l, p0l)

    wstore = DiffCache(zeros(3, Nx, 2))
    fluxstore = DiffCache(zeros(3, Nx))

    prob = ODEProblem(euler1d!, u0, tspan, (; prob, wstore, fluxstore))

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

function get_source!(
    prob::FDProblem{Grid1D,<:Union{EulerStaticGravity,EulerSelfGravity},<:Any,<:Any},
    sourcestore,
    u,
    ϕ,
)
    (; Nx, Δx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        dϕdx = (ϕ[ip] - ϕ[im]) / 2Δx

        sourcestore[2, i] = -u[1, i] * dϕdx
        sourcestore[3, i] = -u[2, i] * dϕdx
    end
end

function euler1d_static_gravity!(du, u, p, t)
    (; prob, wstore, fluxstore, sourcestore, ϕext) = p
    (; Nx, Δx) = prob.grid

    wstore = get_tmp(wstore, u)
    fluxstore = get_tmp(fluxstore, u)
    sourcestore = get_tmp(sourcestore, u)

    reconstruct!(prob, wstore, u, ϕext)
    solve_riemann_problem!(prob, fluxstore, wstore)
    get_source!(prob, sourcestore, u, ϕext)

    for i in 1:Nx
        im = i == 1 ? Nx : i - 1

        @. du[:, i] = -(fluxstore[:, i] - fluxstore[:, im]) / Δx + sourcestore[:, i]
    end
end

function solve(
    prob::FDProblem{Grid1D,EulerStaticGravity,<:Any,<:Any},
    ρ0l,
    v0l,
    P0l,
    ϕext,
    tspan,
)
    (; Nx) = prob.grid

    u0 = setup_initial_state(prob, ρ0l, v0l, P0l)

    wstore = zeros(3, Nx, 2)
    fluxstore = zeros(3, Nx)
    sourcestore = zeros(3, Nx)

    wstore = DiffCache(zeros(3, Nx, 2))
    fluxstore = DiffCache(zeros(3, Nx))
    sourcestore = DiffCache(zeros(3, Nx))

    prob = ODEProblem(
        euler1d_static_gravity!,
        u0,
        tspan,
        (; prob, wstore, fluxstore, sourcestore, ϕext),
    )

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

function euler1d_self_gravity!(du, u, p, t)
    (; prob, wstore, fluxstore, sourcestore, potentialstore) = p
    (; Nx, Δx) = prob.grid

    wstore = get_tmp(wstore, u)
    fluxstore = get_tmp(fluxstore, u)
    sourcestore = get_tmp(sourcestore, u)
    potentialstore = get_tmp(potentialstore, u)

    solve_poisson!(prob, potentialstore, u[1, :])
    reconstruct!(prob, wstore, u, potentialstore)
    solve_riemann_problem!(prob, fluxstore, wstore)
    get_source!(prob, sourcestore, u, potentialstore)

    for i in 1:Nx
        im = i == 1 ? Nx : i - 1

        @. du[:, i] = -(fluxstore[:, i] - fluxstore[:, im]) / Δx + sourcestore[:, i]
    end
end

function solve(prob::FDProblem{Grid1D,EulerSelfGravity,<:Any,<:Any}, ρ0l, v0l, P0l, tspan)
    (; Nx) = prob.grid

    u0 = setup_initial_state(prob, ρ0l, v0l, P0l)

    wstore = DiffCache(zeros(3, Nx, 2))
    fluxstore = DiffCache(zeros(3, Nx))
    sourcestore = DiffCache(zeros(3, Nx))
    potentialstore = DiffCache(zeros(1, Nx))

    prob = ODEProblem(
        euler1d_self_gravity!,
        u0,
        tspan,
        (; prob, wstore, fluxstore, sourcestore, potentialstore),
    )

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