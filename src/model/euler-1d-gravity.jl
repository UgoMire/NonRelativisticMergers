function setup_finite_volume_cache(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any}
)
    (; grid) = prob
    (; Nx) = grid

    wstore = DiffCache(zeros(3, Nx, 2))
    fluxstore = DiffCache(zeros(3, Nx))
    sourcestore = DiffCache(zeros(3, Nx))
    potentialstore = DiffCache(zeros(Nx))

    return (;
        wstore,
        fluxstore,
        sourcestore,
        potentialstore
    )
end

function get_source!(
        prob::FDProblem{
            Grid1D, <:Union{EulerStaticGravity, EulerSelfGravity}, <:Any, <:Any},
        sourcestore,
        u,
        ϕ
)
    (; Nx, Δx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        ρ = u[1, i]
        ρv = u[2, i]

        dϕdx = (ϕ[ip] - ϕ[im]) / 2Δx

        sourcestore[2, i] = -ρ * dϕdx
        sourcestore[3, i] = -ρv * dϕdx
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
        prob::FDProblem{Grid1D, EulerStaticGravity, <:Any, <:Any},
        ρ0l,
        v0l,
        P0l,
        ϕext,
        tspan
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
        (; prob, wstore, fluxstore, sourcestore, ϕext)
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
        progress_steps = 100
    )

    return sol
end

function euler1d_self_gravity!(du, u, p, t)
    (; prob, fv_cache, fft_cache) = p
    (; wstore, fluxstore, sourcestore, potentialstore) = fv_cache
    (; Nx, Δx) = prob.grid

    wstore = get_tmp(wstore, u)
    fluxstore = get_tmp(fluxstore, u)
    sourcestore = get_tmp(sourcestore, u)
    potentialstore = get_tmp(potentialstore, u)

    solve_poisson!(prob, potentialstore, u[1, :], fft_cache)
    reconstruct!(prob, wstore, u, potentialstore)
    solve_riemann_problem!(prob, fluxstore, wstore)
    get_source!(prob, sourcestore, u, potentialstore)

    for j in 1:3, i in 1:Nx
        im = i == 1 ? Nx : i - 1

        du[j, i] = -(fluxstore[j, i] - fluxstore[j, im]) / Δx + sourcestore[j, i]
    end
end

function solve(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any}, ρ0l, v0l, P0l, tspan)
    (; grid) = prob

    u0 = setup_initial_state(prob, ρ0l, v0l, P0l)

    fv_cache = setup_finite_volume_cache(prob)

    fft_cache = setup_fft_cache(grid)

    prob = ODEProblem(
        euler1d_self_gravity!,
        u0,
        tspan,
        (; prob, fv_cache, fft_cache)
    )

    sol = OrdinaryDiffEq.solve(
        prob,
        Tsit5();
        # TRBDF2(; autodiff = AutoFiniteDiff());
        # AutoTsit5(Rosenbrock23());
        # QNDF();
        saveat = range(tspan[1], tspan[2]; length = 100),
        abstol = 1e-8,
        reltol = 1e-8,
        progress = true,
        progress_steps = 100
    )

    return sol
end
