function get_source!(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sourcestore, u, ϕ)
    (; grid) = prob
    (; Nx, Δx, Ny, Δy) = grid

    for i in 1:Nx, j in 1:Ny
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        jp = j == Ny ? 1 : j + 1
        jm = j == 1 ? Ny : j - 1

        dϕdx = (ϕ[ip, j] - ϕ[im, j]) / 2Δx
        dϕdy = (ϕ[i, jp] - ϕ[i, jm]) / 2Δy
    end
end

function euler2d_self_gravity!(du, u, p, t)
    (; prob, wstore, xfluxstore, yfluxstore) = p
    (; xl, yl, Δx, Δy, Nx, Ny) = prob.grid

    ρ = u[1, :, :]

    solve_poisson!(prob, potentialstore, ρ, fft_cache)
    reconstruct!(prob, wstore, u)
    solve_riemann_problem!(prob, xfluxstore, yfluxstore, wstore)
    get_source!(prob, sourcestore, u, potentialstore)

    for j in 1:4, ix in eachindex(xl), iy in eachindex(yl)
        ixm = ix == 1 ? Nx : ix - 1
        iym = iy == 1 ? Ny : iy - 1

        du[j, ix, iy] = -(xfluxstore[j, ix, iy] - xfluxstore[j, ixm, iy]) / Δx -
                        (yfluxstore[j, ix, iy] - yfluxstore[j, ix, iym]) / Δy
    end
end

function solve(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any}, ρ0, vx0, vy0, P0, tspan)
    (; grid) = prob
    (; Nx, Ny) = grid

    u0 = setup_initial_state(prob, ρ0, vx0, vy0, P0)

    wstore = zeros(4, Nx, Ny, 4)
    xfluxstore = zeros(4, Nx, Ny)
    yfluxstore = zeros(4, Nx, Ny)
    potentialstore = zeros(Nx, Ny)

    fft_cache = setup_fft_cache(grid)

    prob = ODEProblem(euler2d!, u0, tspan,
        (; prob, wstore, xfluxstore, yfluxstore, potentialstore, fft_cache))

    sol = OrdinaryDiffEq.solve(
        prob,
        Tsit5();
        # TRBDF2();
        # AutoTsit5(Rosenbrock23());
        # QNDF();
        saveat = range(tspan[1], tspan[2]; length = 100),
        # saveat = [],
        abstol = 1e-8,
        reltol = 1e-8,
        progress = true,
        progress_steps = 100
    )

    return sol
end