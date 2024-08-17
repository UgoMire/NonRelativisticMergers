function setup_finite_volume_cache(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any})
    (; grid) = prob
    (; Nx, Ny) = grid

    wstore = zeros(4, Nx, Ny, 4)
    xfluxstore = zeros(4, Nx, Ny)
    yfluxstore = zeros(4, Nx, Ny)
    sourcestore = zeros(4, Nx, Ny)
    potentialstore = zeros(Nx, Ny)

    return (;
        wstore,
        xfluxstore,
        yfluxstore,
        potentialstore,
        sourcestore
    )
end

function get_source!(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        sourcestore, u, ϕ)
    (; grid) = prob
    (; Nx, Δx, Ny, Δy) = grid

    for ix in 1:Nx, iy in 1:Ny
        ip = ix == Nx ? 1 : ix + 1
        im = ix == 1 ? Nx : ix - 1

        jp = iy == Ny ? 1 : iy + 1
        jm = iy == 1 ? Ny : iy - 1

        ρ = u[1, ix, iy]
        ρvx = u[2, ix, iy]
        ρvy = u[3, ix, iy]

        dϕdx = (ϕ[ip, iy] - ϕ[im, iy]) / 2Δx
        dϕdy = (ϕ[ix, jp] - ϕ[ix, jm]) / 2Δy

        sourcestore[2, ix, iy] = -ρ * dϕdx
        sourcestore[3, ix, iy] = -ρ * dϕdy
        sourcestore[4, ix, iy] = -ρvx * dϕdx - ρvy * dϕdy
    end
end

function euler2d_self_gravity!(du, u, p, t)
    (; prob, fv_cache, fft_cache) = p
    (; xl, yl, Δx, Δy, Nx, Ny) = prob.grid
    (; wstore, xfluxstore, yfluxstore, potentialstore, sourcestore) = fv_cache

    # Solve Poisson equation.
    ρ = @view u[1, :, :]
    solve_poisson!(prob, potentialstore, ρ, fft_cache)

    # Perform the finite volume scheme.
    reconstruct!(prob, wstore, u)
    solve_riemann_problem!(prob, xfluxstore, yfluxstore, wstore)
    get_source!(prob, sourcestore, u, potentialstore)

    for j in 1:4, ix in eachindex(xl), iy in eachindex(yl)
        ixm = ix == 1 ? Nx : ix - 1
        iym = iy == 1 ? Ny : iy - 1

        du[j, ix, iy] = -(xfluxstore[j, ix, iy] - xfluxstore[j, ixm, iy]) / Δx -
                        (yfluxstore[j, ix, iy] - yfluxstore[j, ix, iym]) / Δy +
                        sourcestore[j, ix, iy]
    end
end

function solve(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        ρ0, vx0, vy0, P0, tspan)
    (; grid) = prob

    u0 = setup_initial_state(prob, ρ0, vx0, vy0, P0)
    fv_cache = setup_finite_volume_cache(prob)
    fft_cache = setup_fft_cache(grid, @view u0[1, :, :])

    prob = ODEProblem(euler2d_self_gravity!, u0, tspan,
        (; prob, fv_cache, fft_cache))

    sol = OrdinaryDiffEq.solve(
        prob,
        Tsit5();
        # Vern6();
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