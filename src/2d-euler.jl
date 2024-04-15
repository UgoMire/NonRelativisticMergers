function get_primitive_variables(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, u, ix, iy)
    (; γ) = prob.model

    ρ = u[1, ix, iy]
    vx = u[2, ix, iy] / ρ
    vy = u[3, ix, iy] / ρ
    P = (γ - 1) * (u[4, ix, iy] - 1 / 2 * ρ * (vx^2 + vy^2))

    return ρ, vx, vy, P
end

function Fflux(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ, vx, vy, P)
    (; γ) = prob.model

    E = P / (γ - 1) + 1 / 2 * ρ * (vx^2 + vy^2)

    return (ρ * vx, ρ * vx^2 + P, ρ * vx * vy, (E + P) * vx)
end

function Gflux(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ, vx, vy, P)
    (; γ) = prob.model

    E = P / (γ - 1) + 1 / 2 * ρ * (vx^2 + vy^2)

    return (ρ * vy, ρ * vx * vy, ρ * vy^2 + P, (E + P) * vy)
end

function euler2d!(du, u, p, t)
    (; prob, wstore, xfluxstore, yfluxstore) = p
    (; xl, yl, Δx, Δy, Nx, Ny) = prob.grid

    reconstruct!(prob, wstore, u)

    solve_riemann_problem!(prob, xfluxstore, yfluxstore, wstore)

    for j in 1:4, ix in eachindex(xl), iy in eachindex(yl)
        ixm = ix == 1 ? Nx : ix - 1
        iym = iy == 1 ? Ny : iy - 1

        du[j, ix, iy] =
            -(xfluxstore[j, ix, iy] - xfluxstore[j, ixm, iy]) / Δx -
            (yfluxstore[j, ix, iy] - yfluxstore[j, ix, iym]) / Δy
    end
end

function setup_initial_state(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ0, vx0, vy0, P0)
    (; Nx, Ny) = prob.grid
    (; γ) = prob.model

    u0 = zeros(4, Nx, Ny)

    @. u0[1, :, :] = ρ0
    @. u0[2, :, :] = vx0 * ρ0
    @. u0[3, :, :] = vy0 * ρ0
    @. u0[4, :, :] = P0 / (γ - 1) + 1 / 2 * ρ0 * (vx0^2 + vy0^2)

    return u0
end

function solve(prob::FDProblem{Grid2D,Euler,<:Any,<:Any}, ρ0, vx0, vy0, P0, tspan)
    (; Nx, Ny) = prob.grid

    u0 = setup_initial_state(prob, ρ0, vx0, vy0, P0)

    wstore = zeros(4, Nx, Ny, 4)
    xfluxstore = zeros(4, Nx, Ny)
    yfluxstore = zeros(4, Nx, Ny)

    # wstore = DiffCache(zeros(4, Nx, Ny, 2))
    # xfluxstore = DiffCache(zeros(4, Nx, Ny))
    # yfluxstore = DiffCache(zeros(4, Nx, Ny))

    prob = ODEProblem(euler2d!, u0, tspan, (; prob, wstore, xfluxstore, yfluxstore))

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
        progress_steps = 100,
    )

    return sol
end