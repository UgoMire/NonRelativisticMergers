function solve_riemann_problem!(
    prob::FDProblem{Grid1D,<:Any,<:Any,NaiveRS},
    fluxstore,
    wreconstructed,
)
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1

        Fminux = flux(
            prob,
            wreconstructed[1, i, 2],
            wreconstructed[2, i, 2],
            wreconstructed[3, i, 2],
        )
        Fplus = flux(
            prob,
            wreconstructed[1, ip, 1],
            wreconstructed[2, ip, 1],
            wreconstructed[3, ip, 1],
        )

        @. fluxstore[:, i] = 1 / 2 * (Fminux + Fplus)
    end
end

function solve_riemann_problem!(
    prob::FDProblem{Grid2D,<:Any,<:Any,NaiveRS},
    xfluxstore,
    yfluxstore,
    w,
)
    (; Nx, Ny) = prob.grid

    ρ = @view w[1, :, :, :]
    vx = @view w[2, :, :, :]
    vy = @view w[3, :, :, :]
    P = @view w[4, :, :, :]

    for i in 1:Nx, j in 1:Ny
        ip = i == Nx ? 1 : i + 1
        jp = j == Ny ? 1 : j + 1

        Fminus = Fflux(prob, ρ[i, j, 2], vx[i, j, 2], vy[i, j, 2], P[i, j, 2])
        Fplus = Fflux(prob, ρ[ip, j, 1], vx[ip, j, 1], vy[ip, j, 1], P[ip, j, 1])

        @. xfluxstore[:, i, j] = 1 / 2 * (Fminus + Fplus)
        # @. xfluxstore[:, i, j] = Fplus
        # @. xfluxstore[:, i, j] = Fminus

        Gminus = Gflux(prob, ρ[i, j, 4], vx[i, j, 4], vy[i, j, 4], P[i, j, 4])
        Gplus = Gflux(prob, ρ[i, jp, 3], vx[i, jp, 3], vy[i, jp, 3], P[i, jp, 3])

        @. yfluxstore[:, i, j] = 1 / 2 * (Gminus + Gplus)
        # @. yfluxstore[:, i, j] = Gplus
        # @. yfluxstore[:, i, j] = Gminus
    end
end