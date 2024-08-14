function reconstruct!(prob::FDProblem{Grid2D,<:Any,KT,<:Any}, wstore, u)
    (; Nx, Ny, Δx, Δy) = prob.grid
    (; θ) = prob.reconstructor

    for ix in 1:Nx, iy in 1:Ny
        ixp = ix == Nx ? 1 : ix + 1
        ixm = ix == 1 ? Nx : ix - 1
        iyp = iy == Ny ? 1 : iy + 1
        iym = iy == 1 ? Ny : iy - 1

        # Get the primitive variables.
        ρ, vx, vy, P = get_primitive_variables(prob, u, ix, iy)
        ρ_px, vx_px, vy_px, P_px = get_primitive_variables(prob, u, ixp, iy)
        ρ_mx, vx_mx, vy_mx, P_mx = get_primitive_variables(prob, u, ixm, iy)
        ρ_py, vx_py, vy_py, P_py = get_primitive_variables(prob, u, ix, iyp)
        ρ_my, vx_my, vy_my, P_my = get_primitive_variables(prob, u, ix, iym)

        # Reconstruct the density.
        Dx_ρ = minmod(θ * (ρ - ρ_mx) / Δx, (ρ_px - ρ_mx) / 2Δx, θ * (ρ_px - ρ) / Δx)
        Dy_ρ = minmod(θ * (ρ - ρ_my) / Δy, (ρ_py - ρ_my) / 2Δy, θ * (ρ_py - ρ) / Δy)

        wstore[1, ix, iy, 1] = ρ - Dx_ρ * Δx / 2
        wstore[1, ix, iy, 2] = ρ + Dx_ρ * Δx / 2
        wstore[1, ix, iy, 3] = ρ - Dy_ρ * Δy / 2
        wstore[1, ix, iy, 4] = ρ + Dy_ρ * Δy / 2

        # Reconstruct the x velocity.
        Dx_vx = minmod(θ * (vx - vx_mx) / Δx, (vx_px - vx_mx) / 2Δx, θ * (vx_px - vx) / Δx)
        Dy_vx = minmod(θ * (vx - vx_my) / Δy, (vx_py - vx_my) / 2Δy, θ * (vx_py - vx) / Δy)

        wstore[2, ix, iy, 1] = vx - Dx_vx * Δx / 2
        wstore[2, ix, iy, 2] = vx + Dx_vx * Δx / 2
        wstore[2, ix, iy, 3] = vx - Dy_vx * Δy / 2
        wstore[2, ix, iy, 4] = vx + Dy_vx * Δy / 2

        # Reconstruct the y velocity.
        Dx_vy = minmod(θ * (vy - vy_mx) / Δx, (vy_px - vy_mx) / 2Δx, θ * (vy_px - vy) / Δx)
        Dy_vy = minmod(θ * (vy - vy_my) / Δy, (vy_py - vy_my) / 2Δy, θ * (vy_py - vy) / Δy)

        wstore[3, ix, iy, 1] = vy - Dx_vy * Δx / 2
        wstore[3, ix, iy, 2] = vy + Dx_vy * Δx / 2
        wstore[3, ix, iy, 3] = vy - Dy_vy * Δy / 2
        wstore[3, ix, iy, 4] = vy + Dy_vy * Δy / 2

        # Reconstruct the pressure.
        Dx_P = minmod(θ * (P - P_mx) / Δx, (P_px - P_mx) / 2Δx, θ * (P_px - P) / Δx)
        Dy_P = minmod(θ * (P - P_my) / Δy, (P_py - P_my) / 2Δy, θ * (P_py - P) / Δy)

        wstore[4, ix, iy, 1] = P - Dx_P * Δx / 2
        wstore[4, ix, iy, 2] = P + Dx_P * Δx / 2
        wstore[4, ix, iy, 3] = P - Dy_P * Δy / 2
        wstore[4, ix, iy, 4] = P + Dy_P * Δy / 2
    end
end