function reconstruct!(prob::FDProblem{Grid1D, <:Any, Constant, <:Any}, wstore, u)
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1

        ρ, v, P = get_primitive_variables(prob, u, i)
        ρp, vp, Pp = get_primitive_variables(prob, u, ip)

        # Reconstruct the density.
        wstore[1, i, 1] = ρ
        wstore[1, i, 2] = ρp

        # Reconstruct the velocity.
        wstore[2, i, 1] = v
        wstore[2, i, 2] = vp

        # Reconstruct the pressure.
        wstore[3, i, 1] = P
        wstore[3, i, 2] = Pp
    end
end

function reconstruct!(prob::FDProblem{Grid2D, <:Any, Constant, <:Any}, wstore, u)
    (; Nx, Ny) = prob.grid

    for ix in 1:Nx, iy in 1:Ny
        ixp = ix == Nx ? 1 : ix + 1
        iyp = iy == Ny ? 1 : iy + 1

        # Get the primitive variables.
        ρ, vx, vy, P = get_primitive_variables(prob, u, ix, iy)
        ρ_px, vx_px, vy_px, P_px = get_primitive_variables(prob, u, ixp, iy)
        ρ_py, vx_py, vy_py, P_py = get_primitive_variables(prob, u, ix, iyp)

        # Reconstruct the density.
        wstore[1, ix, iy, 1] = ρ
        wstore[1, ix, iy, 2] = ρ_px
        wstore[1, ix, iy, 3] = ρ
        wstore[1, ix, iy, 4] = ρ_py

        # Reconstruct the x velocity.
        wstore[2, ix, iy, 1] = vx
        wstore[2, ix, iy, 2] = vx_px
        wstore[2, ix, iy, 3] = vx
        wstore[2, ix, iy, 4] = vx_py

        # Reconstruct the y velocity.
        wstore[3, ix, iy, 1] = vy
        wstore[3, ix, iy, 2] = vy_px
        wstore[3, ix, iy, 3] = vy
        wstore[3, ix, iy, 4] = vy_py

        # Reconstruct the pressure.
        wstore[4, ix, iy, 1] = P
        wstore[4, ix, iy, 2] = P_px
        wstore[4, ix, iy, 3] = P
        wstore[4, ix, iy, 4] = P_py
    end
end
