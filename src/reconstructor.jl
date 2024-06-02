function get_primitive_variables_at_boundary(prob, wr, i, j, n)
    if n.x == 1 && n.y == 0
        ip = i == prob.grid.Nx ? 1 : i + 1

        wL = (
            ρ = wr[1, i, j, 2],
            vx = wr[2, i, j, 2],
            vy = wr[3, i, j, 2],
            P = wr[4, i, j, 2],
        )
        wR = (
            ρ = wr[1, ip, j, 2],
            vx = wr[2, ip, j, 2],
            vy = wr[3, ip, j, 2],
            P = wr[4, ip, j, 2],
        )
    elseif n.x == 0 && n.y == 1
        jp = j == prob.grid.Ny ? 1 : j + 1

        wL = (
            ρ = wr[1, i, j, 4],
            vx = wr[2, i, j, 4],
            vy = wr[3, i, j, 4],
            P = wr[4, i, j, 4],
        )
        wR = (
            ρ = wr[1, i, jp, 3],
            vx = wr[2, i, jp, 3],
            vy = wr[3, i, jp, 3],
            P = wr[4, i, jp, 3],
        )
    end

    return (; wL, wR)
end

function reconstruct!(prob::FDProblem{Grid1D,<:Any,Constant,<:Any}, wstore, u)
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

function reconstruct!(prob::FDProblem{Grid2D,<:Any,Constant,<:Any}, wstore, u)
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

# function reconstruct(method::MUSCL, gd, u)
#     κ = method.κ

#     w_reconstruct = zeros(3, gd.Nx, 2)

#     ρl = @view u[1, :]
#     vl = @view u[2, :]
#     pl = @view u[3, :]

#     for i in 1:gd.Nx
#         ip = i == gd.Nx ? 1 : i + 1
#         im = i == 1 ? gd.Nx : i - 1
#         ipp = i == gd.Nx ? 2 : i == gd.Nx - 1 ? 1 : i + 2

#         κ = 1 / 3

#         w_reconstruct[1, i, 1] =
#             ρl[i] + 1 / 4 * ((1 - κ) * (ρl[i] - ρl[im]) + (1 + κ) * (ρl[ip] - ρl[i]))
#         w_reconstruct[1, i, 2] =
#             ρl[ip] - 1 / 4 * ((1 + κ) * (ρl[ip] - ρl[i] + (1 - κ) * (ρl[ipp] - ρl[ip])))
#         w_reconstruct[2, i, 1] =
#             vl[i] + 1 / 4 * ((1 - κ) * (vl[i] - vl[im]) + (1 + κ) * (vl[ip] - vl[i]))
#         w_reconstruct[2, i, 2] =
#             vl[ip] - 1 / 4 * ((1 + κ) * (vl[ip] - vl[i] + (1 - κ) * (vl[ipp] - vl[ip])))
#         w_reconstruct[3, i, 1] =
#             pl[i] + 1 / 4 * ((1 - κ) * (pl[i] - pl[im]) + (1 + κ) * (pl[ip] - pl[i]))
#         w_reconstruct[3, i, 2] =
#             pl[ip] - 1 / 4 * ((1 + κ) * (pl[ip] - pl[i] + (1 - κ) * (pl[ipp] - pl[ip])))
#     end

#     return w_reconstruct
# end

@inline @fastmath function minmod(a1, a2, a3)
    if a1 > 0 && a2 > 0 && a3 > 0
        return min(a1, a2, a3)
    elseif a1 < 0 && a1 < 0 && a3 < 0
        return max(a1, a2, a3)
    else
        return 0
    end
end

@inline @fastmath function kt_reconstruction(prob, q, qp, qm)
    (; θ) = prob.reconstructor
    (; Δx) = prob.grid

    Dq = minmod(θ * (q - qm) / Δx, (qp - qm) / 2Δx, θ * (qp - q) / Δx)

    qL = q - Dq * Δx / 2
    qR = q + Dq * Δx / 2

    return qL, qR
end

@inline @fastmath function reconstruct!(
    prob::FDProblem{GD,MD,RC,RS},
    wstore,
    u,
) where {GD<:Grid1D,MD<:Euler,RC<:KT,RS<:Any}
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        ρ, v, P = get_primitive_variables(prob, u, i)
        ρm, vm, Pm = get_primitive_variables(prob, u, im)
        ρp, vp, Pp = get_primitive_variables(prob, u, ip)

        # Reconstruct the density, x-velocity and pressure.
        wstore[1, i, :] .= kt_reconstruction(prob, ρ, ρp, ρm)
        wstore[2, i, :] .= kt_reconstruction(prob, v, vp, vm)
        wstore[3, i, :] .= kt_reconstruction(prob, P, Pp, Pm)
    end
end

function reconstruct!(
    prob::FDProblem{Grid1D,<:Union{EulerStaticGravity,EulerSelfGravity},KT,<:Any},
    wstore,
    u,
    ϕext,
)
    (; Nx) = prob.grid

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        ρ, v, P = get_primitive_variables(prob, u, i)
        ρm, vm, Pm = get_primitive_variables(prob, u, im)
        ρp, vp, Pp = get_primitive_variables(prob, u, ip)

        # Reconstruct the density and x-velocity.
        wstore[1, i, :] .= kt_reconstruction(prob, ρ, ρp, ρm)
        wstore[2, i, :] .= kt_reconstruction(prob, v, vp, vm)
        # wstore[3, i, :] .= kt_reconstruction(prob, P, Pp, Pm)

        # Reconstruct the pressure with hydrostatic equilibrium preserving scheme.
        ϕ, ϕp, ϕm = ϕext[i], ϕext[ip], ϕext[im]

        P0_left = P + ρ * (ϕ - ϕm) / 2
        P0_right = P - ρ * (ϕp - ϕ) / 2

        P1 = 0
        P1m = Pm - P - (ρm + ρ) / 2 * (ϕ - ϕm)
        P1p = Pp - P + (ρ + ρp) / 2 * (ϕp - ϕ)

        P1_left, P1_right = kt_reconstruction(prob, P1, P1p, P1m)

        wstore[3, i, 1] = P0_left + P1_left
        wstore[3, i, 2] = P0_right + P1_right
    end
end

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