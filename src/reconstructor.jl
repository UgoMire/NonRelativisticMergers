function reconstruct!(prob::FDProblem{<:Any,<:Any,Constant,<:Any}, wstore, u)
    Nx = prob.gd.Nx
    γ = prob.model.γ

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1

        # Reconstruct the density.
        rho = u[1, i]
        rhop = u[1, ip]

        wstore[1, i, 1] = rho
        wstore[1, i, 2] = rhop

        # Reconstruct the velocity.
        v = u[2, i] / u[1, i]
        vp = u[2, ip] / u[1, ip]

        wstore[2, i, 1] = v
        wstore[2, i, 2] = vp

        # Reconstruct the pressure.
        p = (γ - 1) * (u[3, i] - 1 / 2 * u[1, i] * u[2, i]^2)
        pp = (γ - 1) * (u[3, ip] - 1 / 2 * u[1, ip] * u[2, ip]^2)

        wstore[3, i, 1] = p
        wstore[3, i, 2] = pp
    end
end

function reconstruct(method::MUSCL, gd, u)
    κ = method.κ

    w_reconstruct = zeros(3, gd.Nx, 2)

    ρl = @view u[1, :]
    vl = @view u[2, :]
    pl = @view u[3, :]

    for i in 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1
        im = i == 1 ? gd.Nx : i - 1
        ipp = i == gd.Nx ? 2 : i == gd.Nx - 1 ? 1 : i + 2

        κ = 1 / 3

        w_reconstruct[1, i, 1] =
            ρl[i] + 1 / 4 * ((1 - κ) * (ρl[i] - ρl[im]) + (1 + κ) * (ρl[ip] - ρl[i]))
        w_reconstruct[1, i, 2] =
            ρl[ip] - 1 / 4 * ((1 + κ) * (ρl[ip] - ρl[i] + (1 - κ) * (ρl[ipp] - ρl[ip])))
        w_reconstruct[2, i, 1] =
            vl[i] + 1 / 4 * ((1 - κ) * (vl[i] - vl[im]) + (1 + κ) * (vl[ip] - vl[i]))
        w_reconstruct[2, i, 2] =
            vl[ip] - 1 / 4 * ((1 + κ) * (vl[ip] - vl[i] + (1 - κ) * (vl[ipp] - vl[ip])))
        w_reconstruct[3, i, 1] =
            pl[i] + 1 / 4 * ((1 - κ) * (pl[i] - pl[im]) + (1 + κ) * (pl[ip] - pl[i]))
        w_reconstruct[3, i, 2] =
            pl[ip] - 1 / 4 * ((1 + κ) * (pl[ip] - pl[i] + (1 - κ) * (pl[ipp] - pl[ip])))
    end

    return w_reconstruct
end

function minmod(a1, a2, a3)
    if a1 > 0 && a2 > 0 && a3 > 0
        return min(a1, a2, a3)
    elseif a1 < 0 && a1 < 0 && a3 < 0
        return max(a1, a2, a3)
    else
        return 0
    end
end

function reconstruct!(prob::FDProblem{<:Any,<:Any,KT,<:Any}, wstore, u)
    Nx = prob.gd.Nx
    Δx = prob.gd.Δx
    θ = prob.reconstructor.θ
    γ = prob.model.γ

    for i in 1:Nx
        ip = i == Nx ? 1 : i + 1
        im = i == 1 ? Nx : i - 1

        # Reconstruct the density.
        rho = u[1, i]
        rhop = u[1, ip]
        rhom = u[1, im]

        Drho = minmod(θ * (rho - rhom) / Δx, (rhop - rhom) / 2Δx, θ * (rhop - rho) / Δx)

        wstore[1, i, 1] = rho - Drho * Δx / 2
        wstore[1, i, 2] = rho + Drho * Δx / 2

        # Reconstruct the velocity.
        v = u[2, i] / u[1, i]
        vp = u[2, ip] / u[1, ip]
        vm = u[2, im] / u[1, im]

        Dv = minmod(θ * (v - vm) / Δx, (vp - vm) / 2Δx, θ * (vp - v) / Δx)

        wstore[2, i, 1] = v - Dv * Δx / 2
        wstore[2, i, 2] = v + Dv * Δx / 2

        # Reconstruct the pressure.
        p = (γ - 1) * (u[3, i] - 1 / 2 * u[1, i] * u[2, i]^2)
        pp = (γ - 1) * (u[3, ip] - 1 / 2 * u[1, ip] * u[2, ip]^2)
        pm = (γ - 1) * (u[3, im] - 1 / 2 * u[1, im] * u[2, im]^2)

        Dp = minmod(θ * (p - pm) / Δx, (pp - pm) / 2Δx, θ * (pp - p) / Δx)

        wstore[3, i, 1] = p - Dp * Δx / 2
        wstore[3, i, 2] = p + Dp * Δx / 2
    end
end