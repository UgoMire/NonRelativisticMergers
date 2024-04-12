abstract type Reconstruction end

struct Constant <: Reconstruction end
struct MUSCL <: Reconstruction
    κ::Float64

    MUSCL(κ = 1 / 3) = new(κ)
end
struct KT <: Reconstruction
    θ::Float64

    KT(θ = 2) = new(θ)
end

function reconstruct(method::Constant, gd, u)
    w_reconstruct = zeros(3, gd.Nx, 2)

    γ = 5 / 3

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        # Reconstruct the density.
        rho = u[1, i]
        rhop = u[1, ip]

        w_reconstruct[1, i, 1] = rho
        w_reconstruct[1, i, 2] = rhop

        # Reconstruct the velocity.
        v = u[2, i] / u[1, i]
        vp = u[2, ip] / u[1, ip]

        w_reconstruct[2, i, 1] = v
        w_reconstruct[2, i, 2] = vp

        # Reconstruct the pressure.
        p = (γ - 1) * (u[3, i] - 1 / 2 * u[1, i] * u[2, i]^2)
        pp = (γ - 1) * (u[3, ip] - 1 / 2 * u[1, ip] * u[2, ip]^2)

        w_reconstruct[3, i, 1] = p
        w_reconstruct[3, i, 2] = pp
    end

    return w_reconstruct
end

function reconstruct(method::MUSCL, gd, u)
    κ = method.κ

    w_reconstruct = zeros(3, gd.Nx, 2)

    ρl = @view u[1, :]
    vl = @view u[2, :]
    pl = @view u[3, :]

    for i = 1:gd.Nx
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

function reconstruct!(w_store, method::KT, gd, u)
    θ = method.θ
    γ = 5 / 3

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1
        im = i == 1 ? gd.Nx : i - 1

        # Reconstruct the density.
        rho = u[1, i]
        rhop = u[1, ip]
        rhom = u[1, im]

        Drho = minmod(
            θ * (rho - rhom) / gd.Δx,
            (rhop - rhom) / 2gd.Δx,
            θ * (rhop - rho) / gd.Δx,
        )

        w_store[1, i, 1] = rho - Drho * gd.Δx / 2
        w_store[1, i, 2] = rho + Drho * gd.Δx / 2

        # Reconstruct the velocity.
        v = u[2, i] / u[1, i]
        vp = u[2, ip] / u[1, ip]
        vm = u[2, im] / u[1, im]

        Dv = minmod(θ * (v - vm) / gd.Δx, (vp - vm) / 2gd.Δx, θ * (vp - v) / gd.Δx)

        w_store[2, i, 1] = v - Dv * gd.Δx / 2
        w_store[2, i, 2] = v + Dv * gd.Δx / 2

        # Reconstruct the pressure.
        p = (γ - 1) * (u[3, i] - 1 / 2 * u[1, i] * u[2, i]^2)
        pp = (γ - 1) * (u[3, ip] - 1 / 2 * u[1, ip] * u[2, ip]^2)
        pm = (γ - 1) * (u[3, im] - 1 / 2 * u[1, im] * u[2, im]^2)

        Dp = minmod(θ * (p - pm) / gd.Δx, (pp - pm) / 2gd.Δx, θ * (pp - p) / gd.Δx)

        w_store[3, i, 1] = p - Dp * gd.Δx / 2
        w_store[3, i, 2] = p + Dp * gd.Δx / 2
    end
end