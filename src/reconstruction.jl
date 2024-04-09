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

    for i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1

        w_reconstruct[:, i, 1] = u[:, i]
        w_reconstruct[:, i, 2] = u[:, ip]
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

function reconstruct(method::KT, gd, u)
    θ = method.θ
    w_reconstruct = zeros(3, gd.Nx, 2)

    for j = 1:3, i = 1:gd.Nx
        ip = i == gd.Nx ? 1 : i + 1
        im = i == 1 ? gd.Nx : i - 1

        w_reconstruct[j, i, 1] =
            u[j, i] -
            minmod(
                θ * (u[j, i] - u[j, im]) / gd.Δx,
                (u[j, ip] - u[j, im]) / 2gd.Δx,
                θ * (u[j, ip] - u[j, i]) / gd.Δx,
            ) * gd.Δx / 2

        w_reconstruct[j, i, 2] =
            u[j, i] +
            minmod(
                θ * (u[j, i] - u[j, im]) / gd.Δx,
                (u[j, ip] - u[j, im]) / 2gd.Δx,
                θ * (u[j, ip] - u[j, i]) / gd.Δx,
            ) * gd.Δx / 2
    end

    return w_reconstruct
end