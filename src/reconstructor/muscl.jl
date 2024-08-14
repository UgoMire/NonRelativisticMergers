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