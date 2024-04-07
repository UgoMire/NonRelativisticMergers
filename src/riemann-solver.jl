abstract type RiemannSolver end

struct NaiveRS <: RiemannSolver end
struct HLLC <: RiemannSolver end
struct Roe <: RiemannSolver end

function solve_riemann(method::NaiveRS, gd, w_reconstruct)
    Fl = zeros(3, gd.Nx)

    for i = 1:gd.Nx
        Fminus =
            flux(w_reconstruct[1, i, 1], w_reconstruct[2, i, 1], w_reconstruct[3, i, 1])
        Fplus = flux(w_reconstruct[1, i, 2], w_reconstruct[2, i, 2], w_reconstruct[3, i, 2])

        @. Fl[:, i] = 1 / 2 * (Fplus + Fminus)
    end

    return Fl
end