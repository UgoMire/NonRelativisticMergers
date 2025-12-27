import NonRelativisticMergers as NRM

using NonRelativisticMergers
using LinearAlgebra
using FFTW
using GLMakie

# gd = Grid2D(; xmin = -3, xmax = 3, ymin = -3, ymax = 3, Nx = 2^5, Ny = 2^6)
# model = EulerSelfGravity(; γ = 6 / 3, G = 2, ϵ = 2)

# prob = FDProblem(gd, model, KT(), HLLC())

# ρ = [0.1 + 10 * exp(-10 * (x^2 + y^2)) for x in gd.xl, y in gd.yl]
ρ = ρ0

ϕ = NRM.solve_poisson(prob, ρ)

heatmap(ρ)
heatmap(ϕ)

## Define the gradient operator and compute the gradient of the potential.
∇x = Matrix(Tridiagonal(-ones(gd.Nx - 1), zeros(gd.Nx), ones(gd.Nx - 1)))
∇x[1, gd.Nx] = -1
∇x[gd.Nx, 1] = 1
∇x ./= 2gd.Δx

∇y = Matrix(Tridiagonal(-ones(gd.Ny - 1), zeros(gd.Ny), ones(gd.Ny - 1)))
∇y[1, gd.Ny] = -1
∇y[gd.Ny, 1] = 1
∇y ./= 2gd.Δy

∇xϕ = ∇x * ϕ
# ∇yϕ = transpose(∇y * transpose(ϕ))
∇yϕ = ϕ * transpose(∇y)

∇xϕ2 = zeros(gd.Nx, gd.Ny)
for j in 1:(gd.Ny)
    ∇xϕ2[:, j] = ∇x * ϕ[:, j]
end

∇yϕ2 = zeros(gd.Nx, gd.Ny)
for i in 1:(gd.Nx)
    ∇yϕ2[i, :] = ∇y * ϕ[i, :]
end

∇xϕ ≈ ∇xϕ2
∇yϕ ≈ ∇yϕ2

heatmap(∇xϕ)
heatmap(∇yϕ)

## Solve for the pressure.
lhs = -∇x * (ρ .* ∇xϕ) .- (ρ .* ∇yϕ) * transpose(∇y)

kx = 2π * rfftfreq(gd.Nx, 1 / gd.Δx)
ky = 2π * fftfreq(gd.Ny, 1 / gd.Δy)

lhshat = rfft(lhs)

phat = zeros(ComplexF64, gd.Nx ÷ 2 + 1, gd.Ny)

for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
    if kxi == 0 && kyj == 0
        phat[i, j] = 0
    else
        phat[i, j] = -lhshat[i, j] / (kxi^2 + kyj^2)
    end
end

p = irfft(phat, gd.Nx)

## Check the hydrostatic equation.
∇xp = ∇x * p
∇yp = p * transpose(∇y)

norm(∇xp .+ ρ .* ∇xϕ)
norm(∇yp .+ ρ .* ∇yϕ)

##
P0 = p .+ 1