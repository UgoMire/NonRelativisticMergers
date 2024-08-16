"""
    solve_poisson(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        ρ
)

FFT Poisson equation solver in 2D. This implementation is not optimized but
maybe clearer to understand.
"""
function solve_poisson(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        ρ;
        fft_type = :rfft
)
    (; grid, model) = prob
    (; Nx, Ny, Δx, Δy) = grid
    (; G, ϵ) = model

    if fft_type == :fft
        kx = 2π * fftfreq(Nx, 1 / Δx)
        ky = 2π * fftfreq(Ny, 1 / Δy)

        ρhat = fft(ρ)

        uhat = zeros(Complex{Float64}, Nx, Ny)
    elseif fft_type == :rfft
        kx = 2π * rfftfreq(Nx, 1 / Δx)
        ky = 2π * fftfreq(Ny, 1 / Δy)

        ρhat = rfft(ρ)

        uhat = zeros(Complex{Float64}, Nx ÷ 2 + 1, Ny)
    end

    for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
        if kxi == 0 && kyj == 0
            uhat[i, j] = 0
        else
            uhat[i, j] = -4π * G * ρhat[i, j] / (kxi^2 + kyj^2 + ϵ^2)
        end
    end

    if fft_type == :fft
        return real(ifft(uhat))
    elseif fft_type == :rfft
        return irfft(uhat, Nx)
    end
end

function setup_fft_cache(grid::Grid2D)
    (; Nx, Ny) = grid

    return setup_fft_cache(grid, zeros(Nx, Ny))
end

function setup_fft_cache(grid::Grid2D, ρ)
    (; Nx, Δx, Ny, Δy) = grid

    return (;
        kx = 2π * rfftfreq(Nx, 1 / Δx),
        ky = 2π * fftfreq(Ny, 1 / Δy),
        ρ̂_cache = zeros(ComplexF64, Nx ÷ 2 + 1, Ny),
        û_cache = zeros(ComplexF64, Nx ÷ 2 + 1, Ny),
        planned_fft = plan_rfft(ρ),
        planned_ifft = plan_irfft(zeros(ComplexF64, Nx ÷ 2 + 1, Ny), Nx)
    )
end

"""
    solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        cache
)

FFT Poisson equation solver in 1D. This implementation is non-allocating thanks
to the preallocated arrays in `cache`. The cache `cache` can be generated with
`setup_fft_cache`.
"""
function solve_poisson!(
        prob::FDProblem{Grid2D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        (; kx, ky, planned_fft, planned_ifft, ρ̂_cache, û_cache)
)
    (; model) = prob
    (; G, ϵ) = model

    mul!(ρ̂_cache, planned_fft, ρ)

    for (i, kxi) in enumerate(kx), (j, kyj) in enumerate(ky)
        if kxi == 0 && kyj == 0
            û_cache[i, j] = 0
        else
            û_cache[i, j] = -4π * G * ρ̂_cache[i, j] / (kxi^2 + kyj^2 + ϵ^2)
        end
    end

    mul!(potentialstore, planned_ifft, û_cache)
end
