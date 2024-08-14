"""
    solve_poisson(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        ρ
)

FFT Poisson equation solver in 1D. This implementation is not optimized but
maybe clearer to understand.
"""
function solve_poisson(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        ρ;
        fft_type = :fft
)
    (; grid, model) = prob
    (; Nx, Δx) = grid
    (; G, ϵ) = model

    if fft_type == :fft
        kx = 2π * fftfreq(Nx, 1 / Δx)

        ρhat = fft(ρ)

        uhat = zeros(Complex{Float64}, Nx)
    elseif fft_type == :rfft
        kx = 2π * rfftfreq(Nx, 1 / Δx)

        ρhat = rfft(ρ)

        uhat = zeros(Complex{Float64}, Nx ÷ 2 + 1)
    end

    for (i, kxi) in enumerate(kx)
        if kxi == 0
            uhat[i] = 0
        else
            uhat[i] = -4π * G * ρhat[i] / (kxi^2 + ϵ^2)
        end
    end

    if fft_type == :fft
        return real(ifft(uhat))
    elseif fft_type == :rfft
        return real(irfft(uhat, Nx))
    end
end

function setup_fft_cache(grid::Grid1D)
    (; Nx, Δx) = grid

    return (;
        kx = 2π * rfftfreq(Nx, 1 / Δx),
        ρ̂_cache = zeros(ComplexF64, Nx ÷ 2 + 1),
        û_cache = zeros(ComplexF64, Nx ÷ 2 + 1),
        planned_fft = plan_rfft(zeros(Nx)),
        planned_ifft = plan_irfft(zeros(ComplexF64, Nx ÷ 2 + 1), Nx)
    )
end

"""
    solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        cache
)

FFT Poisson equation solver in 1D. This implementation is do not allocate thanks
to the cache.
"""
function solve_poisson!(
        prob::FDProblem{Grid1D, EulerSelfGravity, <:Any, <:Any},
        potentialstore,
        ρ,
        (; kx, planned_fft, planned_ifft, ρ̂_cache, û_cache)
)
    (; model) = prob
    (; G, ϵ) = model

    mul!(ρ̂_cache, planned_fft, ρ)

    for (i, kxi) in enumerate(kx)
        if kxi == 0
            û_cache[i] = 0
        else
            û_cache[i] = -4π * G * ρ̂_cache[i] / (kxi^2 + ϵ^2)
        end
    end

    mul!(potentialstore, planned_ifft, û_cache)
end
