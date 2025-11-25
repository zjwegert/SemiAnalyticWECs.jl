using SemiAnalyticWECs
using Test

@test dispersion_free_surface(0.1, 10, 1) ≈ dispersion_free_surface(0.1+0im, 10, 1) # Real vs Complex input
@test dispersion_free_surface(0.1, 10, sqrt(10)) ≈ dispersion_free_surface(0.1+0im, 10, sqrt(10)) # Real vs Complex input

α = 0.1
h = sqrt(10)
k = dispersion_free_surface(α, 100, h)
f(k) = α + k*tan(k*h)
@test all(@. abs(f(k)) < 1e-10)

α = 0.1+0.1im
h = sqrt(10)
k = dispersion_free_surface(α, 10, h)
f(k) = α + k*tan(k*h)
@test all(@. abs(f(k)) < 1e-10)