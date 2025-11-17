using LinearAlgebra
using Integrals
using Roots
using Roots: Newton

include("DispersionEquations.jl")



function compute_contour_integral(x,ξ,z,ζ,h,K)
  r = sqrt((x-ξ)^2 + (z-ζ)^2)
  r₁ = sqrt((x-ξ)^2 + (z+ζ)^2)
  X = x - ξ

  function F(R)
    μ = R*exp(im*π/4)
    exp(im*π/4)*cos(μ*X)/cosh(μ*h)*( cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )
  end
  I = Float64(real(quadgk(F, big(0.0), big(20000), rtol=1e-8, order = 50)[1]))
  # I = real(quadgk(F, big(0.0), big(40000), rtol=1e-8, order = 50)[1])
  # I = real(quadgk(F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]) # Fails because of overflow times underflow

  k₀ = first(dispersion_free_surface(K,0,h))
  k = k₀/(-im)
  N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
  Im_φ = -π/(k₀*h*N₀²)*cosh(k*(z+h))*cosh(k*(ζ+h))*cos(k*X)

  return log(r/r₁) - 2*I + Im_φ
end

function compute_eigenexpansion(x,ξ,z,ζ,h,K;N=10)
  kₙ = dispersion_free_surface(K,N,h)
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  return -sum(@. π/(kₙ*h*Nₙ²)*cos(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
end

Ic = compute_contour_integral(1,0.5,-1,0.5,1,1)
I_eigen = compute_eigenexpansion(1,0.5,-1,0.5,1,1;N=100)

Ic - I_eigen