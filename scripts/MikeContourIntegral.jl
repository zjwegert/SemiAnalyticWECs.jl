using LinearAlgebra
using QuadGK
using Roots
using Roots: Newton

include("DispersionEquations.jl")

function compute_contour_integral(x,ξ,z,ζ,h,K)
  r = sqrt((x-ξ)^2 + (z-ζ)^2)
  r₁ = sqrt((x-ξ)^2 + (z+ζ)^2)
  X = x - ξ

  function F(R)
    μ = R*exp(im*π/100)
    exp(im*π/100)*cos(μ*X)/cosh(μ*h)*( cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )
  end
  I = (real(quadgk(F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]))

  k₀ = first(dispersion_free_surface(K,0,h))
  k = k₀/(-im)
  N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
  Im_φ = -π/(k₀*h*N₀²)*cosh(k*(z+h))*cosh(k*(ζ+h))*cos(k*X)

  return log(r/r₁) - 2*I + Im_φ
end

function compute_eigenexpansion(x,ξ,z,ζ,h,K;N=10)
  kₙ = dispersion_free_surface(K,N,h)
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  ϕ = -sum(@. π/(kₙ*h*Nₙ²)*cos(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  ∂_z_ϕ = sum(@. π/(Nₙ²)*sin(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  ∂_ζ_ϕ = sum(@. π/(Nₙ²)*cos(kₙ*(z+h))*sin(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  return ϕ, ∂_z_ϕ, ∂_ζ_ϕ
end

using ForwardDiff
Ic = ComplexF64(compute_contour_integral(1,0.5,-1,0.5,1,1))
∂_z_Ic = ComplexF64(ForwardDiff.derivative(z -> compute_contour_integral(1,0.5,z,0.5,1,1), -1))
∂_ζ_Ic = ComplexF64(ForwardDiff.derivative(ζ -> compute_contour_integral(1,0.5,-1,ζ,1,1), 0.5))

I_eigen,∂_z_I_eigen, ∂_ζ_I_eigen = compute_eigenexpansion(1,0.5,-1,0.5,1,1;N=100)

Ic - I_eigen
∂_z_Ic - ∂_z_I_eigen
∂_ζ_Ic - ∂_ζ_I_eigen