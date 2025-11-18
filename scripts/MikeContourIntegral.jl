using LinearAlgebra
using QuadGK
using Roots
using Roots: Newton

include("DispersionEquations.jl")

"""
  Compute Green's function for submerged source using residue calculus.

    ϕ = ... (B.38, Linton 2001)

  where kₙ are roots of dispersion equation for free surface waves.
"""
function compute_contour_integral(x,ξ,z,ζ,h,K)
  r = sqrt((x-ξ)^2 + (z-ζ)^2)
  r₁ = sqrt((x-ξ)^2 + (z+ζ)^2)
  X = x - ξ

  function F(R)
    μ = R*exp(im*π/100)
    exp(im*π/100)*cos(μ*X)/cosh(μ*h)*( cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )
  end
  I = Float64(real(quadgk(F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]))

  k₀ = first(dispersion_free_surface(K,0,h))
  k = k₀/(-im)
  N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
  Im_φ = -π/(k₀*h*N₀²)*cosh(k*(z+h))*cosh(k*(ζ+h))*cos(k*X)

  return log(r/r₁) - 2*I + Im_φ
end

function compute_contour_integral_partials(x,ξ,z,ζ,h,K)
  X = x - ξ

  function ∂²zζ_F(R)
    μ = R*exp(im*π/100)
    exp(im*π/100)*cos(μ*X)/cosh(μ*h)*(μ^2*sinh(μ*(z+h))*sinh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)*μ*cosh(μ*z)*cosh(μ*ζ) )
  end
  ∂I = Float64(real(quadgk(∂²zζ_F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]))
  k₀ = first(dispersion_free_surface(K,0,h))
  k = k₀/(-im)
  N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
  ∂Im_φ = -k^2*π/(k₀*h*N₀²)*sinh(k*(z+h))*sinh(k*(ζ+h))*cos(k*X)

  ∂log = 2*(z-ζ)^2/((x-ξ)^2 + (z-ζ)^2)^2 + 2*(z+ζ)^2/((x-ξ)^2 + (z+ζ)^2)^2 -
    1/((x-ξ)^2 + (z-ζ)^2) - 1/((x-ξ)^2 + (z+ζ)^2)

  return ∂log - 2*∂I + ∂Im_φ
end

"""
  Compute eigenfunction expansion of Green's function for submerged source.

    ϕ = -Σ{ π/(kₙNₙ²h)cos(kₙ(z+h))*cos(kₙ(ζ+h))exp(-kₙ|x-ξ|) }
    ∂²_zζϕ = -Σ{ kₙπ/(Nₙ²h)sin(kₙ(z+h))*sin(kₙ(ζ+h))exp(-kₙ|x-ξ|) }

  where kₙ are roots of dispersion equation for free surface waves.
"""
function compute_eigenexpansion(x,ξ,z,ζ,h,K;N=10)
  kₙ = dispersion_free_surface(K,N,h)
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  ϕ = -sum(@. π/(kₙ*h*Nₙ²)*cos(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  ∂z∂ζ_ϕ = -sum(@. kₙ*π/(h*Nₙ²)*sin(kₙ*(z+h))*sin(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  return ϕ, ∂z∂ζ_ϕ
end

Ic = compute_contour_integral(1.0,0,-1,-1,10,1)
∂Ic = compute_contour_integral_partials(1.0,0,-1,-1,10,1)

I_eigen, ∂I_eigen = compute_eigenexpansion(1.0,0,-1,-1,10,1;N=100)

Ic - I_eigen
abs(∂Ic - ∂I_eigen)