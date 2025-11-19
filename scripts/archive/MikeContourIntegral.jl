using LinearAlgebra
using QuadGK
using Roots
using Roots: Newton

include("DispersionEquations.jl")

### Works, but commiting some mathematical sins!
# """
#   Compute Green's function for submerged source using residue calculus.

#     ϕ = ... (B.38, Linton 2001)

#   where kₙ are roots of dispersion equation for free surface waves.
# """
# function compute_contour_integral(x,ξ,z,ζ,h,K)
#   r = sqrt((x-ξ)^2 + (z-ζ)^2)
#   r₁ = sqrt((x-ξ)^2 + (z+ζ)^2)
#   X = x - ξ

#   function F(R)
#     μ = R*exp(im*π/100)
#     exp(im*π/100)*cos(μ*X)/cosh(μ*h)*( cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )
#   end
#   I = Float64(real(quadgk(F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]))

#   k₀ = -first(dispersion_free_surface(K,0,h))
#   k = k₀/im
#   N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
#   Im_φ = -π/(k₀*h*N₀²)*cosh(k*(z+h))*cosh(k*(ζ+h))*cos(k*X)

#   return log(r/r₁) - 2*I + Im_φ
# end

# function compute_contour_integral_partials(x,ξ,z,ζ,h,K)
#   X = x - ξ

#   function ∂²zζ_F(R)
#     μ = R*exp(im*π/100)
#     exp(im*π/100)*cos(μ*X)/cosh(μ*h)*(μ^2*sinh(μ*(z+h))*sinh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)*μ*cosh(μ*z)*cosh(μ*ζ) )
#   end
#   ∂I = Float64(real(quadgk(∂²zζ_F, big(0.0), big(Inf), rtol=1e-12, order = 21)[1]))
#   k₀ = -first(dispersion_free_surface(K,0,h))
#   k = k₀/im
#   N₀² = 1/2*(1 + sin(2*k₀*h)/(2*k₀*h))
#   ∂Im_φ = -k^2*π/(k₀*h*N₀²)*sinh(k*(z+h))*sinh(k*(ζ+h))*cos(k*X)

#   # ∂log = 2*(z-ζ)^2/((x-ξ)^2 + (z-ζ)^2)^2 + 2*(z+ζ)^2/((x-ξ)^2 + (z+ζ)^2)^2 -
#   #   1/((x-ξ)^2 + (z-ζ)^2) - 1/((x-ξ)^2 + (z+ζ)^2)
#   r₁² = (x-ξ)^2 + (z+ζ)^2
#   ∂log = 2*(z+ζ)^2/r₁²^2 - 1/r₁²

#   # r² = (x-ξ)^2 + (z-ζ)^2
#   # ∂logr = 2*(z-ζ)^2/r²^2 - 1/r²

#   # return ∂log - ∂logr - 2*∂I + ∂Im_φ
#   return ∂log - 2*∂I + ∂Im_φ
# end

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

function compute_eigenexpansion_removed_singularity(x,ξ,z,ζ,h,K;N=10)
  kₙ = dispersion_free_surface(K,N,h)
  kₙ[1] *= -1
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  ∂z∂ζ_ϕ = -sum(@. kₙ*π/(h*Nₙ²)*sin(kₙ*(z+h))*sin(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  r² = (x-ξ)^2 + (z-ζ)^2
  ∂logr = 2*(z-ζ)^2/r²^2 - 1/r²

  ∂z_ϕ = sum(@. π/(h*Nₙ²)*sin(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
  ∂z_logr = (z-ζ)/r²
  return ∂z∂ζ_ϕ - ∂logr, ∂z_ϕ - ∂z_logr
end

# Ic = compute_contour_integral(1.0,0,-1,-1,10,1)
∂Ic = compute_contour_integral_partials(1.0,0,-1,-1,10,1)

# I_eigen, _ = compute_eigenexpansion(1.0,0,-1,-1,10,1;N=100)

∂I_eigen,∂zI_eigen = compute_eigenexpansion_removed_singularity(1,0,-1,-1,10,1;N=100)

∂I_eigen - compute_contour_integral_partials(1,-1,-1,10,1)
∂zI_eigen

# Ic - I_eigen
# abs(∂Ic - ∂I_eigen)

@btime compute_contour_integral_partials(1,-1,-1,10,1);

using ForwardDiff
using BenchmarkTools
@btime ForwardDiff.derivative(z -> compute_contour_integral_partials(1.0,z,-1,10,1), -1.0)

green_submerged_finite_depth_2d

∂I_eigen,∂zI_eigen = compute_eigenexpansion_removed_singularity(2,0,-1,-3,10,1;N=100)

compute_contour_integral_partials(2,-1,-3,10,1)
∂zI_eigen
∂z_regular_greens_submerged_2d(2,-1,-3,10,5)
compute_contour_integral_partials_z(2,-1,-3,10,1)

function compute_contour_integral_partials(x,z,ζ,h,K)
  k₀ = first(dispersion_free_surface(K,0,h))
  N₀² = @. 1/2*(1 - sin(k₀*h)^2/K/h)
  k = k₀/im
  r₁² = x^2 + (z+ζ)^2

  function ∂²zζ_g(μ)
    v = (μ^2*sinh(μ*(z+h))*sinh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)*μ*cosh(μ*z)*cosh(μ*ζ) )/cosh(μ*h)
    if isnan(v) || isinf(v)
      # Use asymptotic expansion for large μ
      return  μ^2*exp(μ*(z+ζ))/(μ-K)+1/2*μ*exp(μ*(z+ζ-2*h));
    else
      return v
    end
  end
  # Compute integrals
  F₁(μ) = exp(im*μ*abs(x))*∂²zζ_g(μ)
  I₁ = quadgk(R -> exp(im*π/4)*F₁(R*exp(im*π/4)), 0, Inf)[1]
  I₁_residue = k*π*im/(h*N₀²)*sinh(k*(z+h))*sinh(k*(ζ+h))*exp(im*k*abs(x));

  F₂(μ) = exp(-im*μ*abs(x))*∂²zζ_g(μ)
  I₂ = quadgk(R -> exp(-im*π/4)*F₂(R*exp(-im*π/4)), 0.0, Inf)[1]
  # Compute derivative of logarithmic terms
  ∂logr1 = 1/r₁² - 2*(z+ζ)^2/r₁²^2

  return - ∂logr1 - I₁ - I₂ - I₁_residue
end

function compute_contour_integral_partials_z(x,z,ζ,h,K)
  k₀ = first(dispersion_free_surface(K,0,h))
  N₀² = @. 1/2*(1 - sin(k₀*h)^2/K/h)
  k = k₀/im
  r₁² = x^2 + (z+ζ)^2

  function ∂z_g(μ)
    v = (μ*sinh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)*cosh(μ*z)*sinh(μ*ζ) )/cosh(μ*h)
    if isnan(v) || isinf(v)
      # Use asymptotic expansion for large μ
      return  μ*exp(μ*(z+ζ))/(μ-K)+1/2*exp(μ*(z+ζ-2*h));
    else
      return v
    end
  end
  # Compute integrals
  F₁(μ) = exp(im*μ*abs(x))*∂z_g(μ)
  I₁ = quadgk(R -> exp(im*π/4)*F₁(R*exp(im*π/4)), 0, Inf)[1]
  I₁_residue = π*im/(h*N₀²)*sinh(k*(z+h))*cosh(k*(ζ+h))*exp(im*k*abs(x));

  F₂(μ) = exp(-im*μ*abs(x))*∂z_g(μ)
  I₂ = quadgk(R -> exp(-im*π/4)*F₂(R*exp(-im*π/4)), 0.0, Inf)[1]
  # Compute derivative of logarithmic terms
  ∂logr1 = (z+ζ)/r₁²

  return - ∂logr1 - I₁ - I₂ - I₁_residue
end