
h = 1.0
α = 1.0
f(μ) = μ*tanh(h*μ) - α

μ0 = first(dispersion_free_surface(α,0,h))*-im

angle(μ0)*180/π

f(μ0)

using FastGaussQuadrature
using LinearAlgebra

function compute_φ(x,ξ,z,ζ,h,K;order=20)
  r = sqrt((x-ξ)^2 + (z-ζ)^2)
  r₁ = sqrt((x-ξ)^2 + (z+ζ)^2)
  X = x - ξ

  function F(r)
    μ = r*(1+im)
    (1+im)*real(exp(im*μ*X))/cosh(μ*h)*( cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )
  end
  R,w=gausslaguerre(order)
  I = (F.(R)⋅w)

  return log(r/r₁) .- 2*I
end

# compute_φ(1,0.5,-1,0.5,1,1;order=10)
# compute_φ(1,0.5,-1,0.5,1,1;order=30)
# compute_φ(1,0.5,-1,0.5,1,1;order=50)
# compute_φ(1,0.5,-1,0.5,1,1;order=60)
# compute_φ(1,0.5,-1,0.5,1,1;order=80)
compute_φ(1,0.5,-1,0.5,1,1;order=100)

function compute_eigenexpansion(x,ξ,z,ζ,h,K;N=10)
  kₙ = dispersion_free_surface(K,N,h)
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  return -sum(@. π/(kₙ*h*Nₙ²)*cos(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x-ξ)))
end

compute_eigenexpansion(1,0.5,-1,0.5,1,1;N=100)