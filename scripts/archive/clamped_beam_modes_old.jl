using Roots
using Roots: Newton
using LinearAlgebra
using Arpack
using FastGaussQuadrature
using ForwardDiff

L = 1.5
N = 10
n = 15
x = -L:L/N:L;

M(μ) = [
     cosh(L*μ)  -sinh(L*μ)    cos(L*μ)  -sin(L*μ)
     cosh(L*μ)   sinh(L*μ)    cos(L*μ)   sin(L*μ)
  -sinh(L*μ) cosh(L*μ)  sin(L*μ) cos(L*μ)
   sinh(L*μ) cosh(L*μ) -sin(L*μ) cos(L*μ)
]

β(μ,(a,b,c,d)) = 1/(2*μ)*(2*(a^2-b^2+c^2+d^2)*L*μ+4*(a*c+b*d)*cosh(L*μ)*sin(L*μ)+(c-d)*(c+d)*sin(2*L*μ)+4*(a*c-b*d)*cos(L*μ)*sinh(L*μ)+(a^2+b^2)*sinh(2*L*μ))

μ = zeros(n)
w = zeros(n,length(x));
d2w = zeros(n,length(x));
α = [zeros(Complex,4) for _ in 1:n];
for i in 1:n # first two roots are zero
  j = i + 1
  g(x) = (-1)^j*tan(x)+tanh(x);
  dg(x) = (-1)^j*sec(x)^2 + sech(x)^2;
  k = floor(j/2)
  approx_root = π*k-(-1)^j*π/4-exp(-2*π*k-π/2)
  μ[i] = find_zero((g,dg),approx_root,Newton())/L
  λ, V = eigen(M(μ[i]));
  P = sortperm(abs.(λ))
  α[i] .= V[:,P[1]]
end

## Equiv
# h(x) = -sech(x)+cos(x);
# dh(x) = sech(x)*tanh(x)-sin(x);
# μ2 = zeros(n)
# for i in 3:n # first two roots are zero
#   j = i - 1
#   k = floor(j/2)
#   approx_root = 2*(π*k-(-1)^j*π/4-exp(-2*π*k-π/2))
#   μ2[i] = find_zero((h,dh),approx_root,Newton())/2/L
# end

α_hat = 1 ./ sqrt.(β.(μ,α)) .* α;

u(x,μ,α) = @. α[1]*cosh(μ*x) + α[2]*sinh(μ*x) + α[3]*cos(μ*x) + α[4]*sin(μ*x);
d2u(x,μ,α) = @. μ^2*(α[1]*cosh(μ*x) + α[2]*sinh(μ*x) - α[3]*cos(μ*x) - α[4]*sin(μ*x));

U = zeros(Complex,n,length(x))
d2U = zeros(Complex,n,length(x));
for i ∈ eachindex(μ)
  U[i,:] = u(x,μ[i],α_hat[i])
  d2U[i,:] = d2u(x,μ[i],α_hat[i])
end

U

# Check orthonormality
xg,wg=gausslegendre(100)
O = zeros(n,n);
for i ∈ 1:n
  for j ∈ 1:n
    O[i,j] = sqrt(abs(wg ⋅ (u(xg,μ[i],α_hat[i]).*u(xg,μ[j],α_hat[j]))))
  end
end
O
# Check derivative
d2U_ad = map(i->map(y->ForwardDiff.derivative(x->ForwardDiff.derivative(x->u(x,μ[i],α_hat[i]),x),y),x),eachindex(μ))
maximum.(abs,eachrow(d2U) - d2U_ad)
# Check PDE
d4U_ad = map(i->map(y->ForwardDiff.derivative(x->ForwardDiff.derivative(x->ForwardDiff.derivative(x->ForwardDiff.derivative(x->u(x,μ[i],α_hat[i]),x),x),x),y),x),eachindex(μ))
maximum.(abs,d4U_ad - eachrow(U).*μ.^4)