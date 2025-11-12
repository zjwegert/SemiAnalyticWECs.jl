# Based on Appendix A. Modes of vibration for a free-free plate (Meylan et al. 2025)

using Roots
using Roots: Newton

L = 1
N = 10
n = 10
x = -L:L/N:L;

μ = zeros(n)
w = zeros(n,length(x));
d2w = zeros(n,length(x));
for i in 3:n # first two roots are zero
  j = i - 1
  g(x) = (-1)^j*tan(x)+tanh(x);
  dg(x) = (-1)^j*sec(x)^2 + sech(x)^2;
  k = floor(j/2)
  approx_root = π*k-(-1)^j*π/4-exp(-2*π*k-π/2)
  μ[i] = find_zero((g,dg),approx_root,Newton())/L
end

for j ∈ eachindex(μ)
    for k ∈ eachindex(x)
        if j-1 == 0               # 1st (rigid) mode
            w[j,k] = 1/sqrt(2*L);
            d2w[j,k] = μ[j]^2/sqrt(2*L);
        elseif j-1 == 1           # 2nd (rigid) mode
            w[j,k] = x[k]/L*sqrt(3/(2*L));
            d2w[j,k] = μ[j]^2*x[k]/L*sqrt(3/(2*L));
        elseif rem(j-1,2) == 0    # symmetric modes
            w[j,k] = 1/sqrt(2*L)*(cos(μ[j]*x[k])/cos(μ[j]*L)+cosh(μ[j]*x[k])/cosh(μ[j]*L));
            d2w[j,k] = μ[j]^2/sqrt(2*L)*(-cos(μ[j]*x[k])/cos(μ[j]*L)+cosh(μ[j]*x[k])/cosh(μ[j]*L));
        else                      # skew-symmetric modes
            w[j,k] =1/sqrt(2*L)*(sin(μ[j]*x[k])/sin(μ[j]*L)+sinh(μ[j]*x[k])/sinh(μ[j]*L));
            d2w[j,k] =μ[j]^2/sqrt(2*L)*(-sin(μ[j]*x[k])/sin(μ[j]*L)+sinh(μ[j]*x[k])/sinh(μ[j]*L));
        end
    end
end