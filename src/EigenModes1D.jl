"""
    eigenmodes_1d(bc_case,N,L,x)

Calculates the first N eigenmodes and eigenvalues of the PDE:
  ∂⁴ₓ u = μ⁴ u,  x ∈ [-L,L]
with boundary conditions specified by `bc_case` at x = -L and x = L
at as tuple. E.g., `bc_case = (:clamped,:clamped)` or equivalently
`bc_case = :clamped`. x is the vector of spatial points where the
modes are evaluated.

Possible boundary conditions are:
- :clamped  => u = 0, ∂ₓu = 0
- :free     => ∂²ₓu = 0, ∂³ₓu = 0
- :simply_supported => u = 0, ∂²ₓu = 0 (at present not implemented)
"""
function eigenmodes_1d(bc_case,N,L,x;debug=false)
  if bc_case == (:clamped,:clamped) || bc_case == :clamped
    return _eigenmodes_clamped_clamped(N,L,x;debug)
  elseif bc_case == (:free,:free) || bc_case == :free
    return _eigenmodes_free_free(N,L,x)
  else
    @error "Boundary condition case $(bc_case), not yet implemented."
  end
end

# free-free

function _eigenmodes_free_free(N,L,x)
  μ = zeros(eltype(x),N)
  U = zeros(eltype(x),N,length(x));
  ∂ₓ²U = zeros(eltype(x),N,length(x));

  # Find solutions of det(M(μᵢ)) = 0
  for i in 3:N # first two roots are zero
    j = i - 1
    g(x) = (-1)^j*tan(x)+tanh(x);
    dg(x) = (-1)^j*sec(x)^2 + sech(x)^2;
    k = floor(j/2)
    approx_root = π*k-(-1)^j*π/4-exp(-2*π*k-π/2)
    μ[i] = find_zero((g,dg),approx_root,Newton())/L
  end

  # We know a nice form for the modes in this case
  for j ∈ eachindex(μ)
    for k ∈ eachindex(x)
      if j-1 == 0               # 1st (rigid) mode
        U[j,k] = 1/sqrt(2*L);
        ∂ₓ²U[j,k] = μ[j]^2/sqrt(2*L);
      elseif j-1 == 1           # 2nd (rigid) mode
        U[j,k] = x[k]/L*sqrt(3/(2*L));
        ∂ₓ²U[j,k] = μ[j]^2*x[k]/L*sqrt(3/(2*L));
      elseif rem(j-1,2) == 0    # symmetric modes
        U[j,k] = 1/sqrt(2*L)*(cos(μ[j]*x[k])/cos(μ[j]*L)+cosh(μ[j]*x[k])/cosh(μ[j]*L));
        ∂ₓ²U[j,k] = μ[j]^2/sqrt(2*L)*(-cos(μ[j]*x[k])/cos(μ[j]*L)+cosh(μ[j]*x[k])/cosh(μ[j]*L));
      else                      # skew-symmetric modes
        U[j,k] =1/sqrt(2*L)*(sin(μ[j]*x[k])/sin(μ[j]*L)+sinh(μ[j]*x[k])/sinh(μ[j]*L));
        ∂ₓ²U[j,k] =μ[j]^2/sqrt(2*L)*(-sin(μ[j]*x[k])/sin(μ[j]*L)+sinh(μ[j]*x[k])/sinh(μ[j]*L));
      end
    end
  end

  return μ,U,∂ₓ²U
end

# clamped-clamped

function _eigenmodes_clamped_clamped(N,L,x;debug)
  u(x,μ,α) = @. α[1]*exp(μ*x) + α[2]*exp(-μ*x) + α[3]*cos(μ*x) + α[4]*sin(μ*x);
  ∂ₓ²u(x,μ,α) = @. μ^2*(α[1]*exp(μ*x) + α[2]*exp(-μ*x) - α[3]*cos(μ*x) - α[4]*sin(μ*x));

  M(μ) = [
    exp(-L*μ)	   exp(L*μ)	       cos(L*μ)	  -sin(L*μ) #   u(−L)=0
    exp(L*μ)	   exp(-L*μ)	     cos(L*μ)	   sin(L*μ) #    u(L)=0
    exp(-L*μ)*μ	-exp(L*μ)*μ	   μ*sin(L*μ)	 μ*cos(L*μ) # ∂ₓu(−L)=0
    exp(L*μ)*μ	-exp(-L*μ)*μ  -μ*sin(L*μ)  μ*cos(L*μ) #  ∂ₓu(L)=0
  ]

  β(μ,(a,b,c,d)) = 1/(2*μ)*(2*(4*a*b+c^2+d^2)*L*μ+4*(b*(c-d)+a*(c+d))*cosh(L*μ)*sin(L*μ)+
    (c-d)*(c+d)*sin(2*L*μ)+4*(a*(c-d)+b*(c+d))*cos(L*μ)*sinh(L*μ)+2*(a^2+b^2)*sinh(2*L*μ))

  μ = zeros(eltype(x),N)
  α = [zeros(eltype(x),4) for _ in 1:N];
  for i in 1:N
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

  U = zeros(eltype(x),N,length(x))
  ∂ₓ²U = zeros(eltype(x),N,length(x));
  α_hat = @. α/sqrt(β(μ,α));
  for i ∈ eachindex(μ)
    U[i,:] = u(x,μ[i],α_hat[i])
    ∂ₓ²U[i,:] = ∂ₓ²u(x,μ[i],α_hat[i])
  end

  if debug
    return μ,U,∂ₓ²U,α_hat
  else
    return μ,U,∂ₓ²U
  end
end