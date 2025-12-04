# Green's function on the surface

"""
    G = matrix_G_surface(α,H,L,n)
Compute G where α is the frequency squared, H is the water depth, L is half the
plate length and 2n+1 is the number of points. The points are evenly space
starting with -L and ending with L (note that because we use Simpson's
rule we require an odd number of points).
"""
function matrix_G_surface(α,H,L,n)
  # Discretisation and Simpson's rule weights
  δx = L/n;
  x = 0:δx:2L;
  w = 4*δx/3*ones(2n+1)
  w[1] = δx/3; w[3:2:2n-1] .*= 1/2; w[2n+1] = δx/3

  ## Compute for j ≠ k
  Gᵥ = greens_surface_2d(α,@view(x[2:2n+1]),H)
  pushfirst!(Gᵥ,0)
  G = zeros(ComplexF64,2*n+1,2*n+1);
  logvalues = zeros(2*n+1,2*n+1);
  for j=1:2*n+1
    for k = 1:j-1
      G[j,k] = w[k]*Gᵥ[abs(j-k)+1];
      logvalues[j,k] = w[k]*log(abs(x[j] - x[k]));
    end
    for k = j+1:2n+1
      G[j,k] = w[k]*Gᵥ[abs(j-k)+1];
      logvalues[j,k] = w[k]*log(abs(x[j] - x[k]));
    end
  end

  ## Compute for j = k
  sumlogvalues = sum(logvalues,dims=2);

  # Find the value of lim R -> 0 G(alpha,R,H) - 1/pi*log(R)
  R = [1e-1,1e-2];
  Gsmall = greens_surface_2d(α,R,H);
  while abs(Gsmall[1] - 1/pi*log(R[1]) - Gsmall[2] + 1/pi*log(R[2])) > 1e-3
    R ./= 10;
    Gsmall = greens_surface_2d(α,R,H);
  end
  lv = Gsmall[2] - 1/pi*log(R[2]);

  G[1,1] = w[1]*lv + 2/pi*(L*log(2*L)-L) - 1/pi*sumlogvalues[1];
  for j=2:2n
      G[j,j] = w[j]*(lv) + 1/pi*((2*L-x[j])*log(2*L-x[j]) - 2*L + x[j]*log(x[j])) - 1/pi*sumlogvalues[j];
  end
  G[2n+1,2n+1] = w[2n+1]*lv + 2/pi*(L*log(2*L)-L) - 1/pi*sumlogvalues[2n+1];

  return G
end

"""
    greens_surface_2d(α,R,H)

Calculates the finite depth Green function for both source and
field point on the free surface in two dimensions. R is the distance between
two points and H is the water depth. The values of R must be ordered from
the smallest to the largest and be positive.

Details can be found on
http://www.wikiwaves.org/index.php/Free-Surface_Green_Function
"""
function greens_surface_2d(α,R,H;ϵ=1e-10,N=10)
  v = zeros(ComplexF64,length(R));
  kₚ = dispersion_free_surface(α,N,H);
  G(r) = k -> exp(-k*r)/(tan(k*H) + H*k*sec(k*H)^2)

  for j = eachindex(R)
    r = R[j]
    v[j] = -sum(G(r),kₚ)
    # Check if need more roots
    last_term = G(r)(kₚ[N+1])
    while abs(last_term / v[j]) > ϵ
        N *= 2;
        kₚ = dispersion_free_surface(α,N,H);
        v[j] = -sum(G(r),kₚ)
        # Check convergence
        last_term = G(r)(kₚ[N+1]);
    end
  end

  return v
end

# Green's function for submerged source

"""
    regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)

Compute the regular Green's function for a submerged source based on Append B.2
of Linton & McIver (2001) (Handbook of mathematical techniques for wave/structure interactions).

In paricular, we compute ϕ - log(r) where r = sqrt(X² + (z-ζ)²).

This can be done either using:
- `method=:residue` : Computes Equation B.38 (Linton & McIver, 2001) using
  contour integration and residue calculus.
- `method=:eigenfunction` : Computes Equation B.43 (Linton & McIver, 2001) using
  eigenfunction expansion method. If this method is used, the number of terms
  in the series is specified by `N`.

Derivatives of ϕ - log(r) can be computed using the following functions:
- `∂z_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)`
- `∂z∂ζ_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)`
"""
function regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)
  if method == :eigenfunction
    return _reg_greens_submerged_2d_eigenfunction_method(x,z,ζ,h,K,N)
  elseif method == :residue
    return _reg_greens_submerged_2d_residue_method(x,z,ζ,h,K)
  else
    error("Unknown method: $method")
  end
end

∂ = ForwardDiff.derivative # Forward AD is quite efficent here, so I think this is cleanest

"""
    ∂z_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)

Compute the z derivative of the regular Green's function for a submerged source based on Append B.2
of Linton & McIver (2001) (Handbook of mathematical techniques for wave/structure interactions).

In paricular, we compute ∂_z(ϕ - log(r)) where r = sqrt(X² + (z-ζ)²) using AD.

This can be done either using:
- `method=:residue` : Computes Equation B.38 (Linton & McIver, 2001) using
  contour integration and residue calculus.
- `method=:eigenfunction` : Computes Equation B.43 (Linton & McIver, 2001) using
  eigenfunction expansion method. If this method is used, the number of terms
  in the series is specified by `N`.
"""
function ∂z_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)
  f(z) = regular_greens_submerged_2d(x,z,ζ,h,K;method,N)
  return ∂(f,z)
end

"""
    ∂z∂ζ_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)

Compute the mixed z-ζ derivative of the regular Green's function for a submerged source based on Append B.2
of Linton & McIver (2001) (Handbook of mathematical techniques for wave/structure interactions).

In paricular, we compute ∂²zζ(ϕ - log(r)) where r = sqrt(X² + (z-ζ)²) using AD.

This can be done either using:
- `method=:residue` : Computes Equation B.38 (Linton & McIver, 2001) using
  contour integration and residue calculus.
- `method=:eigenfunction` : Computes Equation B.43 (Linton & McIver, 2001) using
  eigenfunction expansion method. If this method is used, the number of terms
  in the series is specified by `N`.
"""
function ∂z∂ζ_regular_greens_submerged_2d(x,z,ζ,h,K;method=:residue,N=100)
  g(ζ,z) = regular_greens_submerged_2d(x,z,ζ,h,K;method,N)
  f(ζ) = ∂(z->g(ζ,z),z)
  return ∂(f,ζ)
end

function _reg_greens_submerged_2d_eigenfunction_method(x,z,ζ,h,K,N)
  kₙ = dispersion_free_surface(K,N,h)
  kₙ[1] *= -1
  Nₙ² = @. 1/2*(1 + sin(2*kₙ*h)/(2*kₙ*h))
  ϕ = -sum(@. π/(kₙ*h*Nₙ²)*cos(kₙ*(z+h))*cos(kₙ*(ζ+h))*exp(-kₙ*abs(x)))
  r = sqrt(x^2 + (z-ζ)^2)
  return ϕ - log(r)
end

function _reg_greens_submerged_2d_residue_method(x,z,ζ,h,K)
  k₀ = first(dispersion_free_surface(K,0,h))
  N₀² = @. 1/2*(1 - sin(k₀*h)^2/K/h)
  k = k₀/im
  r₁ = sqrt(x^2 + (z+ζ)^2)

  function g(μ)
    v = (cosh(μ*(z+h))*cosh(μ*(ζ+h))/(μ*sinh(μ*h) - K*cosh(μ*h)) + exp(-μ*h)/μ*sinh(μ*z)*sinh(μ*ζ) )/cosh(μ*h)
    if isnan(v) || isinf(v)
      # Use asymptotic expansion for large μ
      return  exp(μ*(z+ζ))/(μ-K)+1/2/μ*exp(μ*(z+ζ-2*h));
    else
      return v
    end
  end
  # Compute integrals
  F₁(μ) = exp(im*μ*abs(x))*g(μ)
  I₁ = quadgk(R -> exp(im*π/4)*F₁(R*exp(im*π/4)), 0, Inf)[1]
  I₁_residue = π*im/(k*h*N₀²)*cosh(k*(z+h))*cosh(k*(ζ+h))*exp(im*k*abs(x));

  F₂(μ) = exp(-im*μ*abs(x))*g(μ)
  I₂ = quadgk(R -> exp(-im*π/4)*F₂(R*exp(-im*π/4)), 0.0, Inf)[1]

  return - log(r₁) - I₁ - I₂ - I₁_residue
end

# Green's function submerged source (3D)

function _reg_greens_submerged_3d(R,z,ζ,K,H)
  k₀ = first(dispersion_free_surface(K,0,H))
  N₀² = @. 1/2*(1 - sin(k₀*H)^2/K/H)
  k = k₀/im
  r₂ = sqrt(R^2+(z+ζ+2H)^2)

  function g(μ)
    v = (μ + K)*exp(-μ*H)*cosh(μ*(z+H))*cosh(μ*(ζ+H))/(μ*sinh(μ*H) - K*cosh(μ*H))
    if isnan(v) || isinf(v)
      # Use asymptotic expansion for large μ
      return  1/2*(μ + K)*exp(μ*(z + ζ))/(μ - K)
    else
      return v
    end
  end
  # Compute integrals
  F₁(μ) = besselh(0,1,μ*R)*g(μ)
  I₁ = quadgk(R -> exp(im*π/4)*F₁(R*exp(im*π/4)), 0, Inf)[1]
  I₁_residue = π*im/(H*N₀²)*cosh(k*(z+H))*cosh(k*(ζ+H))*besselh(0,1,k*R);

  F₂(μ) = besselh(0,2,μ*R)*g(μ)
  I₂ = quadgk(R -> exp(-im*π/4)*F₂(R*exp(-im*π/4)), 0.0, Inf)[1]

  return 1/r₂ + I₁ + I₂ + I₁_residue
end

function ∂z∂ζ_reg_greens_submerged_3d(R,z,ζ,K,H)
  g(ζ,z) = regular_greens_submerged_3d(R,z,ζ,K,H)
  f(ζ) = ∂(z->g(ζ,z),z)
  return ∂(f,ζ)
end