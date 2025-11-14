using Roots
using Roots: Newton
using LinearAlgebra

# Green's function on the surface

"""
    G = matrix_G_surface(α,n,L,H)
Compute G where α is the frequency squared, H is the water depth, L is half the
plate length and 2N+1 is the number of points. The points are evenly space
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
  Gᵥ = green_surface_finite_depth_2d(α,@view(x[2:2n+1]),H)
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
  Gsmall = green_surface_finite_depth_2d(α,R,H);
  while abs(Gsmall[1] - 1/pi*log(R[1]) - Gsmall[2] + 1/pi*log(R[2])) > 1e-3
    R ./= 10;
    Gsmall = green_surface_finite_depth_2d(α,R,H);
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
    green_surface_finite_depth_2d(α,R,H)

Calculates the finite depth Green function for both source and
field point on the free surface in two dimensions. R is the distance between
two points and H is the water depth. The values of R must be ordered from
the smallest to the largest and be positive.

Details can be found on
http://www.wikiwaves.org/index.php/Free-Surface_Green_Function
"""
function green_surface_finite_depth_2d(α,R,H;ϵ=1e-10,N=10)
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