# Dispersion equation for free surface problem: -α = k*tan(kh)
"""
    dispersion_free_surface(α,N,h)

Calculates the positive imaginary and N first solutions with positive real part
of -α = k*tan(k h) for complex α. It uses three methods:
- (a) homotopy for starting with α = 1
- (b) Guess with linear expansion;
- (c) and a linear expansion.

The first roots is positive imaginary and the next are the first N positive real
ordered from smallest.

Ported from MATLAB code by Prof Mike Meylan.
"""
function dispersion_free_surface(α,N,h=1.0)
  α *= h;
  roots = zeros(ComplexF64,N+1)
  # Treated seperately, no special methods needed
  roots[1] = dispersion_free_homotopy(α,0)
  # Find remaining roots
  if N > 0
    i = 1
    while true
      # (a) Homotopy from previous root
      roots[i+1] = dispersion_free_homotopy(α,i);
      # (b) Check approx root
      if abs(roots[i+1] - (im*i*π + α/(im*i*π))) < 0.01
        while true
          roots[i+1] = dispersion_free_root(α,im*i*π + α/(im*i*π));
          if abs(roots[i+1] - (im*i*π + α/(1im*i*π))) < 1e-8
            # (c) Use linear expansion for remaining roots if converged
            for j ∈ i:N
              roots[j+1] = 1im*j*π + α/(1im*j*π);
            end
            i = N;
            break
          end
          i == N && break
          i += 1
        end
      end
      i == N && break
      i += 1
    end
  end
  roots   .*= -im/h
  roots[1] *= -1
  if iszero(N)
    return roots[1:1]
  end
  return roots
end

function dispersion_free_homotopy(α,N)
  f(α) = z -> z*tanh(z) - α
  df(z) = tanh(z) + z*sech(z)^2
  mroot = N == 0 ? 1+0im : im*N*pi
  mroot = dispersion_free_root(1,mroot)

  step = 0.043;
  c = abs(α) < 1 ? -1 : 1
  αstep = 1:c*step:abs(α);
  for k=2:length(αstep)
    mroot = dispersion_free_root(αstep[k],mroot)
  end
  mroot = dispersion_free_root(abs(α),mroot)

  c = angle(α) > 0 ? 1 : -1
  αstep = 0:c*π/30:angle(α)
  for k=2:length(αstep)
    mroot = dispersion_free_root(abs(α)*exp.(im*αstep[k]),mroot)
  end
  mroot = dispersion_free_root(abs(α)*exp(im*angle(α)),mroot)

  return mroot
end

function dispersion_free_root(α,x₀)
  f(z) = z*tanh(z) - α
  df(z) = tanh(z) + z*sech(z)^2
  find_zero((f,df),x₀,Newton());
end