using Roots
using Roots: Newton

α = 0.5 + im;
N = 5
# h = 2.0

f(α) = z -> z*tanh(z) - α
df(z) = tanh(z) + z*sech(z)^2

# α *= h;

homotopy(α,N)
f(α)(homotopy(α,N))
f(α)(0.0639 +15.6762im)

function homotopy(alpha,N)
  if N == 0;
    mroot = find_zero((f(1),df),1,Newton());
  else
    mroot = find_zero((f(1),df),im*N*pi,Newton());
  end

  step =0.043;
  if abs(alpha) < 1
      alphastep = ([1:-step:abs(alpha);abs(alpha)]);
  else
      alphastep = ([1:step:abs(alpha);abs(alpha)]);
  end

  for k=2:length(alphastep)
          mroot = find_zero((f(alphastep[k]),df),mroot,Newton());
  end

  if angle(alpha) > 0
      alphastep = abs(alpha)*exp.(im*[0:pi/30:angle(alpha);angle(alpha)]);
  else
      alphastep = abs(alpha)*exp.(im*[0:-pi/30:angle(alpha);angle(alpha)]);
  end

  for k=2:length(alphastep)
    mroot = find_zero((f(alphastep[k]),df),mroot,Newton());
  end

  return mroot
end