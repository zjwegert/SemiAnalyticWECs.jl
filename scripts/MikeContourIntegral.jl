
h = 1.0
α = 1.0
f(μ) = μ*tanh(h*μ) - α

μ0 = first(dispersion_free_surface(α,0,h))*-im

angle(μ0)*180/π

f(μ0)