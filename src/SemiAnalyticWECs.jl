module SemiAnalyticWECs

using Roots
using Roots: Newton
using LinearAlgebra
using ForwardDiff
using QuadGK
using SpecialFunctions

include("Materials.jl")
export PVDF_TechMan_material_coefficents
export PZT5H_material_coefficents

include("DispersionEquations.jl")
export dispersion_free_surface

include("EigenModes1D.jl")
export eigenmodes_1d

include("GreensFunctions.jl")
export matrix_G_surface
export greens_surface_2d
export regular_greens_submerged_2d
export ∂ζ_regular_greens_submerged_2d
export ∂z∂ζ_regular_greens_submerged_2d

include("Solvers2D.jl")
export solve_surface_plate_2d
export solve_submerged_plate_2d

"""
    compute_pierson_moskowitz_efficiency(T_p,ϵ,ωs;α=8.1e-3,g=9.81)

Compute the efficiency of a system based on the Pierson-Moskowitz spectrum.

The function takes:
- the peak period `T_p`;
- a vector of energy absorption values `ϵ = 1 - |R|² - |T|²` for each frequency; and
- a vector of frequencies `ωs`.

The optional parameters `α` and `g` are the Phillips constant and gravitational acceleration, respectively.
"""
function compute_pierson_moskowitz_efficiency(T_p,ϵ,Ts;α=8.1e-3,g=9.81)
  ω_p = 2π/T_p;
  S(ω) = α*g^2/ω^5*exp(-5/4*(ω_p/ω)^4);
  ωs = 2π./Ts
  # Reorder
  jj = sortperm(ωs)
  ϵ = ϵ[jj]
  ωs = ωs[jj]
  # Trapezoidal rule for integration
  E = sum(1/2*(S(ωs[i-1]) + S(ωs[i]))*(ωs[i] - ωs[i-1]) for i in 2:length(ωs))
  E_p = sum(1/2*(ϵ[i-1]*S(ωs[i-1]) + ϵ[i]*S(ωs[i]))*(ωs[i] - ωs[i-1]) for i in 2:length(ωs))
  return E_p/E
end

export compute_pierson_moskowitz_efficiency

end
