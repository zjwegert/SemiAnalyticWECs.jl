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
export ∂z_regular_greens_submerged_2d
export ∂z∂ζ_regular_greens_submerged_2d

include("Solvers2D.jl")
export solve_surface_plate_2d
export solve_submerged_plate_2d

end
