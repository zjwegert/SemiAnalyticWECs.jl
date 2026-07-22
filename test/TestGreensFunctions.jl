using SemiAnalyticWECs
using Test 

@test regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ ∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ (-0.041294670059438 - 0.022018584405564im) # Test against matlab code

# FDM
using FiniteDiff
fd_1 = FiniteDiff.finite_difference_derivative(ζ->log(sqrt(5^2 + ((-3)-ζ)^2)) + regular_greens_submerged_2d(5,-3,ζ,10,1.1;method=:residue), 0.0)
fd_2 = FiniteDiff.finite_difference_derivative(ζ->log(sqrt(5^2 + ((-3)-ζ)^2)) + regular_greens_submerged_2d(5,-3,ζ,10,1.1;method=:eigenfunction), 0.0)
ad_1 = ∂ζ_regular_greens_submerged_2d(5,-3,0,10,1.1;method=:residue)
ad_2 = ∂ζ_regular_greens_submerged_2d(5,-3,0,10,1.1;method=:eigenfunction,N=100)

@test fd_1 ≈ ad_1
@test fd_2 ≈ ad_2