using SemiAnalyticWECs
using Test

@test regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂z_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ ∂z_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:eigenfunction,N=100)
@test ∂z_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ (0.085855631870347 - 0.020016903112010im) # Test against matlab code
@test ∂z∂ζ_regular_greens_submerged_2d(5,-2,-3,10,1.1;method=:residue) ≈ (-0.041294670059438 - 0.022018584405564im) # Test against matlab code