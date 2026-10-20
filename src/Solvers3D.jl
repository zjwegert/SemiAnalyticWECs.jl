"""
    solve_submerged_plate_3d(ω::Real, A::Real, χ::Real, H::Real,
      h::Real, ρ::Real, d::Real, vox_trian::Triangulation,
      δx::Real, mode_data::NamedTuple,
      cache::Tuple; g::Real=9.81, ρ_w::Real=1025
    )
Solve the 3D wave-plate interaction problem.

Inputs:
- ω: Frequency
- A: Incoming wave amplitude
- χ: Incoming wave angle
- H: Water depth
- h: Submergence depth
- ρ: Plate density
- d: Plate thickness
- vox_trian: Voxelised triangulation for solving problem
- δx: Voxel side length
- mode_data: Dry mode information computed using `eigenmodes_2d`
- cache: Cache computed via `precompute_cache`. Removing this parameter will
    regenerate the cache for each call.

Optional Inputs:
- g: gravitational acceleration (default 9.81 m/s²)
- ρ_w: water density (default 1025 kg/m³)

Outputs:
- W: Plate displacement
- f: Optical function
- scs: Left-hand side of optical theorem
- ecs: Right-hand side of optical theorem
- P_farfield: Power according to farfield
- P_nearfield: Power according to nearfield
"""
function solve_submerged_plate_3d(ω::Real, A::Real, χ::Real, H::Real,
    h::Real, ρ::Real, d::Real, η11::Real, η12::Real, η22::Real,
    G::Real, C::Real, vox_trian::Triangulation,
    δx::Real, mode_data::NamedTuple,
    cache::Tuple; g::Real=9.81, ρ_w::Real=1025
  )
  ωi, Ui, ∂xxUi, ∂xyUi, ∂yyUi = mode_data
  R0, Xd, Yd, xd, yd, M2 = cache
  inout_idx = vox_trian.tface_to_mface

  # Wave number
  α = ω^2/g;
  k = first(dispersion_free_surface(α,0,H))/im

  # Construct kernel
  M1 = construct_M1_matrix(R0,α,H,h,δx)
  K0 = -1/(4*π)*(M1+M2);
  K = build_kernel(K0, xd, yd, inout_idx);

  # Solve system
  Xd_in = Xd[inout_idx]
  Yd_in = Yd[inout_idx]
  f_vec = @. -im*g*A/(ω*cosh(k*H))*k*sinh(k*(H-h))*exp(im*k*(Xd_in*cos(χ)+Yd_in*sin(χ)));
  rhs = [f_vec im*ω*Ui];
  sol = K\rhs;
  phi_d_jump = sol[:,1];
  phi_ra_jump = sol[:,2:end];

  # dry mode coefficients
  Modes_IP_matrix = Ui'*(Ui*δx^2);
  M_stiffness = ρ*d*Modes_IP_matrix*Diagonal(ωi.^2);
  M_mass = ρ*d*ω^2*Modes_IP_matrix;
  M_added_mass = im*ω*ρ_w*Ui'*(phi_ra_jump*δx^2);
  f_vec = -im*ω*ρ_w*Ui'*(phi_d_jump*δx^2);
  c = (M_stiffness-M_mass+M_added_mass)\f_vec;

  # Construct displacement and Hessian terms
  W = Ui*c
  ∂xxW = ∂xxUi*c
  ∂xyW = ∂xyUi*c
  ∂yyW = ∂yyUi*c

  # Total jump in potential
  sigma = phi_d_jump+phi_ra_jump*c;

  # Conservation of energy (plane wave only)
  r_x = @. sqrt(Xd_in^2+Yd_in^2);
  theta_x = @. atan(Yd_in,Xd_in);
  f(θ) = -k^3*sinh(k*H)*sinh(k*(H-h))/(ω*(2*k*H+sinh(2*k*H)))*transpose(sigma)*(exp.(-im*k*Diagonal(r_x)*cos.(theta_x.-θ))*δx^2);
  scs = 1/π*hquadrature(abs2 ∘ f, χ-π,χ+π)[1]
  ecs = -2real(f(χ))

  # Far-field power
  Cg = ω/(2*k)*(1+2*k*H/sinh(2*k*H))
  P_farfield = ρ_w*g*A^2/k*Cg*(-scs+ecs)

  # Near-field power
  ηij∂ijW = η11*∂xxW + 2*η12*∂xyW + η22*∂yyW
  P_nearfield = ω^2*G/(2(G^2+ω^2*C^2))*(ηij∂ijW'*ηij∂ijW)*δx^2

  return W,f,scs,ecs,P_farfield,P_nearfield,sigma
end

function solve_submerged_plate_3d(ω::Real, A::Real, χ::Real, H::Real,
    h::Real, ρ::Real, d::Real, η11::Real, η12::Real, η22::Real,
    G::Real, C::Real, vox_trian::Triangulation,
    δx::Real, mode_data::NamedTuple; kwargs...
  )
  vox_model = get_background_model(vox_trian)
  cache = precompute_cache(vox_model,δx)
  return solve_submerged_plate_3d(ω, A, χ, H, h, ρ, d,
    η11, η12, η22, G, C, vox_trian,
    δx, mode_data, cache; kwargs...
  )
end

# Precompute data
"""
  precompute_cache(vox_model,δx)

Precompute `cache` for solving submerged 3d problem.
"""
function precompute_cache(vox_model,δx)
  R0, Xd, Yd, xd, yd = construct_R0(vox_model)
  M2 = construct_M2_matrix(R0,Xd,Yd,xd,yd,δx)
  return R0, Xd, Yd, xd, yd, M2
end

# Construct distance matrix to (x[1],y[1]) and mid-point information
function construct_R0(vox_model)
  bg_midpoints = map(mean, get_cell_coordinates(vox_model))
  midpoints = vec(bg_midpoints)
  Xd = getindex.(midpoints,1);
  Yd = getindex.(midpoints,2);
  xd = getindex.(bg_midpoints,1)[:,1]
  yd = getindex.(bg_midpoints,2)[1,:]
  R0 = @. sqrt((Xd-xd[1])^2+(Yd-yd[1])^2); # distances to (x[1],y[1])
  R0 = reshape(R0,length(xd),length(yd));
  return R0, Xd, Yd, xd, yd
end

# Construct reduced-kernel matrix M2 for 3D problem
function construct_M2_matrix(R0,Xd,Yd,xd,yd,δx)
  M2 = zero(R0);
  Threads.@threads for ii ∈ eachindex(M2)
    if iszero(R0[ii])
      M2[ii] = -8*sqrt(2)/δx;
    else
      M2[ii] = r3_int(Xd[ii],Yd[ii],δx,xd[1],yd[1]);
    end
  end

  return M2
end

function r3_int(xp,yp,delta_x,a,b)
  x0 = xp - delta_x/2;
  y0 = yp - delta_x/2;
  x1 = xp + delta_x/2;
  y1 = yp + delta_x/2;
  r((x,y),) = sqrt((x-a)^2 + (y-b)^2);
  f(x) = 1/r(x)^3;
  return first(hcubature(f,(x0,y0),(x1,y1)));
end

# Construct reduced-kernel matrix M1 for 3D problem
function construct_M1_matrix(R0,α,H,h,δx)
  ∂z∂ζG(x) = ∂z∂ζ_regular_greens_submerged_3d(norm(x,2),-h,-h,α,H);
  M1 = zeros(ComplexF64,size(R0));
  Threads.@threads for ii ∈ eachindex(M1)
    if iszero(R0[ii])
      M1[ii] = 4*hcubature(∂z∂ζG,(0,0),(δx/2,δx/2))[1];
    else
      # midpoint
      M1[ii] = ∂z∂ζ_regular_greens_submerged_3d(R0[ii],-h,-h,α,H)*δx^2;
    end
  end

  return M1
end

# Build full kernal matrix K from reduced matrix M for 3D problem
function build_kernel(M::AbstractMatrix{V}, xd::AbstractVector{T}, yd::AbstractVector{T}, inout_idx::AbstractArray{T2}) where {V,T,T2}
    nxg, nyg = length(xd), length(yd)
    ix = repeat(1:nxg, outer=nyg)
    iy = repeat(1:nyg, inner=nxg)
    Nin = length(inout_idx)
    ixv = ix[inout_idx]
    iyv = iy[inout_idx]
    S = Matrix{V}(undef, Nin, Nin)
    Threads.@threads for q in 1:Nin
        ixq, iyq = ixv[q], iyv[q]
        @inbounds @simd for p in 1:Nin
            S[p, q] = M[abs(ixv[p] - ixq) + nxg*abs(iyv[p] - iyq) + 1]
        end
    end
    return S
end