"""
    eigenmodes_2d(bc_case,D,ρ_d,N,model,vox_trian,cache;order=2,μ=order*(order+1),dirichlet_tags=["∂Γ"])

Calculates the first N eigenmodes and eigenvalues of the PDE:
  ∇∇:(D:∇∇u) = ρh*ω²*u on Γ
with boundary conditions `bc_case ∈ (:clamped, :free, :simply_supported)` using
C/DG FEM. Returns modes and second-order derivatives evaluated on the voxelised mesh `vox_trian`.

Homogenous clamped or simply supported boundary conditions are applied on the Dirichlet boundary
specified by `dirichlet_tags`. When these tags only mark a subset of the boundary, the remaining
part of the boundary is treated as free.

Parameters:
- `D`: A symmetric 4-tensor of dimension 2
- `ρ_d`: A function or constant representing the density times thickness
- `N`: The number of eigenmodes to calculate
- `model`: A Gridap `DiscreteModel`. This should have a vector of tags `dirichlet_tags` associated with the Dirichlet boundary.
- `vox_trian`: A Gridap `Triangulation` representing the voxelised mesh.
- `cache`: A precomputed cache of triangulation and operators, obtained from `precompute_eigenmode_cache`.

Optional Parameters:
- `order`: The order of the finite element space (must be ≥2)
- `μ`: The penalty parameter (default is `order*(order+1)`)
- `dirichlet_tags`: A vector of tags for the Dirichlet boundary (default is ["∂Γ"])
"""
function eigenmodes_2d(bc_case,D::SymFourthOrderTensorValue{2,T},ρ_d,N::Int,model,vox_trian,cache;
    order=2,
    μ=order*(order+1),
    dirichlet_tags=["∂Γ"]
  ) where T<:Number
  @assert order ≥ 2 "Order must be ≥ 2"
  @assert N > 0 "Number of eigenmodes must be positive"

  c = bc_case == :clamped ? 1 : 0;
  R,dΓ,dΛ,n_Λ,d∂Γ,n_∂Γ,hₑ = cache.trian_and_spaces
  midpoint_cache = cache.midpoint_cache
  A,M = cache.operators

  a_op, l_op = get_bilinear_forms(D,ρ_d,c,μ,hₑ,dΓ,dΛ,d∂Γ,n_Λ,n_∂Γ)

  assemble_matrix!(a_op,A,R,R)
  assemble_matrix!(l_op,M,R,R)

  # Solve shifted generalised eigenvalue problem
  sigma = -1
  AsigM = A - sigma * M

  F = factorisation(AsigM,T)
  function T_action!(y, x)
      y .= F \ (M * x)
      return y
  end

  T_map = LinearMap{eltype(A)}(T_action!, size(A,1), ismutating=true)
  vals_T, w_vec = eigs(T_map, nev=N, which=:LM)
  ωi = @. sqrt(1/vals_T + sigma)
  _w = eltype(R.vector_type)<:Real ? real(w_vec) : w_vec
  _uhi = [FEFunction(R, (_w[:, i])) for i in 1:N]
  uhi = [FEFunction(R, (_w[:, i]/sqrt(sum(∫(_uhi[i]'*_uhi[i])dΓ)))) for i in 1:N]

  # Evaluate modes and second-order derivatives on voxelised mesh
  Ui, ∂xxUi, ∂xyUi, ∂yyUi = evaluate_at_midpoints(uhi,midpoint_cache)

  return (;ωi, Ui, ∂xxUi, ∂xyUi, ∂yyUi)
end

function eigenmodes_2d(bc_case,D::SymFourthOrderTensorValue{2,T},ρ_d,N::Int,model,vox_trian;
    order=2,
    μ=order*(order+1),
    dirichlet_tags=["∂Γ"]
  ) where T<:Number
  cache = precompute_eigenmode_cache(bc_case,D,ρ_d,model,vox_trian;order,μ,dirichlet_tags)
  return eigenmodes_2d(bc_case,D,ρ_d,N,model,vox_trian,cache;order,μ,dirichlet_tags)
end

function precompute_eigenmode_cache(bc_case,D,ρ_d,model,vox_trian;order=2,μ=order*(order+1),dirichlet_tags=["∂Γ"])
  R,dΓ,dΛ,n_Λ,d∂Γ,n_∂Γ,hₑ = build_trian_and_spaces(bc_case,D,model,order,dirichlet_tags)

  c = bc_case == :clamped ? 1 : 0;
  a_op, l_op = get_bilinear_forms(D,ρ_d,c,μ,hₑ,dΓ,dΛ,d∂Γ,n_Λ,n_∂Γ)
  assem = SparseMatrixAssembler(R,R)
  A = allocate_matrix(a_op,assem,R,R)
  M = allocate_matrix(l_op,assem,R,R)

  Γ = get_triangulation(R)
  pts = map(mean, get_cell_coordinates(vox_trian))
  loc = precompute_point_location(Γ, pts)

  return (;trian_and_spaces=(R,dΓ,dΛ,n_Λ,d∂Γ,n_∂Γ,hₑ), operators=(A,M), midpoint_cache=(loc,pts))
end

function get_bilinear_forms(D,ρ_d,c,μ,hₑ,dΓ,dΛ,d∂Γ,n_Λ,n_∂Γ)
  a_op(w,r) = ∫(((D ⊙ ∇∇(w)) ⊙ ∇∇(r)))dΓ -
      c*∫((((D ⊙ ∇∇(w)) ⋅ n_∂Γ) ⋅ ∇(r)))d∂Γ -
      c*∫((((D ⊙ ∇∇(r)) ⋅ n_∂Γ) ⋅ ∇(w)))d∂Γ +
      c*∫(maximum(abs,D.data)*μ/hₑ*(∇(w) ⊗ n_∂Γ) ⊙ (∇(r) ⊗ n_∂Γ))d∂Γ -
      ∫(mean(D ⊙ ∇∇(w)) ⊙ jump(∇(r) ⊗ n_Λ))dΛ -
      ∫(mean(D ⊙ ∇∇(r)) ⊙ jump(∇(w) ⊗ n_Λ))dΛ +
      ∫(maximum(abs,D.data)*μ/hₑ*jump(∇(w) ⊗ n_Λ) ⊙ jump(∇(r) ⊗ n_Λ))dΛ
  l_op(w,r) = ∫(ρ_d*w*r)dΓ

  return a_op, l_op
end

function build_trian_and_spaces(bc_case,D,model,order,dirichlet_tags)
  Γ = Triangulation(model)
  hₑ = mean(CellField(lazy_map(vol->(vol)^(1/2),get_cell_measure(Γ)),Γ))
  Λ = Skeleton(model)
  dΓ = Measure(Γ,3*order)
  dΛ = Measure(Λ,3*order)
  n_Λ = get_normal_vector(Λ)

  ∂Γ = BoundaryTriangulation(Γ,tags=dirichlet_tags)
  d∂Γ = Measure(∂Γ,3*order)
  n_∂Γ = get_normal_vector(∂Γ)

  vector_type = all(isreal,D) ? Vector{Float64} : Vector{ComplexF64}
  if bc_case ∈ (:clamped,:simply_supported)
    R = TestFESpace(Γ,ReferenceFE(lagrangian,Float64,order),
      dirichlet_tags=dirichlet_tags,vector_type=vector_type)
  else
    R = TestFESpace(Γ,ReferenceFE(lagrangian,Float64,order),
      vector_type=vector_type)
  end

  return R,dΓ,dΛ,n_Λ,d∂Γ,n_∂Γ,hₑ
end

factorisation(A,::Type{<:Complex}) = lu(A)
factorisation(A,::Type{<:Real}) = cholesky(Hermitian(A))

# Evaluation on voxelised mesh
function precompute_point_location(trian::Triangulation, xs::AbstractVector{<:Gridap.Point})
  searchmethod = KDTreeSearch()
  loc_cache = Gridap.CellData._point_to_cell_cache(searchmethod, trian)
  point_to_cell = map(x -> Gridap.CellData._point_to_cell!(loc_cache, x), xs)
  cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell, num_cells(trian))
  cell_to_xs = lazy_map(Broadcasting(Reindex(xs)), cell_to_points)
  return point_to_cell, cell_to_xs, point_to_lpoint
end

function evaluate_at_location(f, loc)
  point_to_cell, cell_to_xs, point_to_lpoint = loc
  cell_to_f = get_array(f)
  cell_to_fxs = lazy_map(Gridap.evaluate, cell_to_f, cell_to_xs)
  point_to_fxs = lazy_map(Reindex(cell_to_fxs), point_to_cell)
  return lazy_map(getindex, point_to_fxs, point_to_lpoint)
end

function evaluate_at_midpoints(uhi,midpoint_cache)
  loc, pts = midpoint_cache
  num_modes = length(uhi)
  num_pts = length(pts)
  T = eltype(get_free_dof_values(first(uhi)))

  uhi_vals = Matrix{T}(undef, num_pts, num_modes)
  ∂xxuhi_vals = Matrix{T}(undef, num_pts, num_modes)
  ∂xyuhi_vals = Matrix{T}(undef, num_pts, num_modes)
  ∂yyuhi_vals = Matrix{T}(undef, num_pts, num_modes)

  u_lazy = evaluate_at_location(uhi[1], loc)
  ∇∇u_lazy = evaluate_at_location(∇∇(uhi[1]), loc)
  u_cache = array_cache(u_lazy)
  ∇∇u_cache = array_cache(∇∇u_lazy)

  for i in 1:num_modes
    u = evaluate_at_location(uhi[i], loc)
    for j in 1:num_pts
      uj = getindex!(u_cache, u, j)
      uhi_vals[j,i] = uj
    end
    ∇∇u = evaluate_at_location(∇∇(uhi[i]), loc)
    for j in 1:num_pts
      ∇∇uj = getindex!(∇∇u_cache, ∇∇u, j)
      ∂xxuhi_vals[j,i] = getindex(∇∇uj, 1, 1)
      ∂xyuhi_vals[j,i] = getindex(∇∇uj, 1, 2)
      ∂yyuhi_vals[j,i] = getindex(∇∇uj, 2, 2)
    end
  end

  return uhi_vals, ∂xxuhi_vals, ∂xyuhi_vals, ∂yyuhi_vals
end