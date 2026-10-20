"""
  disk_2d(; R1=0.25, R2=R1, center=(0.5, 0.5), h=0.02)

Generate a disk-shaped 2D mesh with semi-axes R1 and R2, centered at `center`, and with mesh size `h`.
"""
function disk_2d(;
  R1=0.25,
  R2=R1,
  center=(0.5, 0.5),
  h=0.02
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("disk_model")

  # Build
  gmsh.model.occ.addEllipse(center[1], center[2], 0.0, R1, R2, 1)
  gmsh.model.occ.addCurveLoop([1], 1)
  gmsh.model.occ.addPlaneSurface([1], 1)

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Γ")
  gmsh.model.addPhysicalGroup(1, [1], 2)
  gmsh.model.setPhysicalName(1, 2, "Γ_d")

  # Generate mesh
  gmsh.model.mesh.generate(2)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()

  # bounding box
  bbox_min = center .- (R1, R2)
  bbox_max = center .+ (R1, R2)
  return model, (bbox_min[1], bbox_max[1], bbox_min[2], bbox_max[2])
end

"""
  annulus_2d(; R1=0.25, R2=R1, Rin=0.5 * R1, center=(0.5, 0.5), h=0.02)

Generate an annular 2D mesh with outer semi-axes R1 and R2, inner radius Rin, centered at `center`, and with mesh size `h`.
"""
function annulus_2d(;
  R1=0.25,
  R2=R1,
  Rin=0.5 * R1,
  center=(0.5, 0.5),
  h=0.02
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("annulus_model")

  # Build: outer ellipse minus a concentric inner disk
  gmsh.model.occ.addEllipse(center[1], center[2], 0.0, R1, R2, 1)
  gmsh.model.occ.addCurveLoop([1], 1)
  gmsh.model.occ.addEllipse(center[1], center[2], 0.0, Rin, Rin * R2 / R1, 2)
  gmsh.model.occ.addCurveLoop([2], 2)
  gmsh.model.occ.addPlaneSurface([1, 2], 1)

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Γ")
  gmsh.model.addPhysicalGroup(1, [1,2], 2)
  gmsh.model.setPhysicalName(1, 2, "Γ_d")

  # Generate mesh
  gmsh.model.mesh.generate(2)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()

  # bounding box
  bbox_min = center .- (R1, R2)
  bbox_max = center .+ (R1, R2)
  return model, (bbox_min[1], bbox_max[1], bbox_min[2], bbox_max[2])
end

"""
  annulus_portion_removed_2d(; R1=0.25, R2=R1, Rin=0.5 * R1, center=(0.5, 0.5), h=0.02,
    angle1=-165 * pi / 180, angle2= 165 * pi / 180)

Generate an annular 2D mesh with outer semi-axes R1 and R2, inner radius Rin, centered at `center`, and with mesh size `h`.
This version removes a portion of the annulus between angles `angle1` and `angle2` and labels the boundary of the resulting boundary
as Γ_d.
"""
function annulus_portion_removed_2d(;
  R1=0.25,
  R2=R1,
  Rin=0.5 * R1,
  center=(0.5, 0.5),
  h=0.02,
  angle1=-165 * pi / 180,
  angle2= 165 * pi / 180
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("annulus_model")

  gmsh.model.occ.addEllipse(center[1], center[2], 0.0, R1, R2, 1, angle1, angle2)

  pt_start = 1
  pt_end   = 2
  gmsh.model.occ.addLine(pt_end, pt_start, 3)
  gmsh.model.occ.addCurveLoop([1, 3], 1)
  gmsh.model.occ.addEllipse(center[1], center[2], 0.0, Rin, Rin * R2 / R1, 2)
  gmsh.model.occ.addCurveLoop([2], 2)
  gmsh.model.occ.addPlaneSurface([1, 2], 1)

  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Γ")
  gmsh.model.addPhysicalGroup(1, [3], 2)
  gmsh.model.setPhysicalName(1, 2, "Γ_d")
  gmsh.model.addPhysicalGroup(1, [1, 2], 3)
  gmsh.model.setPhysicalName(1, 3, "∂Γ")


  gmsh.model.mesh.generate(2)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()

  # Bounding box
  bbox_min = center .- (R1, R2)
  bbox_max = center .+ (R1, R2)

  return model, (bbox_min[1], bbox_max[1], bbox_min[2], bbox_max[2])
end

"""
  rectangle_2d(; R1=0.25, R2=R1, center=(0.5, 0.5), h=0.02)

Generate a rectangular 2D mesh with width R1 and height R2, centered at `center`, and with mesh size `h`.
"""
function rectangle_2d(;
  R1=0.25,
  R2=R1,
  center=(0.5, 0.5),
  h=0.02
)
  gmsh.initialize()

  gmsh.option.setNumber("Mesh.SaveAll", 1)
  gmsh.model.add("rectangle_model")

  # Build
  gmsh.model.occ.addRectangle(center[1] - R1, center[2] - R2, 0.0, 2 * R1, 2 * R2, 1)
  # gmsh.model.occ.addCurveLoop([1], 1)
  # gmsh.model.occ.addPlaneSurface([1], 1)

  # Mesh sizes and tags
  gmsh.model.occ.synchronize()
  gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Γ")
  gmsh.model.addPhysicalGroup(1, [1], 2)
  gmsh.model.setPhysicalName(1, 2, "Γ_d_bottom")
  gmsh.model.addPhysicalGroup(1, [2], 3)
  gmsh.model.setPhysicalName(1, 3, "Γ_d_right")
  gmsh.model.addPhysicalGroup(1, [3], 4)
  gmsh.model.setPhysicalName(1, 4, "Γ_d_top")
  gmsh.model.addPhysicalGroup(1, [4], 5)
  gmsh.model.setPhysicalName(1, 5, "Γ_d_left")
  gmsh.model.addPhysicalGroup(0, [1,4], 7)
  gmsh.model.setPhysicalName(0, 7, "Γ_d_left")
  gmsh.model.addPhysicalGroup(0, [2,3], 8)
  gmsh.model.setPhysicalName(0, 8, "Γ_d_right")

  # gmsh.model.occ.synchronize()
  # gmsh.fltk.run()

  # Generate mesh
  gmsh.model.mesh.generate(2)
  model = GmshDiscreteModel(gmsh)
  gmsh.finalize()

  # bounding box
  bbox_min = center .- (R1, R2)
  bbox_max = center .+ (R1, R2)
  return model, (bbox_min[1], bbox_max[1], bbox_min[2], bbox_max[2])
end

function point_in_triangulation(trian, pts::Vector{<:Gridap.Point}; cache=_point_to_cell_cache(KDTreeSearch(), trian))
  for pt in pts
    if !point_in_triangulation(trian, pt; cache=cache)
      return false
    end
  end
  return true
end

function point_in_triangulation(trian, pt; cache=_point_to_cell_cache(KDTreeSearch(), trian))
  try
    Gridap.CellData._point_to_cell!(cache, pt)
    return true
  catch
    return false
  end
end

function points_in_triangulation(trian,pts)
  inout = zeros(Bool, size(pts));
  cache = Gridap.CellData._point_to_cell_cache(KDTreeSearch(), trian)

  for i ∈ Base.OneTo(length(pts))
    pt = pts[i]
    if point_in_triangulation(trian, pt; cache=cache)
      inout[i] = true
    end
  end;
  return inout
end

"""
  generate_voxelised_trian(model, bbox, δx)

Generate a voxelised triangulation of the given `model` within the bounding box `bbox`
with voxel size `δx`. The function returns a new triangulation representing
the voxelised model and the voxel size used.
"""
function generate_voxelised_trian(model, bbox, δx)
  bbox = bbox./δx
  _bbox = (floor(bbox[1]), ceil(bbox[2]), floor(bbox[3]), ceil(bbox[4])).*δx
  nx = round(Int,(_bbox[2]-_bbox[1])/δx)
  ny = round(Int,(_bbox[4]-_bbox[3])/δx)
  vox_model = CartesianDiscreteModel(_bbox,(nx,ny))
  cell_to_coords = get_cell_coordinates(vox_model)
  vox_δx = (cell_to_coords[1][4] - cell_to_coords[1][1])[1]
  vox_δy = (cell_to_coords[1][4] - cell_to_coords[1][1])[2]
  @assert vox_δx ≈ vox_δy "$vox_δx $vox_δy"
  midpoints = map(mean, cell_to_coords);
  inout = points_in_triangulation(Triangulation(model),midpoints)

  return Triangulation(vox_model,inout), vox_δx
end

# This constructs an interpolation of W on Γ by projecting W onto H¹ with a
# weighting parameter α. If α=0.0, this is equivalent to a simple L² projection.
# If α>0.0, the interpolation is smoother.
using ScatteredInterpolation
function interpolate_from_Γ_vox_to_Γ(W,Γ_vox,Γ;α=0.0)
  points = map(mean, get_cell_coordinates(Γ_vox))
  points_mat = hcat(getindex.(points,1), getindex.(points,2))

  intp = ScatteredInterpolation.interpolate(NearestNeighbor(),points_mat',W);
  function IW(x)
    first(ScatteredInterpolation.evaluate(intp, [x[1];x[2]]))
  end
  h = maximum(sqrt,get_cell_measure(Γ))
  V = TestFESpace(get_background_model(Γ), ReferenceFE(lagrangian,Float64,1), vector_type=Vector{ComplexF64})
  dΓ = Measure(Γ,2)
  op = AffineFEOperator((u,v)->∫((α*h)^2*∇(u)⋅∇(v) + u*v)dΓ, v->∫(IW*v)dΓ,V,V)
  return Gridap.solve(op)
end