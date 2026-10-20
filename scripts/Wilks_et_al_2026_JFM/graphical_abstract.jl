using SemiAnalyticWECs
using Gridap, Gridap.TensorValues
using DelimitedFiles
using ScatteredInterpolation
using GridapGmsh, Gmsh
using ProgressMeter

function custom_geometry_2d(;
    lc=0.05,
    petals1=5, petals2=6, petals3=4,
    center1=(0.0, 0.0), center2=(2.5, 1.5), center3=(2.5, -1.5),
    N = 100 # Number of boundary points per flower
)
    gmsh.initialize()

    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.model.add("three_dynamic_flowers")

    # =========================================================================
    # Helper Function: Generates a flower geometry and returns its boundary lines
    # =========================================================================
    function add_flower!(cx, cy, R, amp, num_petals, base_tag, surf_tag)
        pts = Int[]

        # 1. Generate points
        for i in 1:N
            t = (i - 1) * 2 * π / N
            r = R + amp * sin(num_petals * t)

            x = cx + r * cos(t)
            y = cy + r * sin(t)

            pt_tag = base_tag + i
            gmsh.model.occ.addPoint(x, y, 0.0, lc, pt_tag)
            push!(pts, pt_tag)
        end

        lines = Int[]

        # 2. Generate lines
        for i in 1:(N-1)
            l_tag = base_tag + i
            gmsh.model.occ.addLine(pts[i], pts[i+1], l_tag)
            push!(lines, l_tag)
        end
        # Close the loop
        l_tag_last = base_tag + N
        gmsh.model.occ.addLine(pts[N], pts[1], l_tag_last)
        push!(lines, l_tag_last)

        # 3. Curve Loop & Surface
        gmsh.model.occ.addCurveLoop(lines, surf_tag)
        gmsh.model.occ.addPlaneSurface([surf_tag], surf_tag)

        return lines # Return the line tags to create 1D physical boundaries
    end

    # =========================================================================
    # 1. FLOWER 1 (Radius 0.6, Amp 0.2 -> Max extent 0.8)
    # =========================================================================
    lines1 = add_flower!(center1[1], center1[2], 0.6, 0.2, petals1, 1000, 1)
    max_r1 = 0.6 + 0.2

    # =========================================================================
    # 2. FLOWER 2 (Radius 0.6, Amp 0.2 -> Max extent 0.8)
    # =========================================================================
    lines2 = add_flower!(center2[1], center2[2], 0.6, 0.2, petals2, 2000, 2)
    max_r2 = 0.6 + 0.2

    # =========================================================================
    # 3. FLOWER 3 (Radius 0.6, Amp 0.2 -> Max extent 0.8)
    # =========================================================================
    lines3 = add_flower!(center3[1], center3[2], 0.6, 0.2, petals3, 3000, 3)
    max_r3 = 0.6 + 0.2

    # =========================================================================
    # Synchronize, Mesh Sizes, & Physical Groups
    # =========================================================================
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

    # Assign 2D physical groups for the surfaces (Domains)
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "Flower_1")

    gmsh.model.addPhysicalGroup(2, [2], 2)
    gmsh.model.setPhysicalName(2, 2, "Flower_2")

    gmsh.model.addPhysicalGroup(2, [3], 3)
    gmsh.model.setPhysicalName(2, 3, "Flower_3")

    # Assign 1D physical groups for the curves (Boundaries)
    gmsh.model.addPhysicalGroup(1, lines1, 4)
    gmsh.model.setPhysicalName(1, 4, "Flower_1_boundary")

    gmsh.model.addPhysicalGroup(1, lines2, 5)
    gmsh.model.setPhysicalName(1, 5, "Flower_2_boundary")

    gmsh.model.addPhysicalGroup(1, lines3, 6)
    gmsh.model.setPhysicalName(1, 6, "Flower_3_boundary")

    # Generate Mesh
    gmsh.model.mesh.generate(2)

    # Import into Gridap
    model = GmshDiscreteModel(gmsh)
    gmsh.finalize()

    # Calculate dynamic bounding box based on the provided centers and max radii
    bbox_min_x = min(center1[1] - max_r1, center2[1] - max_r2, center3[1] - max_r3)
    bbox_max_x = max(center1[1] + max_r1, center2[1] + max_r2, center3[1] + max_r3)
    bbox_min_y = min(center1[2] - max_r1, center2[2] - max_r2, center3[2] - max_r3)
    bbox_max_y = max(center1[2] + max_r1, center2[2] + max_r2, center3[2] + max_r3)

    return model, (bbox_min_x, bbox_max_x, bbox_min_y, bbox_max_y)
end

output_dir = "results/graphical_abstract_clamped/"
mkpath(output_dir)

# Wave parameters
ω = 10;
ρ_w = 1025;
g = 9.81;
A = 1.0;
χ = 0.0;
H = 2;

# Plate parameters
h = 0.2;
d = 1;
ρ = 1;
num_modes = 90;

# Construct rigidity tensor
B_scale = 5
B11 = B_scale*1;
B22 = B_scale*1;
B12 = B_scale*0.3;
B16 = 0;
B26 = 0;
B66 = B_scale*0.35;
B = [B11 B12 B16;
     B12 B22 B26;
     B16 B26 B66]
D = B |> TensorValue |> from_voigt

# Setup meshes
δx_min = 1/60
δx_vox = 1/120
model, bbox = custom_geometry_2d(;lc=δx_min,petals1=7, petals2=5, petals3=4,center1=(0.0, 0.0), center2=(1.5, 0.75), center3=(1.5, -0.75),N = 200);
bc_case = :clamped; dirichlet_tags = ["Flower_1_boundary","Flower_2_boundary","Flower_3_boundary"];
Γ_vox, δx = generate_voxelised_trian(model,bbox,δx_vox)

writevtk(model,"$output_dir/model")
writevtk(Γ_vox,"$output_dir/voxelised_model")
writevtk(Boundary(model),"$output_dir/model_∂Γ")

# Compute eigenmodes
mode_cache = SemiAnalyticWECs.precompute_eigenmode_cache(bc_case,D,ρ*d,model,Γ_vox;dirichlet_tags);
modes = eigenmodes_2d(bc_case,D,ρ*d,num_modes,model,Γ_vox,mode_cache;dirichlet_tags);

# Solve problem
solve_cache = SemiAnalyticWECs.precompute_cache(get_background_model(Γ_vox),δx)
data = solve_submerged_plate_3d(ω, A, χ, H, h, ρ, d, 0.0, 0.0, 0.0, 1, 1, Γ_vox, δx, modes, solve_cache; g, ρ_w)
W, f, scs, ecs, _, _, σ = data

# Plots
writevtk(Γ_vox,"$output_dir/W", cellfields=["re(W)"=>real.(W),"imag(W)"=>imag.(W)])
writevtk(Triangulation(model),"$output_dir/model")
writedlm("$output_dir/scs_ecs.txt", [scs, ecs],',')

IWh = SemiAnalyticWECs.interpolate_from_Γ_vox_to_Γ(W,Γ_vox,Triangulation(model);α=2.0)
writevtk(Triangulation(model),"$output_dir/W_interpolated",
  cellfields=["re(W)"=>real(IWh),"imag(W)"=>imag(IWh)])

# Free surface
inout_idx = Γ_vox.tface_to_mface
Xd_in, Yd_in = solve_cache[2][inout_idx], solve_cache[3][inout_idx]
model_fs = CartesianDiscreteModel(bbox.*1.25,(100,100))
Γ_fs = Triangulation(model_fs)
midpoints_fs = vec(map(mean,get_cell_coordinates(Γ_fs)))
X = getindex.(midpoints_fs,1)
Y = getindex.(midpoints_fs,2)
R = norm.(midpoints_fs)
α = ω^2/g;
k = first(dispersion_free_surface(α,0,H))/im;

function compute_free_surface_elevation(X,Y,Xd_in,Yd_in,σ,A,χ,k,ω,g,h,α,H,δx)
  cχ, sχ = cos(χ), sin(χ)
  pre = -im*ω/g*δx^2
  η = zeros(ComplexF64,length(X))
  @showprogress Threads.@threads for p ∈ eachindex(X)
    xp, yp = X[p], Y[p]
    acc = zero(ComplexF64)
    @inbounds for ii ∈ eachindex(Xd_in)
      rp = sqrt((Xd_in[ii]-xp)^2+(Yd_in[ii]-yp)^2)
      Gval = -1/4π*SemiAnalyticWECs.∂ζ_greens_submerged_3d(max(rp,eps()),0,-h,α,H)
      acc += Gval*σ[ii]
    end
    η[p] = A*exp(im*k*(xp*cχ+yp*sχ)) + pre*acc
  end
  return η
end

η = compute_free_surface_elevation(X,Y,Xd_in,Yd_in,σ,A,χ,k,ω,g,h,α,H,δx)

function interpolate_same_mesh(W,Γ;α=0.0)
  h = maximum(sqrt,get_cell_measure(Γ))
  V = TestFESpace(get_background_model(Γ), ReferenceFE(lagrangian,Float64,1), vector_type=Vector{ComplexF64})
  dΓ = Measure(Γ,2)
  op = AffineFEOperator((u,v)->∫((α*h)^2*∇(u)⋅∇(v) + u*v)dΓ, v->∫(W*v)dΓ,V,V)
  return solve(op)
end

ηh_0 = interpolate_same_mesh(η,Γ_fs;α=0.0)
ηh_2 = interpolate_same_mesh(η,Γ_fs;α=2.0)

writevtk(Γ_fs,"$output_dir/free_surface",
  cellfields=["re(eta)"=>real.(η),"imag(eta)"=>imag.(η),
    "re(eta0)"=>real(ηh_0),"imag(eta0)"=>imag(ηh_0),
    "re(eta2)"=>real(ηh_2),"imag(eta2)"=>imag(ηh_2)])