using CairoMakie
using SemiAnalyticWECs
using Gridap, Gridap.TensorValues
using DelimitedFiles
using ScatteredInterpolation
using GridapGmsh, Gmsh

function custom_geometry_2d(;
    lc=0.05,
    petals1=5, petals2=6, petals3=4,
    center1=(0.0, 0.0), center2=(2.5, 1.5), center3=(2.5, -1.5),
    N = 100
)
    gmsh.initialize()

    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.model.add("three_dynamic_flowers")

    function add_flower!(cx, cy, R, amp, num_petals, base_tag, surf_tag)
        pts = Int[]

        # Generate points
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

        # Generate lines
        for i in 1:(N-1)
            l_tag = base_tag + i
            gmsh.model.occ.addLine(pts[i], pts[i+1], l_tag)
            push!(lines, l_tag)
        end
        l_tag_last = base_tag + N
        gmsh.model.occ.addLine(pts[N], pts[1], l_tag_last)
        push!(lines, l_tag_last)

        # Curve Loop & Surface
        gmsh.model.occ.addCurveLoop(lines, surf_tag)
        gmsh.model.occ.addPlaneSurface([surf_tag], surf_tag)

        return lines
    end

    lines1 = add_flower!(center1[1], center1[2], 0.6, 0.2, petals1, 1000, 1)
    max_r1 = 0.6 + 0.2

    lines2 = add_flower!(center2[1], center2[2], 0.6, 0.2, petals2, 2000, 2)
    max_r2 = 0.6 + 0.2

    lines3 = add_flower!(center3[1], center3[2], 0.6, 0.2, petals3, 3000, 3)
    max_r3 = 0.6 + 0.2

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "Flower_1")
    gmsh.model.addPhysicalGroup(2, [2], 2)
    gmsh.model.setPhysicalName(2, 2, "Flower_2")
    gmsh.model.addPhysicalGroup(2, [3], 3)
    gmsh.model.setPhysicalName(2, 3, "Flower_3")
    gmsh.model.addPhysicalGroup(1, lines1, 4)
    gmsh.model.setPhysicalName(1, 4, "Flower_1_boundary")
    gmsh.model.addPhysicalGroup(1, lines2, 5)
    gmsh.model.setPhysicalName(1, 5, "Flower_2_boundary")
    gmsh.model.addPhysicalGroup(1, lines3, 6)
    gmsh.model.setPhysicalName(1, 6, "Flower_3_boundary")

    gmsh.model.mesh.generate(2)
    model = GmshDiscreteModel(gmsh)
    gmsh.finalize()

    bbox_min_x = min(center1[1] - max_r1, center2[1] - max_r2, center3[1] - max_r3)
    bbox_max_x = max(center1[1] + max_r1, center2[1] + max_r2, center3[1] + max_r3)
    bbox_min_y = min(center1[2] - max_r1, center2[2] - max_r2, center3[2] - max_r3)
    bbox_max_y = max(center1[2] + max_r1, center2[2] + max_r2, center3[2] + max_r3)

    return model, (bbox_min_x, bbox_max_x, bbox_min_y, bbox_max_y)
end

bc_case = :free;
output_dir = "results/multiple_plate_scattering_$(string(bc_case))/"
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
dirichlet_tags = ["Flower_1_boundary","Flower_2_boundary","Flower_3_boundary"];

Γ_vox, δx = generate_voxelised_trian(model,bbox,δx_vox)

writevtk(model,"$output_dir/model")
writevtk(Γ_vox,"$output_dir/voxelised_model")

# Compute eigenmodes
modes = eigenmodes_2d(bc_case,D,ρ*d,num_modes,model,Γ_vox;dirichlet_tags);

# Solve problem
data = solve_submerged_plate_3d(ω, A, χ, H, h, ρ, d, 0.0, 0.0, 0.0, 1, 1, Γ_vox, δx, modes; g, ρ_w)
W, f, scs, ecs = data

# Plots
theta = range(χ-π,χ+π,2001)
fig = Figure(size = (400, 400));
ax = PolarAxis(fig[1, 1],
  rtickformat = xs -> string.(round.(xs; sigdigits=3)),
  thetaticks=(0:π/3:5π/3,[L"0",L"π/3",L"2π/3",L"π",L"4π/3",L"5π/3"]))
lineobject = lines!(ax, theta, abs2 ∘ f, color = :red);
writedlm("$output_dir/f_of_θ.txt", [theta abs2.(f.(theta))],',')
save("$output_dir/f_of_θ.png",fig)

writevtk(Γ_vox,"$output_dir/W", cellfields=["re(W)"=>real.(W),"imag(W)"=>imag.(W)])
writevtk(Triangulation(model),"$output_dir/model")
writedlm("$output_dir/scs_ecs.txt", [scs, ecs],',')

IWh = SemiAnalyticWECs.interpolate_from_Γ_vox_to_Γ(W,Γ_vox,Triangulation(model);α=0.0)
writevtk(Triangulation(model),"$output_dir/W_interpolated",
  cellfields=["re(W)"=>real(IWh),"imag(W)"=>imag(IWh)])