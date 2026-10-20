using CairoMakie
using SemiAnalyticWECs
using Gridap, Gridap.TensorValues
using DelimitedFiles

isotropic_mat = false
ω = 10;
bc_case = :clamped

output_dir = isotropic_mat ? "results/annulus_isotropic_mixed_omega=$ω/" : "results/annulus_anisotropic_mixed_omega=$ω/"
mkpath(output_dir)

# Wave parameters
ρ_w = 1025;
g = 9.81;
A = 1.0;
χ = 0.0;
H = 2;

# Plate parameters
dirichlet_tags = ["Γ_d"]
h = 0.2;
d = 1;
ρ = 1;
num_modes = 30;

# Construct rigidity tensor
B_scale = 5
if isotropic_mat
    B11 = B_scale*1;
    B22 = B_scale*1;
    B12 = B_scale*0.3;
    B16 = 0;
    B26 = 0;
    B66 = B_scale*0.35;
else
    B11 = B_scale*1;
    B22 = B_scale*1;
    B12 = B_scale*0.9082;
    B16 = B_scale*0.6724;
    B26 = B_scale*0.6724;
    B66 = B_scale*0.9341;
end
B = [B11 B12 B16;
     B12 B22 B26;
     B16 B26 B66]
D = B |> TensorValue |> from_voigt

# Setup meshes
δx_min = 1/200
model, bbox = annulus_portion_removed_2d(;R1=0.5,R2=0.25,h=δx_min,angle1=195 * pi / 180, angle2= 165 * pi / 180)
Γ_vox, δx = generate_voxelised_trian(model,bbox,δx_min)

# Compute eigenmodes
modes = eigenmodes_2d(bc_case,D,ρ*d,num_modes,model,Γ_vox;dirichlet_tags);

# Solve problem
data = solve_submerged_plate_3d(ω, A, χ, H, h, ρ, d, 0, 0, 0, 1, 1, Γ_vox, δx, modes; g, ρ_w)
W, f, scs, ecs = data

# Plots
theta = range(χ-π,χ+π,1001)
fig = Figure(size = (400, 400));
ax = PolarAxis(fig[1, 1],
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