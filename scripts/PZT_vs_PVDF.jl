using SemiAnalyticWECs
using CairoMakie, DataFrames, ProgressMeter, JLD2

# Constants
H = 10;
h = 2;
L = 10;
nu0 = 0.49;
E0 = 3.2e+6;
d0 = 0.01;
rho0 = 1250;
Ts = collect(range(3,9,2000));

## Uncomment to re-run
# R_s = Vector{ComplexF64}[]; T_s = Vector{ComplexF64}[];
# P_FFs = Vector{Float64}[]; P_NFs = Vector{Float64}[];
# mat_names = String[]; bc_names = String[]; problem_names = String[];

# for mat in (PVDF_TechMan_material_coefficents,PZT5H_material_coefficents)
#   for bc in (:simply_supported,:clamped)
#     for solver in (solve_submerged_plate_2d,solve_surface_plate_2d)
#       ## Mat and Bc
#       bc_name = string(bc)
#       mat_name = first(split(string(mat),"_"))
#       prob_name = split(string(solver),"_")[2]
#       println("----- Solving $(prob_name) problem, material = $(mat_name), bc = $(bc_name) -----");
#       C11,e21,K22,rhop = mat();
#       dp = 1.1e-4;
#       G = K22/(2*dp);
#       C = K22/(2*dp);
#       B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
#       η = 1/2*(d0+dp)*e21;
#       Ib = rho0*d0 + 2*rhop*dp;

#       # Solve
#       ω = 2π./Ts;
#       D = @. B + ω*η^2/(im*G+ω*C);
#       data = @showprogress map(1:length(ω)) do i
#         solver(ω[i],D[i],Ib,η,G,C,H,h,L,100,ceil(Int,L/0.05);bc_case=bc)
#       end

#       push!(R_s, getfield.(data,:R)); push!(T_s, getfield.(data,:T));
#       push!(P_FFs, getfield.(data,:P_farfield)); push!(P_NFs, getfield.(data,:P_nearfield));
#       push!(mat_names, mat_name); push!(bc_names, bc_name); push!(problem_names, prob_name);
#     end
#   end
# end

# data = DataFrame("Problem"=>problem_names,"Material"=>mat_names,"BC"=>bc_names,"R"=>R_s,"T"=>T_s,"P_farfield"=>P_FFs,"P_nearfield"=>P_NFs);
# jldsave("$(@__DIR__)/data/pzt_vs_pvdf.jld2";data)

##############################
### Plotting
##############################
data = load("$(@__DIR__)/data/pzt_vs_pvdf.jld2")["data"]

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=4) do
  fig = Figure(size = (1400, 400),figure_padding = (1,95,1,1))
  # Energy (surface)
  ax = Axis(fig[2,1],aspect=2,yscale=log10,xlabel="Period (s)",ylabel=L"1 -|R|^2-|T|^2",xticks=3:9)
  _data = data[data.BC .== "simply_supported" .&& data.Problem .== "surface" .&& data.Material .== "PVDF",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PVDF (Simply supported)")
  _data = data[data.BC .== "clamped" .&& data.Problem .== "surface" .&& data.Material .== "PVDF",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PVDF (Clamped)")
  _data = data[data.BC .== "simply_supported" .&& data.Problem .== "surface" .&& data.Material .== "PZT5H",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PZT5H (Simply supported)",linestyle=:dash)
  _data = data[data.BC .== "clamped" .&& data.Problem .== "surface" .&& data.Material .== "PZT5H",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PZT5H (Clamped)",linestyle=:dash)
  ylims!(ax,0.8e-4,2.5e-1)
  # Energy (submerged)
  ax = Axis(fig[2,2],aspect=2,yscale=log10,xlabel="Period (s)",ylabel=L"1 -|R|^2-|T|^2",xticks=3:9)
  _data = data[data.BC .== "simply_supported" .&& data.Problem .== "submerged" .&& data.Material .== "PVDF",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PVDF (Simply supported)")
  _data = data[data.BC .== "clamped" .&& data.Problem .== "submerged" .&& data.Material .== "PVDF",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PVDF (Clamped)")
  _data = data[data.BC .== "simply_supported" .&& data.Problem .== "submerged" .&& data.Material .== "PZT5H",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PZT5H (Simply supported)",linestyle=:dash)
  _data = data[data.BC .== "clamped" .&& data.Problem .== "submerged" .&& data.Material .== "PZT5H",:]
  lines!(ax,Ts,(1 .-abs.(first(_data.R)).^2 .-abs.(first(_data.T)).^2),label="PZT5H (Clamped)",linestyle=:dash)
  L = Legend(fig[1,1:2],ax,orientation=:horizontal)#,tellwidth=true)
  ylims!(ax,0.8e-4,2.5e-1)
  # L.nbanks = 2
  # resize_to_layout!(fig)
  fig
end

save("$(@__DIR__)/figures/PZT_vs_PVDF.png",fig;dpi=300)