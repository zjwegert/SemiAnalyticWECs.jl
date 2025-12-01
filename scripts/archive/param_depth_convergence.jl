using SemiAnalyticWECs
using CairoMakie, DataFrames, ProgressMeter, JLD2

# Constants
H = 10;
L = 10;
nu0 = 0.49;
E0 = 3.2e+6;
d0 = 0.01;
rho0 = 1250;
Ts = 6;

# Uncomment to re-run
R_s = Vector{ComplexF64}[]; T_s = Vector{ComplexF64}[];
P_FFs = Vector{Float64}[]; P_NFs = Vector{Float64}[];
mat_names = String[]; bc_names = String[]; problem_names = String[];
h_values = Float64[]; δxs = Float64[]; modess = Int[];

for δx in (0.1,0.01,0.001)
  for modes in (20,50,100)
    for h in (1e-4,1e-3,1e-2,1e-1,1,2)
      for mat in (PZT5H_material_coefficents,)
        for bc in (:simply_supported,)
          for solver in (solve_submerged_plate_2d,)
            ## Mat and Bc
            bc_name = string(bc)
            mat_name = first(split(string(mat),"_"))
            prob_name = split(string(solver),"_")[2]
            println("----- Solving $(prob_name) problem, material = $(mat_name), bc = $(bc_name), h=$h -----");
            C11,e21,K22,rhop = mat();
            dp = 1.1e-4;
            G = K22/(2*dp);
            C = K22/(2*dp);
            B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
            η = 1/2*(d0+dp)*e21;
            Ib = rho0*d0 + 2*rhop*dp;

            # Solve
            ω = 2π./Ts;
            D = @. B + ω*η^2/(im*G+ω*C);
            data = @showprogress map(1:length(ω)) do i
              solver(ω[i],D[i],Ib,η,G,C,H,h,L,modes,ceil(Int,L/(δx));bc_case=bc)
            end

            push!(R_s, getfield.(data,:R)); push!(T_s, getfield.(data,:T));
            push!(P_FFs, getfield.(data,:P_farfield)); push!(P_NFs, getfield.(data,:P_nearfield));
            push!(mat_names, mat_name); push!(bc_names, bc_name); push!(problem_names, prob_name);
            push!(h_values, h); push!(δxs, δx); push!(modess, modes);
          end
        end
      end
    end
  end
end

data_fine = DataFrame("Problem"=>problem_names,"Material"=>mat_names,"BC"=>bc_names,"h"=>h_values, "δx"=>δxs, "modes"=>modess, "R"=>R_s,"T"=>T_s,"P_farfield"=>P_FFs,"P_nearfield"=>P_NFs);

data_fine[data_fine.modes .== 100 .&& data_fine.h .== 1e-3,:].P_farfield

data_fine.P_nearfield[1]
# jldsave("$(@__DIR__)/data/param_depth.jld2";data)

##############################
### Plotting
##############################
data = load("$(@__DIR__)/data/param_depth.jld2")["data"]
data[4,:].P_nearfield


data_surface = load("$(@__DIR__)/data/pzt_vs_pvdf.jld2")["data"]
data_surface[6,:].P_nearfield



# fig_h_less_1 = with_theme(theme_latexfonts(),fontsize=24,linewidth=3) do
#   fig = Figure()
#   ax = Axis(fig[2,1],aspect=3,yscale=log10,xlabel="Period (s)",ylabel=L"P~\mathrm{(Wm^{-1})}")
#   for _data in eachrow(data[data.h .< 1,:])
#     lines!(ax,Ts,_data.P_nearfield,
#         label=L"h=%$(_data.h)")
#   end
#   L = Legend(fig[1,1],ax,orientation=:horizontal,tellwidth=true)
#   L.nbanks = 1
#   resize_to_layout!(fig)
#   fig
# end;

# fig_h_geq_1 = with_theme(theme_latexfonts(),fontsize=24,linewidth=3) do
#   fig = Figure()
#   ax = Axis(fig[2,1],aspect=3,yscale=log10,xlabel="Period (s)",ylabel=L"P~\mathrm{(Wm^{-1})}")
#   for _data in eachrow(data[data.h .>= 1,:])
#     lines!(ax,Ts,_data.P_nearfield,
#         label=L"h=%$(_data.h)")
#   end
#   L = Legend(fig[1,1],ax,orientation=:horizontal,tellwidth=true)
#   L.nbanks = 1
#   resize_to_layout!(fig)
#   fig
# end

# save("$(@__DIR__)/figures/depth_less_1.png",fig_h_less_1;dpi=300)
# save("$(@__DIR__)/figures/depth_geq_1.png",fig_h_geq_1;dpi=300)