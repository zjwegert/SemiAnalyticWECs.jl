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

for h in (0.0085,0.0065)#0.01,0.0075,0.005)
  for modes in (100,)
    for δx in (h/10,h/20)
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
jldsave("$(@__DIR__)/data/convergence.jld2";data_fine)

#  h = 0.01, δx=h/10: 1089.0775171154764
#  h = 0.01, δx=h/20: 1088.4415470517333
#  h = 0.0085, δx=h/10: 2278.4019197695816
#  h = 0.0085, δx=h/20: 2282.5638727718756
#  h = 0.0075, δx=h/10: 1408.566322861937
#  h = 0.0075, δx=h/20: 1404.8881220775806
#  h = 0.0065, δx=h/10: 804.1568142349087
#  h = 0.005, δx=h/10: 627.9119153540815

# | P_farfield - P_nearfield | ~ 1e-10

##############################
### Plotting
##############################
data_fine = load("$(@__DIR__)/data/convergence.jld2")["data_fine"]

f = with_theme(theme_latexfonts(),fontsize=24,linewidth=3) do
  fig = Figure()
  ax = Axis(fig[2,1],aspect=2,xlabel=L"h",ylabel=L"P~\mathrm{(Wm^{-1})}")
  ax.xreversed = true
  _d = data_fine[data_fine.δx .== data_fine.h/10,:]
  p = sortperm(_d.h)
  lines!(ax,_d.h[p],first.(_d.P_nearfield)[p],    label=L"\delta_x=h/10")
  _d = data_fine[data_fine.δx .== data_fine.h/20,:]
  p = sortperm(_d.h)
  lines!(ax,_d.h[p],first.(_d.P_nearfield)[p],label=L"\delta_x=h/20",linestyle=:dash)
  L = Legend(fig[1,1],ax,orientation=:horizontal)
  # L.nbanks = 1
  # resize_to_layout!(fig)
  fig
end;

save("$(@__DIR__)/figures/convergence.png",f;dpi=300)