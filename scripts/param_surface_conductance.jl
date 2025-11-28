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
# G_values = Float64[]; G_coeffs = Float64[];

# for coeff in (10 .^ range(-2,2,5))
#   for mat in (PZT5H_material_coefficents,)
#     for bc in (:simply_supported,)
#       for solver in (solve_submerged_plate_2d,)
#         ## Mat and Bc
#         bc_name = string(bc)
#         mat_name = first(split(string(mat),"_"))
#         prob_name = split(string(solver),"_")[2]
#         println("----- Solving $(prob_name) problem, material = $(mat_name), bc = $(bc_name), G_coeff=$coeff -----");
#         C11,e21,K22,rhop = mat();
#         dp = 1.1e-4;
#         G = coeff*K22/(2*dp);
#         C = K22/(2*dp);
#         B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
#         η = 1/2*(d0+dp)*e21;
#         Ib = rho0*d0 + 2*rhop*dp;

#         # Solve
#         ω = 2π./Ts;
#         D = @. B + ω*η^2/(im*G+ω*C);
#         data = @showprogress map(1:length(ω)) do i
#           solver(ω[i],D[i],Ib,η,G,C,H,h,L,100,ceil(Int,L/0.05);bc_case=bc)
#         end

#         push!(R_s, getfield.(data,:R)); push!(T_s, getfield.(data,:T));
#         push!(P_FFs, getfield.(data,:P_farfield)); push!(P_NFs, getfield.(data,:P_nearfield));
#         push!(mat_names, mat_name); push!(bc_names, bc_name); push!(problem_names, prob_name);
#         push!(G_coeffs, coeff); push!(G_values, G);
#       end
#     end
#   end
# end

# data = DataFrame("Problem"=>problem_names,"Material"=>mat_names,"BC"=>bc_names,"G_coeff"=>G_coeffs,"G_value"=>G_values,"R"=>R_s,"T"=>T_s,"P_farfield"=>P_FFs,"P_nearfield"=>P_NFs);
# jldsave("$(@__DIR__)/data/surface_conductance.jld2";data)

##############################
### Plotting
##############################
data = load("$(@__DIR__)/data/surface_conductance.jld2")["data"]

fig = with_theme(theme_latexfonts(),fontsize=24,linewidth=3) do
  fig = Figure()
  ax = Axis(fig[2,1],aspect=3,yscale=log10,xlabel="Period (s)",ylabel=L"P~\mathrm{(Wm^{-1})}")
  for _data in eachrow(data)
    lab = isone(_data.G_coeff) ? L"G=G_0" : L"G=10^{%$(Int(round(log10(_data.G_coeff);sigdigits=1)))}G_0"
    lines!(ax,Ts,_data.P_nearfield,
        label=lab)
  end
  L = Legend(fig[1,1],ax,orientation=:horizontal,tellwidth=true)
  L.nbanks = 1
  resize_to_layout!(fig)
  fig
end

save("$(@__DIR__)/figures/surface_conductance.png",fig;dpi=300)