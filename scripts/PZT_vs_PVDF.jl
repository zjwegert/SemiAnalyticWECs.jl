using SemiAnalyticWECs
using CairoMakie

# Constants
H = 10;
h = 2;
L = 10;
nu0 = 0.49;
E0 = 3.2e+6;
d0 = 0.01;
rho0 = 1250;

## Surface free-free
# PVDF
C11,e21,K22,rhop = PVDF_TechMan_material_coefficents();
dp = 1.1e-4;
G = K22/(2*dp);
C = K22/(2*dp);
B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
η = 1/2*(d0+dp)*e21;
Ib = rho0*d0 + 2*rhop*dp;

# Submerged
Ts = collect(range(3,9,1000));
ω = 2π./Ts;
Tsx = 5;
D = @. B + ω*η^2/(im*G+ω*C);

data = map(1:length(ω)) do i
  println(((i-1)/length(ω)*100), " percent complete");
  solve_submerged_plate_2d(ω[i],D[i],Ib,η,G,C,H,h,L,100,ceil(Int,L/0.05);bc_case=:simply_supported)
end

with_theme(theme_latexfonts(),fontsize=24,linewidth=2.5) do
  fig = Figure()
  ax = Axis(fig[1,1])
  lines!(ax,Ts,getfield.(data,:P_farfield),label="Far-field Power")
  lines!(ax,Ts,getfield.(data,:P_nearfield),label="Near-field Power",linestyle=:dash)
  axislegend(ax)
  fig
end