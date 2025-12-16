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
Ts = 3;

## PZT5H
C11,e21,K22,rhop = PZT5H_material_coefficents()
dp = 1.1e-4;
G = K22/(2*dp);
C = K22/(2*dp);
B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
η = 1/2*(d0+dp)*e21;
Ib = rho0*d0 + 2*rhop*dp;

ω = 2π/Ts;
D = B + ω*η^2/(im*G+ω*C);
data = solve_submerged_plate_2d(ω,D,Ib,η,G,C,H,h,L,100,ceil(Int,L/0.025);bc_case=:simply_supported,return_displacements=true)

eta, w, v = data.η_s, data.w, data.v
XF, x = data.XF, data.x'
# Plotting
omega_t = 0;
xf_range = 75:225

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=3) do
  fig = Figure(size = (1400, 400), figure_padding = (1,99,1,1));
  ax = Axis(fig[1,1],aspect=4,xlabel=L"x~(\mathrm{m})",ylabel=L"z~(\mathrm{m})",xticks=-15:5:15)
  lines!(ax,data.XF[xf_range],vec(real(eta*exp(-im*omega_t)))[xf_range])
  lines!(ax,x,vec(real(w*exp(-im*omega_t))).-h)
  xlims!(ax,-16,16)
  hidexdecorations!(ax;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax2 = Axis(fig[2,1],aspect=4,xlabel=L"x~(\mathrm{m})",ylabel="Voltage (V)",xticks=-15:5:15)
  lines!(ax2,x,vec(real(v*exp(-im*omega_t))))
  xlims!(ax2,-16,16)

  fig
end;

save("$(@__DIR__)/figures/TimeDomain_PZT5H_SimplySupported.png",fig;dpi=300)

## PVDF
C11,e21,K22,rhop = PVDF_TechMan_material_coefficents()
dp = 1.1e-4;
G = K22/(2*dp);
C = K22/(2*dp);
B = d0^3/12*E0/(1-nu0^2) + 2*dp*(dp^2/3+d0*dp/2+d0^2/4)*C11;
η = 1/2*(d0+dp)*e21;
Ib = rho0*d0 + 2*rhop*dp;

ω = 2π/Ts;
D = B + ω*η^2/(im*G+ω*C);
data = solve_submerged_plate_2d(ω,D,Ib,η,G,C,H,h,L,100,ceil(Int,L/0.025);bc_case=:simply_supported,return_displacements=true)

eta, w, v = data.η_s, data.w, data.v
XF, x = data.XF, data.x'
# Plotting
omega_t = 0;
xf_range = 75:225

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=3) do
  fig = Figure(size = (1400, 400), figure_padding = (1,99,1,1));
  ax = Axis(fig[1,1],aspect=4,xlabel=L"x~(\mathrm{m})",ylabel=L"z~(\mathrm{m})",xticks=-15:5:15)
  lines!(ax,data.XF[xf_range],vec(real(eta*exp(-im*omega_t)))[xf_range])
  lines!(ax,x,vec(real(w*exp(-im*omega_t))).-h)
  xlims!(ax,-16,16)
  hidexdecorations!(ax;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax2 = Axis(fig[2,1],aspect=4,xlabel=L"x~(\mathrm{m})",ylabel="Voltage (V)",xticks=-15:5:15)
  lines!(ax2,x,vec(real(v*exp(-im*omega_t))))
  xlims!(ax2,-16,16)

  fig
end

save("$(@__DIR__)/figures/TimeDomain_PVDF_SimplySupported.png",fig;dpi=300)