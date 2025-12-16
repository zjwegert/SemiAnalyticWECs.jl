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

eta, w, ∂ₓ²w, ∂ₓ⁴w, v = data.η_s, data.w, data.∂ₓ²w, data.∂ₓ⁴w, data.v
q = G/(-im*ω)*v
M11 = -B*∂ₓ²w + η*v
p = D*∂ₓ⁴w - ω^2*Ib*w
XF, x = data.XF, data.x'
# Plotting
xf_range = 75:225

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=3) do
  fig = Figure(size = (1400, 800), figure_padding = (1,99,1,1));
  ax11 = Axis(fig[1,1],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Re(w)}",xticks=-15:5:15)
  lines!(ax11,data.XF[xf_range],vec(real(eta))[xf_range])
  lines!(ax11,x,vec(real(w)).-h)
  xlims!(ax11,-16,16)
  hidexdecorations!(ax11;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax21 = Axis(fig[2,1],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Re(v)}",xticks=-15:5:15)
  lines!(ax21,x,vec(real(v)))
  xlims!(ax21,-16,16)
  hidexdecorations!(ax21;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax31 = Axis(fig[3,1],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Re(q)}",xticks=-15:5:15)
  lines!(ax31,x,vec(real(q)))
  xlims!(ax31,-16,16)
  hidexdecorations!(ax31;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax41 = Axis(fig[4,1],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Re(M_{11})}",xticks=-15:5:15)
  lines!(ax41,x,vec(real(M11)))
  xlims!(ax41,-16,16)
  hidexdecorations!(ax41;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax51 = Axis(fig[5,1],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Re(p)}",xticks=-15:5:15)
  lines!(ax51,x,vec(real(p)))
  xlims!(ax51,-16,16)

  ax12 = Axis(fig[1,2],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Im(w)}",xticks=-15:5:15)
  lines!(ax12,data.XF[xf_range],vec(imag(eta))[xf_range])
  lines!(ax12,x,vec(imag(w)).-h)
  xlims!(ax12,-16,16)
  hidexdecorations!(ax12;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax22 = Axis(fig[2,2],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Im(v)}",xticks=-15:5:15)
  lines!(ax22,x,vec(imag(v)))
  xlims!(ax22,-16,16)
  hidexdecorations!(ax22;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax32 = Axis(fig[3,2],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Im(q)}",xticks=-15:5:15)
  lines!(ax32,x,vec(imag(q)))
  xlims!(ax32,-16,16)
  hidexdecorations!(ax32;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax42 = Axis(fig[4,2],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Im(M_{11})}",xticks=-15:5:15)
  lines!(ax42,x,vec(imag(M11)))
  xlims!(ax42,-16,16)
  hidexdecorations!(ax42;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax52 = Axis(fig[5,2],aspect=4,xlabel=L"x",ylabel=L"\mathrm{Im(p)}",xticks=-15:5:15)
  lines!(ax52,x,vec(imag(p)))
  xlims!(ax52,-16,16)

  fig
end

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
  lines!(ax,data.XF[xf_range],vec(real(eta))[xf_range])
  lines!(ax,x,vec(real(w)).-h)
  xlims!(ax,-16,16)
  hidexdecorations!(ax;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax2 = Axis(fig[2,1],aspect=4,xlabel=L"x~(\mathrm{m})",ylabel="Voltage (V)",xticks=-15:5:15)
  lines!(ax2,x,vec(real(v)))
  xlims!(ax2,-16,16)

  fig
end

save("$(@__DIR__)/figures/TimeDomain_PVDF_SimplySupported.png",fig;dpi=300)