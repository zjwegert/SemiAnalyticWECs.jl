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
Ts = 4;

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
xf_range = 91:211

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=3) do
  fig = Figure(size = (1400, 300*3), figure_padding = (0,0,0,0));
  ax11_l = Axis(fig[1,1],width=400,height=300,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(\hat{\zeta})~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax11_r = Axis(fig[1,1],width=400,height=300,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Re(w)-h~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax11_l,data.XF[xf_range],vec(real(eta))[xf_range],color=Makie.wong_colors()[1])
  lines!(ax11_r,x,vec(real(w)).-h,color=Makie.wong_colors()[2])
  xlims!(ax11_r,-12.5,12.5)
  xlims!(ax11_l,-12.5,12.5)
  ylims!(ax11_r,-2.8,1.1)
  ylims!(ax11_l,-2.8,1.1)
  hidespines!(ax11_r)
  hidexdecorations!(ax11_r)
  hidexdecorations!(ax11_l;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax21_l = Axis(fig[2,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(M_{11})~(N)}",xticks=-10:5:10,
    yticks=[-300,0,300],ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax21_r = Axis(fig[2,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Re(p)~(Pa)}",xticks=-10:5:10,
    yticks=[-1000,0,1000],ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax21_l,x,vec(real(M11)),color=Makie.wong_colors()[1])
  lines!(ax21_r,x,vec(real(p)),color=Makie.wong_colors()[2])
  xlims!(ax21_r,-12.5,12.5)
  xlims!(ax21_l,-12.5,12.5)
  ylims!(ax21_r,-1500,1500)
  ylims!(ax21_l,-450,450)
  hidespines!(ax21_r)
  hidexdecorations!(ax21_r)
  hidexdecorations!(ax21_l;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax31 = Axis(fig[3,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(v)~(V)}",xticks=-10:5:10)
  lines!(ax31,x,vec(real(v)))
  xlims!(ax31,-12.5,12.5)
  hidexdecorations!(ax31;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax41 = Axis(fig[4,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(q)~(Cm^{-2})}",xticks=-10:5:10)
  lines!(ax41,x,vec(real(q)))
  xlims!(ax41,-12.5,12.5)


  ax12_l = Axis(fig[1,2],width=400,height=300,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(\hat{\zeta})~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax12_r = Axis(fig[1,2],width=400,height=300,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Im(w)-h~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax12_l,data.XF[xf_range],vec(imag(eta))[xf_range],color=Makie.wong_colors()[1])
  lines!(ax12_r,x,vec(imag(w)).-h,color=Makie.wong_colors()[2])
  xlims!(ax12_r,-12.5,12.5)
  xlims!(ax12_l,-12.5,12.5)
  ylims!(ax12_r,-2.8,1.1)
  ylims!(ax12_l,-2.8,1.1)
  hidespines!(ax12_r)
  hidexdecorations!(ax12_r)
  hidexdecorations!(ax12_l;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax22_l = Axis(fig[2,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(M_{11})~(N)}",xticks=-10:5:10,
    yticks=[-300,0,300],ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax22_r = Axis(fig[2,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Im(p)~(Pa)}",xticks=-10:5:10,
    yticks=[-600,0,600],ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax22_l,x,vec(imag(M11)),color=Makie.wong_colors()[1])
  lines!(ax22_r,x,vec(imag(p)),color=Makie.wong_colors()[2])
  xlims!(ax22_r,-12.5,12.5)
  xlims!(ax22_l,-12.5,12.5)
  ylims!(ax22_r,-800,800)
  ylims!(ax22_l,-400,400)
  hidespines!(ax22_r)
  hidexdecorations!(ax22_r)
  hidexdecorations!(ax22_l;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax32 = Axis(fig[3,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(v)~(V)}",xticks=-10:5:10)
  lines!(ax32,x,vec(imag(v)))
  xlims!(ax32,-12.5,12.5)
  hidexdecorations!(ax32;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax42 = Axis(fig[4,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(q)~(Cm^{-2})}",xticks=-10:5:10)
  lines!(ax42,x,vec(imag(q)))
  xlims!(ax42,-12.5,12.5)

  fig
end

save("$(@__DIR__)/figures/SingleFreq_PZT5H_SimplySupported.png",fig;dpi=300)

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

eta, w, ∂ₓ²w, ∂ₓ⁴w, v = data.η_s, data.w, data.∂ₓ²w, data.∂ₓ⁴w, data.v
q = G/(-im*ω)*v
M11 = -B*∂ₓ²w + η*v
p = D*∂ₓ⁴w - ω^2*Ib*w
XF, x = data.XF, data.x'
# Plotting
xf_range = 91:211

fig = with_theme(theme_latexfonts(),fontsize=28,linewidth=3) do
  fig = Figure(size = (1400, 300*3), figure_padding = (0,0,0,0));
  ax11_l = Axis(fig[1,1],width=400,height=300,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(\hat{\zeta})~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax11_r = Axis(fig[1,1],width=400,height=300,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Re(w)-h~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax11_l,data.XF[xf_range],vec(real(eta))[xf_range],color=Makie.wong_colors()[1])
  lines!(ax11_r,x,vec(real(w)).-h,color=Makie.wong_colors()[2])
  xlims!(ax11_r,-12.5,12.5)
  ylims!(ax11_r,-3,1.1)
  ylims!(ax11_l,-3,1.1)
  hidespines!(ax11_r)
  hidexdecorations!(ax11_r)
  hidexdecorations!(ax11_l;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax21_l = Axis(fig[2,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(M_{11})~(N)}",xticks=-10:5:10,
    yticks=[-100,0,100],ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax21_r = Axis(fig[2,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Re(p)~(Pa)}",xticks=-10:5:10,
    yticks=[-1000,0,1000],ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax21_l,x,vec(real(M11)),color=Makie.wong_colors()[1])
  lines!(ax21_r,x,vec(real(p)),color=Makie.wong_colors()[2])
  xlims!(ax21_r,-12.5,12.5)
  xlims!(ax21_l,-12.5,12.5)
  ylims!(ax21_l,-130,130)
  ylims!(ax21_r,-1300,1300)
  hidespines!(ax21_r)
  hidexdecorations!(ax21_r)
  hidexdecorations!(ax21_l;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax31 = Axis(fig[3,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(v)~(V)}",xticks=-10:5:10)
  lines!(ax31,x,vec(real(v)))
  xlims!(ax31,-12.5,12.5)
  hidexdecorations!(ax31;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax41 = Axis(fig[4,1],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Re(q)~(Cm^{-2})}",xticks=-10:5:10)
  lines!(ax41,x,vec(real(q)))
  xlims!(ax41,-12.5,12.5)


  ax12_l = Axis(fig[1,2],width=400,height=300,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(\hat{\zeta})~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax12_r = Axis(fig[1,2],width=400,height=300,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Im(w)-h~(m)}",xticks=-10:5:10,
    ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax12_l,data.XF[xf_range],vec(imag(eta))[xf_range],color=Makie.wong_colors()[1])
  lines!(ax12_r,x,vec(imag(w)).-h,color=Makie.wong_colors()[2])
  xlims!(ax12_r,-12.5,12.5)
  ylims!(ax12_r,-3,1.1)
  ylims!(ax12_l,-3,1.1)
  hidespines!(ax12_r)
  hidexdecorations!(ax12_r)
  hidexdecorations!(ax12_l;label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)

  ax22_l = Axis(fig[2,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(M_{11})~(N)}",xticks=-10:5:10,
    yticks=[-50,-25,0,25,50],ylabelcolor=Makie.wong_colors()[1],yticklabelcolor=Makie.wong_colors()[1])
  ax22_r = Axis(fig[2,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",yaxisposition=:right,ylabel=L"\mathrm{Im(p)~(Pa)}",xticks=-10:5:10,
    yticks=[-400,-200,0,200,400],ylabelcolor=Makie.wong_colors()[2],yticklabelcolor=Makie.wong_colors()[2])
  lines!(ax22_l,x,vec(imag(M11)),color=Makie.wong_colors()[1])
  lines!(ax22_r,x,vec(imag(p)),color=Makie.wong_colors()[2])
  xlims!(ax22_r,-12.5,12.5)
  xlims!(ax22_l,-12.5,12.5)
  ylims!(ax22_l,-50,70)
  ylims!(ax22_r,-50*8,70*8)
  hidespines!(ax22_r)
  hidexdecorations!(ax22_r)
  hidexdecorations!(ax22_l;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax32 = Axis(fig[3,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(v)~(V)}",xticks=-10:5:10)
  lines!(ax32,x,vec(imag(v)))
  xlims!(ax32,-12.5,12.5)
  hidexdecorations!(ax32;label = true, ticklabels = true, ticks = true, grid = false,
      minorgrid = false, minorticks = false)

  ax42 = Axis(fig[4,2],width=400,height=150,xlabel=L"x~(\mathrm{m})",ylabel=L"\mathrm{Im(q)~(Cm^{-2})}",xticks=-10:5:10)
  lines!(ax42,x,vec(imag(q)))
  xlims!(ax42,-12.5,12.5)

  fig
end

save("$(@__DIR__)/figures/SingleFreq_PVDF_SimplySupported.png",fig;dpi=300)