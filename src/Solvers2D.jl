"""
    solve_surface_plate_2d(
      ω,D,Ib,ηp,Gp,Cp,H,L,N,n;
      g = 9.81,
      ρ_w = 1025,
      return_displacements=false,
      bc_case=:clamped
    )
Solve the 2D wave-plate interaction problem.

Inputs:
- ω: frequency
- D: plate stiffness
- Ib: plate mass
- ηp: piezoelectric coupling
- Gp: surface conductance
- Cp: permittivity
- H: water depth
- L: beam/plate half-length
- N: number of modes
- n: number of calculation points along the plate is 2n+1

Optional Inputs:
- g: gravitational acceleration (default 9.81 m/s²)
- ρ_w: water density (default 1025 kg/m³)
- return_displacements: if true, also return plate displacements and potentials (default false)
- bc_case: boundary condition case for plate eigenmodes (default :clamped)

Outputs:
- R: reflection coefficient
- T: transmission coefficient
- P_farfield: far-field power takeoff
- P_nearfield: near-field power takeoff
- Cg: group velocity
- w (optional): plate displacements
- ϕ (optional): wave potentials
"""
function solve_surface_plate_2d(
  ω,D,Ib,ηp,Gp,Cp,H,L,N,n;
  g = 9.81,
  ρ_w = 1025,
  return_displacements=false,
  bc_case=:clamped
)
  D = D' # Due to choice of exp(-k*x) for incident wave
  np = 2n + 1
  α = ω.^2/g

  # Calculation points
  x = -L:L/n:L;

  # Compute eigenvalues and eigenmodes for plate
  λ,u,∂ₓ²u = eigenmodes_1d(bc_case,N,L,x)

  # Calculate free-surface Green function
  G = matrix_G_surface(α,H,L,n);

  # Solve for the plate radiation potentials
  ϕ_r = (I-α*G)\G*(-im*ω)*transpose(u);

  # Setup Simpson's rule weights
  δx = (2*L)/(np-1);
  diagws = 2*δx/3*ones(np);
  diagws[1] = δx/3; diagws[2:2:np-1] .= 4*δx/3; diagws[end] = δx/3;
  ws = Diagonal(diagws);

  # Damping matrix
  B = ρ_w*transpose(u*ws*ϕ_r);

  # Incident wave potential
  k = first(dispersion_free_surface(α,0,H));
  ϕ_i = exp.(-k*x);

  # Diffraction wave potential
  ϕ_d = (I-α*G)\ϕ_i;

  # Compute ξ
  diff = ρ_w*u*ws*ϕ_d; # calculate diffraction matrix
  f = -im*ω*diff;
  ξ = (ρ_w*g*I + D*Diagonal(λ.^4) - ω.^2*Ib*I + im*ω*B)\f; # calculated α_n

  w = transpose(ξ)*u;
  η_s =  transpose(ξ)*transpose(ϕ_r) + transpose(ϕ_d);

  # Compute reflection and transmission coefficients
  R = -1/(tan(k*H) + k*H*sec(k*H)^2)*transpose(exp.(-k*x)).*(α*η_s - im*ω*w) ⋅ diag(ws)
  T = 1 - 1/(tan(k*H) + k*H*sec(k*H)^2)*transpose(exp.(k*x)).*(α*η_s - im*ω*w) ⋅ diag(ws);
  k = imag(k) # k = k/im;
  Cg = ω/(2*k)*(1+2*k*H/sinh(2*k*H)); # Group velocity

  # Far-field power takeoff
  A = 1;
  P_farfield = 1/2*ρ_w*g*Cg*abs(A)^2*(1 - abs(R)^2 - abs(T)^2);

  # Compute near-field power takeoff
  ∂ₓ²w = transpose(ξ)*∂ₓ²u;
  P_nearfield = Gp*ω^2/2*abs.(ηp*∂ₓ²w/(Gp-im*ω*Cp)).^2 ⋅ diag(ws) * (g/ω)^2; # normalise by (g/ω)^2

  if return_displacements
    v = im*ω*ηp/(Gp-im*ω*Cp)*∂ₓ²w; # voltage
    return (;R, T, P_farfield, P_nearfield, Cg, w, η_s, v)
  else
    return (;R, T, P_farfield, P_nearfield, Cg)
  end
end

solve_surface_plate_2d(ω,D,Ib,ηp,Gp,Cp,H,h,L,N,n;kwargs...) = solve_surface_plate_2d(ω,D,Ib,ηp,Gp,Cp,H,L,N,n;kwargs...)

"""
    solve_submerged_plate_2d(
      ω,D,Ib,ηp,Gp,Cp,H,h,L,N,n;
      g = 9.81,
      ρ_w = 1025,
      return_displacements=false,
      bc_case=:clamped
    )
Solve the 2D wave-plate interaction problem.

Inputs:
- ω: frequency
- D: plate stiffness
- Ib: plate mass
- ηp: piezoelectric coupling
- Gp: surface conductance
- Cp: permittivity
- H: water depth
- h: submergence depth
- L: beam/plate half-length
- N: number of modes
- n: number of calculation points along the plate is 2n+1

Optional Inputs:
- g: gravitational acceleration (default 9.81 m/s²)
- ρ_w: water density (default 1025 kg/m³)
- return_displacements: if true, also return plate displacements and potentials (default false)
- bc_case: boundary condition case for plate eigenmodes (default :clamped)

Outputs:
- R: reflection coefficient
- T: transmission coefficient
- P_farfield: far-field power takeoff
- P_nearfield: near-field power takeoff
- Cg: group velocity
- w (optional): plate displacements
- ϕ (optional): wave potentials
"""
function solve_submerged_plate_2d(
  ω,D,Ib,ηp,Gp,Cp,H,h,L,N,n;
  g = 9.81,
  ρ_w = 1025,
  return_displacements=false,
  bc_case=:clamped
)
  A = 1; # incident wave amplitude
  α = ω.^2/g

  # Calculation points, here we instead use midpoints
  x = collect(range(-L,L,n+1));
  δx = x[2]-x[1];
  x = transpose(x[1:end-1] .+ δx/2); # use midpoints

  # Compute eigenvalues and eigenmodes for plate
  λ,u,∂ₓ²u = eigenmodes_1d(bc_case,N,L,x) # TODO: could be cached
  u = transpose(u);
  ∂ₓ²u = transpose(∂ₓ²u);

  # Computing BIE matrix
  k = first(dispersion_free_surface(α,0,H));
  N₀² = 1/2*(1 - (sin(k*H)^2)/α/H);
  k = k/im;

  G = zeros(ComplexF64,n);
  G[1] = 4/δx+∂z∂ζ_regular_greens_submerged_2d(0,-h,-h,H,α)*δx;
  for ii=2:n
    dist = (ii-1)*δx;
    G[ii] = -1/(dist-δx/2)+1/(dist+δx/2) + ∂z∂ζ_regular_greens_submerged_2d(dist,-h,-h,H,α)*δx;
  end

  K = 1/(2*π)*G[round.(Int,abs.((1:n).-transpose(1:n)).+1)];

  # solving diffraction problem
  f = -transpose(-im*A*g*k/(ω*cosh(k*H))*sinh(k*(H-h))*exp.(im*k*x));
  ϕ_di_jump = (-K)\f; # jump in potential of diffraction problem

  # solving radiation problems
  ϕ_rad_jump = (-K)\(-im*ω*u); # jump in potential of radiation problems

  # Solving coupled elasticity problem from dry modes expansion
  K_stiff = Diagonal(D*λ.^4);
  M_mass = Ib*I;
  M_added = δx*u'*ϕ_rad_jump;

  # Forcing vector
  F = -im*ω*ρ_w*δx*u'*ϕ_di_jump;

  # Coefficients of dry modes expansion
  c = (K_stiff-ω^2*M_mass+im*ω*ρ_w*M_added)\F;

  # Jump in potential of coupled problem
  ϕ_jump = ϕ_di_jump + ϕ_rad_jump*c;

  # Reflection and transmission
  coeff = ω*sinh(k*(H-h))*cosh(k*H)/(2*A*g*H*N₀²);
  R = -first(coeff*exp.(im*k*x)*ϕ_jump*δx);
  T = 1 - first(coeff*exp.(-im*k*x)*ϕ_jump*δx);
  Cg = real(ω/(2*k)*(1+2*k*H/sinh(2*k*H)));
  P_farfield = 1/2*ρ_w*g*abs(A)^2*Cg*(1 - abs(R)^2 - abs(T)^2);

  # Near-field power takeoff
  ∂ₓ²w = ∂ₓ²u*c;
  P_nearfield = Gp*ω^2/2*abs(ηp/(Gp-im*ω*Cp))^2*(abs(first(∂ₓ²w'*∂ₓ²w))*δx);

  if !return_displacements
      return (;R, T, P_farfield, P_nearfield, Cg)
  else
    # Evaluate displacements if plotting
    XF = collect(range(-3*L,3*L,301)); # free surface points
    dzG = zeros(ComplexF64,length(XF),length(x));
    w = u*c; # deflection of plate
    v = im*ω*ηp/(Gp-im*ω*Cp)*∂ₓ²w; # voltage

    # Greens function matrix for free surface
    for ii=axes(dzG,1)
        for jj=axes(dzG,2)
            dzG[ii,jj] = ∂z_regular_greens_submerged_2d(abs(XF[ii]-x[jj]),-h,0,H,α)/2/π;
        end
    end

    # free surface elevation
    η_inc = A*exp.(im*k*XF);
    η_sc = -im*ω/g*dzG*ϕ_jump*δx;
    η_s = η_inc+η_sc;

    return (;R, T, P_farfield, P_nearfield, Cg, w, η_s, v)
  end
end