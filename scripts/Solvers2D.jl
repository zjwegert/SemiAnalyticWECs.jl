include("DispersionEquations.jl")
include("EigenModes1D.jl")
include("GreensFunctions.jl")

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
  ϕ =  transpose(ξ)*transpose(ϕ_r) + transpose(ϕ_d);

  # Compute reflection and transmission coefficients
  R = -1/(tan(k*H) + k*H*sec(k*H)^2)*transpose(exp.(-k*x)).*(α*ϕ - im*ω*w) ⋅ diag(ws)
  T = 1 - 1/(tan(k*H) + k*H*sec(k*H)^2)*transpose(exp.(k*x)).*(α*ϕ - im*ω*w) ⋅ diag(ws);
  k = imag(k) # k = k/im;
  Cg = ω/(2*k)*(1+2*k*H/sinh(2*k*H)); # Group velocity
  A = im*ω/g;

  # Far-field power takeoff
  P_farfield = 1/2*ρ_w*g*Cg*abs(A)^2*(1 - abs(R)^2 - abs(T)^2);

  # Compute near-field power takeoff
  ∂ₓ²w = transpose(ξ)*∂ₓ²u;
  P_nearfield = Gp*ω^2/2*abs.(ηp*∂ₓ²w/(Gp-im*ω*Cp)).^2 ⋅ diag(ws);

  if return_displacements
    return R, T, P_farfield, P_nearfield, Cg, w, ϕ
  else
    return R, T, P_farfield, P_nearfield, Cg
  end
end


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
  α = ω.^2/g

  # Calculation points, here we instead use midpoints
  x = range(-L,L,n+1);
  δx = x[2]-x[1];
  x[1:end-1] .+= δx/2; # use midpoints

  # Compute eigenvalues and eigenmodes for plate
  λ,u,∂ₓ²u = eigenmodes_1d(bc_case,N,L,x)

  # TODO: need following functions
  # - G_int_z_zeta_regular
  # - G_int_z
end