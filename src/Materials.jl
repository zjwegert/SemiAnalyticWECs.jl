"""
    PVDF_TechMan_material_coefficents()

Material coefficients `(C11,e21,K22,rho)` for PVDF in plane-stress configuration.

Values based on:
  (1) - Measurement Specialities Inc. Piezo Film Sensors Technical Manual. Hampton, VA.
  (2) - Kepler RG, Ferroelectric. Pyroelectric, and piezoelectric properties of poly(vinylidene fluoride).
        In: Nalwa HS, editor. Ferroelectric polymers: chemistry, physics and applications. New York: Marcel Dekker; 1995. p. 183–232.
"""
function PVDF_TechMan_material_coefficents()

  #thickness ~ 0.11 mm
  rho = 1780; #kg/m^3) (1)

  d21 = 23e-12;#% C/N (1)
  Ep = 3e+9; # Pa (1)
  nup = 0.392; #% (2)
  epss = 113e-12; # F/m (1)

  # Convert, plane-stress (based on Renzi)
  C11 = Ep/(1-nup^2); # Pa
  e21 = d21*Ep/(1-nup); # C/m^2
  K22 = epss; # F/m

  return (;C11,e21,K22,rho)
end

"""
    PZT5H_material_coefficents(θ)

Material coefficients `(C11,e21,K22,rho)` for PZT5H in plane-stress with poling angle `θ`.

Values based on Yang (2018). An Introduction to the Theory of Piezoelectricity. Doi: 10.1007/978-3-030-03137-4.
"""
function PZT5H_material_coefficents(θ)

  rho = 7500; #(kg/m^3)

  C = [12.6 7.95 8.41 0 0 0
      7.95 12.6 8.41 0 0 0
      8.41 8.41 11.70 0 0 0
      0 0 0 2.30 0 0
      0 0 0 0 2.3 0
      0 0 0 0 0 2.33]*10^10;

  e = [0 0 0 0 17 0
      0 0 0 17 0 0
      -6.5 -6.5 23.3 0 0 0];

  kappa = [1700 0 0
      0 1700 0
      0 0 1470]*8.854*10^-12;

  # Plane-stress formulation
  Sᴱ = inv(C); # Compliance at constant electric field
  d = e*Sᴱ; # Piezoelectric strain coefficients
  Kᵀ = kappa + d*transpose(e); # Permittivity at constant stress

  C11=(8*Sᴱ[2,2])/(Sᴱ[2,2]*Sᴱ[5,5]-sin(θ)^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])+2*cos(θ)^2*sin(θ)^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])-cos(θ)^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5])))
  e31=(8*(cos(θ)^3*(d[3,2]*Sᴱ[1,2]-d[3,1]*Sᴱ[2,2])+cos(θ)*sin(θ)^2*(d[1,5]*Sᴱ[2,2]-d[3,3]*Sᴱ[2,2]+d[3,2]*Sᴱ[2,3])))/(-Sᴱ[2,2]*Sᴱ[5,5]+sin(θ)^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])-2*cos(θ)^2*sin(θ)^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])+cos(θ)^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5])))
  kappa33=1/Sᴱ[2,2]*(-cos(θ)^2*d[3,2]^2+cos(θ)^2*Kᵀ[3,3]*Sᴱ[2,2]+Kᵀ[1,1]*sin(θ)^2*Sᴱ[2,2]+(8*(cos(θ)^3*(d[3,2]*Sᴱ[1,2]-d[3,1]*Sᴱ[2,2])+cos(θ)*sin(θ)^2*(d[1,5]*Sᴱ[2,2]-d[3,3]*Sᴱ[2,2]+d[3,2]*Sᴱ[2,3]))^2)/(-Sᴱ[2,2]*Sᴱ[5,5]+sin(θ)^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])-2*cos(θ)^2*sin(θ)^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])+cos(θ)^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5]))))

  return (;C11,e31,kappa33,rho)
end