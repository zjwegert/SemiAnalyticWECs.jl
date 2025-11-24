# Values based on:
# (1) - Measurement Specialities Inc. Piezo Film Sensors Technical Manual. Hampton, VA.
# (2) - Kepler RG, Ferroelectric. Pyroelectric, and piezoelectric properties of poly(vinylidene fluoride). In: Nalwa HS, editor. Ferroelectric polymers: chemistry, physics and applications. New York: Marcel Dekker; 1995. p. 183–232.
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

  C11=(8*Sᴱ[2,2])/(Sᴱ[2,2]*Sᴱ[5,5]-s^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])+2*c^2*s^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])-c^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5])))
  e31=(8*(c^3*(d[3,2]*Sᴱ[1,2]-d[3,1]*Sᴱ[2,2])+c*s^2*(d[1,5]*Sᴱ[2,2]-d[3,3]*Sᴱ[2,2]+d[3,2]*Sᴱ[2,3])))/(-Sᴱ[2,2]*Sᴱ[5,5]+s^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])-2*c^2*s^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])+c^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5])))
  kappa33=1/Sᴱ[2,2]*(-c^2*d[3,2]^2+c^2*Kᵀ[3,3]*Sᴱ[2,2]+Kᵀ[1,1]*s^2*Sᴱ[2,2]+(8*(c^3*(d[3,2]*Sᴱ[1,2]-d[3,1]*Sᴱ[2,2])+c*s^2*(d[1,5]*Sᴱ[2,2]-d[3,3]*Sᴱ[2,2]+d[3,2]*Sᴱ[2,3]))^2)/(-Sᴱ[2,2]*Sᴱ[5,5]+s^4*(8*Sᴱ[2,3]^2-8*Sᴱ[2,2]*Sᴱ[3,3]+Sᴱ[2,2]*Sᴱ[5,5])-2*c^2*s^2*(8*Sᴱ[1,3]*Sᴱ[2,2]-8*Sᴱ[1,2]*Sᴱ[2,3]+3*Sᴱ[2,2]*Sᴱ[5,5])+c^4*(8*Sᴱ[1,2]^2+Sᴱ[2,2]*(-8*Sᴱ[1,1]+Sᴱ[5,5]))))

  return (;C11,e31,kappa33,rho)
end


# Voigt to Tensor conversions

"""
  Given a material constant given in Voigt notation,
  return a SymFourthOrderTensorValue using ordering from Gridap
"""
function voigt2tensor4(A::Array{M,2}) where M
  if isequal(size(A),(3,3))
    return SymFourthOrderTensorValue(A[1,1], A[3,1], A[2,1],
                                     A[1,3], A[3,3], A[2,3],
                                     A[1,2], A[3,2], A[2,2])
  elseif isequal(size(A),(6,6))
    return SymFourthOrderTensorValue(A[1,1], A[6,1], A[5,1], A[2,1], A[4,1], A[3,1],
                                     A[1,6], A[6,6], A[5,6], A[2,6], A[4,6], A[3,6],
                                     A[1,5], A[6,5], A[5,5], A[2,5], A[4,5], A[3,5],
                                     A[1,2], A[6,2], A[5,2], A[2,2], A[4,2], A[3,2],
                                     A[1,4], A[6,4], A[5,4], A[2,4], A[4,4], A[3,4],
                                     A[1,3], A[6,3], A[5,3], A[2,3], A[4,3], A[3,3])
  else
      @notimplemented
  end
end

"""
  Given a material constant given in Voigt notation,
  return a ThirdOrderTensorValue using ordering from Gridap
"""
function voigt2tensor3(A::Array{M,2}) where M
  if isequal(size(A),(2,3))
    return ThirdOrderTensorValue(A[1,1], A[2,1], A[1,3], A[2,3], A[1,3], A[2,3], A[1,2], A[2,2])
  elseif isequal(size(A),(3,6))
    return ThirdOrderTensorValue(
      A[1,1], A[2,1], A[3,1], A[1,6], A[2,6], A[3,6], A[1,5], A[2,5], A[3,5],
      A[1,6], A[2,6], A[3,6], A[1,2], A[2,2], A[3,2], A[1,4], A[2,4], A[3,4],
      A[1,5], A[2,5], A[3,5], A[1,4], A[2,4], A[3,4], A[1,3], A[2,3], A[3,3])
  else
    @notimplemented
  end
end

"""
  Given a material constant given in Voigt notation,
  return a SymTensorValue using ordering from Gridap
"""
function voigt2tensor2(A::Array{M,2}) where M
  return TensorValue(A)
end