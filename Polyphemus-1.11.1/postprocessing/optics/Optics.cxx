// Copyright (C) 2006-2008, ENPC - INRIA - EDF R&D
// Author(s): Marilyne Tombette
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// Polyphemus is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// Polyphemus is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/


#include <cmath>


template<class T>
T max_(T x, T y)
{
  if (x <= y)
    return y;
  else
    return x;
}


template<class T>
T min_(T x, T y)
{
  if (x <= y)
    return x;
  else
    return y;
}


template<class T>
T sqr(T x)
{
  return x * x;
}


template<class T>
T sqr3(T x)
{
  return x * x * x;
}


template<class T>
T sqrt2(T a)
{
  return sqrt(double(a));
}


template<class T>
T arccos2(T a)
{
  return acos(double(a));
}


template<class T>
T sign2(T a)
{
  if (a != T(0.0))
    return a / T(abs(double(a)));
  else
    return T(1.0);
}


// Square root of a complex number.
template<class T>
void sqrt_for_complex(T a, T b, T& a_out, T& b_out)
{
  T rho_in, phi_in, rho_out, phi_out;
  rho_in = sqrt2(sqr(a) + sqr(b));
  phi_in = sign2(b) * arccos2(T(a / rho_in));
  rho_out = sqrt2(rho_in);
  phi_out = phi_in / T(2.0);
  a_out = rho_out * cos(phi_out);
  b_out = rho_out * sin(phi_out);
}


// Square of a complex number.
template<class T>
void sqr_for_complex(T a, T b, T& a_out, T& b_out)
{
  a_out = sqr(a) - sqr(b);
  b_out = T(2.0) * a * b;
}


// Quotient of two complex number.
template<class T>
void quotient_for_complex(T divisor_real, T divisor_imaginary,
                          T dividend_real, T dividend_imaginary,
                          T& quotient_real, T& quotient_imaginary)
{
  quotient_real = (divisor_real * dividend_real
                   + divisor_imaginary * dividend_imaginary)
    / (sqr(dividend_real) + sqr(dividend_imaginary));

  quotient_imaginary = (divisor_imaginary * dividend_real
                        - divisor_real * dividend_imaginary)
    / (sqr(dividend_real) + sqr(dividend_imaginary));
}


//////////////////////////////////////////////////////////////////////////////


template<class real>
void compute_Hanel_diameter(real dry_radius, real relative_humidity,
                            real& wet_radius)
{
  // epsilon = 0.25 for organics and 0.285 for sulfate.
  const real epsilon = 0.25;
  wet_radius = dry_radius * exp(-epsilon * log(1 - relative_humidity));
}


template<class real>
void compute_Gerber_diameter(real dry_diameter, real relative_humidity,
                             real temperature, real& wet_diameter)
{
  // 'dry_diameter' and 'wet_diameter' unit is \mu m.
  const real C1(0.4989352162271429),
    C2(0.3026183900844475e1),
    C3(0.5372215625062934e-12),
    C4(-0.1371059101078550e1),
    C5(0.3942463621284677e-02);
  const real temperature_ref(298.0);
  real aa, dry_radius, wet_radius;

  dry_radius = (dry_diameter / 2) * 1.e-4;
  aa = C3 * (1.0 + C5 * (temperature_ref - temperature));
  wet_radius = wpower(C1 * wpower(dry_radius, C2) /
                      abs(aa * wpower(dry_radius, C4)
                          - log(relative_humidity))
                      + wpower(dry_radius, real(3.0)), real(1.0 / 3.0));
  wet_diameter = 2.0 * wet_radius * 1.e4;
}


template<class real>
void compute_wet_diameter_from_water_content(real dry_diameter,
                                             real dry_aerosol_concentration,
                                             real water_concentration,
                                             real& wet_diameter)
{
  // 'dry_aerosol_concentration' is the total dry aerosol concentration in
  // aerosols. 'water_concentration' is the total water concentration in
  // aerosols. Units of both concentrations are: \mu g.m^{-3}.

  // Unit of densities is g.m^{-3}.
  const real dry_aerosol_density = 1.4e6;
  const real water_density = 1.0e6;

  const real pi = 3.14159265358979323846264;
  real global_aerosol_concentration, aerosol_number, global_aerosol_density;

  global_aerosol_concentration =
    dry_aerosol_concentration + water_concentration;

  aerosol_number = 6. * dry_aerosol_concentration /
    (pi * wpower(dry_diameter, real(3.0)) * dry_aerosol_density);

  global_aerosol_density = (dry_aerosol_density * dry_aerosol_concentration
                            + water_density * water_concentration)
    / global_aerosol_concentration;

  wet_diameter = wpower(real(6.) * global_aerosol_concentration
                        / (pi * global_aerosol_density * aerosol_number),
                        real(1.0 / 3.0));
}
