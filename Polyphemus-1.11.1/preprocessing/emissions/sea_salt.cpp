// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Marilyne Tombette and Karine Sartelet
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

// This program computes sea salt emissions.


//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////


// Rate of sea-salt aerosol generation.
// Multiply this expression by r**(-3) to have real dF0dr.
template <class T>
T dF0dr(string ss_parameterization, const T& radius, const T& wind_10_m)
{
  T dF0dr;

  if (ss_parameterization == "Monahan")
    {
      // Rate of sea-salt aerosol generation per unit area
      // by indirect mechanism (bubbles), according to Monahan (1986).
      T B = (0.38 - log10(radius)) / 0.65;

      dF0dr = 1.373 * pow(wind_10_m, T(3.41))
        * (1. + 0.057 * pow(radius, T(1.05)))
        * pow(10., 1.19 * exp(-B * B));
    }
  else
    {
      // Rate of sea-salt aerosol generation per unit area
      // by indirect mechanism (bubbles) and direct mechanism (spume)
      //  according to Smith and Harrison (1998).
      T A1 = 0.2 * pow(wind_10_m, T(3.5));
      T A2 = 6.8 * 0.001 * pow(wind_10_m, T(3.0));
      const T f1(1.5), f2(1.0), r01(3.0), r02(30.0);

      dF0dr = pow(radius, T(3.0))
        * (A1 * exp(-f1 * pow(log(radius / r01), T(2.0))) +
           A2 * exp(-f2 * pow(log(radius / r02), T(2.0))));
    }

  return dF0dr;
}

// Computes wet radius of aerosols with Gerber's formula.
template <class T>
T WetRadius(const T& dry_radius,
            const T& relative_humidity,
            const T& temperature)
{
  const T C1(0.7674), C2(3.079), C4(1.424);

  T C3(2.573e-11);
  C3 *= (1.0 + 0.004 * (298.0 - temperature));

  T WetRadius;

  // Conversion: from µm to cm.
  T dry_radius_cm = 1.e-4 * dry_radius;

  WetRadius = pow(C1 * pow(dry_radius_cm, C2)
                  / (C3 * pow(dry_radius_cm, C4) - log(relative_humidity))
                  + pow(dry_radius_cm, T(3.)), T(1. / 3.));

  // Conversion: from cm to µm.
  WetRadius *= 1.e4;

  return WetRadius;
}


// Integrates dF0/dr * r^3 using Simpson method over interval [radius_inf,
// radius_sup].
template <class T>
T SimpsonIntegral(string ss_parameterization,
                  const T& radius_inf,
                  const T& radius_sup,
                  const T& wind_10_m)
{
  T SimpsonIntegral(0.0);

  // Number of iterations for Simpson's method.
  int nb_iteration = 10;

  nb_iteration = 2 * nb_iteration;

  T dradius = (radius_sup - radius_inf) / nb_iteration;
  T sum_unpaired(0.0), sum_paired(0.0);

  for (int i = 1; i < nb_iteration; i++)
    {
      T xr = radius_inf + i * dradius;
      if (i % 2 != 0) // i is odd.
        sum_unpaired += dF0dr(ss_parameterization, xr, wind_10_m);
      else    // i is even.
        sum_paired += dF0dr(ss_parameterization, xr, wind_10_m);
    }

  SimpsonIntegral =  dradius / 3. *
    (dF0dr(ss_parameterization, radius_inf, wind_10_m)
     + 4. * sum_paired
     + 2. * sum_unpaired
     + dF0dr(ss_parameterization, radius_sup, wind_10_m));

  return SimpsonIntegral;
}


//////////
// MAIN //
//////////

int main(int argc, char* argv[])
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("sea_salt.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  Date date_prev(date_beg);
  date_prev.AddDays(-1);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Constants.
  const real pi = 3.14159265358979323846264;

  // Output domain.
  string date_meteo_str;
  int Nt, Nt_sea_salt, Nz, Ny, Nx;
  real Delta_t, Delta_t_sea_salt, Delta_y, Delta_x;
  real t_min, y_min, x_min;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nx", "> 0", Nx);
  configuration.PeekValue("Delta_t", "> 0", Delta_t);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);
  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);

  configuration.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw string("\nThe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  int days = date_beg.GetDaysFrom(date_meteo);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  // Files.
  string directory_out, file_in_wind;

  configuration.SetSection("[paths]");

  configuration.PeekValue("10_m_wind_module_file", file_in_wind);

  configuration.PeekValue("Directory_sea_salt", directory_out);

  // Sea salt.
  real threshold_radius;
  string ss_parameterization;

  configuration.SetSection("[sea_salt]");

  configuration.PeekValue("Parameterization", "Smith | Monahan",
                          ss_parameterization);

  configuration.PeekValue("Threshold_radius", "positive", threshold_radius);
  configuration.PeekValue("Delta_t", "> 0", Delta_t_sea_salt);
  Nt_sea_salt = int(real(Nt) * Delta_t / Delta_t_sea_salt + 0.5);

  // LUC.
  string LUC_file;
  int Nluc, sea_index;

  configuration.SetSection("[LUC]");

  configuration.PeekValue("File", LUC_file);
  configuration.PeekValue("Nb_luc", "positive", Nluc);
  configuration.PeekValue("Sea_index", ">= 0 | < " + to_str(Nluc),
                          sea_index);

  // PM.
  bool sections_computed;
  real dmin(0.01), dmax(10.);
  int Nsections, Nbounds;
  string file_sections("");

  configuration.SetSection("[PM]");

  configuration.PeekValue("Sections_computed", sections_computed);
  if (sections_computed)
    {
      configuration.PeekValue("Diameter_min", "positive", dmin);
      configuration.PeekValue("Diameter_max", "positive | > " + to_str(dmin),
                              dmax);
    }
  else
    configuration.PeekValue("File_sections", file_sections);
  configuration.PeekValue("Nsections", "> 0", Nsections);
  Nbounds = Nsections + 1;


  configuration.SetSection("[sea-salt_composition]");

  // Mass average fraction in sea_salt.
  real fraction_na;
  real fraction_cl;
  real fraction_so4;
  configuration.PeekValue("NA", "positive", fraction_na);
  configuration.PeekValue("CL", "positive", fraction_cl);
  configuration.PeekValue("SO4", "positive", fraction_so4);

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  if (!exists(file_in_wind))
    throw string("10 meter Wind Module file doesn\'t exist");

  bool extra_step = (int(file_size(file_in_wind) / (Nt * Nx * Ny
                                                    * sizeof(float)))
                     > days + 1);

  RegularGrid<real> GridT_meteo(t_min, Delta_t, extra_step ? Nt + 1 : Nt);

  // Output settings.

  // Output grids.
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridT(t_min, Delta_t_sea_salt, Nt_sea_salt);
  RegularGrid<real> GridC(Nluc);
  RegularGrid<real> Gridsection(Nsections);

  // Data.
  Data<real, 3> WindModule_10m_in(GridT_meteo, GridY, GridX);
  Data<real, 3> WindModule_10m(GridT, GridY, GridX);

  Data<real, 3> LUC(GridC, GridY, GridX);

  cout << " done." << endl << endl;


  /////////////////////////
  // METEOROLOGICAL DATA //
  /////////////////////////

  cout << "Generating sea-salt emissions from "
       << date_beg.GetDate("%y-%m-%d %h:%i") << " to "
       << date_end.GetDate("%y-%m-%d %h:%i") << "." << endl;

  cout << "Reading and interpolating meteorological data...";
  cout.flush();

  FormatBinary<real> Meteo;

  Meteo.ReadSteps(file_in_wind, step, WindModule_10m_in);
  LinearInterpolationRegular(WindModule_10m_in, WindModule_10m);

  Meteo.Read(LUC_file, LUC);

  cout << " done." << endl;


  ///////////////
  // EMISSIONS //
  ///////////////

  // Emissions are estimated according to:
  // Monahan, E. C. , Spiel, D. E. , and Davidson, K. L. ,1986,
  // A model of marine aerosol generation via whitecaps and wave disruption,
  // In Monahan, E. C. and MacNiochaill, G. , editors, Oceanic Whitecaps.
  // D. Reidel, Norwell, Mass.

  /*** Dry radius ***/

  cout << "Reading or computing dry radius...";
  cout.flush();

  // Fixed aerosol diameter bounds (µm).
  Data<real, 1> dry_diameters(Nbounds);
  // Fixed aerosol dry radius bounds (µm).
  Data<real, 1> dry_radii(Nbounds);

  // Linear step of discretization.
  real hx = log(dmax / dmin) / Nsections;

  // Diameter and radius at bounds.
  if (sections_computed)
    for (int js = 0; js < Nbounds; js++)
      dry_diameters(js) = dmin * exp(hx * js);  //(µm)

  else
    {
      FormatText Namsections;
      Namsections.Read(file_sections, dry_diameters);
    }

  for (int js = 0; js < Nbounds; js++)
    dry_radii(js) = dry_diameters(js) / 2.; //(µm)

  cout << " done." << endl;

  /*** Computing emission fluxes ***/

  cout << "Computing fluxes of sea salt aerosol for each size bin...";
  cout.flush();

  Data <real, 4> Flux_Salt(Gridsection, GridT, GridY, GridX);
  Data <real, 3> PNA(GridT, GridY, GridX);
  Data <real, 3> PCL(GridT, GridY, GridX);
  Data <real, 3> PSO4(GridT, GridY, GridX);
  Flux_Salt.SetZero();
  PCL.SetZero();
  PNA.SetZero();
  PSO4.SetZero();

  FormatBinary<real> EmissionsFormat;

  real sea_salt_density = 2.165E-6;
  real wind_10_m;
  real c80 = 2.0;

  int h, j, i, js;

  for (h = 0; h < Nt_sea_salt; h++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        {
          if (LUC(sea_index, j, i) != 0.0)
            {
              wind_10_m = WindModule_10m(h, j, i);
              if (wind_10_m < 0)
                wind_10_m = 0;
              for (js = 0; js < Nsections; js++)
                {
                  // Wet radius for RH = 80%.
                  real wet_radius80 = dry_radii(js) * c80;
                  real wet_radius80_sup = dry_radii(js + 1) * c80;
                  if (wet_radius80 > threshold_radius)
                    {
                      // Units: particles.m-2.s-1
                      real integral = SimpsonIntegral(ss_parameterization,
                                                      wet_radius80,
                                                      wet_radius80_sup,
                                                      wind_10_m);
                      // Units: µg.m-2.s-1
                      Flux_Salt(js, h, j, i) = integral
                        / pow(c80, real(3.)) * sea_salt_density * 4. / 3. * pi
                        * LUC(sea_index, j, i);
                    }
                  else if (wet_radius80_sup > threshold_radius)
                    {
                      // Units: particles.m-2.s-1
                      real integral = SimpsonIntegral(ss_parameterization,
                                                      threshold_radius,
                                                      wet_radius80_sup,
                                                      wind_10_m);
                      // Units: µg.m-2.s-1
                      Flux_Salt(js, h, j, i) = integral
                        / pow(c80, real(3.)) * sea_salt_density * 4. / 3. * pi
                        * LUC(sea_index, j, i);
                    }
                }
            }
        }

  cout << " done." << endl;


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  cout << "Writing fluxes of sea salt aerosol for each size bin...";
  cout.flush();

  for (js = 0; js < Nsections; js++)
    {
      PNA.SetZero();
      PCL.SetZero();
      PSO4.SetZero();
      for (h = 0; h < Nt_sea_salt; h++)
        for (i = 0; i < Nx; i++)
          for (j = 0; j < Ny; j++)
            {
              // PNA and PCL are in µg.m-2.s-1.
              PNA(h, j, i) = Flux_Salt(js, h, j, i) * fraction_na;
              PCL(h, j, i) = Flux_Salt(js, h, j, i) * fraction_cl;
              PSO4(h, j, i) = Flux_Salt(js, h, j, i) * fraction_so4;
            }

      if (!PNA.IsZero())
        EmissionsFormat.Append(PNA, directory_out +
                               "PNA_" + to_str(js) + ".bin");
      if (!PCL.IsZero())
        EmissionsFormat.Append(PCL, directory_out +
                               "PHCL_" + to_str(js) + ".bin");
      if (!PSO4.IsZero())
        EmissionsFormat.Append(PSO4, directory_out +
                               "PSO4_" + to_str(js) + ".bin");

    }

  cout << " done." << endl;
  cout << endl;

  END;

  return 0;

}
