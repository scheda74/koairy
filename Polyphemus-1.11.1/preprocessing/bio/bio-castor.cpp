// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
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


// This program is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).
// Base programs: "diagbio.f90" and dependencies.


// This program computes biogenic emissions based on Guenther et al. (1997)
// and Simpson et al. (1999).


// WARNING :
// It is assumed that 't_min' and 'Delta_t' are integers.


//////////////
// INCLUDES //

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////

int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("bio-castor.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  typedef float real;

  cout << "Reading configuration...";
  cout.flush();

  int h, j, i;

  /*** Input domain ***/

  int Nt, Ny, Nx;
  real Delta_t, Delta_y, Delta_x;
  real t_min, y_min, x_min;

  config.SetSection("[domain]");

  config.PeekValue("Ny", "> 0", Ny);
  config.PeekValue("Nx", "> 0", Nx);
  config.PeekValue("Delta_t", "> 0", Delta_t);
  config.PeekValue("Delta_y", "> 0", Delta_y);
  config.PeekValue("Delta_x", "> 0", Delta_x);
  config.PeekValue("y_min", y_min);
  config.PeekValue("x_min", x_min);

  // Dates.
  string date_meteo_str;
  config.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);
  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw string("\nThe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  /*** Files ***/

  string dir_out, surface_temperature_file, wind10m_file, convective_file,
    soil_moisture_file, attenuation_file, land_data_file;

  config.SetSection("[paths]");

  config.PeekValue("SurfaceTemperature", surface_temperature_file);
  config.PeekValue("WindModule_10m", wind10m_file);
  config.PeekValue("SoilMoisture", soil_moisture_file);
  config.PeekValue("Attenuation", attenuation_file);
  config.PeekValue("ConvectiveVelocity", convective_file);

  config.PeekValue("Land_data", land_data_file);
  config.PeekValue("Directory_bio", dir_out);

  /** Options ***/

  config.SetSection("[biogenic]");

  real minimum_wind_velocity;
  config.PeekValue("Minimum_wind_velocity", "positive",
                   minimum_wind_velocity);

  int Nterpenes, Nratios;
  vector<string> terpenes_names;
  vector<real> terpenes_ratios;

  config.Find("Terpenes");
  split(config.GetLine(), terpenes_names);
  Nterpenes = int(terpenes_names.size());

  config.Find("Terpenes_ratios");
  split(config.GetLine(), terpenes_ratios);
  Nratios = int(terpenes_ratios.size());

  if (Nterpenes != Nratios)
    throw string("Error in configuration: Terpenes") +
      " and Terpenes_ratios have not the same number of elements.";

  cout << endl;


  /////////////////////////
  // METEOROLOGICAL DATA //
  /////////////////////////


  cout << "Reading meteorological data...";
  cout.flush();

  // Grids.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);

  // Input data.
  Data<real, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<real, 3> WindModule_10m(GridT, GridY, GridX);
  Data<real, 3> ConvectiveVelocity(GridT, GridY, GridX);
  Data<real, 3> Attenuation(GridT, GridY, GridX);
  Data<real, 3> SoilMoisture(GridT, GridY, GridX);
  Data<real, 3> BioFactor(Ny, Nx, 5);

  // Output data.
  Data<real, 3> Isoprene(GridT, GridY, GridX);
  Data<real, 3> Terpenes(GridT, GridY, GridX);
  Data<real, 3> NO(GridT, GridY, GridX);

  // Reading data.
  FormatText().Read(land_data_file, BioFactor);

  FormatBinary<float> Meteo;
  Meteo.ReadSteps(surface_temperature_file, step, SurfaceTemperature);
  Meteo.ReadSteps(wind10m_file, step, WindModule_10m);
  Meteo.ReadSteps(convective_file, step, ConvectiveVelocity);
  Meteo.ReadSteps(soil_moisture_file, step, SoilMoisture);
  Meteo.ReadSteps(attenuation_file, step, Attenuation);

  cout << " done." << endl;


  ///////////////
  // EMISSIONS //
  ///////////////


  cout << "Computing biogenic emissions...";
  cout.flush();

  // Temporary variables used to compute the zenith angle.
  real sol_time, hour, loc_time, dec, djul;
  real zenith_angle;
  // Cosine of the zenith angle.
  real czenith_angle;
  // Solstice date.
  Date sol_date;
  sol_date.SetYear(date_beg.GetYear());
  sol_date.SetMonth(12);
  sol_date.SetDay(21);
  sol_date.SetHour(12);

  // Current date in the time loop.
  Date current_date(date_beg);

  // Temporary physical variables.
  real w10, wstar, soim, atte, cdnms, vustas, ustas, temp, ct;
  real tsoil, terfac, cr, r0, r1, r2, rt;

  const real pi = 3.14159265358979323846264;
  const real c1 = 95000.0 / 8.314;
  const real c2 = 230000.0 / 8.314;
  const real c3 = 0.09;
  const real ts = 303.0;
  const real tm = 314.0;

  for (h = 0; h < Nt; h++)
    {
      sol_time = current_date.GetSecondsFrom(sol_date) / 3600.;
      hour = real(current_date.GetHour());
      djul = real(current_date.GetNumberOfDays());

      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            // Zenith angle.
            dec = -0.409105176 * cos(7.168657935e-04 * sol_time);
            loc_time = hour - 12.
              + 0.066666666666667 * (x_min + Delta_x * real(i));
            loc_time *= 0.261799387;
            zenith_angle =
              sin((y_min + Delta_y * real(j)) / 180. * pi) * sin(dec)
              + cos((y_min + Delta_y * real(j)) / 180. * pi)
              * cos(dec) * cos(loc_time);
            czenith_angle = max(0.01, zenith_angle);

            w10 = WindModule_10m(h, j, i);
            wstar = ConvectiveVelocity(h, j, i);
            cdnms = .4 / log(10. / 5.e-4);
            vustas = max(w10, minimum_wind_velocity);
            ustas = cdnms * sqrt(vustas * vustas + 1.44 * wstar * wstar);
            temp = SurfaceTemperature(h, j, i);
            soim = SoilMoisture(h, j, i);
            atte = Attenuation(h, j, i);

            tsoil  = 1.03 * (temp - 273.15) + 2.9;
            terfac = exp(c3 * (temp - ts));
            ct = exp(c1 * (1. / ts - 1. / temp))
              / (1. + exp(c2 * (temp - tm) / (ts * temp)));
            if (czenith_angle <= 1.e-2)
              cr = 0.0;
            else
              {
                r0 = 6.89634 * atte;
                r1 = r0 * exp(-0.16 / czenith_angle) * czenith_angle;
                r2 = 0.4 * (r0 - r1) * czenith_angle;
                rt = r1 + r2;
                cr = rt / sqrt(1. + rt * rt);
              }

            if (current_date.GetMonth() >= 11 || current_date.GetMonth() <= 3)
              Isoprene(h, j, i) = 1.e-5;
            else
              Isoprene(h, j, i) = cr * ct * BioFactor(j, i, 0);

            Terpenes(h, j, i) = terfac * BioFactor(j, i, 2)
              + cr * ct * BioFactor(j, i, 1);

            if (current_date.GetMonth() >= 5 && current_date.GetMonth() <= 7)
              NO(h, j, i) = exp(0.071 * tsoil) * BioFactor(j, i, 3);
            else
              NO(h, j, i) = exp(0.071 * tsoil) * BioFactor(j, i, 4);
          }

      current_date.AddHours(int(Delta_t));
    } // Time loop.

  cout << " done." << endl;


  ////////////
  // OUTPUT //
  ////////////


  cout << "Writing output emissions...";
  cout.flush();

  FormatBinary<float> EmissionsFormat;

  EmissionsFormat.Append(Isoprene, dir_out + "Isoprene.bin");
  EmissionsFormat.Append(NO, dir_out + "NO.bin");

  // Terpenes speciation.
  for (int i = 0; i < Nterpenes; i++)
    if (terpenes_ratios[i] > 0.)
      {
        Data<real, 3> Terpenes_tmp(Terpenes);
        Terpenes_tmp.Mlt(terpenes_ratios[i]);
        EmissionsFormat.Append(Terpenes_tmp,
                               dir_out + terpenes_names[i] + ".bin");
      }

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
