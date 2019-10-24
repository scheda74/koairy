// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
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


//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
#include <algorithm>
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

  string configuration_file, sec_config_file, default_name("meteo.cfg");
  Date date_beg;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, default_name);

  Date date_end(date_beg);
  date_end.AddDays(1);

  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  int h, i, j, k;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  // Domain.
  int Nt, Nz, Ny, Nx;
  real Delta_t, Delta_y, Delta_x;
  real t_min, y_min, x_min;
  string vertical_levels, date_str;

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Delta_t", "> 0", Delta_t);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);
  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);
  configuration.PeekValue("Vertical_levels", vertical_levels);
  configuration.PeekValue("Date", date_str);

  Date date_pre(date_str);

  double difference = date_beg.GetSecondsFrom(date_pre);
  if (difference < 0)
    throw string("\nThe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  // Paths.
  string directory_in, directory_rain_in, file_out;

  configuration.SetSection("[paths]");

  configuration.PeekValue("Directory_meteo", directory_in);
  configuration.PeekValue("File_Kz", file_out);

  // Land use.
  string LUC_file;
  int Nc, urban_index;
  configuration.PeekValue("LUC_file", LUC_file);
  if (!exists(LUC_file))
    throw "Unable to open land use cover file \"" + LUC_file + "\".";
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny * Nx);
  configuration.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                          urban_index);

  // Vertical diffusion.
  real Kz_min, Kz_min_urban, Kz_max;
  bool apply_vert;

  configuration.SetSection("[Kz]");

  configuration.PeekValue("Min", "positive", Kz_min);
  configuration.PeekValue("Min_urban", "positive", Kz_min_urban);
  configuration.PeekValue("Max", "positive", Kz_max);
  configuration.PeekValue("Apply_vert", apply_vert);

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Grids.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridC(Nc);

  // Data may be provided on in_interfaces.
  RegularGrid<real> GridZ_interf(Nz + 1);
  RegularGrid<real> GridY_interf(y_min - Delta_y / 2., Delta_y, Ny + 1);
  RegularGrid<real> GridX_interf(x_min - Delta_x / 2., Delta_x, Nx + 1);

  // Reads output altitudes.
  FormatText Heights;
  Heights.Read(vertical_levels, GridZ_interf);
  // Sets values at nodes.
  for (k = 0; k < Nz; k++)
    GridZ(k) = (GridZ_interf(k) + GridZ_interf(k + 1)) / 2.0;


  //////////
  // DATA //
  //////////

  // Input fields.

  Data<real, 3> LUC(GridC, GridY, GridX);

  Data<real, 4> Temperature(GridT, GridZ, GridY, GridX);
  Data<real, 4> PotentialTemperature(GridT, GridZ, GridY, GridX);
  Data<real, 4> Pressure(GridT, GridZ, GridY, GridX);

  Data<real, 3> ConvectiveRain(GridT, GridY, GridX);

  Data<real, 4> MeridionalWind(GridT, GridZ, GridY_interf, GridX);
  Data<real, 4> ZonalWind(GridT, GridZ, GridY, GridX_interf);
  Data<real, 4> VerticalWind(GridT, GridZ_interf, GridY, GridX);

  Data<real, 4> Kz(GridT, GridZ_interf, GridY, GridX);

  cout << " done." << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  FormatBinary<float> PolairFormat;

  cout << "Extracting data...";
  cout.flush();

  PolairFormat.Read(LUC_file, LUC);

  PolairFormat.ReadSteps(directory_in + "Temperature.bin", step, Temperature);
  PolairFormat.ReadSteps(directory_in + "Pressure.bin", step, Pressure);
  PolairFormat.ReadSteps(directory_in + "MeridionalWind.bin", step,
                         MeridionalWind);
  PolairFormat.ReadSteps(directory_in + "ZonalWind.bin", step, ZonalWind);
  PolairFormat.ReadSteps(directory_in + "ConvectiveRain.bin", step,
                         ConvectiveRain);

  cout << " done." << endl;


  /////////////////////////
  // WINDS AND DIFFUSION //
  /////////////////////////

  // Vertical diffusion (Louis formula, 1979).
  cout << "Computing Kz...";
  cout.flush();

  VerticalWind.SetZero();

  ComputePotentialTemperature(Temperature, Pressure, PotentialTemperature);

  ComputeLouisKz(ZonalWind, MeridionalWind, PotentialTemperature, Kz);

  real Kz_max_loc;
  int imax(0);

  for (h = 0; h < Kz.GetLength(0); h++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        {
          Kz_max_loc = 0.0;
          for (k = 0; k < Nz + 1; k++)
            if (Kz(h, k, j, i) >= Kz_max_loc)
              {
                Kz_max_loc = Kz(h, k, j, i);
                imax = k;
              }
          if (ConvectiveRain(h, j, i) > 1. / 6.)
            for (k = imax; k < Nz + 1; k++)
              Kz(h, k, j, i) = Kz_max_loc;
          else
            for (k = imax; k < Nz + 1; k++)
              Kz(h, k, j, i) = (1. - ConvectiveRain(h, j, i) * 6.)
                * Kz(h, k, j, i)
                + 6. * ConvectiveRain(h, j, i) * Kz_max_loc;
        }

  // Thresholds.
  real local_min;
  for (j = 0; j < Ny; j++)
    for (i = 0; i < Nx; i++)
      {
        local_min = Kz_min * (1. - LUC(urban_index, j, i))
          + Kz_min_urban * LUC(urban_index, j, i);
        if (apply_vert)
          {
            for (h = 0; h < Nt; h++)
              for (k = 0; k < Nz + 1; k++)
                if (Kz(h, k, j, i) < local_min)
                  Kz(h, k, j, i) = local_min;
          }
        else
          for (h = 0; h < Nt; h++)
            if (Kz(h, 1, j, i) < local_min)
              Kz(h, 1, j, i) = local_min;
      }

  Kz.ThresholdMax(Kz_max);

  cout << " done." << endl;


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  cout << "Writing data...";
  cout.flush();
  PolairFormat.Append(Kz, file_out);
  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
