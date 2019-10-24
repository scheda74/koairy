// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet and Hervé Njomgang
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


// This program generates LUC data for Polair3D based on USGS data.


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


int main(int argc, char **argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("");

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 default_name);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Units: SI.
  real lon_origin, lat_origin, lon_origin_af, lat_origin_af;
  real lon_upper_left, lat_upper_left, lon_upper_left_af, lat_upper_left_af;

  // Pixel size (LUC, USGS): 1000 meters.
  real step;

  int i, j, k, l;

  // Input files.
  string dir_in, file_in, file_in_af;

  // Output files.
  string dir_out, file_out;

  // Input dimensions.
  int Nx_in_LUC, Ny_in_LUC, Nx_in_af_LUC, Ny_in_af_LUC;

  // Number of LUC categories.
  int Nc;

  // Output dimensions.
  int Nx_out, Ny_out;
  real Delta_y_out, Delta_x_out, x_min_out, y_min_out;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  config.SetSection("[paths]");

  config.PeekValue("Database_luc-usgs", dir_in);
  config.PeekValue("LUC_in_ea", file_in);
  config.PeekValue("LUC_in_af", file_in_af);
  config.PeekValue("Directory_luc-usgs", dir_out);
  config.PeekValue("LUC_out", file_out);

  config.SetSection("[USGS]");

  config.PeekValue("lon_origin_ea", lon_origin);
  config.PeekValue("lat_origin_ea", lat_origin);
  config.PeekValue("lon_origin_af", lon_origin_af);
  config.PeekValue("lat_origin_af", lat_origin_af);
  config.PeekValue("lon_upper_left_ea", lon_upper_left);
  config.PeekValue("lat_upper_left_ea", lat_upper_left);
  config.PeekValue("lon_upper_left_af", lon_upper_left_af);
  config.PeekValue("lat_upper_left_af", lat_upper_left_af);
  config.PeekValue("Step", "> 0", step);
  config.PeekValue("Nx_ea", "> 0", Nx_in_LUC);
  config.PeekValue("Ny_ea", "> 0", Ny_in_LUC);
  config.PeekValue("Nx_af", "> 0", Nx_in_af_LUC);
  config.PeekValue("Ny_af", "> 0", Ny_in_af_LUC);
  config.PeekValue("Nc", "> 0", Nc);

  int sea_index;
  config.PeekValue("Sea_index", ">= 0 | < " + to_str(Nc), sea_index);

  config.SetSection("[domain]");

  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Delta_x", "> 0", Delta_x_out);
  config.PeekValue("Delta_y", "> 0", Delta_y_out);
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("y_min", y_min_out);

  cout << " done." << endl;


  /////////////////
  // ALLOCATIONS //
  /////////////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  // Input grids.
  RegularGrid<int> GridX_in_LUC(Nx_in_LUC), GridY_in_LUC(Ny_in_LUC),
    GridX_in_af_LUC(Nx_in_af_LUC), GridY_in_af_LUC(Ny_in_af_LUC);
  // Input LUC data.
  Data<char, 2, int> LUC_in(GridY_in_LUC, GridX_in_LUC),
    LUC_in_af(GridY_in_af_LUC, GridX_in_af_LUC);

  // Output settings.

  // Output grids.
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);
  RegularGrid<real> GridC(Nc);
  // Output data.
  Data<float, 3, real> LUC_out(GridC, GridY_out, GridX_out);
  Data<float, 3, real> LUC_out_af(GridC, GridY_out, GridX_out);
  Data<float, 2, real> NbCells(GridY_out, GridX_out);
  Data<float, 2, real> NbCells_af(GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////
  // LUC //
  /////////

  cout << "Reading LUC data...";
  cout.flush();
  FormatBinary<char> FormatIn, FormatIn_af;
  FormatIn.Read(dir_in + "/" + file_in, LUC_in);
  FormatIn_af.Read(dir_in + "/" + file_in_af, LUC_in_af);
  cout << " done." << endl;

  // Convertion from Lambert azimuthal equal area to latitude/longitude.
  LaeaToLonlat<real> Convertion(lon_origin, lat_origin);
  LaeaToLonlat<real> Convertion_af(lon_origin_af, lat_origin_af);

  int ik, il;

  cout << "Building LUC data on output grid...";
  cout.flush();

  LUC_out.SetZero();
  NbCells.SetZero();

  real xmin_usgs, xmax_usgs, ymin_usgs, ymax_usgs;
  real lamb_xmin_usgs, lamb_xmax_usgs, lamb_ymin_usgs, lamb_ymax_usgs;
  real xx, yy, lon_usgs, lat_usgs;
  int inx, iny;
  real gridx_min_model, gridy_min_model;
  real gridx_max_model, gridy_max_model;
  real myzero = 0.0;
  real Norm_x, Norm_y;
  real increment;

  for (j = 0; j < Ny_in_LUC; j++)
    for (i = 0; i < Nx_in_LUC; i++)
      if (int(LUC_in(j, i)) != 0)
        {
          xx = i * step + lon_upper_left;
          yy = -j * step + lat_upper_left;
          Convertion(xx, yy, lon_usgs, lat_usgs);
          inx = -1;
          if ((lon_usgs >=  x_min_out - Delta_x_out / 2.0) &&
              (lon_usgs <=  x_min_out + (Nx_out - 0.5) * Delta_x_out))
            inx = int((lon_usgs - x_min_out + Delta_x_out / 2.0) /
                      Delta_x_out);
          iny = -1;
          if ((lat_usgs >=  y_min_out - Delta_y_out / 2.0) &&
              (lat_usgs <=  y_min_out + (Ny_out - 0.5) * Delta_y_out))
            iny = int((lat_usgs - y_min_out + Delta_y_out / 2.0) /
                      Delta_y_out);

          if (inx >= 0 && inx < Nx_out && iny >= 0 && iny < Ny_out)
            {
              lamb_xmin_usgs = (i - 0.5) * step + lon_upper_left;
              lamb_xmax_usgs = (i + 0.5) * step + lon_upper_left;
              lamb_ymax_usgs = -(j - 0.5) * step + lat_upper_left;
              lamb_ymin_usgs = -(j + 0.5) * step + lat_upper_left;
              Convertion(lamb_xmin_usgs, lamb_ymin_usgs, xmin_usgs,
                         ymin_usgs);
              Convertion(lamb_xmax_usgs, lamb_ymax_usgs, xmax_usgs,
                         ymax_usgs);
              gridx_min_model =  x_min_out + (inx - 0.5) * Delta_x_out;
              gridy_min_model =  y_min_out + (iny - 0.5) * Delta_y_out;
              gridx_max_model =  x_min_out + (inx + 0.5) * Delta_x_out;
              gridy_max_model =  y_min_out + (iny + 0.5) * Delta_y_out;
              Norm_x = xmax_usgs - xmin_usgs;
              Norm_y = ymax_usgs - ymin_usgs;
              increment = (min(xmax_usgs, gridx_max_model)
                           - max(gridx_min_model, xmin_usgs))
                * (min(ymax_usgs, gridy_max_model)
                   - max(gridy_min_model, ymin_usgs)) / Norm_x / Norm_y;
              LUC_out(int(LUC_in(j, i)) - 1, iny, inx) += increment;
              NbCells(iny, inx) += increment;
              if ((iny != 0) && (inx != Nx_out - 1))
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out(int(LUC_in(j, i)) - 1, iny - 1, inx + 1)
                    += increment;
                  NbCells(iny - 1, inx + 1) += increment;
                }
              if (inx != Nx_out - 1)
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * (min(ymax_usgs, gridy_max_model)
                       - max(gridy_min_model, ymin_usgs))
                    / Norm_x / Norm_y;
                  LUC_out(int(LUC_in(j, i)) - 1, iny, inx + 1) += increment;
                  NbCells(iny, inx + 1) += increment;
                }
              if ((iny != Ny_out - 1) && (inx != Nx_out - 1))
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y;
                  LUC_out(int(LUC_in(j, i)) - 1, iny + 1, inx + 1)
                    += increment;
                  NbCells(iny + 1, inx + 1) += increment;
                }
              if (inx != 0)
                {
                  increment = max(gridx_min_model - xmin_usgs, myzero)
                    * (min(ymax_usgs, gridy_max_model)
                       - max(gridy_min_model, ymin_usgs))
                    / Norm_x / Norm_y;
                  LUC_out(int(LUC_in(j, i)) - 1, iny, inx - 1) += increment;
                  NbCells(iny, inx - 1) += increment;
                }
              if ((iny != 0) && (inx != 0))
                {
                  increment =
                    max(gridx_min_model - xmin_usgs, myzero)
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out(int(LUC_in(j, i)) - 1, iny - 1, inx - 1)
                    += increment;
                  NbCells(iny - 1, inx - 1) += increment;
                }
              if ((iny != Ny_out - 1) && (inx != 0))
                {
                  increment = max(gridx_min_model - xmin_usgs, myzero)
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out(int(LUC_in(j, i)) - 1, iny + 1, inx - 1)
                    += increment;
                  NbCells(iny + 1, inx - 1) += increment;
                }
              if (iny != Ny_out - 1)
                {
                  increment = (min(xmax_usgs, gridx_max_model)
                               - max(gridx_min_model, xmin_usgs))
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out(int(LUC_in(j, i)) - 1, iny + 1, inx) += increment;
                  NbCells(iny + 1, inx) += increment;
                }
              if (iny != 0)
                {
                  increment = (min(xmax_usgs, gridx_max_model)
                               - max(gridx_min_model, xmin_usgs))
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out(int(LUC_in(j, i)) - 1, iny - 1, inx) += increment;
                  NbCells(iny - 1, inx) += increment;
                }
            }
        }

  for (j = 0; j < Ny_out; j++)
    for (i = 0; i < Nx_out; i++)
      if (NbCells(j, i) != 0.)
        for (k = 0; k < Nc; k++)
          LUC_out(k, j, i) /= NbCells(j, i);

  LUC_out_af.SetZero();
  NbCells_af.SetZero();

  for (j = 0; j < Ny_in_af_LUC; j++)
    for (i = 0; i < Nx_in_af_LUC; i++)
      if (int(LUC_in_af(j, i)) != 0 && int(LUC_in_af(j, i)) < 25)
        {
          xx = i * step + lon_upper_left_af;
          yy = -j * step + lat_upper_left_af;
          Convertion_af(xx, yy, lon_usgs, lat_usgs);
          inx = -1;
          if ((lon_usgs >=  x_min_out - Delta_x_out / 2.0) &&
              (lon_usgs <=  x_min_out + (Nx_out - 0.5) * Delta_x_out))
            inx = int((lon_usgs - x_min_out + Delta_x_out / 2.0)
                      / Delta_x_out);
          iny = -1;
          if ((lat_usgs >=  y_min_out - Delta_y_out / 2.0) &&
              (lat_usgs <=  y_min_out + (Ny_out - 0.5) * Delta_y_out))
            iny = int((lat_usgs - y_min_out + Delta_y_out / 2.0)
                      / Delta_y_out);

          if (inx >= 0 && inx < Nx_out && iny >= 0 && iny < Ny_out)
            {
              lamb_xmin_usgs = (i - 0.5) * step + lon_upper_left_af;
              lamb_xmax_usgs = (i + 0.5) * step + lon_upper_left_af;
              lamb_ymax_usgs = -(j - 0.5) * step + lat_upper_left_af;
              lamb_ymin_usgs = -(j + 0.5) * step + lat_upper_left_af;
              Convertion_af(lamb_xmin_usgs, lamb_ymin_usgs, xmin_usgs,
                            ymin_usgs);
              Convertion_af(lamb_xmax_usgs, lamb_ymax_usgs, xmax_usgs,
                            ymax_usgs);
              gridx_min_model =  x_min_out + (inx - 0.5) * Delta_x_out;
              gridy_min_model =  y_min_out + (iny - 0.5) * Delta_y_out;
              gridx_max_model =  x_min_out + (inx + 0.5) * Delta_x_out;
              gridy_max_model =  y_min_out + (iny + 0.5) * Delta_y_out;
              Norm_x = xmax_usgs - xmin_usgs;
              Norm_y = ymax_usgs - ymin_usgs;
              increment = (min(xmax_usgs, gridx_max_model)
                           - max(gridx_min_model, xmin_usgs))
                * (min(ymax_usgs, gridy_max_model)
                   - max(gridy_min_model, ymin_usgs)) / Norm_x / Norm_y ;
              LUC_out_af(int(LUC_in_af(j, i)) - 1, iny, inx) += increment;
              NbCells_af(iny, inx) += increment;
              if ((iny != 0) && (inx != Nx_out - 1))
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny - 1, inx + 1)
                    += increment;
                  NbCells_af(iny - 1, inx + 1) += increment;
                }
              if (inx != Nx_out - 1)
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * (min(ymax_usgs, gridy_max_model)
                       - max(gridy_min_model, ymin_usgs))
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny, inx + 1)
                    += increment;
                  NbCells_af(iny, inx + 1) += increment;
                }
              if ((iny != Ny_out - 1) && (inx != Nx_out - 1))
                {
                  increment = max(xmax_usgs - gridx_max_model, myzero)
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny + 1, inx + 1)
                    += increment;
                  NbCells_af(iny + 1, inx + 1) += increment;
                }
              if (inx != 0)
                {
                  increment = max(gridx_min_model - xmin_usgs, myzero)
                    * (min(ymax_usgs, gridy_max_model)
                       - max(gridy_min_model, ymin_usgs))
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny, inx - 1)
                    += increment;
                  NbCells_af(iny, inx - 1) += increment;
                }
              if ((iny != 0) && (inx != 0))
                {
                  increment = max(gridx_min_model - xmin_usgs, myzero)
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny - 1, inx - 1)
                    += increment;
                  NbCells_af(iny - 1, inx - 1) += increment;
                }
              if ((iny != Ny_out - 1) && (inx != 0))
                {
                  increment = max(gridx_min_model - xmin_usgs, myzero)
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny + 1, inx - 1)
                    += increment;
                  NbCells_af(iny + 1, inx - 1) += increment;
                }
              if (iny != Ny_out - 1)
                {
                  increment = (min(xmax_usgs, gridx_max_model)
                               - max(gridx_min_model, xmin_usgs))
                    * max(ymax_usgs - gridy_max_model, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1, iny + 1, inx)
                    += increment;
                  NbCells_af(iny + 1, inx) += increment;
                }
              if (iny != 0)
                {
                  increment = (min(xmax_usgs, gridx_max_model)
                               - max(gridx_min_model, xmin_usgs))
                    * max(gridy_min_model - ymin_usgs, myzero)
                    / Norm_x / Norm_y ;
                  LUC_out_af(int(LUC_in_af(j, i)) - 1,
                             iny - 1, inx) += increment;
                  NbCells_af(iny - 1, inx) += increment;
                }
            }
        }

  for (j = 0; j < Ny_out; j++)
    for (i = 0; i < Nx_out; i++)
      if (NbCells_af(j, i) != 0.)
        for (k = 0; k < Nc; k++)
          LUC_out_af(k, j, i) /= NbCells_af(j, i);

  // Merges the two LUC arrays.
  for (j = 0; j < Ny_out; j++)
    for (i = 0; i < Nx_out; i++)
      if ((NbCells(j, i) == 0. || LUC_out(sea_index, j, i) == 1.)
          && NbCells_af(j, i) != 0.)
        for (k = 0; k < Nc; k++)
          LUC_out(k, j, i) = LUC_out_af(k, j, i);

  cout << " done." << endl;
  cout << endl;


  ////////////////////
  // WRITING OUTPUT //
  ////////////////////

  cout << "Writing output data...";
  cout.flush();

  FormatBinary<float> PolairFormat;
  PolairFormat.Write(LUC_out, dir_out + "/" + file_out);

  cout << " done." << endl;


  END;

  return 0;

}
