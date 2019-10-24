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


// This program generates LUC data for Polair3D based on GLCF data.


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
  real lon_origin, lat_origin;

  // Pixel size (LUC, GLCF): 1000 meters.
  real step;

  int i, j, k, l, c;

  // Input files.
  string dir_in, file_in;

  // Output files.
  string dir_out, file_out;

  // Input dimensions.
  int Nx_in_LUC, Ny_in_LUC;

  // Number of LUC categories.
  int Nc;

  // Starting index of LUCs.
  int shift;

  // Output dimensions.
  int Nx_out, Ny_out;
  real Delta_x_out, Delta_y_out, x_min_out, y_min_out;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  config.SetSection("[paths]");

  config.PeekValue("Database_luc-glcf", dir_in);
  config.PeekValue("LUC_in", file_in);
  config.PeekValue("Directory_luc-glcf", dir_out);
  config.PeekValue("LUC_out", file_out);

  config.SetSection("[GLCF]");

  config.PeekValue("Step", "> 0", step);
  config.PeekValue("x_min", lon_origin);
  config.PeekValue("y_min", lat_origin);
  config.PeekValue("Nx", "> 0", Nx_in_LUC);
  config.PeekValue("Ny", "> 0", Ny_in_LUC);
  config.PeekValue("Nc", "> 0", Nc);
  config.PeekValue("Shift", shift);

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
  RegularGrid<real> GridX_in(Nx_in_LUC), GridY_in(Ny_in_LUC);
  RegularGrid<real> GridC(Nc);
  // Output grids.
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out),
    GridY_out(y_min_out, Delta_y_out, Ny_out);

  // Input data.
  Data<unsigned char, 2, real> LUC_in(GridY_in, GridX_in);
  // Output data.
  Data<float, 3, real> LUC_out(GridC, GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////
  // LUC //
  /////////

  cout << "Reading LUC data...";
  cout.flush();
  FormatBinary<unsigned char> FormatIn;
  FormatIn.Read(dir_in + "/" + file_in, LUC_in);
  cout << " done." << endl;


  cout << "Building LUC data on output grid...";
  cout.flush();

  LUC_out.SetZero();

  int pxl_min, pxl_max, line_min, line_max, NbCells;

  for (j = 0 ; j < Ny_out ; j++)
    for (i = 0 ; i < Nx_out ; i++)
      {
        // limits of the output cell for the input data
        pxl_min = int((x_min_out + (i - 0.5) * Delta_x_out
                       - lon_origin - 0.5 * step)
                      / step);
        pxl_max = int((x_min_out + (i + 0.5) * Delta_x_out
                       - lon_origin - 0.5 * step)
                      / step);
        line_min = int((- y_min_out - (j + 0.5) * Delta_y_out
                        + lat_origin - 0.5 * step)
                       / step);
        line_max = int((- y_min_out - (j - 0.5) * Delta_y_out
                        + lat_origin - 0.5 * step)
                       / step);

        NbCells = (line_max - line_min) * (pxl_max - pxl_min);

        if (pxl_max >= Nx_in_LUC || pxl_min < 0
            || line_max >= Ny_in_LUC || line_min < 0)
          {
            cout << endl;
            cout << "(" << i << ", " << j << "): ";
            cout << "This cell exceeds the limits of the input file." << endl;
            return 1;
          }

        for (l = line_min ; l < line_max ; l++)
          for (k = pxl_min ; k < pxl_max ; k++)
            LUC_out(int(LUC_in(l, k) - shift), j, i)++;

        for (c = 0 ; c < Nc ; c++)
          LUC_out(c, j, i) /= real(NbCells);
      }

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
