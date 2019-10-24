// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi
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

// This program extracts land use data from GLCF global data. It
// generates land use description over the EMEP domain, which is
// required to process EMEP emissions.


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

  string configuration_file, sec_config_file;

  parse_argument(argc, argv, configuration_file, sec_config_file, "");


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef float real;

  // Files used.
  string file_in, file_out;

  // Input dimensions.
  int Nx_in, Ny_in;
  real x_min_in, y_max_in, step;

  // Output dimensions.
  int Nx_out,  Ny_out;
  real x_min_out,  y_min_out;

  int i_min, j_min, i, j;


  ////////////////////////
  // CONFIGURATION FILE //
  ////////////////////////


  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  config.SetSection("[paths]");
  config.PeekValue("File_GLCF", file_in);
  config.PeekValue("File_out", file_out);

  config.SetSection("[GLCF]");
  config.PeekValue("x_min", x_min_in);
  config.PeekValue("Nx", "> 0", Nx_in);
  config.PeekValue("y_max", y_max_in);
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Step", "> 0", step);

  config.SetSection("[subdomain]");
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("y_min", y_min_out);
  config.PeekValue("Ny", "> 0", Ny_out);

  cout << " done." << endl;


  /////////////////
  // ALLOCATIONS //
  /////////////////


  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input grids.
  RegularGrid<real> GridX_in(Nx_in), GridY_in(Ny_in);

  // Output grids.
  RegularGrid<real> GridX_out(Nx_out), GridY_out(Ny_out);

  // Input data.
  Data<unsigned char, 2, real> LUC_in(GridY_in, GridX_in);
  // Output data.
  Data<int, 2, real> LUC_out(GridY_out, GridX_out);

  cout << " done." << endl;


  /////////
  // LUC //
  /////////


  cout << "Reading LUC data...";
  cout.flush();

  FormatBinary<unsigned char> FormatIn;
  FormatIn.Read(file_in, LUC_in);
  LUC_in.ReverseData(0);

  cout << " done." << endl;

  cout << "Extraction...";
  cout.flush();

  i_min = int((x_min_out - x_min_in) / step);
  j_min = int((y_min_out - (y_max_in - Ny_in * step)) / step);

  for (j = 0; j < Ny_out; j++)
    for (i = 0; i < Nx_out; i++)
      LUC_out(j, i) = int(LUC_in(j + j_min, i + i_min));

  cout << " done." << endl;


  ////////////////////
  // WRITING OUTPUT //
  ////////////////////


  cout << "Writing output data...";
  cout.flush();

  FormatBinary<int> PolairFormat;
  PolairFormat.Write(LUC_out, file_out);

  cout << " done." << endl;


  END;

  return 0;

}
