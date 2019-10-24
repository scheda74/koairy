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


// This program converts LUC data from one classification to another.


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

  int i, j, k, l;

  // Input files.
  string dir_in, dir_out, file_in, file_out;

  // Input dimensions.
  // Number of categories.
  int Nc_in, Nc_out;

  // Output dimensions.
  int Nx, Ny;
  real Delta_y, Delta_x, y_min, x_min;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  config.SetSection("[domain]");

  config.PeekValue("Nx", "> 0", Nx);
  config.PeekValue("Ny", "> 0", Ny);
  config.PeekValue("Delta_x", "> 0", Delta_x);
  config.PeekValue("Delta_y", "> 0", Delta_y);
  config.PeekValue("x_min", x_min);
  config.PeekValue("y_min", y_min);

  config.SetSection("[paths]");

  config.PeekValue("Database_luc-convert", dir_in);
  config.PeekValue("File_in", file_in);
  config.PeekValue("Directory_luc-convert", dir_out);
  config.PeekValue("File_out", file_out);

  config.SetSection("[dimensions]");

  config.PeekValue("Nc_in", "> 0", Nc_in);
  config.PeekValue("Nc_out", "> 0",  Nc_out);

  cout << " done." << endl;


  /////////////////
  // ALLOCATIONS //
  /////////////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridC_in(Nc_in);
  RegularGrid<real> GridC_out(Nc_out);

  // Input data.
  Data<real, 3, real> LUC_in(GridC_in, GridY, GridX);

  // Output data.
  Data<real, 3, real> LUC_out(GridC_out, GridY, GridX);
  Data<real, 2> Coefficients(Nc_in, Nc_out);
  Data<int, 1> Categories(Nc_in);

  cout << " done." << endl;


  /////////
  // LUC //
  /////////

  cout << "Reading LUC input data...";
  cout.flush();

  FormatBinary<real> InputPolair;
  InputPolair.Read(dir_in + "/" + file_in, LUC_in);

  FormatFormattedText InputFile(string("<e><e ") + to_str<int>(Nc_out) + ">");

  config.Rewind();
  config.SetSection("[coefficients]");
  InputFile.Read(**config.GetCurrent(), "0", Categories);
  config.Rewind();
  config.SetSection("[coefficients]");
  InputFile.Read(**config.GetCurrent(), "1", Coefficients);

  cout << " done." << endl;


  ///////////////////
  // BUILDING DATA //
  ///////////////////

  cout << "Building LUC data on output grid...";
  cout.flush();

  LUC_out.SetZero();
  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      for (k = 0; k < Nc_in; k++)
        for (l = 0; l < Nc_out; l++)
          if (Coefficients(k, l) != 0)
            LUC_out(l, j, i) += Coefficients(k, l) * LUC_in(k, j, i);

  cout << " done." << endl;


  ////////////////////
  // WRITING OUTPUT //
  ////////////////////

  cout << "Writing LUC output data...";
  cout.flush();

  FormatBinary<real> PolairFormat;
  PolairFormat.Write(LUC_out, dir_out + "/" + file_out);

  cout << " done." << endl;

  END;

  return 0;

}
