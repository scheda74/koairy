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


// This program computes the roughness on the output grid on the basis of
// roughness data given by land use coverage categories.


//////////////
// INCLUDES //

#include <vector>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData ;

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

  int i, j, k;

  // Input domain.
  int Ny, Nx, Nc;

  // Files.
  string dir_out, file_out, LUC_file, roughness_data_file;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  config.SetSection("[domain]");

  config.PeekValue("Ny", "> 0", Ny);
  config.PeekValue("Nx", "> 0", Nx);

  config.SetSection("[paths]");

  config.PeekValue("LUC_file", LUC_file);
  Nc = int(file_size(LUC_file)) / sizeof(real) / (Ny * Nx);

  config.PeekValue("Directory_roughness", dir_out);
  config.PeekValue("Roughness_out", file_out);

  config.SetSection("[data]");

  config.PeekValue("Roughness_data_file", roughness_data_file);

  cout << " done." << endl;


  ///////////////
  // ROUGHNESS //
  ///////////////

  cout << "Reading roughness data...";
  cout.flush();
  FormatFormattedText RoughnessData("<e><e><a>");
  Data<real, 1> roughness_data(Nc);
  RoughnessData.Read(roughness_data_file, "1", roughness_data);
  cout << " done." << endl;

  FormatBinary<real> LucFormat ;

  Data<real, 3> LUC(Nc, Ny, Nx);
  Data<real, 2> RoughnessValue(Ny, Nx);

  RoughnessValue.Fill(0.0);
  LucFormat.Read(LUC_file, LUC);

  for (k = 0; k < Nc; k++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        RoughnessValue(j, i) += LUC(k, j, i) * roughness_data(k);


  cout << "Writing roughness binary ...";
  cout.flush();

  FormatBinary<real> RoughnessFormat;
  RoughnessFormat.Write(RoughnessValue, dir_out + "/" + file_out);

  cout << " done." << endl;


  END;

  return 0;

}
