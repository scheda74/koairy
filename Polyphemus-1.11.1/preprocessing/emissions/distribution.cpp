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


// This program generates the vertical distribution of emissions.

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

int main(int argc, char* argv[])
{

  TRY;

  cout << endl;

  typedef float real;

  string configuration_file, sec_config_file, default_name("");

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 default_name);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  // Configuration file.
  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  /*** Input ***/

  int Nz_in, Nz;
  string levels_file_in, levels_file, input_distribution;

  configuration.SetSection("[domain]");
  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Vertical_levels", levels_file);

  configuration.SetSection("[EMEP]");
  configuration.PeekValue("Nz_in", "> 0", Nz_in);
  configuration.PeekValue("Vertical_levels", levels_file_in);
  configuration.PeekValue("Vertical_distribution", input_distribution);

  string output_file;
  configuration.PeekValue("Polair_vertical_distribution", output_file);

  FormatFormattedText levels_format(string("<e ") + to_str(Nz_in + 1) + ">");
  RegularGrid<real> GridZ_in(Nz_in);
  RegularGrid<real> GridZ_interf_in(Nz_in + 1);
  levels_format.Read(levels_file_in, "0", GridZ_interf_in);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridZ_interf(Nz + 1);
  levels_format.SetFormat(string("<e ") + to_str(Nz + 1) + ">");
  levels_format.Read(levels_file, "0", GridZ_interf);

  for (int k = 0; k < Nz_in; k++)
    GridZ_in(k) = (GridZ_interf_in(k) + GridZ_interf_in(k + 1)) / 2.0;
  for (int k = 0; k < Nz; k++)
    GridZ(k) = (GridZ_interf(k) + GridZ_interf(k + 1)) / 2.0;

  /*** EMEP ***/

  int Nsectors(10);

  RegularGrid<real> GridSectors(Nsectors);


  //////////
  // DATA //
  //////////

  Data<real, 1> Ground_part(GridSectors);
  Data<real, 2> vertical_distribution_in(GridSectors, GridZ_in);
  Data<real, 2> vertical_distribution_out(GridSectors, GridZ);


  ///////////////////
  // COMPUTES DATA //
  ///////////////////

  cout << "Computing vertical distribution... ";
  cout.flush();
  // Reads first initial distribution.
  FormatFormattedText distribution_format(string("<e><e ") + to_str(Nz_in) + ">");
  distribution_format.Read(input_distribution, "1", vertical_distribution_in);
  distribution_format.Read(input_distribution, "0", Ground_part);

  ComputeVerticalDistribution(vertical_distribution_in,
                              GridZ_interf_in, GridZ_interf,
                              vertical_distribution_out);
  cout << " done." << endl;


  /////////////////
  // WRITES DATA //
  /////////////////

  cout << "Writing vertical distribution...";
  cout.flush();
  ofstream output_stream(output_file.c_str());
  output_stream.precision(8);
  for (int s = 0; s < Nsectors; s++)
    {
      output_stream << scientific << Ground_part(s);
      for (int k = 0; k < Nz; k++)
        output_stream << scientific << '\t' << vertical_distribution_out(s, k);
      output_stream << endl;
    }
  output_stream.close();
  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
