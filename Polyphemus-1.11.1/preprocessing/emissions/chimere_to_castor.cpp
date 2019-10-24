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

// This program converts Chimere emission files into files in a format
// suitable for Castor (and Polair3D).

#include <iostream>
#include <vector>
#include <list>
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

  string configuration_file, sec_config_file,
    default_name("chimere_to_castor.cfg");
  Date date;

  parse_argument(argc, argv, configuration_file, sec_config_file, date,
                 default_name);


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  int h, i, j, k, s;

  cout << "Reading configuration...";
  cout.flush();

  // Configuration file.
  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  /*** Description ***/

  int Nt = 24;
  int Nx, Ny, Nz;
  string month;

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx);
  configuration.PeekValue("Ny", "> 0", Ny);

  configuration.SetSection("[description]");

  // Number of emission levels.
  configuration.PeekValue("Nz", "> 0", Nz);

  month = to_str(date.GetMonth());

  /*** Path ***/

  string file_in, directory_out;

  configuration.SetSection("[path]");

  configuration.PeekValue("File_in", file_in);
  file_in += fill(month, 2, '0', ostringstream::right);
  configuration.PeekValue("Directory_out", directory_out);

  /*** Species ***/

  int Ns;
  vector<string> species_name;

  configuration.SetSection("[species]");

  // Warning: species must be in the same order as in Chimere file.
  while (!configuration.IsEmpty())
    species_name.push_back(configuration.GetElement());

  Ns = int(species_name.size());

  cout << " done." << endl;


  ////////////////
  // READS DATA //
  ////////////////


  cout << "Reading input emissions...";
  cout.flush();

  // Input emissions. The second dimension (length 3) is for the type of day
  // (weekday, Saturday and Sunday). The last dimension is for hours.
  Data<real, 6> Emission_in(Ns, 3, Ny, Nx, Nz, 24);

  FormatText().Read(file_in, Emission_in);

  cout << " done." << endl;


  ////////////////
  // CONVERSION //
  ////////////////


  cout << "Converting to Castor and Polair3D emissions..." << endl;

  // Output emissions for a given species.
  Data<real, 4> Emission_out(Nt, Nz, Ny, Nx);

  // 0 for weekdays, 1 for Saturday and 2 for Sunday.
  int day_type;

  // For all species.
  for (s = 0; s < Ns; s++)
    {
      cout << "   + " << species_name[s] << endl;
      Date current_date(date);
      current_date.AddHours(-1);
      // Hour per hour.
      for (h = 0; h < Nt; h++)
        {
          day_type = max(current_date.GetWeekDay() - 4, 0);
          for (k = 0; k < Nz; k++)
            for (j = 0; j < Ny; j++)
              for (i = 0; i < Nx; i++)
                Emission_out(h, k, j, i) =
                  Emission_in(s, day_type, j, i, k, current_date.GetHour());
          current_date.AddHours(1);
        }
      FormatBinary<float>().Append(Emission_out,
                                   directory_out + species_name[s] + ".bin");
    }

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
