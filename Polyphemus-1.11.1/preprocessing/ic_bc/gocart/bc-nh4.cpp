// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Edouard Debry
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


// Generates boundary conditions for Polair3D based on Gocart concentrations
// for particulate matter.


//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string main_config_file("../general.cfg"), sec_config_file("");

  if (argc != 2 && argc != 3)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file]";
      mesg += " [secondary config file] \n\n";
      mesg += "Arguments:\n";
      mesg += "  [main configuration file] (optional): general configuration";
      mesg += " file, default is ../general.cfg.\n";
      mesg += "  [secondary configuration file] : NH4 configuration file.\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 2)
    sec_config_file = argv[1];
  else if (argc == 3)
    {
      main_config_file = argv[1];
      sec_config_file = argv[2];
    }

  if (!exists(main_config_file))
    throw string("Unable to find configuration file \"")
      + main_config_file + "\".";


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef float real;
  int j, s, t;

  // Configuration.
  ConfigStreams config(main_config_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  // Input/output directories and files.
  string Directory_out;

  config.SetSection("[paths]");
  config.PeekValue("Directory_bc", Directory_out);

  // Output domain.
  int Nx_out, Ny_out, Nz_out;

  config.SetSection("[domain]");
  config.PeekValue("Nz", "> 0", Nz_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Nx", "> 0", Nx_out);

  // Input species.
  config.NoSection();
  config.FindFromBeginning("[input_species]");

  vector<vector<string> > list_species;

  while (!config.IsEmpty() && config.PeekElement()[0] != '[')
    list_species.push_back(split(config.GetLine(),
                                 config.GetStreams()[0]->GetDelimiters()));

  int Nelec = list_species[0].size() - 1;
  cout << "Number of species for electroneutrality: " << Nelec << endl;

  // Sets species for electroneutrality.
  vector<string> electro_species;
  for (j = 0; j < Nelec; j++)
    electro_species.push_back(list_species[0][j + 1]);

  // Sets stoechiometric coefficient.
  Array<real, 1> stoechio_coefficient(Nelec);
  for (j = 0; j < Nelec; j++)
    stoechio_coefficient(j) = to_num<real>(list_species[1][j + 1]);

  // Sets molar weight of input species.
  Array<real, 1> molar_weight(Nelec);
  for (j = 0; j < Nelec; j++)
    molar_weight(j) = to_num<real>(list_species[2][j + 1]);

  cout << "Species:";
  for (j = 0; j < Nelec; j++)
    cout << "\t" << electro_species[j];
  cout << endl << "Stoechiometric coefficient:";
  for (j = 0; j < Nelec; j++)
    cout << "\t" << stoechio_coefficient(j);
  cout << endl << "Molar weight:";
  for (j = 0; j < Nelec; j++)
    cout << "\t" << molar_weight(j);
  cout << endl;

  // Output species.
  config.SetSection("[output_species]");

  // Name of species.
  string species_name;
  config.PeekValue("Species_name", species_name);

  // Molar weight of species.
  real molar_weight_nh4;
  config.PeekValue("Molar_weight", "positive", molar_weight_nh4);

  // Number of bins.
  int Nbin;
  config.PeekValue("Bin_number", "positive", Nbin);
  cout << "Number of bins: " << Nbin << endl;

  // Number of output time-steps (days) to compute.
  // Guessed with size of given file "*_x.bin"
  int Nt_out;
  string file_x;
  config.PeekValue("File_x", file_x);

  file_x = Directory_out + file_x;

  if (exists(file_x))
    Nt_out = int(file_size(file_x)) / (Nz_out * Ny_out * 2 * sizeof(float));
  else
    throw string("File ") + file_x + string(" not found,\n")
      + string("need this file to guess the number of output time-steps.");

  cout << "Number of days to compute: " << Nt_out << endl;


  // Final coefficient.
  Array<real, 1> coefficient(Nelec);
  for (j = 0; j < Nelec; j++)
    coefficient(j) = stoechio_coefficient(j)
      * molar_weight_nh4 / molar_weight(j);

  // Output settings.
  FormatBinary<real> PolairFormat;

  // Concentrations and temporary arrays.
  Data<real, 3> Conc_out_x(Nz_out, Ny_out, 2);
  Data<real, 3> Conc_out_y(Nz_out, 2, Nx_out);
  Data<real, 3> Conc_out_z(1, Ny_out, Nx_out);

  Data<real, 3> Conc_in_x(Nz_out, Ny_out, 2);
  Data<real, 3> Conc_in_y(Nz_out, 2, Nx_out);
  Data<real, 3> Conc_in_z(1, Ny_out, Nx_out);


  //////////////////////
  // COMPUTES OUTPUTS //
  //////////////////////


  // Loop over sizes and time-steps.
  for (j = 0; j < Nbin; j++)
    {
      // The bin.
      string bin = to_str(j);

      // Output file names.
      string spec_size = species_name + "_" + bin;
      string file_out_x = Directory_out + spec_size + "_x.bin";
      string file_out_y = Directory_out + spec_size + "_y.bin";
      string file_out_z = Directory_out + spec_size + "_z.bin";

      ofstream data_out_x;
      ofstream data_out_y;
      ofstream data_out_z;

      data_out_x.open(file_out_x.c_str());
      data_out_y.open(file_out_y.c_str());
      data_out_z.open(file_out_z.c_str());

      cout << "Computing " << spec_size << "..." << endl;

      // Sets file name of electro species (input files).
      vector<string> file_in_x;
      vector<string> file_in_y;
      vector<string> file_in_z;

      for (s = 0; s < Nelec; s++)
        {
          file_in_x.push_back(Directory_out + electro_species[s]
                              + "_" + bin + "_x.bin");
          file_in_y.push_back(Directory_out + electro_species[s]
                              + "_" + bin + "_y.bin");
          file_in_z.push_back(Directory_out + electro_species[s]
                              + "_" + bin + "_z.bin");
        }

      // Open input files.
      ifstream data_in_x[4];
      ifstream data_in_y[4];
      ifstream data_in_z[4];

      for (s = 0; s < Nelec; s++)
        {
          if (!exists(file_in_x[s]))
            throw string("File ") + file_in_x[s] + string(" not found.");

          if (!exists(file_in_y[s]))
            throw string("File ") + file_in_y[s] + string(" not found.");

          if (!exists(file_in_z[s]))
            throw string("File ") + file_in_z[s] + string(" not found.");

          data_in_x[s].open(file_in_x[s].c_str());
          data_in_y[s].open(file_in_y[s].c_str());
          data_in_z[s].open(file_in_z[s].c_str());
        }

      // Read files day per day.
      for (t = 0; t < Nt_out; t++)
        {
          Conc_out_x.SetZero();
          Conc_out_y.SetZero();
          Conc_out_z.SetZero();

          for (s = 0; s < Nelec; s++)
            {
              // Read input.
              PolairFormat.Read(data_in_x[s], Conc_in_x);
              PolairFormat.Read(data_in_y[s], Conc_in_y);
              PolairFormat.Read(data_in_z[s], Conc_in_z);

              // Computes electroneutrality.
              Conc_out_x.GetArray() = Conc_out_x.GetArray() +
                coefficient(s) * Conc_in_x.GetArray();

              Conc_out_y.GetArray() = Conc_out_y.GetArray() +
                coefficient(s) * Conc_in_y.GetArray();

              Conc_out_z.GetArray() = Conc_out_z.GetArray() +
                coefficient(s) * Conc_in_z.GetArray();
            }

          PolairFormat.Write(Conc_out_x, data_out_x);
          PolairFormat.Write(Conc_out_y, data_out_y);
          PolairFormat.Write(Conc_out_z, data_out_z);
        }

      cout << "Min mean max concentration in file " + spec_size + "_x.bin :"
           << Conc_out_x.GetMin() << "\t" << Conc_out_x.Mean()
           << "\t" << Conc_out_x.GetMax() << endl;
      cout << "Min mean max concentration in file " + spec_size + "_y.bin :"
           << Conc_out_y.GetMin() << "\t" << Conc_out_y.Mean()
           << "\t" << Conc_out_y.GetMax() << endl;
      cout << "Min mean max concentration in file " + spec_size + "_z.bin :"
           << Conc_out_z.GetMin() << "\t" << Conc_out_z.Mean()
           << "\t" << Conc_out_z.GetMax() << endl;

      // Check clipping.
      if (Conc_out_x.GetMin() < 0.0 ||
          Conc_out_y.GetMin() < 0.0 ||
          Conc_out_z.GetMin() < 0.0)
        throw string("Negative concentrations, ")
                                + string("possible error in electroneutrality.");

      // Close all files.
      data_out_x.close();
      data_out_y.close();
      data_out_z.close();

      for (s = 0; s < Nelec; s++)
        {
          data_in_x[s].close();
          data_in_y[s].close();
          data_in_z[s].close();
        }
    }

  END;

  return 0;
}
