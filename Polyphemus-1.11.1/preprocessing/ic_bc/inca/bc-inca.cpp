// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Denis Qu�lo
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

using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_NETCDF

#include "AtmoData.hxx"
using namespace AtmoData;


// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string main_config_file("bc-inca.cfg"), sec_config_file("");

  if (argc != 3 && argc != 4)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file] [secondary config file] [Inca file]\n";
      mesg += string("  ") + argv[0] + " [main configuration file] [Inca file]\n";
      mesg += string("  ") + argv[0] + " [Inca file]\n\n";
      mesg += "Arguments:\n";
      mesg += "  [main configuration file] (optional): main configuration file. Default: bc-inca.cfg\n";
      mesg += "  [secondary configuration file] (optional): secondary configuration file.\n";
      mesg += "  [Inca file]: output file from Inca.\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 4)
    sec_config_file = argv[2];
  if (argc > 2)
    main_config_file = argv[1];
  string inca_file = argv[argc - 1];

  if (!exists(main_config_file))
    throw string("Unable to find configuration file \"")
      + main_config_file + "\".";


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Configuration.
  ConfigStreams config(main_config_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  // Input domain.
  int Nz_in, Ny_in, Nx_in;
  real Delta_y_in, Delta_x_in;
  real y_min_in, x_min_in;
  string species_file;

  config.SetSection("[bc_input_domain]");
  config.PeekValue("Nz", "> 0", Nz_in);
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Nx", "> 0", Nx_in);
  config.PeekValue("Delta_y", "> 0", Delta_y_in);
  config.PeekValue("Delta_x", "> 0", Delta_x_in);
  config.PeekValue("y_min", y_min_in);
  config.PeekValue("x_min", x_min_in);

  config.PeekValue("Species", species_file);

  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;
  string hybrid_coefficient_file;

  config.SetSection("[domain]");
  config.PeekValue("Nz", "> 0", Nz_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("Delta_y", "> 0", Delta_y_out);
  config.PeekValue("Delta_x", "> 0", Delta_x_out);
  config.PeekValue("y_min", y_min_out);
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("Vertical_levels", hybrid_coefficient_file);

  // Adds the cells around the domain.
  Nx_out += 2;
  Ny_out += 2;
  y_min_out -= Delta_y_out;
  x_min_out -= Delta_x_out;

  // Input/output directories and files.
  string Directory_out, file_out;

  config.SetSection("[bc_files]");
  config.PeekValue("Nt", "> 0", Nt_out);
  config.PeekValue("Directory_bc", Directory_out);


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  // Input grids.
  RegularGrid<real> GridY_in(y_min_in, Delta_y_in, Ny_in);
  RegularGrid<real> GridX_in(x_min_in, Delta_x_in, Nx_in);
  // Pressure levels.
  RegularGrid<real> GridZ_in(Nz_in);

  // Output settings.
  FormatBinary<float> ctm;

  // Output grids.
  RegularGrid<real> GridZ_out(Nz_out);
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  // Output vertical grid.
  FormatFormattedText InputLevels("<e><e>");
  Array<real, 1> alpha(Nz_out);
  Array<real, 1> beta(Nz_out);
  InputLevels.Read(hybrid_coefficient_file, "0", alpha);
  InputLevels.Read(hybrid_coefficient_file, "1", beta);

  // Reads the species.
  vector<string> species_name;
  ExtStream species(species_file);
  if (!species.is_open())
    throw string("Unable to open file \"") + species_file + "\".";
  while (!species.IsEmpty())
    species_name.push_back(species.GetElement());

  int Ns = species_name.size();

  cout << " done" << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  cout << "Reads file...";
  cout.flush();

  Data<real, 4> inca_data(Nz_in, Ny_in, Nx_in, 5 + Ns);
  FormatText().Read(inca_file, inca_data);

  cout << " done" << endl;


  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 2> SurfacePressure(GridY_in, GridX_in);
  Data<real, 3> Pressure(GridZ_in, GridY_in, GridX_in);
  Data<real, 3> Conc_in(GridZ_in, GridY_in, GridX_in);

  // Output fields.
  Data<real, 2> SurfacePressure_h_out(GridY_out, GridX_out);
  Data<real, 3> Pressure_h_out(GridZ_in, GridY_out, GridX_out);
  Data<real, 3> PressureCastor(GridZ_out, GridY_out, GridX_out);
  Data<real, 3> PressureCastor_interf(GridZ_out, GridY_out, GridX_out);
  Data<real, 3> Conc_h_out(GridZ_in, GridY_out, GridX_out);
  Data<real, 3> Conc_out(GridZ_out, GridY_out, GridX_out);
  Data<real, 3> Conc_out_x(Nz_out, Ny_out - 2, 2);
  Data<real, 3> Conc_out_y(Nz_out, 2, Nx_out - 2);
  Data<real, 2> Conc_out_z(Ny_out - 2, Nx_out - 2);


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  cout << "Input data processing...";
  cout.flush();

  inca_data.ReverseData(1);

  Pressure.SubData(inca_data, Range::all(), Range::all(), Range::all(), 4);
  SurfacePressure.SubData(Pressure, 0, Range::all(), Range::all());
  SurfacePressure.Mlt(1.008);

  // Indices and weights for horizontal interpolation.
  Array<int, 1> Index_i_in(Nx_out);
  Array<int, 1> Index_j_in(Ny_out);
  Array<real, 1> Weight_i_in(Nx_out);
  Array<real, 1> Weight_j_in(Ny_out);

  InterpolationRegularInput(GridX_in, GridX_out, Index_i_in, Weight_i_in);
  InterpolationRegularInput(GridY_in, GridY_out, Index_j_in, Weight_j_in);

  // Horizontal interpolation.
  HorizontalInterpolation(SurfacePressure(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in,
                          SurfacePressure_h_out());
  HorizontalInterpolation(Pressure(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in,
                          Pressure_h_out());

  for (int k = 0; k < Nz_out; k++)
    for (int j = 0; j < Ny_out; j++)
      for (int i = 0; i < Nx_out; i++)
        PressureCastor(k, j, i) = alpha(k) * 1.e5
          + beta(k) * SurfacePressure_h_out(j, i);

  for (int j = 0; j < Ny_out; j++)
    for (int i = 0; i < Nx_out; i++)
      PressureCastor_interf(0, j, i) = 0.5
        * (SurfacePressure_h_out(j, i) + PressureCastor(0, j, i));

  for (int k = 1; k < Nz_out; k++)
    for (int j = 0; j < Ny_out; j++)
      for (int i = 0; i < Nx_out; i++)
        PressureCastor_interf(k, j, i) = 0.5
          * (PressureCastor(k - 1, j, i) + PressureCastor(k, j, i));

  cout << " done" << endl;


  /////////////////////////////
  // BOUNDARY CONCENTRATIONS //
  /////////////////////////////

  cout << "Species: " << endl;

  for (int l = 0; l < Ns; l++)
    {
      cout << "\t" << species_name[l] << " ...";
      cout.flush();

      Conc_in.SubData(inca_data, Range::all(), Range::all(), Range::all(),
                      5 + l);

      HorizontalInterpolation(Conc_in(), Index_i_in, Index_j_in,
                              Weight_i_in, Weight_j_in,
                              Conc_h_out());
      for (int j = 0; j < Ny_out; j++)
        for (int i = 0; i < Nx_out; i++)
          {
            int k_in = 0;
            for (int k = 0; k < Nz_out; k++)
              {
                while (PressureCastor_interf(k, j, i)
                       <= Pressure_h_out(k_in, j, i))
                  k_in++;
                Conc_out(k, j, i) = Conc_h_out(k_in, j, i);
              }
          }

      // Beware that a Blitz range is including: [a, b] and not [a, b[.
      Conc_out_x.SubData(Conc_out,
                         Range::all(), // Nz_out
                         Range(1, Ny_out - 2), // Ny_out - 2
                         Range(0, Nx_out - 1, Nx_out - 1)); // 2

      Conc_out_y.SubData(Conc_out,
                         Range::all(), // Nz_out
                         Range(0, Ny_out - 1, Ny_out - 1),
                         Range(1, Nx_out - 2));

      Conc_out_z.SubData(Conc_out,
                         Nz_out - 1,
                         Range(1, Ny_out - 2),  // Ny_out - 2
                         Range(1, Nx_out - 2)); // Nx_out - 2

      for (int h = 0; h < Nt_out; h++)
        {
          ctm.Append(Conc_out_x,
                     Directory_out + species_name[l] + "_x.bin");
          ctm.Append(Conc_out_y,
                     Directory_out + species_name[l] + "_y.bin");
          ctm.Append(Conc_out_z,
                     Directory_out + species_name[l] + "_z.bin");
        }

      cout << " done" << endl;
    }

  cout << endl;

  END;

  return 0;

}
