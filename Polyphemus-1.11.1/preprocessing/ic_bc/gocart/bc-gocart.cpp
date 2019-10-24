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
#include <fstream>
#include <vector>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "modal_distribution.hxx"

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string main_config_file("../general.cfg"), gocart_config_file("");
  string gocart_file("");
  Date date_gocart(0);
  int number_of_days = 0;

  if (argc != 5 && argc != 6)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file]";
      mesg += "  [Gocart config file] [Gocart file] [Gocart date] [Number of days] \n\n";
      mesg += "Arguments:\n";
      mesg += "  [main config file] (optional): general configuration";
      mesg += " file, default is " + main_config_file + ".\n";
      mesg += "  [Gocart config file] : Gocart species configuration file.\n";
      mesg += "  [Gocart file]: output file from Gocart.\n";
      mesg += "  [Gocart date]: date of Gocart file, YYYYMM.\n";
      mesg += "  [Number of days]: number of days to compute.\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 6)
    {
      main_config_file = argv[1];
      gocart_config_file = argv[2];
      gocart_file = argv[3];
      date_gocart.SetDate(convert<int>(argv[4]) * 100 + 1);
      number_of_days = convert<int>(argv[5]);
    }
  else if (argc == 5)
    {
      gocart_config_file = argv[1];
      gocart_file = argv[2];
      date_gocart.SetDate(convert<int>(argv[3]) * 100 + 1);
      number_of_days = convert<int>(argv[4]);
    }

  if (!exists(main_config_file))
    throw string("Unable to find configuration file \"")
      + main_config_file + "\".";


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef float real;

  int i, j;

  // Configuration.
  ConfigStreams config(main_config_file);
  if (exists(gocart_config_file))
    config.AddFile(gocart_config_file);

  // Input/output directories and files.
  string file_temp;
  string Directory_out;

  config.SetSection("[paths]");
  config.PeekValue("Directory_bc", Directory_out);

  // Input domain.
  int Nt_in, Nz_in, Ny_in, Nx_in, Ns_in;
  real Delta_y_in, Delta_x_in, scale_height;
  real surface_pressure, top_pressure;
  real y_min_in, x_min_in;
  string sigma_levels;

  config.SetSection("[bc_input_domain]");
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Nx", "> 0", Nx_in);
  config.PeekValue("Delta_y", "> 0", Delta_y_in);
  config.PeekValue("Delta_x", "> 0", Delta_x_in);
  config.PeekValue("y_min", y_min_in);
  config.PeekValue("x_min", x_min_in);

  config.PeekValue("Nz", "> 0", Nz_in);
  config.PeekValue("Sigma_levels", sigma_levels);
  config.PeekValue("Scale_height", "positive", scale_height);
  config.PeekValue("Surface_pressure", "positive", surface_pressure);
  config.PeekValue("Top_pressure", "positive", top_pressure);

  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out, Ns_out;
  int t_min_out;
  real Delta_y_out, Delta_x_out, Delta_t_meteo;
  real y_min_out, x_min_out;
  string vertical_levels;

  config.SetSection("[domain]");
  config.PeekValue("Nz", "> 0", Nz_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("Delta_t", "> 0", Delta_t_meteo);
  config.PeekValue("Delta_y", "> 0", Delta_y_out);
  config.PeekValue("Delta_x", "> 0", Delta_x_out);
  config.PeekValue("y_min", y_min_out);
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("Vertical_levels", vertical_levels);

  // Dates.
  string date_meteo_str;
  config.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  if (date_gocart.GetYear() != date_meteo.GetYear())
    throw string("Gocart and meteorological years mismatch: ")
      + to_str(date_gocart.GetYear()) + " and " + to_str(date_meteo.GetYear())
      + " respectively.";

  cout << "Gocart month: " << date_gocart.GetMonth() << endl;

  if (date_gocart.GetMonth() > date_meteo.GetMonth())
    t_min_out = 0;
  else
    t_min_out = date_meteo.GetDay() - 1;

  // Size distribution.
  int Nbin;
  real Dmin, Dmax, Ntot[3], Dmean[3], Sigma[3];

  config.SetSection("[size_distribution]");
  config.PeekValue("Bin_number", "positive", Nbin);
  config.PeekValue("Diameter_min", "positive", Dmin);
  config.PeekValue("Diameter_max", "positive | > " + to_str(Dmin), Dmax);
  config.PeekValue("Ntot_nuclei", "positive", Ntot[0]);
  config.PeekValue("Dmean_nuclei", "positive", Dmean[0]);
  config.PeekValue("Sigma_nuclei", "positive", Sigma[0]);
  config.PeekValue("Ntot_accumulation", "positive", Ntot[1]);
  config.PeekValue("Dmean_accumulation", "positive", Dmean[1]);
  config.PeekValue("Sigma_accumulation", "positive", Sigma[1]);
  config.PeekValue("Ntot_coarse", "positive", Ntot[2]);
  config.PeekValue("Dmean_coarse", "positive", Dmean[2]);
  config.PeekValue("Sigma_coarse", "positive", Sigma[2]);

  // Gocart input species.
  config.NoSection();
  config.FindFromBeginning("[input_species]");

  vector<vector<string> > data_size;

  while (!config.IsEmpty() && config.PeekElement()[0] != '[')
    data_size.push_back(split(config.GetLine(),
                              config.GetStreams()[0]->GetDelimiters()));

  Ns_in = data_size.size();

  Array<real, 2> size_range(Ns_in, 2);

  for (i = 0; i < Ns_in; i++)
    for (j = 0; j < 2; j++)
      size_range(i, j) = to_num<real>(data_size[i][j + 1]);

  // Computes size partition coefficient.
  Array<real, 2> size_distribution(Ns_in, Nbin);

  BinDist discretization(Nbin, Dmin, Dmax);

  vector<string> data_bin;

  for (j = 0; j < Nbin; j++)
    data_bin.push_back(to_str(j));

  for (i = 0; i < Ns_in; i++)
    {
      // Creates the aerosol modal density.
      ModalAerosol trimodal(size_range(i, 0), size_range(i, 1),
                            Ntot, Dmean, Sigma);

      Array<real, 1> vol_bin(Nbin);

      real vol_tot = 0.0;
      for (j = 0; j < Nbin; j++)
        {
          vol_bin(j) = trimodal.VolQuantity(discretization.GetLowBound(j),
                                            discretization.GetUpBound(j));
          vol_tot += vol_bin(j);
        }

      real size_tot = 0.0;
      for (j = 0; j < Nbin; j++)
        {
          size_distribution(i, j) = vol_bin(j) / vol_tot;
          size_tot += size_distribution(i, j);
        }

      // Checks if coefficients sum is equal to one.
      if (abs(size_tot - 1.0) > 1.e-6)
        throw string("ERROR: sum of partition coefficients is ")
          + to_str(size_tot) + string(" instead of 1.");
    }

  // Polair3D output species.
  config.FindFromBeginning("[output_species]");

  vector<vector<string> > data_species;

  while (!config.IsEmpty() && config.PeekElement()[0] != '[')
    data_species.push_back(split(config.GetLine(),
                                 config.GetStreams()[0]->GetDelimiters()));

  if (int(data_species[0].size() - 1) != Ns_in)
    throw to_str(data_species[0].size() - 1)
      + to_str(" lines found in [output_species] ")
      + to_str(" section: inconsistant with number of input Gocart species (")
      + to_str(Ns_in) + string(" species).");

  Ns_out = data_species.size();

  Array<real, 2> species_distribution(Ns_out, Ns_in);

  for (i = 0; i < Ns_out; i++)
    for (j = 0; j < Ns_in; j++)
      species_distribution(i, j) = to_num<real>(data_species[i][j + 1]);

  // Computes the number of days.
  Nt_in = int(file_size(gocart_file)) / sizeof(float)
    / (3 + Nx_in * Ny_in * Nz_in) / Ns_in;
  if (Nt_in < 28 || Nt_in > 31)
    throw to_str(Nt_in) + string(" days found. Please check the number ")
      + string("of available species in Gocart file (Ns).");

  cout << "Number of days in Gocart file: " << Nt_in << endl;

  cout << "Number of days given in arguments: "
       << number_of_days << endl;

  Date date_meteo_end(date_meteo);
  date_meteo_end.AddDays(number_of_days - 1);

  if (date_meteo_end.GetYear() == date_gocart.GetYear())
    {
      if (date_gocart.GetMonth() < date_meteo_end.GetMonth())
        Nt_out = Nt_in;
      else
        Nt_out = date_meteo_end.GetDay();
    }
  else
    Nt_out = Nt_in;

  Nt_out -= t_min_out;
  cout << "Number of days to compute: " << Nt_out << endl;

  ///////////
  // GRIDS //
  ///////////


  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  // Input grids.
  RegularGrid<real> GridS_in(Ns_in), GridT_in(Nt_in), GridZ_in(Nz_in),
    GridY_in(y_min_in, Delta_y_in, Ny_in),
    GridX_in(x_min_in, Delta_x_in, Nx_in);

  // Reads values.
  Data<real, 1> SigmaLevels(Nz_in);
  FormatFormattedText GocartLevels("<e " + to_str(Nz_in) + ">");
  GocartLevels.Read(sigma_levels, "0", SigmaLevels);

  for (int k = 0; k < Nz_in; k++)
    {
      real sigma = SigmaLevels(k);
      real pressure = sigma * (surface_pressure - top_pressure)
        + top_pressure;
      GridZ_in(k) = scale_height * log(surface_pressure / pressure);
    }

  // Output settings.

  // Output grids.
  RegularGrid<real> GridT_out(t_min_out, 1., Nt_out), GridZ_out(Nz_out),
    GridY_out(y_min_out, Delta_y_out, Ny_out),
    GridX_out(x_min_out, Delta_x_out, Nx_out);

  // Checks if time output grid is inside time output grid.
  if (GridT_out(0) < GridT_in(0) ||
      GridT_out(Nt_out - 1) > GridT_in(Nt_in - 1))
    throw string("Warning: output time grid exceeds input time output\n")
      + string("GridT_in min, max: ") + to_str(GridT_in(0)) + string("\t")
      + to_str(GridT_in(Nt_in - 1)) + string("\n")
      + string("GridT_out min, max: ") + to_str(GridT_out(0)) + string("\t")
      +  to_str(GridT_out(Nt_out - 1));

  // For boundary conditions.

  // All interfaces along z.
  RegularGrid<real> GridZ_all_interf_out(Nz_out + 1);

  // The interface where boundary concentration are required.
  RegularGrid<real> GridZ_interf_out(1);

  // Boundary layers along x and y.
  RegularGrid<real> GridY_interf_out(y_min_out - Delta_y_out / 2.,
                                     (Ny_out + 1) * Delta_y_out, 2);
  RegularGrid<real> GridX_interf_out(x_min_out - Delta_x_out / 2.,
                                     (Nx_out + 1) * Delta_x_out, 2);
  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels, GridZ_all_interf_out);

  // Sets values at nodes.
  for (int k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_all_interf_out(k) + GridZ_all_interf_out(k + 1)) / 2.;

  // Sets the boundary layer altitude.
  GridZ_interf_out(0) = 1.5 * GridZ_all_interf_out(Nz_out - 1)
    - 0.5 * GridZ_all_interf_out(Nz_out - 2);

  FormatBinary<real> PolairFormat;


  //////////
  // DATA //
  //////////


  // Input fields.
  Data<real, 5> Conc_in(GridS_in, GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Conc_sum_in(GridT_in, GridZ_in, GridY_in, GridX_in);

  // Concentrations and temporary arrays used to sum concentrations.
  Data<real, 4> Conc_out_x(GridT_out, GridZ_out, GridY_out, GridX_interf_out);
  Data<real, 4> Conc_out_y(GridT_out, GridZ_out, GridY_interf_out, GridX_out);
  Data<real, 4> Conc_out_z(GridT_out, GridZ_interf_out, GridY_out, GridX_out);
  // Previous contributions.
  Data<real, 4> Conc_prv_x(GridT_out, GridZ_out, GridY_out, GridX_interf_out);
  Data<real, 4> Conc_prv_y(GridT_out, GridZ_out, GridY_interf_out, GridX_out);
  Data<real, 4> Conc_prv_z(GridT_out, GridZ_interf_out, GridY_out, GridX_out);

  cout << " done" << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////


  cout << "Reading input concentrations...";
  cout.flush();

  FormatBinary<float> Gocart;

  ifstream gocart_stream(gocart_file.c_str());

  for (int h = 0; h < Nt_in; h++)
    for (int s = 0; s < Ns_in; s++)
      {
        // Skips dates and species number.
        gocart_stream.seekg(4 * 3, ios_base::cur);

        Array<real, 3> Conc_in_3D =
          Conc_in.GetArray()(s, h, Range::all(), Range::all(), Range::all());

        Gocart.Read(gocart_stream, Conc_in_3D);
      }

  // Converts from big endian to little endian.
  swap(Conc_in.GetArray());

  // Converts from g/m^3 to µg/m^3.
  Conc_in.GetArray() *= 1.0e06;

  cout << " done" << endl;


  //////////////////////
  // COMPUTES OUTPUTS //
  //////////////////////

  cout << endl << "Writes output in directory: "
       << Directory_out << endl << endl;

  // Gets minimum file size.
  int big_number = int(1.e9);
  int size_x_min = big_number;
  int size_y_min = big_number;
  int size_z_min = big_number;
  for (i = 0; i < Ns_out; i++)
    for (j = 0; j < Nbin; j++)
      {
        string spec_size = data_species[i][0] + "_" + data_bin[j];
        string file_x = Directory_out + spec_size + "_x.bin";
        string file_y = Directory_out + spec_size + "_y.bin";
        string file_z = Directory_out + spec_size + "_z.bin";

        if (exists(file_x))
          size_x_min = min(size_x_min, int(file_size(file_x)));
        else
          size_x_min = 0;

        if (exists(file_y))
          size_y_min = min(size_y_min, int(file_size(file_y)));
        else
          size_y_min = 0;

        if (exists(file_z))
          size_z_min = min(size_z_min, int(file_size(file_z)));
        else
          size_z_min = 0;
      }

  if (size_x_min > 0)
    cout << "Files *_x exist, next boundary conditions "
         << "appended to existing file." << endl;
  else
    cout << "Creates files *_x." << endl;

  if (size_y_min > 0)
    cout << "Files *_y exist, next boundary conditions "
         << "appended to existing file." << endl;
  else
    cout << "Creates files *_y." << endl;

  if (size_z_min > 0)
    cout << "Files *_z exist, next boundary conditions "
         << "appended to existing file." << endl;
  else
    cout << "Creates files *_z." << endl;

  cout << endl;

  // Loop over Polair3D species.
  for (i = 0; i < Ns_out; i++)
    {
      // Loop over sizes.
      for (j = 0; j < Nbin; j++)
        {
          // Part of file name common to all files.
          string spec_size = data_species[i][0] + "_" + data_bin[j];

          cout << "Computing concentrations " << spec_size << " ..." << endl;
          cout.flush();

          // Set name of each file for one Polair3D species.
          string file_x = Directory_out + spec_size + "_x.bin";
          string file_y = Directory_out + spec_size + "_y.bin";
          string file_z = Directory_out + spec_size + "_z.bin";

          // Zero initialization.
          Conc_out_x.SetZero();
          Conc_out_y.SetZero();
          Conc_out_z.SetZero();

          // Reads previous contribution if necessary, else sets to zero.
          if (exists(file_x) && int(file_size(file_x)) > size_x_min)
            {
              cout << "A contribution already exists for file "
                   <<  spec_size + "_x.bin" << endl;
              // a contribution has already been added.
              ifstream prev_data_x;
              prev_data_x.open(file_x.c_str());
              PolairFormat.Read(prev_data_x, Conc_prv_x);
              prev_data_x.seekg(size_x_min, ios_base::beg);
              prev_data_x.close();
            }
          else
            // no contribution.
            Conc_prv_x.SetZero();

          // Same for y.
          if (exists(file_y) && int(file_size(file_y)) > size_y_min)
            {
              cout << "A contribution already exists for file "
                   << spec_size + "_y.bin" << endl;
              ifstream prev_data_y;
              prev_data_y.open(file_y.c_str());
              PolairFormat.Read(prev_data_y, Conc_prv_y);
              prev_data_y.seekg(size_y_min, ios_base::beg);
              prev_data_y.close();
            }
          else
            Conc_prv_y.SetZero();

          // Same for z.
          if (exists(file_z) && int(file_size(file_z)) > size_z_min)
            {
              cout << "A contribution already exists for file "
                   << spec_size + "_z.bin" << endl;
              ifstream prev_data_z;
              prev_data_z.open(file_z.c_str());
              PolairFormat.Read(prev_data_z, Conc_prv_z);
              prev_data_z.seekg(size_z_min, ios_base::beg);
              prev_data_z.close();
            }
          else
            Conc_prv_z.SetZero();

          // Distributes input concentrations.
          Conc_sum_in.SetZero();
          for (int s = 0; s < Ns_in; s++)
            if (species_distribution(i, s) > 0.0
                && size_distribution(s, j) > 0.0)
              Conc_sum_in.GetArray() = Conc_sum_in.GetArray()
                + species_distribution(i, s) * size_distribution(s, j)
                * Conc_in.GetArray()(s, Range::all(), Range::all(),
                                     Range::all(), Range::all());

          // Checks clipping.
          if (Conc_sum_in.GetMin() < 0.0)
            throw string("Negative concentrations")
              + string("error in size or species coefficients.");

          // Interpolates from Gocart grid to Polair3D grid.
          LinearInterpolationRegular(Conc_sum_in, Conc_out_x);
          LinearInterpolationRegular(Conc_sum_in, Conc_out_y);
          LinearInterpolationRegular(Conc_sum_in, Conc_out_z);

          // Checks clipping.
          if (Conc_out_x.GetMin() < 0.0 ||
              Conc_out_y.GetMin() < 0.0 ||
              Conc_out_z.GetMin() < 0.0)
            throw  string("Negative concentrations, error in interpolation.");

          // Adds new contribution to old one.
          Conc_out_x.GetArray() = Conc_out_x.GetArray()
                                    + Conc_prv_x.GetArray();
          Conc_out_y.GetArray() = Conc_out_y.GetArray()
                                    + Conc_prv_y.GetArray();
          Conc_out_z.GetArray() = Conc_out_z.GetArray()
                                    + Conc_prv_z.GetArray();

          // Opens file for writing, created if not exists.
          ofstream prev_data_x;
          ofstream prev_data_y;
          ofstream prev_data_z;

          prev_data_x.open(file_x.c_str(), ofstream::app);
          prev_data_y.open(file_y.c_str(), ofstream::app);
          prev_data_z.open(file_z.c_str(), ofstream::app);

          // Makes sure pointer is at the right place.
          prev_data_x.seekp(size_x_min, ios_base::beg);
          prev_data_y.seekp(size_y_min, ios_base::beg);
          prev_data_z.seekp(size_z_min, ios_base::beg);

          cout << "Min mean max concentration in file " << spec_size + "_x.bin"
               << ": " << Conc_out_x.GetMin() << "\t" << Conc_out_x.Mean()
               << "\t" << Conc_out_x.GetMax() << endl;
          cout << "Min mean max concentration in file " << spec_size + "_y.bin"
               << ": " << Conc_out_y.GetMin() << "\t" << Conc_out_y.Mean()
               << "\t" << Conc_out_y.GetMax() << endl;
          cout << "Min mean max concentration in file " << spec_size + "_z.bin"
               << ": " << Conc_out_z.GetMin() << "\t" << Conc_out_z.Mean()
               << "\t" << Conc_out_z.GetMax() << endl;

          // Writes final concentration, overwrites if already one.
          PolairFormat.Write(Conc_out_x, prev_data_x);
          PolairFormat.Write(Conc_out_y, prev_data_y);
          PolairFormat.Write(Conc_out_z, prev_data_z);

          // Closes files of given Polair3D species.
          prev_data_x.close();
          prev_data_y.close();
          prev_data_z.close();

          cout << " done" << endl;
        }
    }

  cout << endl;

  END;

  return 0;
}
