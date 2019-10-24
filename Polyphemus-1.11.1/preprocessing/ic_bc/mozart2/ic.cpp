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


// Generates initial conditions for Polair based on Mozart simulations downloaded
// from the NCAR community data portal.
//
// NCAR community data portal: http://cdp.ucar.edu/home/home.htm
// Mozart homepage:
// http://www.mpimet.mpg.de/en/wissenschaft/modelle/mozart.html

// Usage:
//
// [program] [main config file] [secondary config file]


//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_NETCDF

#include "SeldonData.hxx"
using namespace SeldonData;

#include "AtmoData.hxx"
using namespace AtmoData;

// INCLUDES //
//////////////


// List of species associated with a ratio.
class SpeciesInfo
{
public:
  vector<string> species_in;
  vector<float> ratio;
};


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string main_config_file("ic.cfg"), sec_config_file("");

  if (argc != 2 && argc != 3)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file] [secondary config file]\n";
      mesg += string("  ") + argv[0] + " [main configuration file]\n";
      mesg += "Arguments:\n";
      mesg += "  [main configuration file] (optional): main configuration file. Default: ic.cfg\n";
      mesg += "  [secondary configuration file] (optional): secondary configuration file.\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 3)
    sec_config_file = argv[2];
  if (argc > 1)
    main_config_file = argv[1];

  if (!exists(main_config_file))
    throw string("Unable to find configuration file \"")
      + main_config_file + "\".";


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Temporary values.
  int h, i, j, k, l;

  // Configuration.
  ConfigStreams config(main_config_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);

  // Constants.
  const real R = 8.314;

  // Input domain.
  int Nt_in, Nz_in, Ny_in, Nx_in, Nx_in_tmp;
  real Delta_t_in;
  string mozart_directory;
  string date_str;

  config.SetSection("[ic_input_domain]");
  config.PeekValue("Nt", "> 0", Nt_in);
  config.PeekValue("Delta_t", "> 0", Delta_t_in);
  config.PeekValue("Nz", "> 0", Nz_in);
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Nx", "> 0", Nx_in_tmp);
  // Grid along X is doubled compared to the related Mozart grid so that
  // interpolation along X could be performed between -360 degrees and
  // +357.1875 degrees.
  Nx_in = 2 * Nx_in_tmp;

  config.PeekValue("Database_ic", mozart_directory);

  // Date.
  config.PeekValue("Date_ic", date_str);
  Date date(date_str);
  cout << "Date: " << date.GetDate("%y-%m-%d") << endl << endl;

  int file_index = 40 + (date.GetNumberOfDays() + 6) / 10;
  string mozart_file = mozart_directory + string("h00")
    + to_str(file_index) + ".nc";


  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;
  string vertical_levels;

  config.SetSection("[domain]");
  Nt_out = 1;
  config.PeekValue("Nz", "> 0", Nz_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("Delta_y", "> 0", Delta_y_out);
  config.PeekValue("Delta_x", "> 0", Delta_x_out);
  config.PeekValue("y_min", y_min_out);
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("Vertical_levels", vertical_levels);

  // Input/output directories and files.
  string species_file, molecular_weights;
  string Directory_out, file_out;

  config.SetSection("[ic_files]");
  config.PeekValue("Directory_ic", Directory_out);
  config.PeekValue("Species", species_file);
  config.PeekValue("Molecular_weights", molecular_weights);


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  // Input grids.
  // The suffix '_in_tmp' for grids indicates it will contain the initial
  // Mozart grid.

  RegularGrid<real> GridT_in(Nt_in);
  RegularGrid<double> GridT_in_double(Nt_in);
  // Vertical levels depend on t, z, y and x.
  GeneralGrid<real, 4> GridZ_in(shape(Nt_in, Nz_in, Ny_in, Nx_in),
                                1, shape(0, 1, 2, 3));
  GeneralGrid<real, 4> GridZ_in_tmp(shape(Nt_in, Nz_in, Ny_in, Nx_in_tmp),
                                    1, shape(0, 1, 2, 3));
  RegularGrid<real> GridY_in(Ny_in);
  RegularGrid<real> GridX_in(Nx_in), GridX_in_tmp(Nx_in_tmp);
  // Pressure levels.
  RegularGrid<real> GridP_in(Nz_in);
  // Grids are shared.
  GridZ_in.SetVariable(1);
  GridZ_in.SetDuplicate(false);
  GridZ_in_tmp.SetVariable(1);
  GridZ_in_tmp.SetDuplicate(false);

  // Reads values.
  FormatNetCDF<float> Mozart;
  Mozart.Read(mozart_file, "lat", GridY_in);
  Mozart.Read(mozart_file, "lon", GridX_in_tmp);
  Mozart.Read(mozart_file, "lev", GridP_in);
  Mozart.Read(mozart_file, "time", GridT_in_double);

  // Modifies the grid to allow interpolation between -360 degrees and
  // +357.1875 degrees.
  for (i = 0; i < Nx_in_tmp; i++)
    GridX_in(i) = GridX_in_tmp(i) - 360.;
  for (i = Nx_in_tmp; i < Nx_in; i++)
    GridX_in(i) = GridX_in_tmp(i - Nx_in_tmp);

  // Output settings.

  FormatBinary<float> Polair;

  // Output grids.
  RegularGrid<real> GridT_out(Nt_out);
  GridT_out(0) = real(365 * 9 + 3 + date.GetNumberOfDays());
  RegularGrid<real> GridZ_out(Nz_out);
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  // All interfaces along z.
  RegularGrid<real> GridZ_all_interf_out(Nz_out + 1);
  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels, GridZ_all_interf_out);
  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_all_interf_out(k) + GridZ_all_interf_out(k + 1)) / 2.;
  // Sets time steps.
  for (i = 0; i < GridT_in.GetLength(); i++)
    GridT_in(i) = GridT_in_double(i);

  // Checks whether the needed date is available.
  if (GridT_out(0) < GridT_in(0) || GridT_out(0) > GridT_in(Nt_in - 1))
    throw string("Date ") + date.GetDate("%y-%m-%d") + " not found.";


  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 4> Temperature(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Temperature_tmp(GridT_in, GridZ_in_tmp, GridY_in,
                                GridX_in_tmp);

  // Fixed pressure levels.
  Data<real, 1> alpha(Nz_in), beta(Nz_in);
  Data<real, 3> SurfacePressure(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SurfacePressure_tmp(GridT_in, GridY_in, GridX_in_tmp);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Conc_in(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Conc_in_tmp(GridT_in, GridZ_in_tmp, GridY_in, GridX_in_tmp);

  // Output fields.
  Data<real, 4> Pressure_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> Temperature_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  // Concentrations and temporary arrays used to sum concentrations.
  Data<real, 4> Conc_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> Conc_out_tmp(GridT_out, GridZ_out, GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  cout << "Extracting temperature...";
  cout.flush();
  Mozart.Read(mozart_file, "T", Temperature_tmp);

  // From the temperature read in the Mozart file, builds a temperature field
  // that allows interpolation along X in an extended range between
  // -360 degrees and +357.1875 degrees.
  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in; j++)
        {
          for (i = 0; i < Nx_in_tmp; i++)
            Temperature(h, k, j, i) = Temperature_tmp(h, k, j, i);
          for (i = Nx_in_tmp; i < Nx_in; i++)
            Temperature(h, k, j, i) = Temperature_tmp(h, k, j, i - Nx_in_tmp);
        }

  // Note: data is provided from the top level to the bottom level
  // (altitude). The interpolation function requires coordinates
  // to be sorted in increasing order. So, data must be reversed...
  Temperature.ReverseData(1);

  cout << " done." << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  cout << "Computing level heights in meter and interpolating pressure and temperature...";
  cout.flush();

  // Hybrid coefficients.
  Mozart.Read(mozart_file, "hyam", alpha);
  Mozart.Read(mozart_file, "hybm", beta);
  Mozart.Read(mozart_file, "PS", SurfacePressure_tmp);

  // From the surface pressure read in the Mozart file, builds a surface
  // pressure field that allows interpolation along X in an extended range
  // between -360 degrees and +357.1875 degrees.
  for (h = 0; h < Nt_in; h++)
    for (j = 0; j < Ny_in; j++)
      {
        for (i = 0; i < Nx_in_tmp; i++)
          SurfacePressure(h, j, i) = SurfacePressure_tmp(h, j, i);
        for (i = Nx_in_tmp; i < Nx_in; i++)
          SurfacePressure(h, j, i) = SurfacePressure_tmp(h, j, i - Nx_in_tmp);
      }

  ComputePressure(alpha, beta, SurfacePressure, Pressure, real(100000.));
  Pressure.ReverseData(1);

  ComputeHeight(SurfacePressure, Pressure, Temperature, GridZ_in);

  // To Polair grid...
  LinearInterpolationOneGeneral(Pressure, Pressure_out, 1);
  LinearInterpolationOneGeneral(Temperature, Temperature_out, 1);

  cout << " done." << endl;
  cout << endl;


  ////////////////
  // MOZART MAP //
  ////////////////

  cout << "Reading, interpolating and writing concentrations..." << endl;
  cout << "Maps input and output species...";
  cout.flush();

  // File "species.txt" format (one line):
  // [Mozart species] [Polair species] {[ratio] {[Polair species] [ratio] {...}}}
  ifstream species(species_file.c_str());

  if (!species.is_open())
    throw string("Unable to open file \"") + species_file + "\".";

  string line, species_in, species_out;

  string delim = " \t\n|";

  // To all output (Polair) species, we associate a list of input (Mozart)
  // species and their ratios (i.e. the proportion of the input species to
  // be put in the output species).
  map<string, SpeciesInfo> species_map;
  map<string, SpeciesInfo>::iterator pos;

  // For all lines.
  while (getline(species, line))
    {

      // 'species_data' will be filled with all words on the line.
      vector<string> species_data;

      string::size_type beg_index, end_index;
      beg_index = line.find_first_not_of(delim);

      // For each "word".
      while (beg_index != string::npos)
        {
          end_index = line.find_first_of(delim, beg_index);
          if (end_index == string::npos)
            end_index = line.length();

          // Put the word in 'species_data'.
          species_data.push_back(line.substr(beg_index, end_index - beg_index));

          beg_index = line.find_first_not_of(delim, end_index);
        }

      // Number of words on the line (Mozart species and Polair species and ratios).
      int size = species_data.size();

      // If the input species is associated with any output species.
      if (size > 1)
        {
          // Number of output species.
          int nb = max((size - 1) / 2, 1);
          // For all output species, update its list of input species and ratios.
          for (i = 0; i < nb; i++)
            // If the output species is not in the map yet.
            if ((pos = species_map.find(species_data[2 * i + 1]))
                == species_map.end())
              {
                SpeciesInfo info;
                info.species_in.push_back(species_data[0]);
                if (size % 2 == 1) // If ratios are specified.
                  info.ratio.push_back(to_num<float>(species_data[2 * (i + 1)]));
                else  // No ratio: assumed to be 1.
                  info.ratio.push_back(1.);
                species_map.insert(make_pair(species_data[2 * i + 1], info));
              }
            else  // the output is in the map.
              {
                pos->second.species_in.push_back(species_data[0]);
                if (size % 2 == 1) // If ratios are specified.
                  pos->second.ratio.push_back(to_num<float>(species_data[2 * (i + 1)]));
                else  // No ratio: assumed to be 1.
                  pos->second.ratio.push_back(1.);
              }
        }

    }  // Loop over 'species_file' lines.

  species.close();


  ////////////////////////////
  // OUTPUT SPECIES WEIGHTS //
  ////////////////////////////

  // Molecular weights of output species.
  map<string, float> weights_map;
  map<string, float>::iterator weights_pos;

  ifstream weights(molecular_weights.c_str());

  if (!weights.is_open())
    throw string("Unable to open file \"") + molecular_weights + "\".";

  float weight;

  // For all lines.
  while (getline(weights, line))
    {
      istringstream sline(line);
      sline >> species_out;
      // If there is still species.
      if (sline.good())
        {
          sline >> weight;
          weights_map.insert(make_pair(species_out, weight));
        }
    }  // Loop over 'molecular_weights' lines.

  weights.close();

  cout << " done." << endl;


  /////////////////////////////
  // BOUNDARY CONCENTRATIONS //
  /////////////////////////////

  cout << "Species: " << endl;

  // For all output species (in the map 'species_map').
  for (pos = species_map.begin(); pos != species_map.end(); pos++)
    {
      // pos->first: output species name.
      // pos->second: input species names and ratios.

      // Displays the output species name.
      cout << "Output: " << pos->first << "; input:";
      cout.flush();

      Conc_out.SetZero();

      // For all input species.
      for (l = 0; l < int(pos->second.species_in.size()); l++)
        {
          // Input species name.
          cout << " " << pos->second.species_in[l];
          cout.flush();
          // Reads input species concentrations.
          Mozart.Read(mozart_file, pos->second.species_in[l] + "_VMR_avrg", Conc_in_tmp);

          Conc_in_tmp.ReverseData(1);

          for (h = 0; h < Nt_in; h++)
            for (k = 0; k < Nz_in; k++)
              for (j = 0; j < Ny_in; j++)
                {
                  for (i = 0; i < Nx_in_tmp; i++)
                    Conc_in(h, k, j, i) = Conc_in_tmp(h, k, j, i);
                  for (i = Nx_in_tmp; i < Nx_in; i++)
                    Conc_in(h, k, j, i) = Conc_in_tmp(h, k, j, i - Nx_in_tmp);
                }

          // Interpolation.
          LinearInterpolationOneGeneral(Conc_in, Conc_out_tmp, 1);
          Conc_out_tmp.ThresholdMin(0.);

          // Adds concentrations with the weighting coefficient 'ratio'.
          for (h = 0; h < Nt_out; h++)
            for (k = 0; k < Nz_out; k++)
              for (j = 0; j < Ny_out; j++)
                for (i = 0; i < Nx_out; i++)
                  Conc_out(h, k, j, i) += pos->second.ratio[l]
                    * Conc_out_tmp(h, k, j, i);
        }

      float weight_loc;
      if (weights_map.find(pos->first) == weights_map.end())
        throw string("Unable to find the output species ")
          + pos->first + ".";
      else
        weight_loc = weights_map.find(pos->first)->second;

      // Conversion from ppt (?) to micrograms per cubic meter.
      for (h = 0; h < Nt_out; h++)
        for (k = 0; k < Nz_out; k++)
          for (j = 0; j < Ny_out; j++)
            for (i = 0; i < Nx_out; i++)
              Conc_out(h, k, j, i) *= Pressure_out(h, k, j, i)
                * weight_loc / R / Temperature_out(h, k, j, i) * 1000000.;
      cout << endl;
      cout << "   (min : " << Conc_out.GetMin() << ")";
      cout << "   (mean: " << Conc_out.Mean() << ")";
      cout << "   (max: " << Conc_out.GetMax() << ")";
      cout << endl;

      Polair.Write(Conc_out, Directory_out + pos->first + ".bin");
    }

  cout << " done." << endl;
  cout << endl;

  END;

  return 0;

}
