// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Victor Winiarek
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

// Generates boundary conditions for Polair3D based on Mozart simulations
// downloaded from the NCAR community data portal.
//
// NCAR community data portal: http://cdp.ucar.edu/home/home.htm
// Mozart homepage:
// http://www.mpimet.mpg.de/en/wissenschaft/modelle/mozart.html


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

  string main_config_file("bc-dates.cfg"), sec_config_file("");

  if (argc != 3 && argc != 4 && argc != 5)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file]" +
        " [secondary config file] [First date] [Second date]\n";
      mesg += string("  ") + argv[0] + " [main configuration file]" +
        " [First date] [Second date]\n";
      mesg += string("  ") + argv[0] + " [First Date] [Second date]\n\n";
      mesg += "Arguments:\n";
      mesg += string("  [main configuration file] (optional): main") +
        " configuration file. Default: bc-dates.cfg\n";
      mesg += string("  [secondary configuration file] (optional):") +
        " secondary configuration file.\n";
      mesg += "  [First date]: first date of the simulation.\n";
      mesg += "  [Second date]: end date of the simulation\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 5)
    sec_config_file = argv[2];
  if (argc > 3)
    main_config_file = argv[1];
  string date_beg_str = argv[argc - 2];
  string date_end_str = argv[argc - 1];

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

  config.SetSection("[bc_input_domain]");
  config.PeekValue("Nt", "> 0", Nt_in);
  config.PeekValue("Delta_t", "> 0", Delta_t_in);
  config.PeekValue("Nz", "> 0", Nz_in);
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Nx", "> 0", Nx_in_tmp);
  // Grid along X is doubled compared to the related Mozart grid so that
  // interpolation along X could be performed between -360 degrees and
  // +357.1875 degrees.
  Nx_in = 2 * Nx_in_tmp;

  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;
  string vertical_levels;

  config.SetSection("[domain]");
  Nt_out = Nt_in;
  config.PeekValue("Nz", "> 0", Nz_out);
  config.PeekValue("Ny", "> 0", Ny_out);
  config.PeekValue("Nx", "> 0", Nx_out);
  config.PeekValue("Delta_y", "> 0", Delta_y_out);
  config.PeekValue("Delta_x", "> 0", Delta_x_out);
  config.PeekValue("y_min", y_min_out);
  config.PeekValue("x_min", x_min_out);
  config.PeekValue("Vertical_levels", vertical_levels);

  // Input/output directories and files.
  string species_file, molecular_weights, mozart_directory;
  string Directory_out, file_out;

  config.SetSection("[bc_files]");
  config.PeekValue("Directory_bc", Directory_out);
  config.PeekValue("Species", species_file);
  config.PeekValue("Molecular_weights", molecular_weights);
  config.PeekValue("Directory_mozart", mozart_directory);

  //////////////////////////////
  //      OUTPUT GRIDS        //

  // Output settings.

  FormatBinary<float> Polair;

  // Output grids.
  RegularGrid<real> GridT_out(Nt_out);
  RegularGrid<real> GridT_out_half(Nt_out / 2);
  RegularGrid<real> GridZ_out(Nz_out);
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

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

  //      OUTPUT GRIDS        //
  //////////////////////////////


  //////////////////////////////
  //MOZART FILES CORRESPONDING//
  // TO FIRST AND END DATE    //
  //////////////////////////////

  Date date_beg(date_beg_str);
  Date date_end(date_end_str);

  string mesg31  = "";

  if (date_beg.GetMonth() == 12 && date_beg.GetDay() == 31)
    {
      date_beg.AddYears(1);
      date_beg.SetMonth(1);
      date_beg.SetDay(1);

      mesg31 += " WARNING : When the first or the end date is Dec. 31, it";
      mesg31 += " is considered as being in the next year.";
    }

  if (date_end.GetMonth() == 12 && date_end.GetDay() == 31)
    {
      date_end.AddYears(1);
      date_end.SetMonth(1);
      date_end.SetDay(1);

      if (mesg31 == "")
        {
          mesg31 += " WARNING : When the first or the end date is Dec. 31, ";
          mesg31 += " it is considered as being in the next year.";
        }
    }

  int first_file_index = 40 + (date_beg.GetNumberOfDays() + 6) / 10;
  int last_file_index = 40 + (date_end.GetNumberOfDays() + 6) / 10;
  int file_index_tmp = first_file_index;
  int last_file_index_tmp = last_file_index;
  int NumberOfYears = date_end.GetYear() - date_beg.GetYear();

  // the last file index is artificially increased in order
  // that the preprocessing can run through year changes.
  last_file_index_tmp += NumberOfYears * 37;

  cout << "\n" << "INFORMATION : \n" << "     ||     \n";
  cout << "     ||     " << "First mozart file : h00";
  cout << to_str(first_file_index) << ".nc";
  cout << "\n";

  ////////////////////////////////////////
  // First day of the first mozart file //
  // (needed in the file polair3d-data.cfg)

  Date init_date = date_beg;

  if (first_file_index == 40)
    {
      init_date.AddYears(-1);
      init_date.SetMonth(12);
      init_date.SetDay(31);
    }
  else
    {
      int DayNumber = init_date.GetNumberOfDays();
      while ((DayNumber + 6) % 10 != 0)
        {
          init_date.AddDays(-1);
          DayNumber = init_date.GetNumberOfDays();
        }
    }

  cout << "     ||     " << "   The first day of the first Mozart file is : ";
  cout << to_str(init_date) << "\n" << "     ||     \n";

  // first day of the first mozart file //
  ////////////////////////////////////////

  cout << "     ||     " << "The year changes " << to_str(NumberOfYears);
  cout << " time(s)\n";

  cout << "     ||     " << mesg31 << "\n" << "     ||     \n";

  cout << "     ||     " << "Last Mozart file : h00";
  cout << to_str(last_file_index) << ".nc\n" << "     ||     \n";


  ///////////////////////////
  // MOZART MAP OF SPECIES //
  ///////////////////////////

  cout << "Reading, interpolating and writing concentrations..."
       << endl;
  cout << "Maps input and output species...";
  cout.flush();

  // File "species.txt" format (one line):
  // [Mozart species] [Polair species] {[ratio] {[Polair species]
  //                                                    [ratio] {...}}}
  ifstream species(species_file.c_str());

  if (!species.is_open())
    throw string("Unable to open file \"") + species_file + "\".";

  string line, species_in, species_out;

  string delim = " \t\n|";

  // To all output (Polair) species, we associate a list of input
  //(Mozart) species and their ratios (i.e. the proportion of the
  //input species to be put in the output species).
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
          species_data.push_back(line.substr(beg_index,
                                             end_index - beg_index));

          beg_index = line.find_first_not_of(delim, end_index);
        }

      // Number of words on the line (Mozart species and Polair
      //                                        species and ratios).

      int size = species_data.size();

      // If the input species is associated with any output species
      if (size > 1)
        {
          // Number of output species.
          int nb = max((size - 1) / 2, 1);
          // For all output species, update its list of input
          //                                species and ratios.
          for (i = 0; i < nb; i++)
            // If the output species is not in the map yet.
            if ((pos = species_map.find(species_data[2 * i + 1]))
                == species_map.end())
              {
                SpeciesInfo info;
                info.species_in.push_back(species_data[0]);
                if (size % 2 == 1) // If ratios are specified.
                  {
                    float spd = to_num<float>(species_data[2 * (i + 1)]);
                    info.ratio.push_back(spd);
                  }
                else  // No ratio: assumed to be 1.
                  info.ratio.push_back(1.);
                species_map.insert(make_pair(species_data[2 * i + 1], info));
              }
            else  // the output is in the map.
              {
                pos->second.species_in.push_back(species_data[0]);
                if (size % 2 == 1) // If ratios are specified.
                  {
                    float spd = to_num<float>(species_data[2 * (i + 1)]);
                    pos->second.ratio.push_back(spd);
                  }
                else  // No ratio: assumed to be 1.
                  pos->second.ratio.push_back(1.);
              }
        }

    }  // Loop over 'species_file' lines.

  species.close();


  ////////////////////////
  //     MAIN LOOP      //
  ////////////////////////

  // the main loop runs over the period if it is less than a year
  // and over a year if the period is longer than a year.

  int last_file_index_max = min(last_file_index_tmp, first_file_index + 36);

  for (file_index_tmp = first_file_index;
       file_index_tmp <= last_file_index_max; file_index_tmp++)
    {
      int file_index = file_index_tmp;

      while (file_index > 76)
        file_index -= 37;

      string mozart_file = mozart_directory + string("h00") +
        to_str(file_index) + ".nc";

      cout << endl;
      cout << "Mozart file : " << mozart_file << endl;


      //////////////////
      // MOZART GRIDS //
      //////////////////

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
      RegularGrid<real> GridT_out_half(Nt_out / 2);
      RegularGrid<real> GridZ_out(Nz_out);
      RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
      RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

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
      for (k = 0; k < Nz_out; k++)
        GridZ_out(k) = (GridZ_all_interf_out(k) +
                        GridZ_all_interf_out(k + 1)) / 2.;
      // Sets the boundary layer altitude.
      GridZ_interf_out(0) = 1.5 * GridZ_all_interf_out(Nz_out - 1)
        - 0.5 * GridZ_all_interf_out(Nz_out - 2);
      // Sets time steps.
      for (i = 0; i < GridT_out.GetLength(); i++)
        GridT_out(i) = GridT_in(i) = GridT_in_double(i);


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
      Data<real, 4> Conc_in_tmp(GridT_in, GridZ_in_tmp, GridY_in,
                                GridX_in_tmp);

      // Output fields.
      Data<real, 4> Pressure_out(GridT_out, GridZ_out, GridY_out,
                                 GridX_out);
      Data<real, 4> Temperature_out(GridT_out, GridZ_out, GridY_out,
                                    GridX_out);
      // Concentrations and temporary arrays used to sum concentrations.
      Data<real, 4> Conc_out_x(GridT_out, GridZ_out, GridY_out,
                               GridX_interf_out);
      Data<real, 4> Conc_out_x_tmp(GridT_out, GridZ_out,
                                   GridY_out, GridX_interf_out);
      Data<real, 4> Conc_out_y(GridT_out, GridZ_out, GridY_interf_out,
                               GridX_out);
      Data<real, 4> Conc_out_y_tmp(GridT_out, GridZ_out,
                                   GridY_interf_out, GridX_out);
      Data<real, 4> Conc_out_z(GridT_out, GridZ_interf_out, GridY_out,
                               GridX_out);
      Data<real, 4> Conc_out_z_tmp(GridT_out, GridZ_interf_out,
                                   GridY_out, GridX_out);

      cout << " done" << endl;
      cout << endl;


      /////////////////
      // READS INPUT //
      /////////////////

      cout << "Extracting temperature...";
      cout.flush();
      Mozart.Read(mozart_file, "T", Temperature_tmp);
      // From the temperature read in the Mozart file, builds a temperature
      // field that allows interpolation along X in an extended range between
      // -360 degrees and +357.1875 degrees.
      for (h = 0; h < Nt_in; h++)
        for (k = 0; k < Nz_in; k++)
          for (j = 0; j < Ny_in; j++)
            {
              for (i = 0; i < Nx_in_tmp; i++)
                Temperature(h, k, j, i) = Temperature_tmp(h, k, j, i);
              for (i = Nx_in_tmp; i < Nx_in; i++)
                Temperature(h, k, j, i) = Temperature_tmp(h, k, j,
                                                          i - Nx_in_tmp);
            }

      // Note: data is provided from the top level to the bottom level
      // (altitude). The interpolation function requires coordinates
      // to be sorted in increasing order. So, data must be reversed...
      Temperature.ReverseData(1);

      cout << " done" << endl;


      ///////////////////////////
      // INPUT DATA PROCESSING //
      ///////////////////////////

      cout << string("Computing level heights in meter") +
        " interpolating pressure and temperature...";
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
              SurfacePressure(h, j, i) = SurfacePressure_tmp(h, j,
                                                             i - Nx_in_tmp);
          }

      ComputePressure(alpha, beta, SurfacePressure, Pressure,
                      real(100000.));
      Pressure.ReverseData(1);

      ComputeHeight(SurfacePressure, Pressure, Temperature, GridZ_in);

      // To Polair grid...
      LinearInterpolationOneGeneral(Pressure, Pressure_out, 1);
      LinearInterpolationOneGeneral(Temperature, Temperature_out, 1);

      cout << " done" << endl;
      cout << endl;


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

      cout << " done" << endl;


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

          Conc_out_x.SetZero();
          Conc_out_y.SetZero();
          Conc_out_z.SetZero();

          // For all input species.
          for (l = 0; l < int(pos->second.species_in.size()); l++)
            {
              // Input species name.
              cout << " " << pos->second.species_in[l];
              cout.flush();
              // Reads input species concentrations.
              Mozart.Read(mozart_file, pos->second.species_in[l]
                          + "_VMR_avrg", Conc_in_tmp);

              Conc_in_tmp.ReverseData(1);

              for (h = 0; h < Nt_in; h++)
                for (k = 0; k < Nz_in; k++)
                  for (j = 0; j < Ny_in; j++)
                    {
                      for (i = 0; i < Nx_in_tmp; i++)
                        Conc_in(h, k, j, i) = Conc_in_tmp(h, k, j, i);
                      for (i = Nx_in_tmp; i < Nx_in; i++)
                        Conc_in(h, k, j, i) = Conc_in_tmp(h, k, j,
                                                          i - Nx_in_tmp);
                    }

              // Interpolation.
              LinearInterpolationOneGeneral(Conc_in, Conc_out_x_tmp, 1);
              Conc_out_x_tmp.ThresholdMin(0.);
              LinearInterpolationOneGeneral(Conc_in, Conc_out_y_tmp, 1);
              Conc_out_y_tmp.ThresholdMin(0.);
              LinearInterpolationOneGeneral(Conc_in, Conc_out_z_tmp, 1);
              Conc_out_z_tmp.ThresholdMin(0.);

              // Adds concentrations with the weighting coefficient 'ratio'.
              for (h = 0; h < Nt_out; h++)
                {
                  for (k = 0; k < Nz_out; k++)
                    for (j = 0; j < Ny_out; j++)
                      for (i = 0; i < 2; i++)
                        Conc_out_x(h, k, j, i) += pos->second.ratio[l]
                          * Conc_out_x_tmp(h, k, j, i);
                  for (k = 0; k < Nz_out; k++)
                    for (j = 0; j < 2; j++)
                      for (i = 0; i < Nx_out; i++)
                        Conc_out_y(h, k, j, i) += pos->second.ratio[l]
                          * Conc_out_y_tmp(h, k, j, i);
                  for (k = 0; k < 1; k++)
                    for (j = 0; j < Ny_out; j++)
                      for (i = 0; i < Nx_out; i++)
                        Conc_out_z(h, k, j, i) += pos->second.ratio[l]
                          * Conc_out_z_tmp(h, k, j, i);
                }
            }

          // Conversion from ppt (?) to micrograms per cubic meter.
          for (h = 0; h < Nt_out; h++)
            {
              for (k = 0; k < Nz_out; k++)
                for (j = 0; j < Ny_out; j++)
                  {
                    Conc_out_x(h, k, j, 0) *= Pressure_out(h, k, j, 0)
                      * weights_map.find(pos->first)->second
                      / R / Temperature_out(h, k, j, 0) * 1000000.;
                    Conc_out_x(h, k, j, 1) *= Pressure_out(h, k, j, Nx_out - 1)
                      * weights_map.find(pos->first)->second
                      / R / Temperature_out(h, k, j, Nx_out - 1) * 1000000.;
                  }
              for (k = 0; k < Nz_out; k++)
                for (i = 0; i < Nx_out; i++)
                  {
                    Conc_out_y(h, k, 0, i) *= Pressure_out(h, k, 0, i)
                      * weights_map.find(pos->first)->second
                      / R / Temperature_out(h, k, 0, i) * 1000000.;
                    Conc_out_y(h, k, 1, i) *= Pressure_out(h, k, Ny_out - 1, i)
                      * weights_map.find(pos->first)->second
                      / R / Temperature_out(h, k, Ny_out - 1, i) * 1000000.;
                  }
              for (j = 0; j < Ny_out; j++)
                for (i = 0; i < Nx_out; i++)
                  Conc_out_z(h, 0, j, i) *= Pressure_out(h, Nz_out - 1, j, i)
                    * weights_map.find(pos->first)->second
                    / R / Temperature_out(h, Nz_out - 1, j, i) * 1000000.;
            }

          // for the file h0040.nc, only the second half of
          // its period is taken.

          if (file_index == 40)
            {
              Data<real, 4> Conc_out_x_trunc(GridT_out_half, GridZ_out,
                                             GridY_out, GridX_interf_out);
              Data<real, 4> Conc_out_y_trunc(GridT_out_half, GridZ_out,
                                             GridY_interf_out, GridX_out);
              Data<real, 4> Conc_out_z_trunc(GridT_out_half, GridZ_interf_out,
                                             GridY_out, GridX_out);

              for (h = 0; h < Nt_out / 2; h++)
                {
                  for (k = 0; k < Nz_out; k++)
                    for (j = 0; j < Ny_out; j++)
                      for (i = 0; i < 2; i++)
                        Conc_out_x_trunc(h, k, j, i) = Conc_out_x(h + Nt_out / 2,
                                                                  k, j, i);
                  for (k = 0; k < Nz_out; k++)
                    for (j = 0; j < 2; j++)
                      for (i = 0; i < Nx_out; i++)
                        Conc_out_y_trunc(h, k, j, i) = Conc_out_y(h + Nt_out / 2,
                                                                  k, j, i);
                  for (k = 0; k < 1; k++)
                    for (j = 0; j < Ny_out; j++)
                      for (i = 0; i < Nx_out; i++)
                        Conc_out_z_trunc(h, k, j, i) = Conc_out_z(h + Nt_out / 2,
                                                                  k, j, i);
                }

              cout << endl;
              cout << "   (min x: " << Conc_out_x_trunc.GetMin() << ")";
              cout << "   (min y: " << Conc_out_y_trunc.GetMin() << ")";
              cout << "   (min z: " << Conc_out_z_trunc.GetMin() << ")";
              cout << endl;
              cout << "   (mean x: " << Conc_out_x_trunc.Mean() << ")";
              cout << "   (mean y: " << Conc_out_y_trunc.Mean() << ")";
              cout << "   (mean z: " << Conc_out_z_trunc.Mean() << ")";
              cout << endl;
              cout << "   (max x: " << Conc_out_x_trunc.GetMax() << ")";
              cout << "   (max y: " << Conc_out_y_trunc.GetMax() << ")";
              cout << "   (max z: " << Conc_out_z_trunc.GetMax() << ")";
              cout << endl;

              Polair.Append(Conc_out_x_trunc,
                            Directory_out + pos->first + "_x.bin");
              Polair.Append(Conc_out_y_trunc,
                            Directory_out + pos->first + "_y.bin");
              Polair.Append(Conc_out_z_trunc,
                            Directory_out + pos->first + "_z.bin");

            }
          else
            {
              cout << endl;
              cout << "   (min x: " << Conc_out_x.GetMin() << ")";
              cout << "   (min y: " << Conc_out_y.GetMin() << ")";
              cout << "   (min z: " << Conc_out_z.GetMin() << ")";
              cout << endl;
              cout << "   (mean x: " << Conc_out_x.Mean() << ")";
              cout << "   (mean y: " << Conc_out_y.Mean() << ")";
              cout << "   (mean z: " << Conc_out_z.Mean() << ")";
              cout << endl;
              cout << "   (max x: " << Conc_out_x.GetMax() << ")";
              cout << "   (max y: " << Conc_out_y.GetMax() << ")";
              cout << "   (max z: " << Conc_out_z.GetMax() << ")";
              cout << endl;

              Polair.Append(Conc_out_x,
                            Directory_out + pos->first + "_x.bin");
              Polair.Append(Conc_out_y,
                            Directory_out + pos->first + "_y.bin");
              Polair.Append(Conc_out_z,
                            Directory_out + pos->first + "_z.bin");
            }
        }

      cout << " done" << endl;
      cout << endl;

    }

  // END OF THE MAIN LOOP //
  //////////////////////////

  //////////////////////////////////////////
  // IF THE PERIOD IS LONGER THAN A YEAR  //
  //////////////////////////////////////////

  if (last_file_index_tmp > first_file_index + 36)

    {

      cout << " INFORMATION :\n";
      cout << "     ||     \n";
      cout << "     ||      The period is longer than a year.\n";
      cout << "     ||     \n";
      cout << "     ||      At this point, the boundary conditions have been";
      cout << " generated over a year. These\n";
      cout << "     ||      results (.bin) will be used to write the";
      cout << " boundary conditions for the rest of the period.\n";

      ////////////////////////////////////////////
      //      concatenation of full year(s)     //
      // if the period is longer than two years //

      while (last_file_index_tmp > first_file_index + 72)
        {

          cout << " --- Writing of boundary conditions over a complete";
          cout << " year --- \n" << endl;

          for (pos = species_map.begin(); pos != species_map.end(); pos++)
            {
              int Nt_year = 365;

              RegularGrid<real> GridT_year(Nt_year);

              Data<real, 4> Conc_out_x_year(GridT_year, GridZ_out,
                                            GridY_out, GridX_interf_out);
              Data<real, 4> Conc_out_y_year(GridT_year, GridZ_out,
                                            GridY_interf_out, GridX_out);
              Data<real, 4> Conc_out_z_year(GridT_year, GridZ_interf_out,
                                            GridY_out, GridX_out);

              cout << "Output : " << pos->first;

              Conc_out_x_year.SetZero();
              Conc_out_y_year.SetZero();
              Conc_out_z_year.SetZero();

              cout << "..........Reading";

              Polair.ReadSteps(Directory_out + pos->first + "_x.bin", 0,
                               Conc_out_x_year);
              Polair.ReadSteps(Directory_out + pos->first + "_y.bin", 0,
                               Conc_out_y_year);
              Polair.ReadSteps(Directory_out + pos->first + "_z.bin", 0,
                               Conc_out_z_year);

              cout << "..........Writing";

              Polair.Append(Conc_out_x_year, Directory_out + pos->first +
                            "_x.bin");
              Polair.Append(Conc_out_y_year, Directory_out + pos->first +
                            "_y.bin");
              Polair.Append(Conc_out_z_year, Directory_out + pos->first +
                            "_z.bin");

              cout << "..........Done." << endl;
            }

          cout << " Done for all species.\n" << endl;
          last_file_index_tmp -= 37;
        }

      //      concatenation of full year(s)     //
      // if the period is longer than two years //
      ////////////////////////////////////////////

      ///////////////////////////////////////////////
      //    concatenation of the outputs for the   //
      //      end of the period (non-full year)    //

      if (last_file_index_tmp > first_file_index + 36)
        {

          cout << " --- Writing of boundary conditions over the end of the";
          cout << " period ---\n" << endl;

          for (pos = species_map.begin(); pos != species_map.end(); pos++)
            {
              int Nt_res = (last_file_index_tmp - first_file_index - 36) * 10;
              RegularGrid<real> GridT_res(Nt_res);

              Data<real, 4> Conc_out_x_res(GridT_res, GridZ_out,
                                           GridY_out, GridX_interf_out);
              Data<real, 4> Conc_out_y_res(GridT_res, GridZ_out,
                                           GridY_interf_out, GridX_out);
              Data<real, 4> Conc_out_z_res(GridT_res, GridZ_interf_out,
                                           GridY_out, GridX_out);

              cout << "Output : " << pos->first;

              Conc_out_x_res.SetZero();
              Conc_out_y_res.SetZero();
              Conc_out_z_res.SetZero();

              cout << "..........Reading";

              Polair.ReadSteps(Directory_out + pos->first + "_x.bin", 0,
                               Conc_out_x_res);
              Polair.ReadSteps(Directory_out + pos->first + "_y.bin", 0,
                               Conc_out_y_res);
              Polair.ReadSteps(Directory_out + pos->first + "_z.bin", 0,
                               Conc_out_z_res);

              cout << "..........Writing";

              Polair.Append(Conc_out_x_res, Directory_out + pos->first +
                            "_x.bin");
              Polair.Append(Conc_out_y_res, Directory_out + pos->first +
                            "_y.bin");
              Polair.Append(Conc_out_z_res, Directory_out + pos->first +
                            "_z.bin");

              cout << "..........Done." << endl;
            }

          cout << " Done for all species.\n" << endl;
        }

      //    concatenation of the outputs for the    //
      //      end of the period (non-full year)     //
      ////////////////////////////////////////////////

    }

  // END of output writing for periods longer than a year //
  //////////////////////////////////////////////////////////

  END;

  return 0;

}
