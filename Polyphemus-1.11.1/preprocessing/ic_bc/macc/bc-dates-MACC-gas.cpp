// Modified version : Valentin RAFFORT, 15 May 2014

//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <string>
#include <stdlib.h>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_NETCDF

#include "Common.cxx"
using namespace Polyphemus;

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
  int Nt_in, Nz_in, Ny_in, Nx_in, Nz_in_tmp, Nlevel_in;
  real Delta_t_in;

  config.SetSection("[bc_input_domain]");
  config.PeekValue("Delta_t", "> 0", Delta_t_in);
  config.PeekValue("Nz", "> 0", Nz_in_tmp);
  config.PeekValue("Ny", "> 0", Ny_in);
  config.PeekValue("Nx", "> 0", Nx_in);
  config.PeekValue("Nlevel", "> 0", Nlevel_in);

  Nz_in = Nz_in_tmp; // +1 ;

  int Days_In_Month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; //30-31
  int Days_In_Month_Biss[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; //30-31

  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;
  string vertical_levels;

  config.SetSection("[domain]");
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
  //MOZART FILES CORRESPONDING//
  // TO FIRST AND END DATE    //
  //////////////////////////////

  Date date_beg(date_beg_str);
  Date date_end(date_end_str);

  string mesg31  = "";

  int first_file_index = date_beg.GetYear() * 10000 + date_beg.GetMonth() * 100 + date_beg.GetDay();
  int last_file_index = date_end.GetYear() * 10000 + date_end.GetMonth() * 100 + date_end.GetDay();
  int file_index_tmp = first_file_index;
  int last_file_index_tmp = last_file_index;
  int NumberOfYears = date_end.GetYear() - date_beg.GetYear();
  int NumberOfMonths = 0;
  int NumberOfDays = 0;
  int current_year = date_beg.GetYear();
  int current_month = date_beg.GetMonth();
  int biss = 0;
  int n_m = 0;

  if (current_year % 400 == 0)
    {
      biss = 1;
    }
  else if (current_year % 4 == 0 and current_year % 100 != 0)
    {
      biss = 1;
    }
  else
    {
      biss = 0;
    }
  cout << "Biss: " << biss << endl;

  if (NumberOfYears == 0)
    {
      NumberOfMonths = date_end.GetMonth() - date_beg.GetMonth() + 1;
    }
  else
    {
      NumberOfMonths = 12 - date_beg.GetMonth() + 1 + (NumberOfYears - 1) * 12 + date_end.GetMonth();
    }


  if (NumberOfMonths == 0)
    {
      NumberOfDays = date_end.GetDay() - date_beg.GetDay() + 1;
    }
  else if (biss == 1)
    {
      NumberOfDays = Days_In_Month_Biss[current_month - 1] - date_beg.GetDay();
      for (n_m = current_month;
           n_m < NumberOfMonths; n_m ++)
        {
          NumberOfDays += Days_In_Month_Biss[n_m];
        }
      NumberOfDays += date_end.GetDay() + 1;
    }
  else if (biss == 0)
    {
      NumberOfDays = Days_In_Month[current_month - 1] - date_beg.GetDay();
      for (n_m = current_month;
           n_m < NumberOfMonths; n_m ++)
        {
          NumberOfDays += Days_In_Month[n_m];
        }
      NumberOfDays += date_end.GetDay() + 1;
    }


  cout << "\n" << "INFORMATION : \n" << "     ||     \n";
  cout << "     ||     " << "First MACC file : GRG-";
  cout << to_str(first_file_index) << ".nc";
  cout << "\n";

  ////////////////////////////////////////
  // First day of the first mozart file //
  // (needed in the file polair3d-data.cfg)

  Date init_date = date_beg;

  cout << "     ||     " << "   The first day of the first MACC file is : ";
  cout << to_str(init_date) << "\n" << "     ||     \n";

  // first day of the first mozart file //
  ////////////////////////////////////////

  cout << "     ||     " << "The year changes " << to_str(NumberOfYears);
  cout << " time(s)\n";
  cout << "     ||     " << "Number of days " << to_str(NumberOfDays) << "\n";
  cout << "     ||     " << mesg31 << "\n" << "     ||     \n";
  cout << "     ||     " << "Last MACC file : GRG-";
  cout << to_str(last_file_index) << ".nc\n" << "     ||     \n";

  // Building the array with indices of needed Mozart files
  int *index_mozart_files = new int[NumberOfDays];
  int year_beg = date_beg.GetYear();
  int month_beg = date_beg.GetMonth();
  int day_beg = date_beg.GetDay();

  int first_day = year_beg * 10000 + month_beg * 100 + day_beg;
  int first_month = year_beg * 100 + month_beg;
  int day_tmp = day_beg;
  int month_tmp = month_beg;
  int year_tmp = year_beg;

  index_mozart_files[0] = first_day;
  for (int i = 1 ; i < NumberOfDays; i++)
    {
      if (biss == 1)
        {
          if (day_tmp == Days_In_Month_Biss[month_tmp - 1])
            {
              year_tmp = (index_mozart_files[i - 1] / 10000);
              month_tmp = (index_mozart_files[i - 1] / 100) - (index_mozart_files[i - 1] / 10000) * 100 + 1;
              day_tmp = 1;
              index_mozart_files[i] = year_tmp * 10000 + month_tmp * 100 + day_tmp;
            }
          else
            {
              year_tmp = (index_mozart_files[i - 1] / 10000);
              month_tmp = (index_mozart_files[i - 1] / 100) - (index_mozart_files[i - 1] / 10000) * 100;
              day_tmp = index_mozart_files[i - 1] - (index_mozart_files[i - 1] / 100) * 100 + 1;
              index_mozart_files[i] = year_tmp * 10000 + month_tmp * 100 + day_tmp;
            }
        }
      else
        {
          if (day_tmp == Days_In_Month[month_tmp - 1])
            {
              year_tmp = (index_mozart_files[i - 1] / 10000);
              month_tmp = (index_mozart_files[i - 1] / 100) - (index_mozart_files[i - 1] / 10000) * 100 + 1;
              day_tmp = 1;
              index_mozart_files[i] = year_tmp * 10000 + month_tmp * 100 + day_tmp;
            }
          else
            {
              year_tmp = (index_mozart_files[i - 1] / 10000);
              month_tmp = (index_mozart_files[i - 1] / 100) - (index_mozart_files[i - 1] / 10000) * 100;
              day_tmp = index_mozart_files[i - 1] - (index_mozart_files[i - 1] / 100) * 100 + 1;
              index_mozart_files[i] = year_tmp * 10000 + month_tmp * 100 + day_tmp;
            }
        }
    }

  for (int i = 1 ; i < NumberOfDays; i++)
    {
      cout << "index : " << index_mozart_files[i] << endl;
    }

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

  int last_file_index_max = last_file_index_tmp;
  //Adap the main loop
  for (file_index_tmp = 0;
       file_index_tmp < NumberOfDays; file_index_tmp++)
    {
      int file_index = index_mozart_files[file_index_tmp];
      string mozart_file = mozart_directory + string("GRG_") +
        to_str(file_index) + "_EU_AQ.nc";
      cout << endl;
      cout << "MACC file : " << mozart_file << endl;
      string temperature_file = mozart_directory + string("MACC_") +
        to_str(file_index) + "_temp_EU_AQ.nc";
      string pressure_file = mozart_directory + string("PS_") +
        to_str(file_index) + "_EU_AQ.nc";
      cout << endl;
      cout << "file index : " << file_index << endl;
      cout << "NbFile index : " << file_index_tmp << endl;
      cout << "Number of Days: " << NumberOfDays << endl;

      day_tmp = index_mozart_files[file_index_tmp] -
        (index_mozart_files[file_index_tmp] / 100) * 100;
      month_tmp = (index_mozart_files[file_index_tmp] / 100) -
        (index_mozart_files[file_index_tmp] / 10000) * 100;
      year_tmp = index_mozart_files[file_index_tmp] / 10000;
      Nt_in = 8;


      cout << "Nt_in = " + to_str(Nt_in) << endl;
      Nt_out = Nt_in;
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
      GeneralGrid<real, 4> GridZ_in_tmp(shape(Nt_in, Nz_in_tmp, Ny_in, Nx_in),
                                        1, shape(0, 1, 2, 3));
      RegularGrid<real> GridY_in(Ny_in);
      RegularGrid<real> GridX_in(Nx_in);
      RegularGrid<real> GridZ_surf(1);

      // Grids are shared.
      GridZ_in.SetVariable(1);
      GridZ_in.SetDuplicate(false);
      GridZ_in_tmp.SetVariable(1);
      GridZ_in_tmp.SetDuplicate(false);

      // Reads values.
      FormatNetCDF<float> Mozart;
      Mozart.Read(mozart_file, "lat", GridY_in);
      Mozart.Read(mozart_file, "lon", GridX_in);
      Mozart.Read(mozart_file, "time", GridT_in_double);

      // Modifies the grid to allow interpolation between -360 degrees and
      // +357.1875 degrees.
      cout << "Conversion from <0 to 360> to <-180 to 180> ... ";
      cout.flush();
      for (i = 0; i < Nx_in; i++)
        {
          if (GridX_in(i) >= 180 && GridX_in(i) <= 360)
            {
              GridX_in(i) = GridX_in(i) - 360;
            }
        }
      cout << " done" << endl;

      // Output settings.

      FormatBinary<float> Polairx;
      FormatBinary<float> Polairy;
      FormatBinary<float> Polairz;


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
        {
          GridT_out(i) = GridT_in(i) = GridT_in_double(i);
        }

      //////////
      // DATA //
      //////////

      // Input fields.
      // Fixed pressure levels.
      Data<real, 1> alpha(Nlevel_in), beta(Nlevel_in);
      Data<real, 4> SurfacePressure(GridT_in, GridZ_surf, GridY_in, GridX_in);
      Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
      Data<real, 4> Temperature(GridT_in, GridZ_in, GridY_in, GridX_in);
      Data<real, 4> Temperature_tmp(GridT_in, GridZ_in_tmp, GridY_in, GridX_in);
      Data<real, 4> Pressure_tmp(GridT_in, GridZ_in_tmp, GridY_in, GridX_in);
      Data<real, 4> Conc_in(GridT_in, GridZ_in, GridY_in, GridX_in);
      Data<real, 4> Conc_in_tmp(GridT_in, GridZ_in_tmp, GridY_in, GridX_in);


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

      ///////////////////////////
      // INPUT DATA PROCESSING //
      ///////////////////////////

      cout << string("Computing level heights in meter") +
        " interpolating pressure and temperature...";
      cout.flush();

      // Hybrid coefficients.
      Mozart.Read(pressure_file, "hyam", alpha);
      Mozart.Read(pressure_file, "hybm", beta);
      Mozart.Read(pressure_file, "var152", SurfacePressure);
      Mozart.Read(temperature_file, "var130", Temperature_tmp);
      alpha.ReverseData(0);
      beta.ReverseData(0);
      Temperature_tmp.ReverseData(1);

      Compute4DPressure(alpha, beta, SurfacePressure, Pressure_tmp, real(100000.));
      Compute4DHeight(SurfacePressure, Pressure_tmp, Temperature_tmp, GridZ_in_tmp);


      // To Polair grid...
      cout << "New GridZ_in, Temperature and Pressure" << endl;
      for (h = 0; h < Nt_in; h++)
        {
          for (j = 0; j < Ny_in; j++)
            {
              for (i = 0; i < Nx_in; i++)
                {
                  GridZ_in.Value(h, 0, j, i) = GridZ_out(0);
                  Temperature(h, 0, j, i) = Temperature_tmp(h, 0, j, i);
                  Pressure(h, 0, j, i) = Pressure_tmp(h, 0, j, i);
                  for (k = 1; k < Nz_in; k++)
                    {
                      GridZ_in.Value(h, k, j, i) = GridZ_in_tmp.Value(h, k - 1, j, i);
                      Temperature(h, k, j, i) = Temperature_tmp(h, k - 1, j, i);
                      Pressure(h, k, j, i) = Pressure_tmp(h, k - 1, j, i);
                    }
                }
            }
        }
      cout << "linear interpolation on Polair Grid ..." << endl;
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
          cout << "Output: " << pos->first << "; input:" << "; weights :" << weights_map.find(pos->first)->second ;
          cout.flush(); //changer ici

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
              Mozart.Read(mozart_file, pos->second.species_in[l], Conc_in_tmp);

              Conc_in_tmp.ReverseData(1);
              for (h = 0; h < Nt_in; h++)
                {
                  for (j = 0; j < Ny_in; j++)
                    {
                      for (i = 0; i < Nx_in; i++)
                        {
                          Conc_in(h, 0, j, i) = Conc_in_tmp(h, 0, j, i);
                          for (k = 1; k < Nz_in; k++)
                            {
                              Conc_in(h, k, j, i) = Conc_in_tmp(h, k - 1, j, i);
                            }
                        }
                    }
                }

              // Interpolation.
              LinearInterpolationOneGeneral(Conc_in, Conc_out_x_tmp, 1);
              Conc_out_x_tmp.ThresholdMin(0.);
              LinearInterpolationOneGeneral(Conc_in, Conc_out_y_tmp, 1);
              Conc_out_y_tmp.ThresholdMin(0.);
              LinearInterpolationOneGeneral(Conc_in, Conc_out_z_tmp, 1);
              Conc_out_z_tmp.ThresholdMin(0.);

              // Adds concentrations with the weighting coefficient 'ratio'
              for (h = 0; h < Nt_out; h++)
                {
                  for (k = 0; k < Nz_out; k++)
                    {
                      for (j = 0; j < Ny_out; j++)
                        {
                          for (i = 0; i < 2; i++)
                            {
                              Conc_out_x(h, k, j, i) += pos->second.ratio[l]
                                * Conc_out_x_tmp(h, k, j, i);
                            }//i
                        }//j
                    }//k
                  for (k = 0; k < Nz_out; k++)
                    {
                      for (j = 0; j < 2; j++)
                        {
                          for (i = 0; i < Nx_out; i++)
                            {
                              Conc_out_y(h, k, j, i) += pos->second.ratio[l]
                                * Conc_out_y_tmp(h, k, j, i);
                            }//i
                        }//j
                    }//k
                  for (k = 0; k < 1; k++)
                    {
                      for (j = 0; j < Ny_out; j++)
                        {
                          for (i = 0; i < Nx_out; i++)
                            {
                              Conc_out_z(h, k, j, i) += pos->second.ratio[l]
                                * Conc_out_z_tmp(h, k, j, i);
                            }//i
                        }//j
                    }//k
                }//h
            }//l

          // Conversion to µg/m3.
          for (h = 0; h < Nt_out; h++)
            {
              for (k = 0; k < Nz_out; k++)
                {
                  for (j = 0; j < Ny_out; j++)
                    {
                      Conc_out_x(h, k, j, 0) *= Pressure_out(h, k, j, 0)
                        * weights_map.find(pos->first)->second
                        / R / Temperature_out(h, k, j, 0) * 1000000.;
                      Conc_out_x(h, k, j, 1) *= Pressure_out(h, k, j, Nx_out - 1)
                        * weights_map.find(pos->first)->second
                        / R / Temperature_out(h, k, j, Nx_out - 1) * 1000000.;
                    }//j
                }//k
              for (k = 0; k < Nz_out; k++)
                {
                  for (i = 0; i < Nx_out; i++)
                    {
                      Conc_out_y(h, k, 0, i) *= Pressure_out(h, k, 0, i)
                        * weights_map.find(pos->first)->second
                        / R / Temperature_out(h, k, 0, i) * 1000000.;
                      Conc_out_y(h, k, 1, i) *= Pressure_out(h, k, Ny_out - 1, i)
                        * weights_map.find(pos->first)->second
                        / R / Temperature_out(h, k, Ny_out - 1, i) * 1000000.;
                    }//i
                }//k
              for (j = 0; j < Ny_out; j++)
                {
                  for (i = 0; i < Nx_out; i++)
                    {
                      Conc_out_z(h, 0, j, i) *= Pressure_out(h, Nz_out - 1, j, i)
                        * weights_map.find(pos->first)->second
                        / R / Temperature_out(h, Nz_out - 1, j, i) * 1000000.;
                    }//i
                }//j
            }//h


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

          Polairx.Append(Conc_out_x,
                         Directory_out + pos->first + "_x.bin");
          Polairy.Append(Conc_out_y,
                         Directory_out + pos->first + "_y.bin");
          Polairz.Append(Conc_out_z,
                         Directory_out + pos->first + "_z.bin");
        }

      cout << " done" << endl;
      cout << endl;

    }

  // END OF THE MAIN LOOP //
  //////////////////////////

  END;

  return 0;

}
