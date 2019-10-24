// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet and Vincent Picavet
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

#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

#include "fastJX.hxx"


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("MM5-meteo.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);

  Date date_prev = date_beg;
  date_prev.AddDays(-1);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  int h, i, j, k;

  // Constants.
  const real pi = 3.14159265358979323846264;
  const real r = 287.;
  const real cp = 1005.;
  const real g = 9.81;

  //// OUTPUT INFORMATION. ////

  // Output dimensions.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  // Output steps.
  real Delta_t_out, Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;
  string vertical_levels;
  // Output directory.
  string directory_out;

  //// INPUT INFORMATION. ////

  // Input dimensions for meteological data.
  int Nt_in, Nz_in, Ny_in, Nx_in;
  // Input-data steps.
  real Delta_t_in, Delta_y_in, Delta_x_in;
  real y_min_in, x_min_in;
  // Input files.
  string file_in, file_in_prev = "";
  // Input domain's projection type.
  // (=1 for Lambert-CC, =2 for Stereographic,
  //  =3 for Mercator, =6 for lat/lon (not implemented yet).
  int projection_type;

  //// PARAMETERIZATION. ////

  // Accumulated rain.
  bool prev_accumulated_rain;

  // Photolysis.
  int photolysis_tabulation;
  bool with_meteo, with_photolysis, ice_cloud;
  int attenuation_type;
  vector<string> photolysis_reaction_list;
  int Nr_photolysis, Nwavelength, Nlegendre;

  // Kz thresholds.
  real Kz_min, Kz_min_urban;
  bool apply_vert;
  real Kz_max;

  // Clouds.
  real min_height;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  //// Domain. ////

  cout << "  + Reading domain information...";
  cout.flush();

  // Input domain.
  configuration.SetSection("[MM5]");

  configuration.PeekValue("Nt", "> 0", Nt_in);
  configuration.PeekValue("Nz", "> 0", Nz_in);
  configuration.PeekValue("Ny", "> 0", Ny_in);
  configuration.PeekValue("Nx", "> 0", Nx_in);

  configuration.PeekValue("Delta_t", "> 0", Delta_t_in);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_in);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_in);

  configuration.PeekValue("y_min", y_min_in);
  configuration.PeekValue("x_min", x_min_in);

  configuration.PeekValue("projection_type", "= 1 2 3", projection_type);

  // Output domain.
  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx_out);
  configuration.PeekValue("Ny", "> 0", Ny_out);
  configuration.PeekValue("Nz", "> 0", Nz_out);
  configuration.PeekValue("Vertical_levels", vertical_levels);

  configuration.PeekValue("Delta_t", "> 0", Delta_t_out);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_out);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_out);

  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);

  Nt_out = compute_Nt(date_beg, date_end, Delta_t_out);

  date_end.AddSeconds(- Delta_t_out * 3600);

  // Paths.
  configuration.SetSection("[paths]");
  configuration.PeekValue("Directory_meteo", directory_out);
  configuration.PeekValue("Database_MM5-meteo", file_in);

  cout << endl;

  //// Input Files. ////

  cout << "  + Reading input files paths and names...";
  cout.flush();

  // Processes the input file name to prepare date expansion.
  file_in = find_replace(file_in, "&D", "%y-%m-%d");
  file_in = find_replace(file_in, "&y", "%y");
  file_in = find_replace(file_in, "&m", "%m");
  file_in = find_replace(file_in, "&d", "%d");
  file_in_prev = date_prev.GetDate(file_in);
  file_in = date_beg.GetDate(file_in);

  if (!exists(file_in))
    throw string("Unable to find MM5 file \"") + file_in + "\".";

  cout << endl;

  //// Parameterizations. ////

  // Land use.
  string LUC_file;
  int Nc, urban_index;
  configuration.PeekValue("LUC_file", LUC_file);
  if (!exists(LUC_file))
    throw "Unable to open land use cover file \"" + LUC_file + "\".";
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny_out * Nx_out);
  configuration.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                          urban_index);

  cout << "  + Reading physical parameterizations...";
  cout.flush();

  // Meteo.
  configuration.SetSection("[meteo]");
  configuration.PeekValue("Compute_Meteo", with_meteo);

  // Accumulated rain.
  configuration.SetSection("[accumulated_rain]");
  configuration.PeekValue("Prev_accumulated_rain", prev_accumulated_rain);

  // Photolysis rates.
  configuration.SetSection("[photolysis_rates]");
  configuration.PeekValue("Compute_Photolysis_Data", with_photolysis);
  string directory_attenuation, directory_photolysis, fastJ_parameter_files;
  Array<char, 2> photolysis_specie_name;
  if (with_photolysis)
    {
      configuration.PeekValue("Directory_attenuation", directory_attenuation);
      configuration.PeekValue("Directory_photolysis_rates", directory_photolysis);
      configuration.PeekValue("FastJ_parameter_files", fastJ_parameter_files);
      configuration.PeekValue("Photolysis_tabulation_option", photolysis_tabulation);
      configuration.PeekValue("Attenuation_Type", "= 1 2", attenuation_type);
      configuration.PeekValue("Ice_cloud", ice_cloud);
      configuration.Find("Species");
      split(configuration.GetLine(), photolysis_reaction_list);
      Nr_photolysis = int(photolysis_reaction_list.size());
      photolysis_specie_name.resize(Nr_photolysis, 10);

      photolysis_specie_name = ' ';
      for (int i = 0; i < Nr_photolysis; i++)
        {
          int Nr_photolysis_reaction = photolysis_reaction_list[i].size();
          for (int j = 0; j < Nr_photolysis_reaction; j++)
            photolysis_specie_name(i, j) = photolysis_reaction_list[i][j];
        }
    }

  // Kz.
  configuration.SetSection("[Kz]");
  configuration.PeekValue("Min", "positive", Kz_min);
  configuration.PeekValue("Min_urban", "positive", Kz_min_urban);
  configuration.PeekValue("Max", "positive", Kz_max);
  configuration.PeekValue("Apply_vert", apply_vert);

  // Clouds.
  configuration.SetSection("[clouds]");
  configuration.PeekValue("Min_height", "positive", min_height);

  cout << endl;

  Date date_beg_meteo = read_date_MM5(file_in);
  Date date_end_meteo = date_beg_meteo;
  date_end_meteo.AddSeconds(Nt_in * Delta_t_in * 3600);

  if (date_beg_meteo > date_beg || date_end_meteo < date_end)
    throw string("The MM5 file \"") + file_in
      + "\" does not contain data for all dates.";

  real t_min_in, t_min_out, differences;
  t_min_in = real(date_beg_meteo.GetHour())
    + real(date_beg_meteo.GetMinutes()) / 60.
    + real(date_beg_meteo.GetSeconds()) / 3600.;
  differences = real(date_beg.GetSecondsFrom(date_beg_meteo));
  t_min_out = t_min_in + differences / 3600.;

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for grids..." << endl;
  cout.flush();

  // Input settings.

  // Input grids.
  RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
  // Land use categories.
  RegularGrid<real> GridC(Nc);
  // Grid for Z in height. Z depends on Z, Y, X.
  GeneralGrid<real, 3> GridZ_in(shape(Nz_in, Ny_in - 1, Nx_in - 1),
                                1, shape(1, 2, 3));

  // Vertical levels are shared.
  GridZ_in.SetVariable(1);
  GridZ_in.SetDuplicate(false);

  // Grid for Z in height (interfaces). Z depends on Z, Y, X.
  GeneralGrid<real, 3> GridZ_interf_in(shape(Nz_in + 1, Ny_in - 1, Nx_in - 1),
                                       1, shape(1, 2, 3));


  // Vertical levels for dot grids are not the same.
  GeneralGrid<real, 3> GridZ_Dot_in(shape(Nz_in, Ny_in, Nx_in),
                                    1, shape(1, 2, 3));
  GridZ_Dot_in.SetVariable(1);
  GridZ_Dot_in.SetDuplicate(false);


  // MM5 Input Grid for Z in height. Z depends on Z, X, Y
  GeneralGrid<real, 3> GridZ_MM5_in(shape(Nz_in, Nx_in - 1, Ny_in - 1),
                                    1, shape(1, 2, 3));
  // Vertical levels are shared.
  GridZ_MM5_in.SetVariable(1);
  GridZ_MM5_in.SetDuplicate(false);

  // Vertical levels for dot grids are not the same.
  GeneralGrid<real, 3> GridZ_MM5_Dot_in(shape(Nz_in, Nx_in, Ny_in),
                                        1, shape(1, 2, 3));
  GridZ_MM5_Dot_in.SetVariable(1);
  GridZ_MM5_Dot_in.SetDuplicate(false);

  RegularGrid<real> GridY_in(y_min_in, Delta_y_in, Ny_in - 1);
  RegularGrid<real> GridX_in(x_min_in, Delta_x_in, Nx_in - 1);

  // Interfaces.
  RegularGrid<real> GridX_interf_in(x_min_in - Delta_x_in / 2.,
                                    Delta_x_in, Nx_in);
  RegularGrid<real> GridY_interf_in(y_min_in - Delta_y_in / 2.,
                                    Delta_y_in, Ny_in);

  // Output settings.
  cout << "  + Output grids..." << endl;
  cout.flush();
  // Output grids.
  RegularGrid<real> GridT_out(t_min_out, Delta_t_out, Nt_out);
  RegularGrid<real> GridZ_out(Nz_out);
  // Data may be provided on interfaces.
  RegularGrid<real> GridZ_interf_out(Nz_out + 1);

  // Latlon
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);
  RegularGrid<real> GridY_interf_out(y_min_out - Delta_y_out / 2.,
                                     Delta_y_out, Ny_out + 1);
  RegularGrid<real> GridX_interf_out(x_min_out - Delta_x_out / 2.,
                                     Delta_x_out, Nx_out + 1);

  // Lcc.
  GeneralGrid<real, 2> GridX_2D_out(shape(Ny_out, Nx_out), 2, shape(1, 2));
  GeneralGrid<real, 2> GridY_2D_out(shape(Ny_out, Nx_out), 1, shape(1, 2));
  GeneralGrid<real, 2> GridX_3D_out(shape(Ny_out, Nx_out), 3, shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_out(shape(Ny_out, Nx_out), 2, shape(2, 3));
  GeneralGrid<real, 2> GridX_3D_Gen_out(shape(Ny_out + 1, Nx_out),
                                        3, shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_Gen_out(shape(Ny_out, Nx_out + 1),
                                        2, shape(2, 3));
  GeneralGrid<real, 2> GridX_3D_interf_out(shape(Ny_out, Nx_out + 1),
                                           3, shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_interf_out(shape(Ny_out + 1, Nx_out),
                                           2, shape(2, 3));

  GridX_2D_out.SetVariable(2);
  GridX_2D_out.SetDuplicate(false);
  GridY_2D_out.SetVariable(1);
  GridY_2D_out.SetDuplicate(false);
  GridX_3D_out.SetVariable(3);
  GridX_3D_out.SetDuplicate(false);
  GridY_3D_out.SetVariable(2);
  GridY_3D_out.SetDuplicate(false);
  GridX_3D_Gen_out.SetVariable(3);
  GridX_3D_Gen_out.SetDuplicate(false);
  GridY_3D_Gen_out.SetVariable(2);
  GridY_3D_Gen_out.SetDuplicate(false);
  GridX_3D_interf_out.SetVariable(3);
  GridX_3D_interf_out.SetDuplicate(false);
  GridY_3D_interf_out.SetVariable(2);
  GridY_3D_interf_out.SetDuplicate(false);

  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels, GridZ_interf_out);
  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_interf_out(k) + GridZ_interf_out(k + 1)) / 2.0;

  // Reads land use data.
  Data<real, 3> LUC(GridC, GridY_out, GridX_out);
  FormatBinary<float>().Read(LUC_file, LUC);

  // Photolysis rates Grid.
  Nwavelength = 4;
  Nlegendre = 8;
  RegularGrid<real> GridP(Nr_photolysis);
  RegularGrid<real> GridWavelength(Nwavelength);
  RegularGrid<real> GridLegendre(Nlegendre);

  cout << " done." << endl;


  /////////////////
  // OUTPUT DATA //
  /////////////////

  cout << "Memory allocation for output data fields...";
  cout.flush();

  // 3D Output fields.
  // Cross fields.
  Data<real, 4> PotentialTemperature_out(GridT_out, GridZ_out,
                                         GridY_out, GridX_out);
  Data<real, 4> Temperature_out(GridT_out, GridZ_out,
                                GridY_3D_out, GridX_3D_out);
  Data<real, 4> Richardson_out(GridT_out, GridZ_out,
                               GridY_3D_out, GridX_3D_out);
  Data<real, 4> Pressure_out(GridT_out, GridZ_out,
                             GridY_3D_out, GridX_3D_out);
  Data<real, 4> SpecificHumidity_out(GridT_out, GridZ_out,
                                     GridY_3D_out, GridX_3D_out);
  Data<real, 4> LiquidWaterContent_out(GridT_out, GridZ_out,
                                       GridY_3D_out, GridX_3D_out);
  // Photolysis cross fields.
  Data<real, 4> Attenuation_out(GridT_out, GridZ_out,
                                GridY_3D_out, GridX_3D_out);
  Data<real, 4> IceWaterContent_out(GridT_out, GridZ_out,
                                    GridY_3D_out, GridX_3D_out);
  Data<real, 4> WindModule_out(GridT_out, GridZ_out,
                               GridY_out, GridX_out);
  Data<real, 4> CloudFraction_out(GridT_out, GridZ_out,
                                  GridY_3D_out, GridX_3D_out);
  Data<real, 4> LiquidWaterExtinction_out(GridT_out, GridZ_out,
                                          GridY_3D_out, GridX_3D_out);
  Data<real, 4> IceWaterExtinction_out(GridT_out, GridZ_out,
                                       GridY_3D_out, GridX_3D_out);
  // Dot fields.
  // Fields defined on X- or Y-interfaces.
  Data<real, 4> MeridionalWind_out(GridT_out, GridZ_out,
                                   GridY_3D_interf_out, GridX_3D_Gen_out);
  Data<real, 4> ZonalWind_out(GridT_out, GridZ_out,
                              GridY_3D_Gen_out, GridX_3D_interf_out);
  // Fields defined on Z_interface.
  Data<real, 4> Kz_out(GridT_out, GridZ_interf_out,
                       GridY_out, GridX_out);

  // 2D Output fields.
  Data<real, 3> SurfacePressure_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfaceTemperature_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SkinTemperature_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfaceRichardson_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> BoundaryHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> CloudBaseHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FrictionModule_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SolarRadiation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> ConvectiveRain_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> Rain_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SoilWater_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SensibleHeat_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> Evaporation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FirstLevelWindModule_out(GridT_out, GridY_2D_out,
                                         GridX_2D_out);
  Data<real, 3> WindModule10_out(GridT_out, GridY_2D_out, GridX_2D_out);

  // Post interpolation computed fields.
  Data<real, 3> PARdb_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PARdiff_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PAR_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SurfacePotentialTemperature_out(GridT_out,
                                                GridY_out, GridX_out);
  // Photolysis.
  Data<real, 4> OpticalDepth_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> IceOpticalDepth_out(GridT_out, GridZ_out, GridY_out,
                                    GridX_out);
  Data<real, 5> PhotolysisRate_out(GridP, GridT_out, GridZ_out,
                                   GridY_out, GridX_out);

  cout << " done." << endl;


  ////////////////////////////////
  // SIGMA TO HEIGHT CONVERSION //
  ////////////////////////////////

  cout << "Conversion from sigma levels to heights..." << endl;
  cout << "  + Computing sigma-levels and terrain elevation data..." << endl;
  cout.flush();

  FormatMM5 InputMeteo;

  // Data to convert z from sigma to heights.
  Data<real, 2> SigmaH(Nt_in, Nz_in);
  // Terrain is 2D Field read from MM5 file (-> XY order).
  Data<real, 3> Terrain(GridT_in, GridX_in, GridY_in);
  Data<real, 3> Terrain_Dot(GridT_in, GridY_interf_in, GridX_interf_in);

  // Reads data to convert Z from sigma levels to heights.
  InputMeteo.ReadWholeField(file_in, "SIGMAH", SigmaH);
  InputMeteo.ReadWholeField(file_in, "TERRAIN", Terrain);
  Terrain.SwitchDimensions(shape(0, 2, 1), GridT_in, GridY_in, GridX_in);

  // Converts Z from sigma levels to heights at cross points.

  // Big header data and general constants.
  Array<int, 2> BHI;
  Array<float, 2> BHR;
  Array<string, 2> BHIC;
  Array<string, 2> BHRC;
  real p00, ptop, Ts0, A, Tiso, Piso, Ziso, aux;
  real ps0, refpress;

  InputMeteo.ReadBigHeader(file_in, BHI, BHR, BHIC, BHRC);

  p00 = BHR(4, 1);
  ptop = BHR(1, 1);
  Ts0 = BHR(4, 2);
  A = BHR(4, 3);
  Tiso = BHR(4, 4);
  Piso = p00 * exp((Tiso - Ts0) / A);
  aux = log(Piso / p00);
  Ziso = - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g;


  // At cross points.
  cout << "  + Computing vertical levels in meter at cross points...";
  cout.flush();

  for (k = 0; k < SigmaH.GetLength(1); k++)
    for (j = 0; j < Ny_in - 1; j++)
      for (i = 0; i < Nx_in - 1; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_in.Value(0, SigmaH.GetLength(1) - k - 1, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_in.Value(0, SigmaH.GetLength(1) - k - 1, j, i)
                = Ziso - r * Tiso * aux / g - Terrain(0, j, i);
            }
        }
  cout << endl;

  // Interfaces.
  for (j = 0; j < Ny_in - 1; j++)
    for (i = 0; i < Nx_in - 1; i++)
      GridZ_interf_in.Value(0, 0, j, i) = 0.;
  for (k = 1; k < SigmaH.GetLength(1); k++)
    for (j = 0; j < Ny_in - 1; j++)
      for (i = 0; i < Nx_in - 1; i++)
        GridZ_interf_in.Value(0, k, j, i) = 2 * GridZ_in.Value(0, k - 1, j, i)
          - GridZ_interf_in.Value(0, k - 1, j, i);

  // Interpolates Terrain at dot points.
  LinearInterpolationRegular(Terrain, Terrain_Dot);

  // Converts z from sigma levels to heights at dot points.
  for (k = 0; k < SigmaH.GetLength(1); k++)
    for (j = 0; j < Ny_in; j++)
      for (i = 0; i < Nx_in; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain_Dot(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_Dot_in.Value(0, SigmaH.GetLength(1) - k - 1, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain_Dot(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_Dot_in.Value(0, SigmaH.GetLength(1) - k - 1, j, i)
                = Ziso - r * Tiso * aux / g - Terrain_Dot(0, j, i);
            }
        }

  cout << " done." << endl;



  ///////////////////////////////////////////////////////
  // Converts from lat/lon output grids to MM5 indices //
  ///////////////////////////////////////////////////////


  cout << "Converting from latlon to MM5 indices...";
  cout.flush();

  if (projection_type == 1)
    // Converts from latitude/longitude to MM5 indices in Lambert conformal
    // conic projection.
    {
      LonlatToMM5LccInd<float> Lcc(BHI(0, 5), BHI(0, 4), BHR(0, 10), BHR(0, 9),
                                   BHR(0, 1), BHR(0, 2), BHR(0, 4), BHR(0, 5),
                                   BHR(0, 0), BHI(0, 19));

      // 3D output grids.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          Lcc(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
              GridX_3D_out.Value(0, 0, j, i), GridY_3D_out.Value(0, 0, j, i));

      // 2D output grids.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            GridX_2D_out.Value(0, j, i) = GridX_3D_out.Value(0, 0, j, i);
            GridY_2D_out.Value(0, j, i) = GridY_3D_out.Value(0, 0, j, i);
          }

      // 3D output grids on interfaces.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out + 1; i++)
          Lcc(x_min_out + Delta_x_out * (real(i) - 0.5),
              y_min_out + Delta_y_out * j,
              GridX_3D_interf_out.Value(0, 0, j, i),
              GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          Lcc(x_min_out + Delta_x_out * i,
              y_min_out + Delta_y_out * (real(j) - 0.5),
              GridX_3D_Gen_out.Value(0, 0, j, i),
              GridY_3D_interf_out.Value(0, 0, j, i));
    }
  else if (projection_type == 2)
    // Converts from latitude/longitude to MM5 indices in Mercator projection.
    {
      // Mercator projection.
      LonlatToMM5MercInd<float> Merc(BHI(0, 5), BHI(0, 4), BHR(0, 10), BHR(0, 9),
                                     BHR(0, 1), BHR(0, 2), BHR(0, 4),
                                     BHR(0, 0), BHI(0, 19));

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          Merc(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
               GridX_3D_out.Value(0, 0, j, i),
               GridY_3D_out.Value(0, 0, j, i));

      // 2D output grids.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            GridX_2D_out.Value(0, j, i) = GridX_3D_out.Value(0, 0, j, i);
            GridY_2D_out.Value(0, j, i) = GridY_3D_out.Value(0, 0, j, i);
          }

      // 3D output grids on interfaces.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out + 1; i++)
          Merc(x_min_out + Delta_x_out * (real(i) - 0.5),
               y_min_out + Delta_y_out * j,
               GridX_3D_interf_out.Value(0, 0, j, i),
               GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          Merc(x_min_out + Delta_x_out * i,
               y_min_out + Delta_y_out * (real(j) - 0.5),
               GridX_3D_Gen_out.Value(0, 0, j, i),
               GridY_3D_interf_out.Value(0, 0, j, i));
    }
  else if (projection_type == 3)
    // Converts from latitude/longitude to MM5 indices in stereographic
    // projection.
    {
      LonlatToMM5StereInd<float> Stereo(BHI(0, 5), BHI(0, 4), BHR(0, 10),
                                        BHR(0, 9), BHR(0, 1), BHR(0, 2),
                                        BHR(0, 4), BHR(0, 0), BHI(0, 19));

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          Stereo(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
                 GridX_3D_out.Value(0, 0, j, i),
                 GridY_3D_out.Value(0, 0, j, i));

      // 2D output grids.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            GridX_2D_out.Value(0, j, i) = GridX_3D_out.Value(0, 0, j, i);
            GridY_2D_out.Value(0, j, i) = GridY_3D_out.Value(0, 0, j, i);
          }

      // 3D output grids on interfaces.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out + 1; i++)
          Stereo(x_min_out + Delta_x_out * (real(i) - 0.5),
                 y_min_out + Delta_y_out * j,
                 GridX_3D_interf_out.Value(0, 0, j, i),
                 GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          Stereo(x_min_out + Delta_x_out * i,
                 y_min_out + Delta_y_out * (real(j) - 0.5),
                 GridX_3D_Gen_out.Value(0, 0, j, i),
                 GridY_3D_interf_out.Value(0, 0, j, i));
    }

  cout << " done." << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  cout << "Applying transformation to read fields...";
  cout.flush();

  // MM5 3D fields.
  Data<real, 4> PressurePerturbation(GridT_in, GridZ_MM5_in,
                                     GridX_in, GridY_in);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> ReferencePressure(GridT_in, GridX_in, GridY_in);

  InputMeteo.ReadWholeField(file_in, "PSTARCRS", ReferencePressure);
  ReferencePressure.SwitchDimensions(shape(0, 2, 1),
                                     GridT_in, GridY_in, GridX_in);

  InputMeteo.ReadWholeField(file_in, "PP", PressurePerturbation);
  PressurePerturbation.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                        GridY_in, GridX_in);
  PressurePerturbation.ReverseData(1);

  Data<real, 4> Temperature(GridT_in, GridZ_MM5_in, GridX_in, GridY_in);

  InputMeteo.ReadWholeField(file_in, "T", Temperature);
  Temperature.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                               GridY_in, GridX_in);
  Temperature.ReverseData(1);

  Data<real, 3> SkinTemperature(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "GROUND T", SkinTemperature);
  SkinTemperature.SwitchDimensions(shape(0, 2, 1),
                                   GridT_in, GridY_in, GridX_in);

  cout << " done." << endl;

  // Computes pressure.
  // See MM5 documentation, section 11.4.1.
  cout << "Computing pressure...";
  cout.flush();

  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in - 1; j++)
        for (i = 0; i < Nx_in - 1; i++)
          Pressure(h, k, j, i) = ReferencePressure(h, j, i)
            * SigmaH(0, Nz_in - k - 1) + ptop
            + PressurePerturbation(h, k, j, i);

  PressurePerturbation.Resize();
  ReferencePressure.Resize();

  cout << " done." << endl;

  // Computes surface pressure.
  // See MM5 documentation, section 7.3.
  cout << "Computing surface pressure...";
  cout.flush();

  Data<real, 3> SurfacePressure(GridT_in, GridY_in, GridX_in);
  for (h = 0; h < Nt_in; h++)
    for (j = 0; j < Ny_in - 1; j++)
      for (i = 0; i < Nx_in - 1; i++)
        SurfacePressure(h, j, i) = Pressure(h, 0, j, i) +
          (Pressure(h, 0, j, i) - Pressure(h, 1, j, i))
          * GridZ_in.Value(0, 0, j, i)
          / (GridZ_in.Value(0, 1, j, i) - GridZ_in.Value(0, 0, j, i));

  // Frees memory.
  Terrain.Resize();
  Terrain_Dot.Resize();

  cout << " done." << endl;

  // Interpolations.
  cout << "Interpolations..." << endl;
  cout.flush();

  // Pressure.
  Data<real, 4> Pressure_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);

  LinearInterpolationDimension(Pressure, Pressure_tmp, 1);
  LinearInterpolationRegularToGeneral(Pressure_tmp, Pressure_out);

  LinearInterpolationRegularToGeneral(SurfacePressure, SurfacePressure_out);

  // Temperature.
  Data<real, 4> Temperature_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);

  LinearInterpolationDimension(Temperature, Temperature_tmp, 1);
  LinearInterpolationRegularToGeneral(Temperature_tmp, Temperature_out);


  ////////////
  // CLOUDS //
  ////////////

  cout << "  + Computing relative humidity and critical relative humidity...";
  cout.flush();

  // Reads specific humidity.
  Data<real, 4> SpecificHumidity(GridT_in, GridZ_MM5_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "Q", SpecificHumidity);
  SpecificHumidity.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                    GridY_in, GridX_in);
  SpecificHumidity.ReverseData(1);
  SpecificHumidity.Threshold(0., 1.);

  // Relative humidity.
  Data<real, 4> RelativeHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  ComputeRelativeHumidity(SpecificHumidity, Temperature, Pressure,
                          RelativeHumidity);
  RelativeHumidity.ThresholdMax(1.);
  Temperature.Resize();

  // Critical relative humidity.
  Data<real, 4> CRH(GridT_in, GridZ_in, GridY_in, GridX_in);
  ComputeCriticalRelativeHumidity(SurfacePressure, Pressure, CRH);

  cout << " done." << endl;

  cout << "  + Computing cloud profile...";
  cout.flush();

  Data<real, 3> BoundaryHeight(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "PBL HGT", BoundaryHeight);
  BoundaryHeight.SwitchDimensions(shape(0, 2, 1),
                                  GridT_in, GridY_in, GridX_in);
  Data<real, 4> CloudFraction(GridT_in, GridZ_in, GridY_in, GridX_in);
  ComputeCloudFraction(BoundaryHeight, RelativeHumidity, CRH, CloudFraction);
  CRH.Resize();

  Data<int, 4> LowIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> MediumIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> HighIndices(Nt_in, Ny_in, Nx_in, 2);

  Data<real, 3> LowCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness(GridT_in, GridY_in, GridX_in);

  ComputeCloudiness(CloudFraction, Pressure, GridZ_interf_in, LowIndices,
                    MediumIndices, HighIndices, LowCloudiness,
                    MediumCloudiness, HighCloudiness);
  Pressure.Resize();

  Data<real, 3> CloudBaseHeight(GridT_in, GridY_in, GridX_in);
  ComputeCloudBaseHeight(LowIndices, MediumIndices, HighIndices,
                         GridZ_interf_in, CloudBaseHeight);

  // Reads liquid water content.
  Data<real, 4> LiquidWaterContent(GridT_in, GridZ_MM5_in,
                                   GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "CLW", LiquidWaterContent);
  LiquidWaterContent.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                      GridY_in, GridX_in);
  LiquidWaterContent.ReverseData(1);
  LiquidWaterContent.ThresholdMin(0.);

  cout << " done." << endl;

  // Fields used by meteo only
  if (with_meteo)
    {
      // Winds.
      // Winds are given at dot points. Z at dots are differents from
      // Z at cross.
      Data<real, 4> MeridionalWind(GridT_in, GridZ_MM5_Dot_in,
                                   GridX_interf_in, GridY_interf_in);
      Data<real, 4> ZonalWind(GridT_in, GridZ_MM5_Dot_in,
                              GridX_interf_in, GridY_interf_in);

      // Wind rotation coefficients.
      // It requires longitudes and latitudes at cross points.
      Data<real, 3> Longitude(Nt_in, Nx_in - 1, Ny_in - 1);
      Data<real, 3> Latitude(Nt_in, Nx_in - 1, Ny_in - 1);
      InputMeteo.ReadWholeField(file_in, "LONGICRS", Longitude);
      InputMeteo.ReadWholeField(file_in, "LATITCRS", Latitude);

      Data<real, 2> Cosine(Nx_in, Ny_in);
      Data<real, 2> Sine(Nx_in, Ny_in);
      real delta_x, delta_y, delta_x2, delta_y2;
      for (i = 0; i < Nx_in - 2; i++)
        for (j = 0; j < Ny_in - 2; j++)
          {
            delta_x = (Longitude(0, i + 1, j) - Longitude(0, i, j))
              * cos(pi / 180. * Latitude(0, i, j));
            delta_y = Latitude(0, i + 1, j) - Latitude(0, i, j);
            delta_x2 = delta_x * delta_x;
            delta_y2 = delta_y * delta_y;
            Cosine(i + 1, j + 1) = delta_x / sqrt(delta_x2 + delta_y2);
            Sine(i + 1, j + 1) = delta_y / sqrt(delta_x2 + delta_y2);
          }
      // Extrapolates on boundaries.
      for (i = 1; i < Nx_in - 1; i++)
        {
          Cosine(i, 0) = 2. * Cosine(i, 1) - Cosine(i, 2);
          Sine(i, 0) = 2. * Sine(i, 1) - Sine(i, 2);
          Cosine(i, Ny_in - 1) = 2. * Cosine(i, Ny_in - 2)
            - Cosine(i, Ny_in - 3);
          Sine(i, Ny_in - 1) = 2. * Sine(i, Ny_in - 2)
            - Sine(i, Ny_in - 3);
        }
      for (j = 0; j < Ny_in; j++)
        {
          Cosine(0, j) = 2. * Cosine(1, j) - Cosine(2, j);
          Sine(0, j) = 2. * Sine(1, j) - Sine(2, j);
          Cosine(Nx_in - 1, j) = 2. * Cosine(Nx_in - 2, j)
            - Cosine(Nx_in - 3, j);
          Sine(Nx_in - 1, j) = 2. * Sine(Nx_in - 2, j)
            - Sine(Nx_in - 3, j);
        }

      InputMeteo.ReadWholeField(file_in, "V", MeridionalWind);
      InputMeteo.ReadWholeField(file_in, "U", ZonalWind);

      // Wind rotation.
      real meridional_wind, zonal_wind;
      for (h = 0; h < Nt_in; h++)
        for (k = 0; k < Nz_in; k++)
          for (i = 0; i < Nx_in; i++)
            for (j = 0; j < Ny_in; j++)
              {
                meridional_wind = MeridionalWind(h, k, i, j);
                zonal_wind = ZonalWind(h, k, i, j);
                ZonalWind(h, k, i, j) = Cosine(i, j) * zonal_wind
                  - Sine(i, j) * meridional_wind;
                MeridionalWind(h, k, i, j) = Sine(i, j) * zonal_wind
                  + Cosine(i, j) * meridional_wind;
              }

      MeridionalWind.SwitchDimensions(shape(0, 1, 3, 2), GridT_in,
                                      GridZ_Dot_in, GridY_interf_in,
                                      GridX_interf_in);
      MeridionalWind.ReverseData(1);

      ZonalWind.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_Dot_in,
                                 GridY_interf_in, GridX_interf_in);
      ZonalWind.ReverseData(1);

      Data<real, 4> MeridionalWind_tmp(GridT_in, GridZ_out,
                                       GridY_interf_in, GridX_interf_in);
      Data<real, 4> ZonalWind_tmp(GridT_in, GridZ_out,
                                  GridY_interf_in, GridX_interf_in);

      LinearInterpolationDimension(MeridionalWind, MeridionalWind_tmp, 1);
      LinearInterpolationRegularToGeneral(MeridionalWind_tmp,
                                          MeridionalWind_out);

      LinearInterpolationDimension(ZonalWind, ZonalWind_tmp, 1);
      LinearInterpolationRegularToGeneral(ZonalWind_tmp, ZonalWind_out);

      cout << " done." << endl;


      ////////////////////////
      // RICHARDSON NUMBERS //
      ////////////////////////

      cout << "Computing Richardson numbers..." << endl;
      cout.flush();

      Data<real, 4> PotentialTemperature(GridT_in, GridZ_out,
                                         GridY_in, GridX_in);
      Data<real, 3> SurfacePotentialTemperature(GridT_in,
                                                GridY_in, GridX_in);

      ComputePotentialTemperature(Temperature_tmp, Pressure_tmp,
                                  PotentialTemperature);
      ComputePotentialTemperature(SkinTemperature, SurfacePressure,
                                  SurfacePotentialTemperature);
      cout << "  + compute potential temp" << endl;


      Data<real, 4> Richardson(GridT_in, GridZ_out, GridY_in, GridX_in);
      Data<real, 3> SurfaceRichardson(GridT_in, GridY_in, GridX_in);

      Data<real, 4> MeridionalWind_cross(GridT_in, GridZ_out,
                                         GridY_in, GridX_in);
      Data<real, 4> ZonalWind_cross(GridT_in, GridZ_out,
                                    GridY_in, GridX_in);
      Data<real, 3> FirstLevelWindModule(GridT_in,
                                         GridY_in, GridX_in);

      // Interpolates ZonalWind and Meridional wind at cross points.
      LinearInterpolationRegular(ZonalWind_tmp, ZonalWind_cross);
      LinearInterpolationRegular(MeridionalWind_tmp, MeridionalWind_cross);
      cout << "  + LinearInt" << endl;

      ZonalWind_tmp.Resize();
      MeridionalWind_tmp.Resize();

      // Richardson number.
      ComputeRichardson(ZonalWind_cross, MeridionalWind_cross,
                        PotentialTemperature, Richardson);
      ComputeModule(ZonalWind_cross, MeridionalWind_cross,
                    FirstLevelWindModule);
      ComputeRichardson(FirstLevelWindModule, SurfacePotentialTemperature,
                        PotentialTemperature, SurfaceRichardson);

      PotentialTemperature.Resize();
      SurfacePotentialTemperature.Resize();

      LinearInterpolationRegularToGeneral(FirstLevelWindModule,
                                          FirstLevelWindModule_out);
      FirstLevelWindModule.Resize();

      ZonalWind_cross.Resize();
      MeridionalWind_cross.Resize();
      ZonalWind.Resize();
      MeridionalWind.Resize();

      LinearInterpolationRegularToGeneral(SurfaceRichardson,
                                          SurfaceRichardson_out);
      SurfaceRichardson.Resize();

      LinearInterpolationRegularToGeneral(Richardson, Richardson_out);
      Richardson.Resize();

      cout << " done." << endl;

    }

  SurfacePressure.Resize();
  Temperature_tmp.Resize();
  Pressure_tmp.Resize();

  // Fields used by photolysis only
  Data<real, 4> Attenuation(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> LiquidWaterExtinction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> IceWaterExtinction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> IceWaterContent(GridT_in, GridZ_MM5_in, GridX_in, GridY_in);

  if (with_photolysis)
    {
      if (ice_cloud &&
          (photolysis_tabulation == 2 || photolysis_tabulation == 3))
        InputMeteo.ReadWholeField(file_in, "ICE", IceWaterContent);
      else
        IceWaterContent.Fill(0.);
      IceWaterContent.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                       GridY_in, GridX_in);
      IceWaterContent.ReverseData(1);
      IceWaterContent.ThresholdMin(0.);

      if (photolysis_tabulation == 1)
        {
          cout << "  + Computing attenuation...";
          cout.flush();

          if (attenuation_type == 1)
            ComputeAttenuation_LWC(LiquidWaterContent,
                                   LowIndices, MediumIndices, HighIndices,
                                   MediumCloudiness, HighCloudiness,
                                   date_beg, Delta_t_in, Attenuation);
          else if (attenuation_type == 2)
            ComputeAttenuation_ESQUIF(MediumCloudiness, HighCloudiness,
                                      RelativeHumidity, Attenuation);
          cout << " done." << endl;
        }

      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          cout << "Computing Extinction...";
          cout.flush();

          ComputeExtinction(LiquidWaterContent, IceWaterContent,
                            LiquidWaterExtinction, IceWaterExtinction);

          cout << " done." << endl;
        }
    }

  // Frees memory.
  LowIndices.Resize();
  MediumIndices.Resize();
  HighIndices.Resize();

  LowCloudiness.Resize();
  MediumCloudiness.Resize();
  HighCloudiness.Resize();

  ///////////////////////////
  // LINEAR INTERPOLATIONS //
  ///////////////////////////

  cout << "Linear interpolations..." << endl;

  // Interpolates 3D fields with data at cross points.

  // Common Fields.

  // LiquidWaterContent.
  Data<real, 4> LiquidWaterContent_tmp(GridT_in, GridZ_out,
                                       GridY_in, GridX_in);
  LinearInterpolationDimension(LiquidWaterContent,
                               LiquidWaterContent_tmp, 1);
  LinearInterpolationRegularToGeneral(LiquidWaterContent_tmp,
                                      LiquidWaterContent_out);
  LiquidWaterContent_tmp.Resize();
  LiquidWaterContent.Resize();
  cout << "  + Liquid Water content" << endl;

  // Interpolates 2D Field (with data at cross points)

  // SpecificHumidity.
  Data<real, 4> SpecificHumidity_tmp(GridT_in, GridZ_out,
                                     GridY_in, GridX_in);
  LinearInterpolationDimension(SpecificHumidity,
                               SpecificHumidity_tmp, 1);
  LinearInterpolationRegularToGeneral(SpecificHumidity_tmp,
                                      SpecificHumidity_out);
  SpecificHumidity_tmp.Resize();
  SpecificHumidity.Resize();
  cout << "  + SpecificHumidity" << endl;

  if (with_meteo)
    {
      // CloudBaseHeight.
      LinearInterpolationRegularToGeneral(CloudBaseHeight, CloudBaseHeight_out);
      CloudBaseHeight_out.ThresholdMin(min_height);
      CloudBaseHeight.Resize();
      cout << "  + CloudBaseHeight" << endl;

      // SkinTemperature.
      LinearInterpolationRegularToGeneral(SkinTemperature,
                                          SkinTemperature_out);
      SkinTemperature.Resize();
      cout << "  + SkinTemperature" << endl;

      // SensibleHeat.
      Data<real, 3> SensibleHeat(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "SHFLUX", SensibleHeat);
      SensibleHeat.SwitchDimensions(shape(0, 2, 1),
                                    GridT_in, GridY_in, GridX_in);
      LinearInterpolationRegularToGeneral(SensibleHeat, SensibleHeat_out);

      // Divides by (\rho_{air} * cp).
      for (h = 0; h < Nt_out; h++)
        for (j = 0; j < Ny_out; j++)
          for (i = 0; i < Nx_out; i++)
            SensibleHeat_out(h, j, i) *= r * SkinTemperature_out(h, j, i)
              / (SurfacePressure_out(h, j, i) * cp);
      SensibleHeat.Resize();
      cout << "  + SensibleHeat" << endl;

      // Evaporation.
      Data<real, 3> LatentHeat(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "LHFLUX", LatentHeat);
      LatentHeat.SwitchDimensions(shape(0, 2, 1),
                                  GridT_in, GridY_in, GridX_in);
      LinearInterpolationRegularToGeneral(LatentHeat, Evaporation_out);
      // Divides by \rho_{water} * L.
      Evaporation_out.Mlt(1. / 2.5e9);
      LatentHeat.Resize();
      cout << "  + Evaporation" << endl;

      // SurfaceTemperature.
      Data<real, 3> SurfaceTemperature(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "T2", SurfaceTemperature);
      SurfaceTemperature.SwitchDimensions(shape(0, 2, 1),
                                          GridT_in, GridY_in, GridX_in);
      LinearInterpolationRegularToGeneral(SurfaceTemperature,
                                          SurfaceTemperature_out);
      SurfaceTemperature.Resize();
      cout << "  + SurfaceTemperature" << endl;

      // SoilWater.
      Data<real, 3> SoilWater(GridT_in, GridX_in, GridY_in);
      // TODO: Should check that ISOIL=2 before reading "SOIL W 1" and maybe
      // issue a warning if the data is not available.
      // see http://www2.mmm.ucar.edu/mm5/On-Line-Tutorial/graph/extras/avail_fields.html
      // There is the same problem for other MM5 fields.
      SoilWater.SetZero();
      InputMeteo.ReadWholeField(file_in,
                                "SOIL W 1", SoilWater);
      SoilWater.SwitchDimensions(shape(0, 2, 1),
                                 GridT_in, GridY_in, GridX_in);
      SoilWater.Threshold(0., 1.);
      LinearInterpolationRegularToGeneral(SoilWater, SoilWater_out);
      SoilWater.Resize();
      cout << "  + SoilWater" << endl;

      // SolarRadiation.
      Data<real, 3> SolarRadiation(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "SWDOWN", SolarRadiation);
      SolarRadiation.SwitchDimensions(shape(0, 2, 1),
                                      GridT_in, GridY_in, GridX_in);
      LinearInterpolationRegularToGeneral(SolarRadiation,
                                          SolarRadiation_out);
      SolarRadiation.Resize();
      cout << "  + SolarRadiation" << endl;

      // ConvectiveRain.
      Data<real, 3> ConvectiveRain(GridT_in, GridX_in, GridY_in);
      Data<real, 3> ConvectiveRain_prev(GridT_in, GridX_in, GridY_in);

      InputMeteo.ReadWholeField(file_in,
                                "RAIN CON", ConvectiveRain);
      ConvectiveRain.SwitchDimensions(shape(0, 2, 1),
                                      GridT_in, GridY_in, GridX_in);

      ConvectiveRain.ThresholdMin(0.);

      if (exists(file_in_prev) && prev_accumulated_rain)
        InputMeteo.ReadWholeField(file_in_prev,
                                  "RAIN CON", ConvectiveRain_prev);
      else
        ConvectiveRain_prev.Fill(0.);
      ConvectiveRain_prev.SwitchDimensions(shape(0, 2, 1),
                                           GridT_in, GridY_in, GridX_in);
      ConvectiveRain_prev.ThresholdMin(0.);

      Decumulate(ConvectiveRain, Nt_in, 0);
      for (j = 0; j < Ny_in - 1; j++)
        for (i = 0; i < Nx_in - 1; i++)
          ConvectiveRain(0, j, i) -= ConvectiveRain_prev(Nt_in - 1, j, i);
      // Converts cm/step to mm/h.
      ConvectiveRain.Mlt(10. / Delta_t_in);

      // Rain.
      // First reads NonConvectiveRain.
      Data<real, 3> Rain(GridT_in, GridX_in, GridY_in);
      Data<real, 3> Rain_prev(GridT_in, GridX_in, GridY_in);

      InputMeteo.ReadWholeField(file_in, "RAIN NON", Rain);
      Rain.SwitchDimensions(shape(0, 2, 1), GridT_in, GridY_in, GridX_in);
      Rain.ThresholdMin(0.);

      if (exists(file_in_prev) && prev_accumulated_rain)
        InputMeteo.ReadWholeField(file_in_prev,
                                  "RAIN NON", Rain_prev);
      else
        Rain_prev.Fill(0.);
      Rain_prev.SwitchDimensions(shape(0, 2, 1), GridT_in,
                                 GridY_in, GridX_in);
      Rain_prev.ThresholdMin(0.);

      Decumulate(Rain, Nt_in, 0);
      for (j = 0; j < Ny_in - 1; j++)
        for (i = 0; i < Nx_in - 1; i++)
          Rain(0, j, i) -= Rain_prev(Nt_in - 1, j, i);
      // Converts cm/step to mm/h.
      Rain.Mlt(10. / Delta_t_in);

      // Adds convective rain to NonConvectiveRain.
      for (h = 0; h < Nt_in; h++)
        for (j = 0; j < Ny_in - 1; j++)
          for (i = 0; i < Nx_in - 1; i++)
            Rain(h, j, i) = Rain(h, j, i) + ConvectiveRain(h, j, i);

      // Interpolates to Rain_out.
      LinearInterpolationRegularToGeneral(Rain, Rain_out);
      Rain_out.ThresholdMin(0.);
      Rain.Resize();

      // Interpolates and frees ConvectiveRain.
      LinearInterpolationRegularToGeneral(ConvectiveRain, ConvectiveRain_out);
      ConvectiveRain_out.ThresholdMin(0.);
      ConvectiveRain.Resize();
      cout << "  + Rain" << endl;

      // FrictionModule.
      Data<real, 3> FrictionModule(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "UST", FrictionModule);
      FrictionModule.SwitchDimensions(shape(0, 2, 1),
                                      GridT_in, GridY_in, GridX_in);
      LinearInterpolationRegularToGeneral(FrictionModule, FrictionModule_out);
      FrictionModule.Resize();
      cout << "  + FrictionModule" << endl;

      // BoundaryHeight.
      LinearInterpolationRegularToGeneral(BoundaryHeight, BoundaryHeight_out);
      BoundaryHeight.Resize();
      BoundaryHeight_out.ThresholdMin(GridZ_interf_out(1));
      cout << "  + BoundaryHeight" << endl;

      // Wind Module at 10m
      Data<real, 3> U10(GridT_in, GridX_in, GridY_in);
      InputMeteo.ReadWholeField(file_in, "U10", U10);
      U10.SwitchDimensions(shape(0, 2, 1),
                           GridT_in, GridY_in, GridX_in);
      Data<real, 3> V10(GridT_in, GridX_in, GridY_in);
      Data<real, 3> WindModule10(GridT_in, GridY_in, GridX_in);
      InputMeteo.ReadWholeField(file_in, "V10", V10);
      V10.SwitchDimensions(shape(0, 2, 1),
                           GridT_in, GridY_in, GridX_in);
      ComputeModule(U10, V10, WindModule10);
      LinearInterpolationRegularToGeneral(WindModule10, WindModule10_out);

      U10.Resize();
      V10.Resize();
      WindModule10.Resize();
      cout << "  + WindModule10" << endl;
    }

  // Fields used by photolysis only
  if (with_photolysis)
    {
      // CloudFraction.
      Data<real, 4> CloudFraction_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);
      LinearInterpolationDimension(CloudFraction, CloudFraction_tmp, 1);
      LinearInterpolationRegularToGeneral(CloudFraction_tmp, CloudFraction_out);
      CloudFraction_tmp.Resize();
      CloudFraction.Resize();
      cout << "  + Cloud Fraction" << endl;

      // IceWaterContent.
      if (ice_cloud && photolysis_tabulation != 1)
        {
          Data<real, 4> IceWaterContent_tmp(GridT_in, GridZ_out,
                                            GridY_in, GridX_in);

          LinearInterpolationDimension(IceWaterContent,
                                       IceWaterContent_tmp, 1);
          LinearInterpolationRegularToGeneral(IceWaterContent_tmp,
                                              IceWaterContent_out);

          IceWaterContent_tmp.Resize();
          IceWaterContent.Resize();

          cout << "  + Ice Water content" << endl;
        }

      // Attenuation.
      if (photolysis_tabulation == 1)
        {
          Data<real, 4> Attenuation_tmp(GridT_in, GridZ_out,
                                        GridY_in, GridX_in);

          LinearInterpolationDimension(Attenuation, Attenuation_tmp, 1);
          LinearInterpolationRegularToGeneral(Attenuation_tmp,
                                              Attenuation_out);
          Attenuation_out.Threshold(0., 2.);
          Attenuation_tmp.Resize();
          Attenuation.Resize();
          cout << "  + Attenuation" << endl;
        }

      // Extinction.
      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          Data<real, 4> LiquidWaterExtinction_tmp(GridT_in, GridZ_out,
                                                  GridY_in, GridX_in);
          LinearInterpolationDimension(LiquidWaterExtinction,
                                       LiquidWaterExtinction_tmp, 1);
          LinearInterpolationRegularToGeneral(LiquidWaterExtinction_tmp,
                                              LiquidWaterExtinction_out);
          LiquidWaterExtinction_tmp.Resize();
          LiquidWaterExtinction.Resize();

          cout << "  + Cloud Extinction" << endl;
        }

      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          Data<real, 4> IceWaterExtinction_tmp(GridT_in, GridZ_out,
                                               GridY_in, GridX_in);
          LinearInterpolationDimension(IceWaterExtinction,
                                       IceWaterExtinction_tmp, 1);
          LinearInterpolationRegularToGeneral(IceWaterExtinction_tmp,
                                              IceWaterExtinction_out);
          IceWaterExtinction_tmp.Resize();
          IceWaterExtinction.Resize();

          cout << "  + Ice Cloud Extinction" << endl;
        }
    } // with_photolysis

  cout << " done." << endl;

  // Resizes grids with latlon coord.
  Pressure_out.ResizeGrid(GridT_out, GridZ_out, GridY_out, GridX_out);
  Temperature_out.ResizeGrid(GridT_out, GridZ_out, GridY_out, GridX_out);
  SpecificHumidity_out.ResizeGrid(GridT_out, GridZ_out, GridY_out, GridX_out);
  LiquidWaterContent_out.ResizeGrid(GridT_out, GridZ_out,
                                    GridY_out, GridX_out);
  BoundaryHeight_out.ResizeGrid(GridT_out, GridY_out, GridX_out);

  if (with_meteo)
    {
      Richardson_out.ResizeGrid(GridT_out, GridZ_out, GridY_out, GridX_out);
      MeridionalWind_out.ResizeGrid(GridT_out, GridZ_out,
                                    GridY_interf_out, GridX_out);
      ZonalWind_out.ResizeGrid(GridT_out, GridZ_out,
                               GridY_out, GridX_interf_out);
      SurfaceTemperature_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      SurfacePressure_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      SkinTemperature_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      SoilWater_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      SolarRadiation_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      ConvectiveRain_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      Rain_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      FrictionModule_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      SurfaceRichardson_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      FirstLevelWindModule_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      WindModule10_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
    }

  if (with_photolysis)
    {
      CloudFraction_out.ResizeGrid(GridT_out, GridZ_out,
                                   GridY_out, GridX_out);
      if (photolysis_tabulation == 1)
        Attenuation_out.ResizeGrid(GridT_out, GridZ_out,
                                   GridY_out, GridX_out);
      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          LiquidWaterExtinction_out.ResizeGrid(GridT_out, GridZ_out,
                                               GridY_out, GridX_out);
          IceWaterExtinction_out.ResizeGrid(GridT_out, GridZ_out,
                                            GridY_out, GridX_out);
        }
    }

  if (with_photolysis)
    {
      ///////////////////
      // OPTICAL DEPTH //
      ///////////////////

      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          cout << "Computing OD...";
          cout.flush();

          real dz, lext, iext, cf;

          for (h = 0; h < Nt_out; h++)
            for (j = 0; j < Ny_out; j++)
              for (i = 0; i < Nx_out; i++)
                for (k = Nz_out - 1; k >= 0; k--)
                  {
                    cf = CloudFraction_out(h, k, j, i);
                    dz = GridZ_interf_out(k + 1) - GridZ_interf_out(k);
                    lext = LiquidWaterExtinction_out(h, k, j, i);
                    iext = IceWaterExtinction_out(h, k, j, i);
                    OpticalDepth_out(h, k, j, i) = lext * dz * pow(cf, 1.5);
                    if (ice_cloud)
                      IceOpticalDepth_out(h, k, j, i) = iext * dz *
                        pow(cf, 1.5);
                  }

          if (!ice_cloud)
            IceOpticalDepth_out.Fill(0.);

          cout << " done." << endl;
        }

      //////////////////////
      // PHOTOLYSIS RATES //
      //////////////////////

      if (photolysis_tabulation == 2)
        {
          int Day, Year;
          real  Hour, Minutes, FloatHour;

          Data<real, 3> OpticalDepthAerosol(GridZ_out, GridY_out, GridX_out);
          Data<real, 4> SingleScatteringAlbedo(GridWavelength, GridZ_out,
                                               GridY_out, GridZ_out);
          Data<real, 4> MeanExtinctionEfficiencyFactor(GridWavelength,
                                                       GridZ_out, GridY_out,
                                                       GridX_out);
          Data<real, 5> PhaseFunction(GridWavelength,
                                      GridZ_out, GridY_out,
                                      GridX_out, GridLegendre);
          Data<real, 4> PhotolysisRate_TimeStep(GridP, GridZ_out,
                                                GridY_out, GridX_out);
          Data<real, 3> OpticalDepth_TimeStep(GridZ_out,
                                              GridY_out, GridX_out);
          Data<real, 3> IceOpticalDepth_TimeStep(GridZ_out,
                                                 GridY_out, GridX_out);
          Data<real, 2> SurfacePressure_TimeStep(GridY_out, GridX_out);

          Day = date_beg.GetOrdinalDay();
          Year = date_beg.GetYear();
          Hour = date_beg.GetHour();
          Minutes = date_beg.GetMinutes();

          FloatHour = Hour + Minutes / 60.;

          OpticalDepthAerosol.Fill(0.0);
          SingleScatteringAlbedo.Fill(1.0);
          MeanExtinctionEfficiencyFactor.Fill(0.1);
          PhaseFunction.Fill(1.0);

          cout << "FastJX...";

          const char* DirectoryParameter = fastJ_parameter_files.c_str();
          int DirectoryParameter_len = strlen(DirectoryParameter);

          for (i = 0; i < Nt_out - 1; i++)
            {
              OpticalDepth_TimeStep.SubData(OpticalDepth_out, i,
                                            Range::all(),
                                            Range::all(),
                                            Range::all());
              IceOpticalDepth_TimeStep.SubData(IceOpticalDepth_out, i,
                                               Range::all(),
                                               Range::all(),
                                               Range::all());
              SurfacePressure_TimeStep.SubData(SurfacePressure_out, i,
                                               Range::all(),
                                               Range::all());

              _fastjx(&Nx_out, &Ny_out, &Nz_out, &x_min_out, &Delta_x_out,
                      &y_min_out, &Delta_y_out,
                      GridZ_interf_out.GetArray().data(),
                      OpticalDepth_TimeStep.GetData(),
                      IceOpticalDepth_TimeStep.GetData(),
                      SurfacePressure_TimeStep.GetData(),
                      OpticalDepthAerosol.GetData(),
                      SingleScatteringAlbedo.GetData(),
                      MeanExtinctionEfficiencyFactor.GetData(),
                      PhaseFunction.GetData(),
                      PhotolysisRate_TimeStep.GetData(),
                      photolysis_specie_name.dataFirst(), &Nr_photolysis,
                      &Year, &Day, &FloatHour,
                      DirectoryParameter, &DirectoryParameter_len);

              for (int h = 0; h < Nz_out; h++)
                for (int j = 0; j < Ny_out; j++)
                  for (int k = 0; k < Nx_out; k++)
                    for (int w = 0; w < Nr_photolysis ; w++)
                      PhotolysisRate_out(w, i, h, j, k) =
                        PhotolysisRate_TimeStep(w, h, j, k);

              FloatHour = FloatHour + Delta_t_out;
            }
          cout << " done." << endl;
        }
    }

  if (with_meteo)
    {
      /////////////////////////
      // WINDS AND DIFFUSION //
      /////////////////////////

      // Vertical diffusion (Louis formula, 1979).
      cout << "Computing Kz (Louis formula)...";
      cout.flush();

      ComputePotentialTemperature(Temperature_out, Pressure_out,
                                  PotentialTemperature_out);

      ComputeLouisKz(ZonalWind_out, MeridionalWind_out,
                     PotentialTemperature_out, Kz_out);

      real Kz_max_loc;
      int imax(0);

      for (h = 0; h < Nt_out; h++)
        for (j = 0; j < Ny_out; j++)
          for (i = 0; i < Nx_out; i++)
            {
              Kz_max_loc = 0.0;
              for (k = 0; k < Nz_out + 1; k++)
                if (Kz_out(h, k, j, i) >= Kz_max_loc)
                  {
                    Kz_max_loc = Kz_out(h, k, j, i);
                    imax = k;
                  }
              if (ConvectiveRain_out(h, j, i) > 1. / 6.)
                for (k = imax; k < Nz_out + 1; k++)
                  Kz_out(h, k, j, i) = Kz_max_loc;
              else
                for (k = imax; k < Nz_out + 1; k++)
                  Kz_out(h, k, j, i) = (1. - ConvectiveRain_out(h, j, i) * 6.)
                    * Kz_out(h, k, j, i)
                    + 6. * ConvectiveRain_out(h, j, i) * Kz_max_loc;
            }

      // Kz_min.
      real local_min;
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            local_min = Kz_min * (1. - LUC(urban_index, j, i))
              + Kz_min_urban * LUC(urban_index, j, i);
            if (apply_vert)
              {
                for (h = 0; h < Nt_out; h++)
                  for (k = 0; k < Nz_out + 1; k++)
                    if (Kz_out(h, k, j, i) < local_min)
                      Kz_out(h, k, j, i) = local_min;
              }
            else
              for (h = 0; h < Nt_out; h++)
                if (Kz_out(h, 1, j, i) < local_min)
                  Kz_out(h, 1, j, i) = local_min;
          }

      // Kz_max.
      Kz_out.ThresholdMax(Kz_max);

      cout << " done." << endl;


      /////////////////////////////////////////
      // PHOTOSYNTHETICALLY ACTIVE RADIATION //
      /////////////////////////////////////////

      // For biogenic emissions.
      cout << "Computing PAR...";
      cout.flush();
      real zenith_angle;
      real optical_thickness, visible_beam, visible_diff,
        water_abs, beam, diff, visible_rad, IR_rad,
        ratio, ratio0;
      real fvb, fvd, PARdb, PARdiff;
      real watt2umol(4.6);

      Date current_date(date_beg);
      for (h = 0; h < Nt_out; h++)
        {
          for (j = 0; j < Ny_out; j++)
            for (i = 0; i < Nx_out; i++)
              {
                real ut = real(GridT_out(h) - int(GridT_out(h) / 24) * 24.);
                zenith_angle = ZenithAngle(GridX_out(i), GridY_out(j),
                                           current_date, ut) / 180. * pi;

                if (zenith_angle >= 1.51844
                    || SolarRadiation_out(h, j, i) <= 0.)
                  {
                    PARdb_out(h, j, i) = 0.;
                    PARdiff_out(h, j, i) = 0.;
                    PAR_out(h, j, i) = 0.;
                  }
                else
                  {
                    optical_thickness = Pressure_out(h, 0, j, i) / 101325.
                      / cos(zenith_angle);
                    visible_beam = 600. * exp(-.185 * optical_thickness)
                      * cos(zenith_angle);
                    visible_diff = 0.42 * (600. - visible_beam)
                      * cos(zenith_angle);
                    water_abs = 1320. * .077 * pow(2.
                                                   * optical_thickness, 0.3);
                    beam = (720. * exp(-0.06 * optical_thickness)
                            - water_abs) * cos(zenith_angle);
                    diff = 0.65 * (720. - water_abs - beam)
                      * cos(zenith_angle);

                    visible_rad = visible_beam + visible_diff;
                    IR_rad  = beam + diff;
                    ratio  = visible_rad / (visible_rad + IR_rad);
                    ratio0 = SolarRadiation_out(h, j, i)
                      / (visible_rad + IR_rad);

                    if (ratio0 >= 0.89)
                      fvb = visible_beam / visible_rad * 0.941124;
                    else if (ratio0 <= 0.21)
                      fvb = visible_beam / visible_rad * 9.55e-3;
                    else
                      fvb = visible_beam / visible_rad
                        * (1. - pow((0.9 - ratio0) / 0.7, 0.666667));
                    fvd = 1. - fvb;

                    PARdb  = SolarRadiation_out(h, j, i) * ratio * fvb
                      * watt2umol;
                    PARdiff = SolarRadiation_out(h, j, i) * ratio * fvd
                      * watt2umol;

                    PARdb_out(h, j, i) = PARdb;
                    PARdiff_out(h, j, i) = PARdiff;
                    PAR_out(h, j, i) = PARdb + PARdiff;
                  }
              }
          current_date.AddSeconds(int(Delta_t_out * 3600));
        }

      cout << " done." << endl;
    }


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  FormatBinary<float> OutputMeteo;

  cout << "Writing data...";
  cout.flush();

  cout << "  (directory_out=" << to_str(directory_out) << ")" << endl;

  if (with_meteo)
    {
      OutputMeteo.Append(Pressure_out, directory_out + "Pressure.bin");
      OutputMeteo.Append(SurfacePressure_out, directory_out
                         + "SurfacePressure.bin");
      OutputMeteo.Append(Temperature_out, directory_out + "Temperature.bin");
      OutputMeteo.Append(SurfaceTemperature_out,
                         directory_out + "SurfaceTemperature.bin");
      OutputMeteo.Append(SkinTemperature_out, directory_out
                         + "SkinTemperature.bin");
      OutputMeteo.Append(Richardson_out, directory_out + "Richardson.bin");
      OutputMeteo.Append(SurfaceRichardson_out,
                         directory_out + "SurfaceRichardson.bin");
      OutputMeteo.Append(SpecificHumidity_out, directory_out
                         + "SpecificHumidity.bin");
      OutputMeteo.Append(LiquidWaterContent_out, directory_out
                         + "LiquidWaterContent.bin");
      OutputMeteo.Append(CloudBaseHeight_out, directory_out + "CloudBaseHeight.bin");
      OutputMeteo.Append(SolarRadiation_out, directory_out
                         + "SolarRadiation.bin");
      OutputMeteo.Append(Rain_out, directory_out + "Rain.bin");
      OutputMeteo.Append(PARdb_out, directory_out + "PARdb.bin");
      OutputMeteo.Append(PARdiff_out, directory_out + "PARdiff.bin");
      OutputMeteo.Append(PAR_out, directory_out + "PAR.bin");
      OutputMeteo.Append(ZonalWind_out, directory_out + "ZonalWind.bin");
      OutputMeteo.Append(MeridionalWind_out, directory_out
                         + "MeridionalWind.bin");
      ComputeModule(MeridionalWind_out, ZonalWind_out, WindModule_out);
      OutputMeteo.Append(WindModule_out, directory_out + "WindModule.bin");
      OutputMeteo.Append(FrictionModule_out, directory_out
                         + "FrictionModule.bin");
      OutputMeteo.Append(BoundaryHeight_out, directory_out
                         + "BoundaryHeight.bin");
      OutputMeteo.Append(Kz_out, directory_out + "Kz_Louis.bin");
      OutputMeteo.Append(SoilWater_out, directory_out + "SoilWater.bin");
      OutputMeteo.Append(SensibleHeat_out, directory_out
                         + "SensibleHeat.bin");
      OutputMeteo.Append(Evaporation_out, directory_out
                         + "Evaporation.bin");
      OutputMeteo.Append(PotentialTemperature_out, directory_out
                         + "PT.bin");
      OutputMeteo.Append(FirstLevelWindModule_out,
                         directory_out + "FirstLevelWindModule.bin");
      OutputMeteo.Append(WindModule10_out,
                         directory_out + "WindModule10.bin");
    }
  if (with_photolysis)
    {
      OutputMeteo.Append(CloudFraction_out, directory_out + "CloudFraction.bin");
      if (photolysis_tabulation == 1)
        OutputMeteo.Append(Attenuation_out, directory_attenuation
                           + "Attenuation.bin");
      if (photolysis_tabulation == 2)
        {
          Data<real, 4> PhotolysisRate4D(GridT_out, GridZ_out, GridY_out,
                                         GridX_out);
          for (i = 0; i < Nr_photolysis; i++)
            {
              PhotolysisRate4D.SubData(PhotolysisRate_out, i, Range::all(),
                                       Range::all(), Range::all(),
                                       Range::all());
              OutputMeteo.Append(PhotolysisRate4D, directory_photolysis
                                 + photolysis_reaction_list[i] + ".bin");
            }
        }
      if (photolysis_tabulation == 3)
        {
          OutputMeteo.Append(OpticalDepth_out, directory_out +
                             "CloudOpticalDepth.bin");
          OutputMeteo.Append(IceOpticalDepth_out, directory_out +
                             "IceOpticalDepth.bin");
        }
      if (ice_cloud && (photolysis_tabulation == 2 || photolysis_tabulation == 3))
        OutputMeteo.Append(IceWaterContent_out, directory_out +
                           "IceWaterContent.bin");
    }


  cout << endl << endl;
  cout << "||||||||||||||||||||||||||||||||||||||||||" << endl;
  cout << "|   Successful completion of MM5-meteo   |" << endl;
  cout << "||||||||||||||||||||||||||||||||||||||||||" << endl << endl;

  END;

  return 0;
}

