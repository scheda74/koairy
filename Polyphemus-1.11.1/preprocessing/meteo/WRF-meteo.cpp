// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Victor Winiarek and Vivien Mallet
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
#define SELDONDATA_DEBUG_CHECK_DIMENSIONS
#define SELDONDATA_WITH_NETCDF
#define COMMON_WITH_NETCDF

#include "SeldonData.hxx"
using namespace SeldonData;

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

#include "fastJX.hxx"



int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("WRF-meteo.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);


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
  int Nt_in, Nz_in, Ny_in, Nx_in, Nsoil_in;
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
  int Nr_photolysis = 0, Nwavelength = 0, Nlegendre = 0;

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

  // Input domain.
  // These dimensions are directly read in the WRF file.
  cout << "    + Reading WRF file's information...";
  cout.flush();

  FormatNetCDF<float> InputFile;

  //// Paths. ////

  configuration.SetSection("[paths]");
  configuration.PeekValue("Directory_meteo", directory_out);
  configuration.PeekValue("Database_WRF-meteo", file_in);

  cout << endl;

  //// Input Files. ////

  cout << "  + Reading input files paths and names...";
  cout.flush();

  // Processes the input file name to prepare date expansion.
  file_in = find_replace(file_in, "&D", "%y-%m-%d");
  file_in = find_replace(file_in, "&y", "%y");
  file_in = find_replace(file_in, "&m", "%m");
  file_in = find_replace(file_in, "&d", "%d");
  //file_in_prev = date_prev.GetDate(file_in);
  file_in = date_beg.GetDate(file_in);
  if (!exists(file_in))
    throw string("Unable to find WRF file \"") + file_in + "\".";

  // Land use.
  string LUC_file;
  configuration.PeekValue("LUC_file", LUC_file);
  if (!exists(LUC_file))
    throw "Unable to open land use cover file \"" + LUC_file + "\".";
  int Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny_out * Nx_out);
  int urban_index;
  configuration.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                          urban_index);

  cout << endl;

  //// Parameterizations. ////

  // Read the size of Times variable.
  InputFile.ReadDimension(file_in, "Times", 0, Nt_in);
  // Read the number of soil layers.
  InputFile.ReadDimension(file_in, "ZS", 1, Nsoil_in);
  // Read in the netcdf file's global attributes.
  InputFile.ReadAttribute(file_in, "BOTTOM-TOP_PATCH_END_UNSTAG", Nz_in);
  InputFile.ReadAttribute(file_in, "SOUTH-NORTH_PATCH_END_UNSTAG", Ny_in);
  InputFile.ReadAttribute(file_in, "WEST-EAST_PATCH_END_UNSTAG", Nx_in);

  // First date in the WRF file.
  Date date_beg_meteo = read_date_WRF(file_in);
  // WRF file's timestep (in hours).
  Delta_t_in = read_delta_t_WRF(file_in);

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

  Date date_end_meteo = date_beg_meteo;
  date_end_meteo.AddSeconds(Nt_in * Delta_t_in * 3600);

  if (date_beg_meteo > date_beg || date_end_meteo < date_end)
    throw string("The WRF file \"") + file_in
      + "\" does not contain data for all dates.";

  real t_min_in, t_min_out, differences;
  t_min_in = real(date_beg_meteo.GetHour())
    + real(date_beg_meteo.GetMinutes()) / 60.
    + real(date_beg_meteo.GetSeconds()) / 3600.;
  differences = real(date_beg.GetSecondsFrom(date_beg_meteo));
  t_min_out = t_min_in + differences / 3600.;

  InputFile.ReadAttribute(file_in, "MAP_PROJ", projection_type);
  // the case "map_proj = 6" (lat/lon) is not implemented yet.
  real xmin_in, ymin_in, deltax_in, deltay_in;
  if (projection_type == 6)
    {
      //   throw "lon/lat projection in WRF's use is not implemented yet.";
      float lon_0, lon_1, lat_0, lat_1;
      Data<real, 3> Longitude(Nt_in, Ny_in, Nx_in);
      InputFile.Read(file_in, "XLONG", Longitude);
      lon_0 = Longitude(0, 0, 0);
      lon_1 = Longitude(0, 0, Nx_in - 1);

      cout << lon_0 << " lon_0 " << lon_1 << " lon_1" << endl;
      cout << endl;
      Data<real, 3> Latitude(Nt_in, Ny_in, Nx_in);
      Data<real, 1> BottomLatitude(Nx_in);
      Data<real, 1> UpLatitude(Nx_in);
      InputFile.Read(file_in, "XLAT", Latitude);
      for (i = 0; i < Nx_in; i++)
        {
          BottomLatitude(i) = Latitude(0, 0, i);
          UpLatitude(i) = Latitude(0, Ny_in - 1, i);
        }
      lat_0 = BottomLatitude.GetMax();
      lat_1 = UpLatitude.GetMax();
      cout << lat_0 << " lat_0 " << lat_1 << " lat_1" << endl;
      xmin_in = lon_0;
      ymin_in = lat_0;
      deltax_in =  Longitude(0, 0, 1) - Longitude(0, 0, 0);
      deltay_in =  Latitude(0, 1, 0) - Latitude(0, 0, 0);
      cout << deltax_in << " deltax_in - deltay_in " << deltay_in << endl;
    }
  else
    {
      xmin_in = 0.5;
      ymin_in = 0.5;
      deltax_in = 1.;
      deltay_in = 1.;
    }

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for grids..." << endl;
  cout.flush();

  // Input settings.
  cout << "  + Input grids...";

  // Input grids.
  RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
  // Land use categories.
  RegularGrid<real> GridC(Nc);
  // Soil layers.
  RegularGrid<real> GridSoil_in(Nsoil_in);

  // Horizontal grids.
  RegularGrid<real> GridX_in(xmin_in, deltax_in, Nx_in);
  RegularGrid<real> GridY_in(ymin_in, deltay_in, Ny_in);
  RegularGrid<real> GridX_interf_in(xmin_in - deltax_in / 2., deltax_in, Nx_in + 1);
  RegularGrid<real> GridY_interf_in(ymin_in - deltay_in / 2., deltay_in, Ny_in + 1);

  // Grid for Z in height. Z depends on Z, Y, X.
  // These are the height levels (in m), calculated from the sigma levels
  // which are directly read in the WRF file.
  GeneralGrid<real, 3> GridZ_in(shape(Nz_in, Ny_in, Nx_in),
                                1, shape(1, 2, 3));
  GeneralGrid<real, 3> GridZ_interfX_in(shape(Nz_in, Ny_in, Nx_in + 1),
                                        1, shape(1, 2, 3));
  GeneralGrid<real, 3> GridZ_interfY_in(shape(Nz_in, Ny_in + 1, Nx_in),
                                        1, shape(1, 2, 3));
  GeneralGrid<real, 3> GridZ_Dot_in(shape(Nz_in, Ny_in + 1, Nx_in + 1),
                                    1, shape(1, 2, 3));
  GeneralGrid<real, 3> GridZ_interfZ_in(shape(Nz_in + 1, Ny_in, Nx_in),
                                        1, shape(1, 2, 3));

  // Vertical levels are shared.
  GridZ_in.SetDuplicate(false);
  GridZ_interfX_in.SetDuplicate(false);
  GridZ_interfY_in.SetDuplicate(false);
  GridZ_interfZ_in.SetDuplicate(false);

  cout << endl;


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


  // Horizontal grids in WRF indices.
  GeneralGrid<real, 2> GridX_2D_out(shape(Ny_out, Nx_out), 2, shape(1, 2));
  GeneralGrid<real, 2> GridY_2D_out(shape(Ny_out, Nx_out), 1, shape(1, 2));
  GeneralGrid<real, 2> GridX_3D_out(shape(Ny_out, Nx_out), 3, shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_out(shape(Ny_out, Nx_out), 2, shape(2, 3));
  GeneralGrid<real, 2> GridX_3D_Gen_out(shape(Ny_out + 1, Nx_out), 3,
                                        shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_Gen_out(shape(Ny_out, Nx_out + 1), 2,
                                        shape(2, 3));
  GeneralGrid<real, 2> GridX_3D_interf_out(shape(Ny_out, Nx_out + 1), 3,
                                           shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_interf_out(shape(Ny_out + 1, Nx_out), 2,
                                           shape(2, 3));

  GridX_2D_out.SetDuplicate(false);
  GridY_2D_out.SetDuplicate(false);
  GridX_3D_out.SetDuplicate(false);
  GridY_3D_out.SetDuplicate(false);
  GridX_3D_Gen_out.SetDuplicate(false);
  GridY_3D_Gen_out.SetDuplicate(false);
  GridX_3D_interf_out.SetDuplicate(false);
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
                                         GridY_3D_out, GridX_3D_out);
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
                               GridY_3D_out, GridX_3D_out);
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
  // Fields defined on Z_interface
  Data<real, 4> Kz_out(GridT_out, GridZ_interf_out,
                       GridY_3D_out, GridX_3D_out);

  // 2D Output fields.
  Data<real, 3> SurfacePressure_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfaceTemperature_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SkinTemperature_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfaceRichardson_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> BoundaryHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> CloudBaseHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> CloudTopHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FrictionModule_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SolarRadiation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> ConvectiveRain_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> Rain_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SoilWater_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SnowHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> CanopyWater_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SensibleHeat_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> Evaporation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FirstLevelWindModule_out(GridT_out, GridY_2D_out,
                                         GridX_2D_out);
  Data<real, 3> WindModule10_out(GridT_out, GridY_2D_out, GridX_2D_out);

  // Post interpolation computed fields.
  Data<real, 3> PARdb_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> PARdiff_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> PAR_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfacePotentialTemperature_out(GridT_out,
                                                GridY_2D_out, GridX_2D_out);
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

  // Data to convert z from sigma to heights.
  Data<real, 2> SigmaH(Nt_in, Nz_in);
  // Terrain is 2D Field read from WRF file.
  Data<real, 3> Terrain(GridT_in, GridY_in, GridX_in);
  // Terrain at interfaces are interpolated from Terrain at cross points.
  Data<real, 3> Terrain_interfX(GridT_in, GridY_in, GridX_interf_in);
  Data<real, 3> Terrain_interfY(GridT_in, GridY_interf_in, GridX_in);
  // Terrain at dot points
  Data<real, 3> Terrain_Dot(GridT_in, GridY_interf_in, GridX_interf_in);

  // Reads data to convert Z from sigma levels to heights.
  InputFile.Read(file_in, "ZNU", SigmaH);
  InputFile.Read(file_in, "HGT", Terrain);

  //Interpolations
  LinearInterpolationRegular(Terrain, Terrain_interfX);
  LinearInterpolationRegular(Terrain, Terrain_interfY);
  LinearInterpolationRegular(Terrain, Terrain_Dot);

  cout << endl;

  // Converts Z from sigma levels to heights at cross points.

  // General constants.
  // See WRF documentations.

  cout << "  + Computing physical constants...";
  cout.flush();

  // Some are directly read in the WRF file.
  real Tiso, Ts0, p00, A, ptop;
  // Some others are calculated from the previous ones.
  real Piso, Ziso, aux, ps0, refpress;
  // Data read in the WRF file to obtain the values of constants.
  Data<real, 1> PTOP(Nt_in);
  Data<real, 1> TISO(Nt_in);
  Data<real, 1> LAPSE(Nt_in);
  Data<real, 1> P00(Nt_in);
  Data<real, 1> TS0(Nt_in);

  InputFile.Read(file_in, "P_TOP", PTOP);
  ptop = PTOP(0);
  InputFile.Read(file_in, "TISO", TISO);
  Tiso = TISO(0);
  InputFile.Read(file_in, "P00", P00);
  p00 = P00(0);
  InputFile.Read(file_in, "T00", TS0);
  Ts0 = TS0(0);
  InputFile.Read(file_in, "TLP", LAPSE);
  A = LAPSE(0);

  // Default Values if no value is found in the WRF file.
  if (A == 0.)
    A = 50.;
  if (p00 == 0.)
    p00 = 100000.;
  if (Tiso == 0.)
    Ts0 = 290.;
  Piso = p00 * exp((Tiso - Ts0) / A);
  aux = log(Piso / p00);
  Ziso = - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g;


  // At cross points.
  cout << "  + Computing vertical levels in meter at cross points...";
  cout.flush();

  for (k = 0; k < Nz_in; k++)
    for (j = 0; j < Ny_in; j++)
      for (i = 0; i < Nx_in; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_in.Value(0, k, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_in.Value(0, k, j, i)
                = Ziso - r * Tiso * aux / g - Terrain(0, j, i);
            }
        }
  cout << endl;

  // At X interfaces.
  cout << "  + Computing vertical levels in meter at x-interface...";
  cout.flush();

  for (k = 0; k < Nz_in; k++)
    for (j = 0; j < Ny_in; j++)
      for (i = 0; i < Nx_in + 1; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain_interfX(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_interfX_in.Value(0, k, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain_interfX(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_interfX_in.Value(0, k, j, i)
                = Ziso - r * Tiso * aux / g - Terrain_interfX(0, j, i);
            }
        }
  cout << endl;

  // At Y interfaces.
  cout << "  + Computing vertical levels in meter at y-interfaces...";
  cout.flush();

  for (k = 0; k < Nz_in; k++)
    for (j = 0; j < Ny_in + 1; j++)
      for (i = 0; i < Nx_in; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain_interfY(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_interfY_in.Value(0, k, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain_interfY(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_interfY_in.Value(0, k, j, i)
                = Ziso - r * Tiso * aux / g - Terrain_interfY(0, j, i);
            }
        }
  cout << endl;

  // At dot points.
  cout << "  + Computing vertical levels in meter at dot points...";
  cout.flush();

  for (k = 0; k < Nz_in; k++)
    for (j = 0; j < Ny_in + 1; j++)
      for (i = 0; i < Nx_in + 1; i++)
        {
          ps0 = p00 * exp(-Ts0 / A + sqrt(Ts0 * Ts0 / (A * A)
                                          - 2. * g * Terrain_Dot(0, j, i)
                                          / (A * r))) - ptop;
          refpress = SigmaH(0, k) * ps0 + ptop;
          if (refpress >= Piso)
            {
              aux = log(refpress / p00);
              GridZ_Dot_in.Value(0, k, j, i) =
                - r * A * aux * aux / (2.0 * g) - r * Ts0 * aux / g
                - Terrain_Dot(0, j, i);
            }
          else
            {
              aux = log(refpress / Piso);
              GridZ_Dot_in.Value(0, k, j, i)
                = Ziso - r * Tiso * aux / g - Terrain_Dot(0, j, i);
            }
        }
  cout << " done." << endl;

  // At Z interfaces.
  cout << "  + Computing vertical interfaces in meter at cross points...";
  cout.flush();

  for (j = 0; j < Ny_in; j++)
    for (i = 0; i < Nx_in; i++)
      GridZ_interfZ_in.Value(0, 0, j, i) = 0.;
  for (k = 1; k < Nz_in + 1; k++)
    for (j = 0; j < Ny_in; j++)
      for (i = 0; i < Nx_in; i++)
        GridZ_interfZ_in.Value(0, k, j, i) =
          2 * GridZ_in.Value(0, k - 1, j, i)
          - GridZ_interfZ_in.Value(0, k - 1, j, i);
  cout << endl;

  // Frees memory.
  Terrain.Resize();
  Terrain_interfX.Resize();
  Terrain_interfY.Resize();
  Terrain_Dot.Resize();


  ///////////////////////////////////////////////////////
  // Converts from lat/lon output grids to WRF indices //
  ///////////////////////////////////////////////////////


  cout << "Converting from latlon to WRF indices";
  cout.flush();

  float lat_ref, lon_ref, lon_cen;
  float true_lat1, true_lat2;
  float grid_x_dist, grid_y_dist;

  InputFile.ReadAttribute(file_in, "CEN_LAT", lat_ref);
  InputFile.ReadAttribute(file_in, "STAND_LON", lon_ref);
  InputFile.ReadAttribute(file_in, "CEN_LON", lon_cen);
  InputFile.ReadAttribute(file_in, "TRUELAT1", true_lat1);
  InputFile.ReadAttribute(file_in, "TRUELAT2", true_lat2);
  InputFile.ReadAttribute(file_in, "DX", grid_x_dist);
  InputFile.ReadAttribute(file_in, "DY", grid_y_dist);

  if (projection_type == 1)
    // Converts from latitude/longitude to WRF indices in Lambert conformal
    // conic projection.
    {
      cout << " (Lambert conformal conic projection)...";
      cout.flush();
      // Projection used to determinate the offset between the center (lon_cen, lat_ref)
      // of the domain and its reference (lon_ref, lat_ref).
      LonlatToWRFLccInd<float> LccWRFOffset(0, 0, lon_ref, lat_ref, 0, 0,
                                            true_lat1, true_lat2,
                                            grid_x_dist, grid_y_dist);

      // index (i, j) of the center of the domain relative to the reference.
      float i_cen, j_cen;
      LccWRFOffset(lon_cen, lat_ref, i_cen, j_cen);

      // Declares attributes of the projection.
      LonlatToWRFLccInd<float> LccWRF(Nx_in, Ny_in, lon_ref, lat_ref,
                                      i_cen, j_cen, true_lat1, true_lat2,
                                      grid_x_dist, grid_y_dist);

      int isinthegridy, isinthegridx;
      // 3D output grids.
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            LccWRF(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
                   GridX_3D_out.Value(0, 0, j, i),
                   GridY_3D_out.Value(0, 0, j, i));
            isinthegridx = 0;
            isinthegridy = 0;
            if (GridY_3D_out.Value(0, 0, j, i) >= GridY_in(0))
              if (GridY_3D_out.Value(0, 0, j, i) <= GridY_in(Ny_in - 1))
                isinthegridy = 1;
            if (GridX_3D_out.Value(0, 0, j, i) >= GridX_in(0))
              if (GridX_3D_out.Value(0, 0, j, i) <= GridX_in(Nx_in - 1))
                isinthegridx = 1;
            if (isinthegridx == 0 || isinthegridy == 0)
              cout << " coordinates " <<  x_min_out + Delta_x_out * i << " " <<
                y_min_out + Delta_y_out * j << " not found " <<  endl;
          }

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
          LccWRF(x_min_out + Delta_x_out * (real(i) - 0.5),
                 y_min_out + Delta_y_out * j,
                 GridX_3D_interf_out.Value(0, 0, j, i),
                 GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          LccWRF(x_min_out + Delta_x_out * i,
                 y_min_out + Delta_y_out * (real(j) - 0.5),
                 GridX_3D_Gen_out.Value(0, 0, j, i),
                 GridY_3D_interf_out.Value(0, 0, j, i));
    }
  else if (projection_type == 2)
    // Converts from latitude/longitude to WRF indices in stereographic
    // projection.
    {
      cout << " (Polar stereographic projection)...";
      cout.flush();

      LonlatToWRFStereInd<float> StereoWRF(Nx_in, Ny_in, lon_cen,
                                           lat_ref, true_lat1,
                                           grid_x_dist, grid_y_dist);

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          StereoWRF(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
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
          StereoWRF(x_min_out + Delta_x_out * (real(i) - 0.5),
                    y_min_out + Delta_y_out * j,
                    GridX_3D_interf_out.Value(0, 0, j, i),
                    GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          StereoWRF(x_min_out + Delta_x_out * i,
                    y_min_out + Delta_y_out * (real(j) - 0.5),
                    GridX_3D_Gen_out.Value(0, 0, j, i),
                    GridY_3D_interf_out.Value(0, 0, j, i));
    }
  else if (projection_type == 3)
    // Converts from latitude/longitude to WRF indices in Mercator projection.
    {
      cout << " (Mercator projection)...";
      cout.flush();
      // Mercator projection.
      LonlatToWRFMercInd<float> MercWRF(Nx_in, Ny_in, lon_cen, lat_ref,
                                        true_lat1, grid_x_dist, grid_y_dist);

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          MercWRF(x_min_out + Delta_x_out * i, y_min_out + Delta_y_out * j,
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
          MercWRF(x_min_out + Delta_x_out * (real(i) - 0.5),
                  y_min_out + Delta_y_out * j,
                  GridX_3D_interf_out.Value(0, 0, j, i),
                  GridY_3D_Gen_out.Value(0, 0, j, i));
      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          MercWRF(x_min_out + Delta_x_out * i,
                  y_min_out + Delta_y_out * (real(j) - 0.5),
                  GridX_3D_Gen_out.Value(0, 0, j, i),
                  GridY_3D_interf_out.Value(0, 0, j, i));
    }
  else if (projection_type == 6)
    // Converts from latitude/longitude to  latitude/longitude.
    {
      cout << " (lat-lon projection)...";
      cout.flush();

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            GridX_3D_out.Value(0, 0, j, i) = x_min_out + Delta_x_out * i;
            GridY_3D_out.Value(0, 0, j, i) = y_min_out + Delta_y_out * j;
          }

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
          {
            GridX_3D_interf_out.Value(0, 0, j, i) = x_min_out + Delta_x_out * (real(i) - 0.5);
            GridY_3D_Gen_out.Value(0, 0, j, i) = y_min_out + Delta_y_out * j;
          }

      for (j = 0; j < Ny_out + 1; j++)
        for (i = 0; i < Nx_out; i++)
          {
            GridX_3D_Gen_out.Value(0, 0, j, i) = x_min_out + Delta_x_out * i;
            GridY_3D_interf_out.Value(0, 0, j, i) = y_min_out + Delta_y_out * (real(j) - 0.5);
          }

      if (GridX_in(0) > GridX_2D_out.Value(0, 0, 0))
        throw("Values of Grid X out are too low");
      if (GridY_in(0) > GridY_2D_out.Value(0, 0, 0))
        throw("Values of Grid Y out are too low");
      if (GridX_in(Nx_in - 1) < GridX_2D_out.Value(0, Ny_out - 1, Nx_out - 1))
        throw("Values of Grid X out are too high");
      if (GridY_in(Ny_in - 1) < GridY_2D_out.Value(0, Ny_out - 1, Nx_out - 1))
        throw("Values of Grid Y out are too high");
    }

  cout << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  // Common fields to photolysis and meteo computation

  //// Pressure. ////

  // See WRF documentation, Section 5.
  cout << "Computing Pressure...";
  cout.flush();

  Data<real, 4> PressurePerturbation(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> BaseStatePressure(GridT_in, GridZ_in, GridY_in, GridX_in);

  InputFile.Read(file_in, "PB", BaseStatePressure);
  InputFile.Read(file_in, "P", PressurePerturbation);

  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          Pressure(h, k, j, i) = BaseStatePressure(h, k, j, i)
            + PressurePerturbation(h, k, j, i);

  // Frees memory.
  PressurePerturbation.Resize();
  BaseStatePressure.Resize();

  // Interpolations.
  Data<real, 4> Pressure_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);

  LinearInterpolationDimension(Pressure, Pressure_tmp, 1);
  LinearInterpolationRegularToGeneral(Pressure_tmp, Pressure_out);

  cout << endl;

  //// Surface Pressure. ////

  cout << "Computing Surface Pressure...";
  cout.flush();

  Data<real, 3> SurfacePressure(GridT_in, GridY_in, GridX_in);
  InputFile.Read(file_in, "PSFC", SurfacePressure);

  LinearInterpolationRegularToGeneral(SurfacePressure, SurfacePressure_out);
  cout << endl;

  //// Potential Temperature. ////

  cout << "Computing Potential Temperature ...";
  cout.flush();

  // Some physics constants.
  const real Theta0(300.);

  // Declaring data.
  Data<real, 4> PotentialTemperaturePerturbation(GridT_in, GridZ_in,
                                                 GridY_in, GridX_in);
  Data<real, 4> PotentialTemperature(GridT_in, GridZ_in,
                                     GridY_in, GridX_in);

  // Reading fields.
  InputFile.Read(file_in, "T", PotentialTemperaturePerturbation);

  // Applying transformation.
  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          PotentialTemperature(h, k, j, i) =
            PotentialTemperaturePerturbation(h, k, j, i) + Theta0;

  // Frees memory.
  PotentialTemperaturePerturbation.Resize();

  // Interpolations.
  Data<real, 4> PotentialTemperature_tmp(GridT_in, GridZ_out,
                                         GridY_in, GridX_in);

  LinearInterpolationDimension(PotentialTemperature,
                               PotentialTemperature_tmp, 1);
  LinearInterpolationRegularToGeneral(PotentialTemperature_tmp,
                                      PotentialTemperature_out);

  cout << endl;

  //// Temperature. ////

  cout << "Computing Temperature...";
  cout.flush();

  Data<real, 4> Temperature(GridT_in, GridZ_in, GridY_in, GridX_in);

  ComputeTemperature(PotentialTemperature, Pressure, Temperature);

  Data<real, 4> Temperature_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);

  LinearInterpolationDimension(Temperature, Temperature_tmp, 1);
  LinearInterpolationRegularToGeneral(Temperature_tmp, Temperature_out);
  Temperature_tmp.Resize();

  cout << endl;

  //// Clouds. ////

  cout << " Computing relative humidity and critical relative humidity...";

  // Reads specific humidity.
  Data<real, 4> SpecificHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  InputFile.Read(file_in, "QVAPOR", SpecificHumidity);
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
  cout << endl;

  cout << "  + Computing cloud profile...";
  cout.flush();

  Data<real, 3> BoundaryHeight(GridT_in, GridY_in, GridX_in);
  InputFile.Read(file_in, "PBLH", BoundaryHeight);

  Data<real, 4> CloudFraction(GridT_in, GridZ_in, GridY_in, GridX_in);
  ComputeCloudFraction(BoundaryHeight, RelativeHumidity, CRH, CloudFraction);

  Data<int, 4> LowIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> MediumIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> HighIndices(Nt_in, Ny_in, Nx_in, 2);

  Data<real, 3> LowCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness(GridT_in, GridY_in, GridX_in);

  ComputeCloudiness(CloudFraction, Pressure, GridZ_interfZ_in, LowIndices,
                    MediumIndices, HighIndices, LowCloudiness,
                    MediumCloudiness, HighCloudiness);
  Pressure.Resize();

  Data<real, 3> CloudBaseHeight(GridT_in, GridY_in, GridX_in);
  ComputeCloudBaseHeight(LowIndices, MediumIndices, HighIndices,
                         GridZ_interfZ_in, CloudBaseHeight);
  Data<real, 3> CloudTopHeight(GridT_in, GridY_in, GridX_in);
  ComputeCloudTopHeight(LowIndices, MediumIndices, HighIndices,
                        GridZ_interfZ_in, CloudTopHeight);

  cout << endl;

  // Reads liquid water content.
  Data<real, 4> LiquidWaterContent(GridT_in, GridZ_in,
                                   GridY_in, GridX_in);
  InputFile.Read(file_in, "QCLOUD", LiquidWaterContent);

  LiquidWaterContent.ThresholdMin(0.);


  // Fields used by meteo only

  Data<real, 3> SkinTemperature(GridT_in, GridY_in, GridX_in);
  if (with_meteo)
    {
      //// Winds. ////

      cout << "Computing Winds..." << endl;

      // Winds are given at x and y interfaces. Z at interfaces are differents
      //from Z at cross.
      Data<real, 4> ZonalWind(GridT_in, GridZ_interfX_in,
                              GridY_in, GridX_interf_in);
      Data<real, 4> MeridionalWind(GridT_in, GridZ_interfY_in,
                                   GridY_interf_in, GridX_in);

      // Wind rotation coefficients at interfaces.
      // It requires longitudes and latitudes at interfaces.
      // Rotations are calculated at interfaces. To do so, ZonalWind
      // is interpolated at x-interface and MeridionalWind at y-interface.
      cout << "  + Computing local Cosine and Sine fields...";
      cout.flush();

      Data<real, 3> Longitude_interfX(Nt_in, Ny_in, Nx_in + 1);
      Data<real, 3> Latitude_interfX(Nt_in, Ny_in, Nx_in + 1);
      Data<real, 3> Longitude_interfY(Nt_in, Ny_in + 1, Nx_in);
      Data<real, 3> Latitude_interfY(Nt_in, Ny_in + 1, Nx_in);

      InputFile.Read(file_in, "XLONG_U", Longitude_interfX);
      InputFile.Read(file_in, "XLAT_U", Latitude_interfX);
      InputFile.Read(file_in, "XLONG_V", Longitude_interfY);
      InputFile.Read(file_in, "XLAT_V", Latitude_interfY);

      Data<real, 2> Cosine_interfX(Ny_in, Nx_in + 1);
      Data<real, 2> Sine_interfX(Ny_in, Nx_in + 1);
      Data<real, 2> Cosine_interfY(Ny_in + 1, Nx_in);
      Data<real, 2> Sine_interfY(Ny_in + 1, Nx_in);
      real delta_x, delta_y, delta_x2, delta_y2;

      for (j = 1; j < Ny_in + 1; j++)
        {
          for (i = 1; i < Nx_in; i++)
            {
              // Cosine and Sine at x-interface.
              delta_x = (Longitude_interfY(0, j - 1, i)
                         - Longitude_interfY(0, j - 1, i - 1))
                * cos(pi / 180. * Latitude_interfY(0, j - 1, i - 1));
              delta_y = Latitude_interfY(0, j - 1, i)
                - Latitude_interfY(0, j - 1, i - 1);
              delta_x2 = delta_x * delta_x;
              delta_y2 = delta_y * delta_y;
              Cosine_interfX(j - 1, i) = delta_x / sqrt(delta_x2 + delta_y2);
              Sine_interfX(j - 1, i) = delta_y / sqrt(delta_x2 + delta_y2);

              // At y-interface.
              delta_x = (Longitude_interfX(0, j - 1, i + 1)
                         - Longitude_interfX(0, j - 1, i))
                * cos(pi / 180. * Latitude_interfX(0, j - 1, i));
              delta_y = Latitude_interfX(0, j - 1, i + 1)
                - Latitude_interfX(0, j - 1, i);
              delta_x2 = delta_x * delta_x;
              delta_y2 = delta_y * delta_y;
              Cosine_interfY(j, i) = delta_x / sqrt(delta_x2 + delta_y2);
              Sine_interfY(j, i) = delta_y / sqrt(delta_x2 + delta_y2);
            }
          // At y-interface for i = 0.
          delta_x = (Longitude_interfX(0, j - 1, 1)
                     - Longitude_interfX(0, j - 1, 0))
            * cos(pi / 180. * Latitude_interfX(0, j - 1, 0));
          delta_y = Latitude_interfX(0, j - 1, 1)
            - Latitude_interfX(0, j - 1, 0);
          delta_x2 = delta_x * delta_x;
          delta_y2 = delta_y * delta_y;
          Cosine_interfY(j, 0) = delta_x / sqrt(delta_x2 + delta_y2);
          Sine_interfY(j, 0) = delta_y / sqrt(delta_x2 + delta_y2);
        }

      // Extrapolates on boundaries.
      for (i = 0; i < Nx_in; i++)
        {
          Cosine_interfY(0, i) = 2. * Cosine_interfY(1, i) -
            Cosine_interfY(2, i);
          Sine_interfY(0, i) = 2. * Sine_interfY(1, i) -
            Sine_interfY(2, i);
        }
      for (j = 0; j < Ny_in; j++)
        {
          Cosine_interfX(j, 0) = 2. * Cosine_interfX(j, 1)
            - Cosine_interfX(j, 2);
          Sine_interfX(j, 0) = 2. * Sine_interfX(j, 1)
            - Sine_interfX(j, 2);
          Cosine_interfX(j, Nx_in) = 2. * Cosine_interfX(j, Nx_in - 1)
            - Cosine_interfX(j, Nx_in - 2);
          Sine_interfX(j, Nx_in) = 2. * Sine_interfX(j, Nx_in - 1)
            - Sine_interfX(j, Nx_in - 2);
        }

      cout << endl;
      cout << "  + Computing winds rotation...";
      cout.flush();

      InputFile.Read(file_in, "V", MeridionalWind);
      InputFile.Read(file_in, "U", ZonalWind);

      // Vertical interpolation of winds so that they are defined only
      // on regular grids.
      Data<real, 4> ZonalWind_tmp(GridT_in, GridZ_out,
                                  GridY_in, GridX_interf_in);
      Data<real, 4> MeridionalWind_tmp(GridT_in, GridZ_out,
                                       GridY_interf_in, GridX_in);
      LinearInterpolationDimension(ZonalWind, ZonalWind_tmp, 1);
      LinearInterpolationDimension(MeridionalWind, MeridionalWind_tmp, 1);

      // Interpolations of winds at interfaces to perform rotations.
      Data<real, 4> ZonalWind_interfY(GridT_in, GridZ_out,
                                      GridY_interf_in, GridX_in);
      Data<real, 4> MeridionalWind_interfX(GridT_in, GridZ_out,
                                           GridY_in, GridX_interf_in);
      LinearInterpolationRegular(ZonalWind_tmp, ZonalWind_interfY);
      LinearInterpolationRegular(MeridionalWind_tmp, MeridionalWind_interfX);

      // Wind rotations.
      real meridional_wind, zonal_wind;
      for (h = 0; h < Nt_in; h++)
        for (k = 0; k < Nz_out; k++)
          {
            for (j = 0; j < Ny_in + 1; j++)
              for (i = 0; i < Nx_in; i++)
                {
                  meridional_wind = MeridionalWind_tmp(h, k, j, i);
                  zonal_wind = ZonalWind_interfY(h, k, j, i);
                  MeridionalWind_tmp(h, k, j, i) = Sine_interfY(j, i)
                    * zonal_wind + Cosine_interfY(j, i) * meridional_wind;
                }
            for (j = 0; j < Ny_in; j++)
              for (i = 0; i < Nx_in + 1; i++)
                {
                  meridional_wind = MeridionalWind_interfX(h, k, j, i);
                  zonal_wind = ZonalWind_tmp(h, k, j, i);
                  ZonalWind_tmp(h, k, j, i) = Cosine_interfX(j, i)
                    * zonal_wind - Sine_interfX(j, i) * meridional_wind;
                }
          }

      // Frees memory.
      MeridionalWind_interfX.Resize();
      ZonalWind_interfY.Resize();
      Sine_interfX.Resize();
      Sine_interfY.Resize();
      Cosine_interfX.Resize();
      Cosine_interfY.Resize();
      Longitude_interfX.Resize();
      Longitude_interfY.Resize();
      Latitude_interfX.Resize();
      Latitude_interfY.Resize();

      // Interpolations at, respectively, y- and x- interfaces.
      LinearInterpolationRegularToGeneral(ZonalWind_tmp, ZonalWind_out);
      LinearInterpolationRegularToGeneral(MeridionalWind_tmp,
                                          MeridionalWind_out);

      // Frees memory.
      MeridionalWind.Resize();
      ZonalWind.Resize();

      cout << endl;


      ////////////////////////
      // RICHARDSON NUMBERS //
      ////////////////////////

      cout << "Computing Richardson numbers...";
      cout.flush();

      Data<real, 3> SurfacePotentialTemperature(GridT_in, GridY_in, GridX_in);
      //    Data<real, 3> SkinTemperature(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "TSK", SkinTemperature);
      cout << "" << endl;

      ComputePotentialTemperature(SkinTemperature, SurfacePressure,
                                  SurfacePotentialTemperature);

      Data<real, 4> Richardson_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);
      Data<real, 3> SurfaceRichardson(GridT_in, GridY_in, GridX_in);

      Data<real, 4> MeridionalWind_cross(GridT_in, GridZ_out,
                                         GridY_in, GridX_in);
      Data<real, 4> ZonalWind_cross(GridT_in, GridZ_out, GridY_in, GridX_in);
      Data<real, 3> FirstLevelWindModule(GridT_in, GridY_in, GridX_in);

      // Interpolates ZonalWind and Meridional wind at cross points.
      LinearInterpolationRegular(ZonalWind_tmp, ZonalWind_cross);
      LinearInterpolationRegular(MeridionalWind_tmp, MeridionalWind_cross);


      // Frees memory.
      ZonalWind_tmp.Resize();
      MeridionalWind_tmp.Resize();

      // Richardson number.
      ComputeRichardson(ZonalWind_cross, MeridionalWind_cross,
                        PotentialTemperature_tmp, Richardson_tmp);
      ComputeModule(ZonalWind_cross, MeridionalWind_cross,
                    FirstLevelWindModule);
      ComputeRichardson(FirstLevelWindModule, SurfacePotentialTemperature,
                        PotentialTemperature, SurfaceRichardson);

      // Frees memory.
      PotentialTemperature.Resize();
      PotentialTemperature_tmp.Resize();
      SurfacePotentialTemperature.Resize();
      Pressure_tmp.Resize();

      LinearInterpolationRegularToGeneral(FirstLevelWindModule,
                                          FirstLevelWindModule_out);
      FirstLevelWindModule.Resize();
      ZonalWind_cross.Resize();
      MeridionalWind_cross.Resize();

      LinearInterpolationRegularToGeneral(SurfaceRichardson,
                                          SurfaceRichardson_out);
      SurfaceRichardson.Resize();

      LinearInterpolationRegularToGeneral(Richardson_tmp, Richardson_out);
      Richardson_tmp.Resize();

      cout << endl;
    }


  // Fields used by photolysis only
  Data<real, 4> Attenuation(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> LiquidWaterExtinction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> IceWaterExtinction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> IceWaterContent(GridT_in, GridZ_in, GridY_in, GridX_in);

  if (with_photolysis)
    {
      if (ice_cloud &&
          (photolysis_tabulation == 2 || photolysis_tabulation == 3))
        InputFile.Read(file_in, "QICE", IceWaterContent);
      else
        IceWaterContent.Fill(0.);
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

  // Interpolates 2D and 3D fields with data.

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

  Data<real, 4> SpecificHumidity_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);

  LinearInterpolationDimension(SpecificHumidity, SpecificHumidity_tmp, 1);
  LinearInterpolationRegularToGeneral(SpecificHumidity_tmp,
                                      SpecificHumidity_out);

  SpecificHumidity_tmp.Resize();
  SpecificHumidity.Resize();

  cout << "  + SpecificHumidity" << endl;

  // Test.

  Data<real, 4> RelativeHumidity_tmp(GridT_in, GridZ_out,
                                     GridY_in, GridX_in);
  Data<real, 4> RelativeHumidity_out(GridT_out, GridZ_out,
                                     GridY_out, GridX_out);

  LinearInterpolationDimension(RelativeHumidity, RelativeHumidity_tmp, 1);
  LinearInterpolationRegularToGeneral(RelativeHumidity_tmp,
                                      RelativeHumidity_out);
  RelativeHumidity_tmp.Resize();
  RelativeHumidity.Resize();

  cout << "  + RelativeHumidity" << endl;

  Data<real, 4> CRH_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);
  Data<real, 4> CRH_out(GridT_out, GridZ_out, GridY_out, GridX_out);

  LinearInterpolationDimension(CRH, CRH_tmp, 1);
  LinearInterpolationRegularToGeneral(CRH_tmp,
                                      CRH_out);

  CRH_tmp.Resize();
  CRH.Resize();

  cout << "  + CRH" << endl;


  if (with_meteo)
    {
      // ClouBasedHeight and CloudTopHeight.
      LinearInterpolationRegularToGeneral(CloudBaseHeight, CloudBaseHeight_out);
      LinearInterpolationRegularToGeneral(CloudTopHeight, CloudTopHeight_out);

      CloudBaseHeight_out.ThresholdMin(min_height);
      CloudBaseHeight.Resize();

      CloudTopHeight_out.ThresholdMin(min_height);
      CloudTopHeight.Resize();
      cout << "  + CloudBaseHeight CloudTopHeight" << endl;

      // SkinTemperature.
      LinearInterpolationRegularToGeneral(SkinTemperature,
                                          SkinTemperature_out);
      SkinTemperature.Resize();
      cout << "  + SkinTemperature" << endl;

      // SensibleHeat.
      Data<real, 3> SensibleHeat(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "HFX", SensibleHeat);
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
      Data<real, 3> LatentHeat(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "LH", LatentHeat);
      LinearInterpolationRegularToGeneral(LatentHeat, Evaporation_out);
      // Divides by \rho_{water} * L.
      Evaporation_out.Mlt(1. / 2.5e9);
      LatentHeat.Resize();
      cout << "  + Evaporation" << endl;

      // SurfaceTemperature.
      Data<real, 3> SurfaceTemperature(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "T2", SurfaceTemperature);
      LinearInterpolationRegularToGeneral(SurfaceTemperature,
                                          SurfaceTemperature_out);
      SurfaceTemperature.Resize();
      cout << "  + SurfaceTemperature" << endl;

      // SoilWater.
      Data<real, 4> SoilWater3D(GridT_in, GridSoil_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "SH2O", SoilWater3D);
      Data<real, 3> SoilWater(GridT_in, GridY_in, GridX_in);
      for (k = 0; k < Nt_in; k++)
        for (j = 0; j < Ny_in; j++)
          for (i = 0; i < Nx_in; i++)
            SoilWater(k, j, i) = SoilWater3D(k, Nsoil_in - 1, j, i);
      SoilWater3D.Resize();
      SoilWater.Threshold(0., 1.);
      LinearInterpolationRegularToGeneral(SoilWater, SoilWater_out);
      SoilWater.Resize();
      cout << "  + SoilWater" << endl;

      //// SnowHeight. ////

      Data<real, 3> SnowHeight(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "SNOWH", SnowHeight);

      LinearInterpolationRegularToGeneral(SnowHeight, SnowHeight_out);

      SnowHeight.Resize();

      cout << "  + SnowHeight" << endl;

      //// CanopyWater. ////

      Data<real, 3> CanopyWater(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "CANWAT", CanopyWater);

      LinearInterpolationRegularToGeneral(CanopyWater, CanopyWater_out);

      CanopyWater.Resize();

      cout << "  + CanopyWater" << endl;

      //// SolarRadiation. ////

      Data<real, 3> SolarRadiation(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "SWDOWN", SolarRadiation);

      LinearInterpolationRegularToGeneral(SolarRadiation, SolarRadiation_out);

      SolarRadiation.Resize();
      cout << "  + SolarRadiation" << endl;

      // ConvectiveRain.
      Data<real, 3> ConvectiveRain(GridT_in, GridY_in, GridX_in);
      Data<real, 3> ConvectiveRain_prev(GridT_in, GridY_in, GridX_in);

      InputFile.Read(file_in, "RAINC", ConvectiveRain);

      ConvectiveRain.ThresholdMin(0.);

      if (exists(file_in_prev) && prev_accumulated_rain)
        InputFile.Read(file_in_prev, "RAINC", ConvectiveRain_prev);
      else
        ConvectiveRain_prev.Fill(0.);

      ConvectiveRain_prev.ThresholdMin(0.);

      Decumulate(ConvectiveRain, Nt_in, 0);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          ConvectiveRain(0, j, i) -= ConvectiveRain_prev(Nt_in - 1, j, i);
      // Converts mm/step to mm/h.
      ConvectiveRain.Mlt(1. / Delta_t_in);

      //// Rain. ////

      // First reads NonConvectiveRain.
      Data<real, 3> Rain(GridT_in, GridY_in, GridX_in);
      Data<real, 3> Rain_prev(GridT_in, GridY_in, GridX_in);

      InputFile.Read(file_in, "RAINNC", Rain);
      Rain.ThresholdMin(0.);

      if (exists(file_in_prev) && prev_accumulated_rain)
        InputFile.Read(file_in_prev, "RAINNC", Rain_prev);
      else
        Rain_prev.Fill(0.);

      Rain_prev.ThresholdMin(0.);

      Decumulate(Rain, Nt_in, 0);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          Rain(0, j, i) -= Rain_prev(Nt_in - 1, j, i);
      // Converts mm/step to mm/h.
      Rain.Mlt(1. / Delta_t_in);

      // Adds convective rain to NonConvectiveRain.
      for (h = 0; h < Nt_in; h++)
        for (j = 0; j < Ny_in; j++)
          for (i = 0; i < Nx_in; i++)
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

      //// FrictionModule. ////

      Data<real, 3> FrictionModule(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "UST", FrictionModule);

      LinearInterpolationRegularToGeneral(FrictionModule, FrictionModule_out);
      FrictionModule.Resize();
      cout << "  + FrictionModule" << endl;

      // BoundaryHeight.
      LinearInterpolationRegularToGeneral(BoundaryHeight, BoundaryHeight_out);
      BoundaryHeight.Resize();
      BoundaryHeight_out.ThresholdMin(GridZ_interf_out(1));
      cout << "  + BoundaryHeight" << endl;

      //// Wind Module at 10m ////

      Data<real, 3> U10(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "U10", U10);

      Data<real, 3> V10(GridT_in, GridY_in, GridX_in);
      InputFile.Read(file_in, "V10", V10);

      Data<real, 3> WindModule10(GridT_in, GridY_in, GridX_in);

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
      // CloudFraction
      Data<real, 4> CloudFraction_tmp(GridT_in, GridZ_out, GridY_in, GridX_in);
      LinearInterpolationDimension(CloudFraction, CloudFraction_tmp, 1);
      LinearInterpolationRegularToGeneral(CloudFraction_tmp, CloudFraction_out);
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
      SnowHeight_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
      CanopyWater_out.ResizeGrid(GridT_out, GridY_out, GridX_out);
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

          cout << "CALL to FastJX" << endl;

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

      Date current_date = date_beg;
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

  cout << "directory_out=" << to_str(directory_out) << endl;

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
      OutputMeteo.Append(CloudTopHeight_out,
                         directory_out + "CloudTopHeight.bin");
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
      OutputMeteo.Append(SnowHeight_out, directory_out + "SnowHeight.bin");
      OutputMeteo.Append(CanopyWater_out, directory_out
                         + "CanopyWater.bin");
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
              OutputMeteo.Append(PhotolysisRate4D, directory_photolysis +
                                 photolysis_reaction_list[i] + ".bin");
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
  cout << "|   Successful completion of WRF-meteo   |" << endl;
  cout << "||||||||||||||||||||||||||||||||||||||||||" << endl << endl;

  END;

  return 0;
}

