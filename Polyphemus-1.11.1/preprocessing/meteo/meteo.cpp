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


//////////////
// INCLUDES //

#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_GRIB

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

#include "fastJX.hxx"


template <class T>
void exp_(T& x)
{
  x = exp(x);
}

template <class T>
void sqrt_(T& x)
{
  x = sqrt(x);
}

template <class T>
void abs_(T& x)
{
  x = abs(x);
}


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("meteo.cfg");
  Date date_beg;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, default_name);

  Date date_end = date_beg;
  date_end.AddDays(1);

  Date date_prev = date_beg;
  date_prev.AddDays(-1);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  int h, i, j, k;

  // Constants.
  const real pi = 3.14159265358979323846264;
  const real P0 = 101325.;
  const real r = 287.;
  const real cp = 1005.;
  const real g = 9.81;


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
  int Nt_in, Nz_in, Ny_in, Nx_in;
  real Delta_t_in, Delta_y_in, Delta_x_in;
  real t_min_in, y_min_in, x_min_in;

  configuration.SetSection("[ECMWF]");

  configuration.PeekValue("Nt", "> 0", Nt_in);
  configuration.PeekValue("Nz", "> 0", Nz_in);
  configuration.PeekValue("Ny", "> 0", Ny_in);
  configuration.PeekValue("Nx", "> 0", Nx_in);

  configuration.PeekValue("Delta_t", "> 0", Delta_t_in);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_in);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_in);

  configuration.PeekValue("t_min", t_min_in);
  configuration.PeekValue("y_min", y_min_in);
  configuration.PeekValue("x_min", x_min_in);

  // Output domain.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_t_out, Delta_y_out, Delta_x_out;
  real t_min_out, y_min_out, x_min_out;
  string vertical_levels;

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx_out);
  configuration.PeekValue("Ny", "> 0", Ny_out);
  configuration.PeekValue("Nz", "> 0", Nz_out);

  configuration.PeekValue("Delta_t", "> 0", Delta_t_out);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_out);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_out);

  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);
  configuration.PeekValue("Vertical_levels", vertical_levels);

  Nt_out = compute_Nt(date_beg, date_end, Delta_t_out);
  t_min_out = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    +  real(date_beg.GetSeconds()) / 3600.;

  // Files.
  string directory_in, directory_out, roughness_in_file;
  configuration.SetSection("[paths]");

  configuration.PeekValue("Database_meteo", directory_in);
  configuration.PeekValue("Directory_meteo", directory_out);
  configuration.PeekValue("Roughness_in", roughness_in_file);

  string file_in = directory_in + date_beg.GetDate("ECMWF-%y%m%d.grb");
  string file_in_prev = directory_in + date_prev.GetDate("ECMWF-%y%m%d.grb");


  // Files.
  bool richardson_with_roughness, with_meteo;

  configuration.SetSection("[meteo]");

  configuration.PeekValue("Compute_Meteo", with_meteo);

  configuration.PeekValue("Richardson_with_roughness",
                          richardson_with_roughness);

  // Accumulated data
  int accumulated_time, accumulated_index;

  configuration.SetSection("[accumulated_data]");

  configuration.PeekValue("Accumulated_time", "positive", accumulated_time);
  configuration.PeekValue("Accumulated_index", "positive", accumulated_index);


  // Photolysis rates
  int photolysis_tabulation, attenuation_type;;
  bool ice_cloud, with_photolysis;
  vector<string> photolysis_reaction_list;
  int Nr_photolysis, Nwavelength, Nlegendre;
  string directory_attenuation, directory_photolysis, fastJ_parameter_files;
  Array<char, 2> photolysis_specie_name(Nr_photolysis, 10);

  configuration.SetSection("[photolysis_rates]");
  configuration.PeekValue("Compute_Photolysis_Data", with_photolysis);
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

      photolysis_specie_name = ' ';
      for (int i = 0; i < Nr_photolysis; i++)
        for (int j = 0; j < (int) photolysis_reaction_list[i].size(); j++)
          photolysis_specie_name(i, j) = photolysis_reaction_list[i][j];
    }

  // Clouds
  real min_height;
  int   option_critical_relative_humidity;
  configuration.SetSection("[clouds]");
  configuration.PeekValue("Critical_relative_humidity", "= 1 2",
                          option_critical_relative_humidity);
  configuration.PeekValue("Min_height", "> 0", min_height);


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input settings.

  // Input grids.
  RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
  // Vertical levels depend on t, z, y and x.
  GeneralGrid<real, 4> GridZ_in(shape(Nt_in, Nz_in, Ny_in, Nx_in),
                                1, shape(0, 1, 2, 3));
  GeneralGrid<real, 4> GridZ_interf_in(shape(Nt_in, Nz_in + 1, Ny_in, Nx_in),
                                       1, shape(0, 1, 2, 3));
  RegularGrid<real> GridY_in(y_min_in, Delta_y_in, Ny_in);
  RegularGrid<real> GridX_in(x_min_in, Delta_x_in, Nx_in);
  // Vertical levels are shared.
  GridZ_in.SetVariable(1);
  GridZ_in.SetDuplicate(false);

  // Output settings.

  // Output grids.
  RegularGrid<real> GridT_out(t_min_out, Delta_t_out, Nt_out);
  RegularGrid<real> GridZ_out(Nz_out);
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  // Data may be provided on interfaces.
  RegularGrid<real> GridZ_interf_out(Nz_out + 1);
  RegularGrid<real> GridY_interf_out(y_min_out - Delta_y_out / 2.,
                                     Delta_y_out, Ny_out + 1);
  RegularGrid<real> GridX_interf_out(x_min_out - Delta_x_out / 2.,
                                     Delta_x_out, Nx_out + 1);

  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels, GridZ_interf_out);
  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_interf_out(k) + GridZ_interf_out(k + 1)) / 2.0;

  // Photolysis rates Grid
  Nwavelength = 4;
  Nlegendre = 8;
  RegularGrid<real> GridP(Nr_photolysis);
  RegularGrid<real> GridWavelength(Nwavelength);
  RegularGrid<real> GridLegendre(Nlegendre);

  cout << " done." << endl;

  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 2> Roughness(GridY_in, GridX_in);
  Data<real, 4> Temperature(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> SurfaceTemperature(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SkinTemperature(GridT_in, GridY_in, GridX_in);
  Data<real, 4> PotentialTemperature(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> SurfacePotentialTemperature(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SurfaceRichardson(GridT_in, GridY_in, GridX_in);
  Data<real, 4> Richardson(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Pressure_interf(GridT_in, GridZ_interf_in,
                                GridY_in, GridX_in);
  Data<real, 3> SurfacePressure(GridT_in, GridY_in, GridX_in);
  Data<real, 4> SpecificHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> RelativeHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> LowCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 4> MeridionalWind(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> ZonalWind(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> FirstLevelWindModule(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SolarRadiation(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SolarRadiation_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> BoundaryHeight(GridT_in, GridY_in, GridX_in);
  Data<real, 4> CloudFraction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> CRH(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> CloudHeight(GridT_in, GridY_in, GridX_in);
  Data<real, 3> U_star(GridT_in, GridY_in, GridX_in);
  Data<real, 3> V_star(GridT_in, GridY_in, GridX_in);
  Data<real, 3> U_star_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> V_star_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SoilWater(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SensibleHeat(GridT_in, GridY_in, GridX_in);
  Data<real, 3> SensibleHeat_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> Evaporation(GridT_in, GridY_in, GridX_in);
  Data<real, 3> Evaporation_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> ConvectiveRain(GridT_in, GridY_in, GridX_in);
  Data<real, 3> ConvectiveRain_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> LargeScaleRain(GridT_in, GridY_in, GridX_in);
  Data<real, 3> LargeScaleRain_prev(GridT_in, GridY_in, GridX_in);

  // Output fields.
  Data<real, 4> Temperature_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 3> SurfaceTemperature_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SkinTemperature_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SurfaceRichardson_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> Richardson_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> Pressure_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 3> SurfacePressure_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> SpecificHumidity_out(GridT_out, GridZ_out, GridY_out,
                                     GridX_out);
  Data<real, 4> LiquidWaterContent_out(GridT_out, GridZ_out, GridY_out,
                                       GridX_out);
  Data<real, 4> IceWaterContent_out(GridT_out, GridZ_out, GridY_out,
                                    GridX_out);
  Data<real, 4> Attenuation_out(GridT_out, GridZ_out, GridY_out,
                                GridX_out);
  Data<real, 4> LiquidWaterExtinction_out(GridT_out, GridZ_out,
                                          GridY_out, GridX_out);
  Data<real, 4> IceWaterExtinction_out(GridT_out, GridZ_out,
                                       GridY_out, GridX_out);
  Data<real, 4> MeridionalWind_out(GridT_out, GridZ_out,
                                   GridY_interf_out, GridX_out);
  Data<real, 4> ZonalWind_out(GridT_out, GridZ_out,
                              GridY_out, GridX_interf_out);
  Data<real, 3> BoundaryHeight_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> U_star_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> V_star_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> WindModule_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 3> FrictionModule_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SolarRadiation_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PARdb_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PARdiff_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PAR_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SoilWater_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SensibleHeat_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> Evaporation_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> FirstLevelWindModule_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> CloudHeight_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> ConvectiveRain_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> LargeScaleRain_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> MediumCloudiness_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> HighCloudiness_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> CloudFraction_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> OpticalDepth_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> IceOpticalDepth_out(GridT_out, GridZ_out,
                                    GridY_out, GridX_out);
  Data<real, 5> PhotolysisRate_out(GridP, GridT_out, GridZ_out,
                                   GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  FormatBinary<float> InputBinaryMeteo;
  FormatGrib InputMeteo;

  cout << "Extracting data...";
  cout.flush();

  // Common fields to photolysis and meteo computation
  InputMeteo.Read(file_in, 167, SurfaceTemperature);
  SurfaceTemperature.ReverseData(1);

  InputMeteo.Read(file_in, 152, SurfacePressure);
  SurfacePressure.ReverseData(1);

  InputMeteo.Read(file_in, 130, Temperature);
  Temperature.ReverseData(1);
  Temperature.ReverseData(2);

  InputMeteo.Read(file_in, 133, SpecificHumidity);
  SpecificHumidity.ReverseData(1);
  SpecificHumidity.ReverseData(2);

  InputMeteo.Read(file_in, 187, MediumCloudiness);
  MediumCloudiness.ReverseData(1);

  InputMeteo.Read(file_in, 188, HighCloudiness);
  HighCloudiness.ReverseData(1);

  InputMeteo.Read(file_in, 143, ConvectiveRain);
  ConvectiveRain.ReverseData(1);

  InputMeteo.Read(file_in, 159, BoundaryHeight);
  BoundaryHeight.ReverseData(1);


  // Fields used by meteo only
  if (with_meteo)
    {
      InputMeteo.Read(file_in, 235, SkinTemperature);
      SkinTemperature.ReverseData(1);

      InputMeteo.Read(file_in, 132, MeridionalWind);
      MeridionalWind.ReverseData(1);
      MeridionalWind.ReverseData(2);

      InputMeteo.Read(file_in, 131, ZonalWind);
      ZonalWind.ReverseData(1);
      ZonalWind.ReverseData(2);

      InputMeteo.Read(file_in, 180, U_star);
      U_star.ReverseData(1);

      InputMeteo.Read(file_in, 181, V_star);
      V_star.ReverseData(1);

      InputMeteo.Read(file_in_prev, 180, U_star_prev);
      U_star_prev.ReverseData(1);

      InputMeteo.Read(file_in_prev, 181, V_star_prev);
      V_star_prev.ReverseData(1);

      InputMeteo.Read(file_in, 169, SolarRadiation);
      SolarRadiation.ReverseData(1);

      InputMeteo.Read(file_in_prev, 169, SolarRadiation_prev);
      SolarRadiation_prev.ReverseData(1);

      InputMeteo.Read(file_in, 39, SoilWater);
      SoilWater.ReverseData(1);

      InputMeteo.Read(file_in, 146, SensibleHeat);
      SensibleHeat.ReverseData(1);

      InputMeteo.Read(file_in_prev, 146, SensibleHeat_prev);
      SensibleHeat_prev.ReverseData(1);

      InputMeteo.Read(file_in, 182, Evaporation);
      Evaporation.ReverseData(1);

      InputMeteo.Read(file_in_prev, 182, Evaporation_prev);
      Evaporation_prev.ReverseData(1);

      if (richardson_with_roughness)
        InputBinaryMeteo.Read(roughness_in_file, Roughness);

      InputMeteo.Read(file_in_prev, 143, ConvectiveRain_prev);
      ConvectiveRain_prev.ReverseData(1);

      InputMeteo.Read(file_in, 142, LargeScaleRain);
      LargeScaleRain.ReverseData(1);

      InputMeteo.Read(file_in_prev, 142, LargeScaleRain_prev);
      LargeScaleRain_prev.ReverseData(1);
    }

  // Transformations.

  SurfacePressure.Apply(exp_);

  SpecificHumidity.Threshold(0., 1.);

  // Transformations usefull for meteo only
  if (with_meteo)
    {
      ConvectiveRain.ThresholdMin(0.);
      ConvectiveRain_prev.ThresholdMin(0.);
      LargeScaleRain.ThresholdMin(0.);
      LargeScaleRain_prev.ThresholdMin(0.);

      SoilWater.Threshold(0., 1.);
      SolarRadiation.ThresholdMin(0.);
      SolarRadiation_prev.ThresholdMin(0.);

      Decumulate(U_star, accumulated_time, accumulated_index);
      Decumulate(V_star, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          {
            U_star(0, j, i) -= U_star_prev(Nt_in - 2, j, i);
            V_star(0, j, i) -= V_star_prev(Nt_in - 2, j, i);
          }
      U_star.Mlt(1. / (Delta_t_in * 3600.));
      V_star.Mlt(1. / (Delta_t_in * 3600.));

      Decumulate(ConvectiveRain, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          ConvectiveRain(0, j, i) -= ConvectiveRain_prev(Nt_in - 2, j, i);
      // To mm/h.
      ConvectiveRain.Mlt(1000. / Delta_t_in);

      Decumulate(LargeScaleRain, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          LargeScaleRain(0, j, i) -= LargeScaleRain_prev(Nt_in - 2, j, i);
      // To mm/h.
      LargeScaleRain.Mlt(1000. / Delta_t_in);

      // Stress to friction velocity.
      // Divides by \rho_{air}.
      for (h = 0; h < Nt_in; h++)
        for (j = 0; j < Ny_in; j++)
          for (i = 0; i < Nx_in; i++)
            {
              U_star(h, j, i) *= r * SkinTemperature(h, j, i)
                / SurfacePressure(h, j, i);
              V_star(h, j, i) *= r * SkinTemperature(h, j, i)
                / SurfacePressure(h, j, i);
            }

      Decumulate(SolarRadiation, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          SolarRadiation(0, j, i) -= SolarRadiation_prev(Nt_in - 2, j, i);
      SolarRadiation.Mlt(1. / (Delta_t_in * 3600.));
      SolarRadiation.ThresholdMin(0.);

      Decumulate(Evaporation, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          Evaporation(0, j, i) -= Evaporation_prev(Nt_in - 2, j, i);
      // Multiplies by \rho_{water} / \rho_{air}.
      for (h = 0; h < Nt_in; h++)
        for (j = 0; j < Ny_in; j++)
          for (i = 0; i < Nx_in; i++)
            Evaporation(h, j, i) *= 1000.* r * SkinTemperature(h, j, i)
              / SurfacePressure(h, j, i);
      // -1. to change the direction of the flux.
      Evaporation.Mlt(-1. / (Delta_t_in * 3600.));

      Decumulate(SensibleHeat, accumulated_time, accumulated_index);
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          SensibleHeat(0, j, i) -= SensibleHeat_prev(Nt_in - 2, j, i);
      // Divides by (\rho_{air} * cp).
      for (h = 0; h < Nt_in; h++)
        for (j = 0; j < Ny_in; j++)
          for (i = 0; i < Nx_in; i++)
            SensibleHeat(h, j, i) *= r * SkinTemperature(h, j, i)
              / (SurfacePressure(h, j, i) * cp);
      // -1. to change the direction of the flux.
      SensibleHeat.Mlt(-1. / (Delta_t_in * 3600.));
    }

  cout << " done." << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  // Computes level heights with pressure levels.
  cout << "Computing level heights in meter...";
  cout.flush();

  ExtStream file_hybrid_coefficients(directory_in +
                                     "hybrid_coefficients.dat");
  int valid_lines_count(0);
  while (file_hybrid_coefficients.GetLine().length() > 0)
    valid_lines_count++;
  file_hybrid_coefficients.Close();

  Data<real, 1> alpha(valid_lines_count), beta(valid_lines_count);

  FormatFormattedText InputCoefficients("<e><e><e><e 2>");
  InputCoefficients.Read(directory_in + "hybrid_coefficients.dat",
                         "1", alpha);
  InputCoefficients.Read(directory_in + "hybrid_coefficients.dat", "2", beta);

  alpha.ReverseData(0);
  alpha.Mlt(1. / P0);
  beta.ReverseData(0);

  ComputePressure(alpha, beta, SurfacePressure, Pressure_interf);

  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          Pressure(h, k, j, i) = (Pressure_interf(h, k + 1, j, i)
                                  + Pressure_interf(h, k, j, i)) / 2.0;

  ComputeInterfHeight(Pressure_interf, Temperature, GridZ_interf_in);
  ComputeMiddleHeight(Pressure_interf, Temperature, GridZ_interf_in,
                      GridZ_in);

  cout << " done." << endl;

  /////////////
  // CLOUDS  //
  /////////////

  cout << "Computing relative humidity and critical relative humidity...";
  cout.flush();

  ComputeRelativeHumidity(SpecificHumidity, Temperature, Pressure,
                          RelativeHumidity);
  RelativeHumidity.ThresholdMax(1.);
  if (option_critical_relative_humidity == 1)
    ComputeCriticalRelativeHumidity(SurfacePressure, Pressure, CRH);
  else if (option_critical_relative_humidity == 2)
    ComputeCriticalRelativeHumidity(Pressure, CRH);

  cout << " done." << endl;

  cout << "Computing cloud profile...";
  cout.flush();

  ComputeCloudFraction(BoundaryHeight, RelativeHumidity, CRH,
                       CloudFraction);

  Data<int, 4> LowIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> MediumIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> HighIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> MediumIndices_out(Nt_out, Ny_out, Nx_out, 2);

  ComputeCloudiness(CloudFraction, Pressure, GridZ_interf_in, LowIndices,
                    MediumIndices, HighIndices, LowCloudiness_diag,
                    MediumCloudiness_diag, HighCloudiness_diag);
  ComputeCloudBaseHeight(LowIndices, MediumIndices, HighIndices, GridZ_interf_in,
                         CloudHeight);

  // Reads liquid water content.
  Data<real, 4> LiquidWaterContent(GridT_in, GridZ_in, GridY_in, GridX_in);
  InputMeteo.Read(file_in, 246, LiquidWaterContent);
  LiquidWaterContent.ReverseData(1);
  LiquidWaterContent.ReverseData(2);
  LiquidWaterContent.ThresholdMin(0.);

  // data processing usefull for meteo only
  if (with_meteo)
    {

      ///////////////////////
      // RICHARDSON NUMBER //
      ///////////////////////

      cout << "Computing Richardson number...";
      cout.flush();

      // 3D Richardson number.
      ComputePotentialTemperature(Temperature, Pressure,
                                  PotentialTemperature);
      ComputeRichardson(ZonalWind, MeridionalWind,
                        PotentialTemperature, Richardson);

      // Surface Richardson number.
      ComputePotentialTemperature(SkinTemperature, SurfacePressure,
                                  SurfacePotentialTemperature);
      ComputeModule(ZonalWind, MeridionalWind, FirstLevelWindModule);
      if (richardson_with_roughness)
        // No more than 5m, knowing that the first ECMWF level is at ~9m.
        Roughness.ThresholdMax(5.0);
      if (richardson_with_roughness)
        ComputeRichardson(Roughness, FirstLevelWindModule,
                          SurfacePotentialTemperature,
                          PotentialTemperature, SurfaceRichardson);
      else
        ComputeRichardson(FirstLevelWindModule, SurfacePotentialTemperature,
                          PotentialTemperature, SurfaceRichardson);

      cout << " done." << endl;
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
        InputMeteo.Read(file_in, 247, IceWaterContent);
      else
        IceWaterContent.Fill(0.);
      IceWaterContent.ReverseData(1);
      IceWaterContent.ReverseData(2);
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

      cout << " done." << endl;
      cout << endl;
    }

  ///////////////////////////
  // LINEAR INTERPOLATIONS //
  ///////////////////////////

  cout << "Linear interpolations...";
  cout.flush();

  LinearInterpolationOneGeneral(CloudFraction, CloudFraction_out, 1);
  LinearInterpolationRegular(SurfacePressure, SurfacePressure_out);

  if (with_meteo)
    {
      LinearInterpolationOneGeneral(Pressure, Pressure_out, 1);
      LinearInterpolationOneGeneral(Temperature, Temperature_out, 1);
      LinearInterpolationOneGeneral(SpecificHumidity,
                                    SpecificHumidity_out, 1);
      LinearInterpolationOneGeneral(LiquidWaterContent,
                                    LiquidWaterContent_out, 1);

      LinearInterpolationOneGeneral(MeridionalWind, MeridionalWind_out, 1);
      LinearInterpolationOneGeneral(ZonalWind, ZonalWind_out, 1);
      LinearInterpolationRegular(SolarRadiation, SolarRadiation_out);
      LinearInterpolationRegular(SurfaceRichardson, SurfaceRichardson_out);
      LinearInterpolationOneGeneral(Richardson, Richardson_out, 1);
      LinearInterpolationRegular(BoundaryHeight, BoundaryHeight_out);
      BoundaryHeight_out.ThresholdMin(GridZ_interf_out(1));
      LinearInterpolationRegular(U_star, U_star_out);
      LinearInterpolationRegular(V_star, V_star_out);
      LinearInterpolationRegular(SurfaceTemperature, SurfaceTemperature_out);
      LinearInterpolationRegular(SkinTemperature, SkinTemperature_out);
      LinearInterpolationRegular(SoilWater, SoilWater_out);
      LinearInterpolationRegular(Evaporation, Evaporation_out);
      LinearInterpolationRegular(SensibleHeat, SensibleHeat_out);
      LinearInterpolationRegular(FirstLevelWindModule,
                                 FirstLevelWindModule_out);
      LinearInterpolationRegular(CloudHeight, CloudHeight_out);
      LinearInterpolationRegular(ConvectiveRain, ConvectiveRain_out);
      LinearInterpolationRegular(LargeScaleRain, LargeScaleRain_out);
      LinearInterpolationRegular(MediumCloudiness, MediumCloudiness_out);
      LinearInterpolationRegular(HighCloudiness, HighCloudiness_out);
      LinearInterpolationRegular(MediumIndices, MediumIndices_out);
      SpecificHumidity_out.Threshold(0., 1.);
      LiquidWaterContent_out.ThresholdMin(0.);
      IceWaterContent_out.ThresholdMin(0.);
      SoilWater_out.Threshold(0., 1.);
      SolarRadiation_out.ThresholdMin(0.);
      CloudHeight_out.ThresholdMin(min_height);
      ConvectiveRain_out.ThresholdMin(0.);
      LargeScaleRain_out.ThresholdMin(0.);
    }

  if (with_photolysis)
    {
      if (photolysis_tabulation == 1)
        LinearInterpolationOneGeneral(Attenuation, Attenuation_out, 1);
      if (photolysis_tabulation == 2 || photolysis_tabulation == 3)
        {
          LinearInterpolationOneGeneral(LiquidWaterExtinction,
                                        LiquidWaterExtinction_out, 1);
          LinearInterpolationOneGeneral(IceWaterExtinction,
                                        IceWaterExtinction_out, 1);
        }
      if (ice_cloud && (photolysis_tabulation == 2 || photolysis_tabulation == 3))
        LinearInterpolationOneGeneral(IceWaterContent,
                                      IceWaterContent_out, 1);
    }

  cout << " done." << endl;

  ///////////////////////////
  // OUTPUT DATA PROCESSING //
  ///////////////////////////


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
                      IceOpticalDepth_out(h, k, j, i) = iext * dz * pow(cf, 1.5);
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
                                               GridY_out, GridX_out);
          Data<real, 4> MeanExtinctionEfficiencyFactor(GridWavelength,
                                                       GridZ_out, GridY_out,
                                                       GridX_out);
          Data<real, 5> PhaseFunction(GridWavelength, GridZ_out, GridY_out,
                                      GridX_out, GridLegendre);
          Data<real, 4> PhotolysisRate_TimeStep(GridP, GridZ_out,
                                                GridY_out, GridX_out);
          Data<real, 3> OpticalDepth_TimeStep(GridZ_out, GridY_out,
                                              GridX_out);
          Data<real, 3> IceOpticalDepth_TimeStep(GridZ_out, GridY_out,
                                                 GridX_out);
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

      cout.flush();

      ComputeModule(ZonalWind_out, MeridionalWind_out, WindModule_out);

      U_star_out.Apply(::abs_);
      V_star_out.Apply(::abs_);

      U_star_out.Apply(sqrt_);
      V_star_out.Apply(sqrt_);

      ComputeModule(U_star_out, V_star_out, FrictionModule_out);

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
                real ut = real(current_date.GetHour()) +
                  current_date.GetMinutes() / 60.
                  + current_date.GetSeconds() / 3600.;
                zenith_angle = ZenithAngle(GridX_out(i), GridY_out(j),
                                           current_date, ut) / 180. * pi;

                if (zenith_angle >= 1.51844 ||
                    SolarRadiation_out(h, j, i) <= 0.)
                  {
                    PAR_out(h, j, i) = 0.;
                    PARdb_out(h, j, i) = 0.;
                    PARdiff_out(h, j, i) = 0.;
                  }
                else
                  {
                    optical_thickness = Pressure_out(h, 0, j, i)
                      / 101325. / cos(zenith_angle);
                    visible_beam = 600. * exp(-.185 * optical_thickness)
                      * cos(zenith_angle);
                    visible_diff = 0.42 * (600. - visible_beam)
                      * cos(zenith_angle);
                    water_abs = 1320. * .077 *
                      pow(2. * optical_thickness, 0.3);
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

                    PARdb  = SolarRadiation_out(h, j, i) * ratio
                      * fvb * watt2umol;
                    PARdiff = SolarRadiation_out(h, j, i) * ratio
                      * fvd * watt2umol;

                    PARdb_out(h, j, i) = PARdb;
                    PARdiff_out(h, j, i) = PARdiff;
                    PAR_out(h, j, i) = PARdb + PARdiff;
                  }
              }
          current_date.AddSeconds(int(Delta_t_out * 3600));
        }

      cout << " done." << endl;
    }

  cout << " done." << endl;


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  FormatBinary<float> OutputMeteo;

  cout << "Writing data...";
  cout.flush();

  // meteo data
  if (with_meteo)
    {
      OutputMeteo.Append(Pressure_out, directory_out + "Pressure.bin");
      OutputMeteo.Append(SurfacePressure_out, directory_out +
                         "SurfacePressure.bin");
      OutputMeteo.Append(Temperature_out, directory_out + "Temperature.bin");
      OutputMeteo.Append(SurfaceTemperature_out,
                         directory_out + "SurfaceTemperature.bin");
      OutputMeteo.Append(SkinTemperature_out, directory_out +
                         "SkinTemperature.bin");
      OutputMeteo.Append(Richardson_out, directory_out + "Richardson.bin");
      OutputMeteo.Append(SurfaceRichardson_out,
                         directory_out + "SurfaceRichardson.bin");
      OutputMeteo.Append(SpecificHumidity_out, directory_out +
                         "SpecificHumidity.bin");
      OutputMeteo.Append(LiquidWaterContent_out, directory_out +
                         "LiquidWaterContent.bin");
      OutputMeteo.Append(SolarRadiation_out, directory_out +
                         "SolarRadiation.bin");
      OutputMeteo.Append(PARdb_out, directory_out + "PARdb.bin");
      OutputMeteo.Append(PARdiff_out, directory_out + "PARdiff.bin");
      OutputMeteo.Append(PAR_out, directory_out + "PAR.bin");
      OutputMeteo.Append(ZonalWind_out, directory_out + "ZonalWind.bin");
      OutputMeteo.Append(MeridionalWind_out, directory_out +
                         "MeridionalWind.bin");
      OutputMeteo.Append(WindModule_out, directory_out + "WindModule.bin");
      OutputMeteo.Append(FrictionModule_out, directory_out +
                         "FrictionModule.bin");
      OutputMeteo.Append(BoundaryHeight_out, directory_out +
                         "BoundaryHeight.bin");
      OutputMeteo.Append(SoilWater_out, directory_out + "SoilWater.bin");
      OutputMeteo.Append(Evaporation_out, directory_out + "Evaporation.bin");
      OutputMeteo.Append(SensibleHeat_out,
                         directory_out + "SensibleHeat.bin");
      OutputMeteo.Append(FirstLevelWindModule_out, directory_out
                         + "FirstLevelWindModule.bin");
      OutputMeteo.Append(CloudHeight_out, directory_out + "CloudHeight.bin");
      OutputMeteo.Append(ConvectiveRain_out,
                         directory_out + "ConvectiveRain.bin");
      LargeScaleRain_out.GetArray() = LargeScaleRain_out.GetArray()
        + ConvectiveRain_out.GetArray();
      OutputMeteo.Append(LargeScaleRain_out, directory_out + "Rain.bin");
    }

  // photolysis data
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
              PhotolysisRate4D.SubData(PhotolysisRate_out, i,
                                       Range::all(), Range::all(), Range::all(),
                                       Range::all());
              OutputMeteo.Append(PhotolysisRate4D,
                                 directory_photolysis
                                 + photolysis_reaction_list[i] + ".bin");
            }
        }
      if (photolysis_tabulation == 3)
        {
          OutputMeteo.Append(OpticalDepth_out,
                             directory_out + "CloudOpticalDepth.bin");
          OutputMeteo.Append(IceOpticalDepth_out,
                             directory_out + "IceOpticalDepth.bin");
        }
      if (ice_cloud && (photolysis_tabulation == 2 || photolysis_tabulation == 3))
        OutputMeteo.Append(IceWaterContent_out, directory_out +
                           "IceWaterContent.bin");
    }

  cout << " done." << endl;

  cout << endl << endl;
  cout << "||||||||||||||||||||||||||||||||||||||||||" << endl;
  cout << "|     Successful completion of meteo     |" << endl;
  cout << "||||||||||||||||||||||||||||||||||||||||||" << endl << endl;

  END;

  return 0;
}

