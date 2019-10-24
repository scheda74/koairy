// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet
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

#include "Diagnosis.cxx"

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("water-plume.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef float real;

  int h, k, j, i;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////


  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (sec_config_file != "")
    configuration.AddFile(sec_config_file);

  /*** Input domain ***/

  configuration.SetSection("[domain]");

  // Domain and meteorological time discretization.
  string date_meteo_string;
  Date date_meteo;
  int Nt, Nz, Ny, Nx;
  real Delta_t, Delta_y, Delta_x;
  real t_min, y_min, x_min;

  configuration.PeekValue("Date", date_meteo_string);
  date_meteo = date_meteo_string;
  configuration.PeekValue("Nx", Nx);
  configuration.PeekValue("Ny", Ny);
  configuration.PeekValue("Nz", Nz);
  configuration.PeekValue("Delta_t", Delta_t);
  configuration.PeekValue("Delta_y", Delta_y);
  configuration.PeekValue("Delta_x", Delta_x);
  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  /*** Simulation ***/

  configuration.SetSection("[simulation]");

  // Simulation time discretization.
  real Delta_t_sim;
  string date_sim_string;

  configuration.PeekValue("Date", date_sim_string);
  configuration.PeekValue("Delta_t", Delta_t_sim);

  Date date_sim = date_sim_string;
  int Nt_sim = compute_Nt(date_beg, date_end, Delta_t_sim);

  // Plume total water content.
  string plume_water_file;
  configuration.PeekValue("PlumeWater", plume_water_file);
  real factor;
  configuration.PeekValue("Factor", factor);


  /*** Meteorological files ***/

  configuration.SetSection("[meteo]");

  string temperature_file, pressure_file, specific_humidity_file,
    liquid_water_content_file;

  // Meteorological files.
  configuration.PeekValue("Temperature", temperature_file);
  configuration.PeekValue("Pressure", pressure_file);
  configuration.PeekValue("SpecificHumidity", specific_humidity_file);
  configuration.PeekValue("LiquidWaterContent", liquid_water_content_file);


  /*** Parameters ***/

  configuration.SetSection("[parameters]");

  // Initial conditions (at source point).
  // Liquid water potential temperature in K at the source.
  real source_temperature;
  // Total water content at the source in kg/kg.
  real source_water_content;

  configuration.PeekValue("source_temperature", source_temperature);
  configuration.PeekValue("source_water_content", source_water_content);

  /*** Output ***/

  configuration.SetSection("[output]");

  string plume_liquid_file, unit, option;

  configuration.PeekValue("Option", option);
  configuration.PeekValue("Unit", unit);
  configuration.PeekValue("LiquidWaterContent", plume_liquid_file);

  if (option != "plume" && option != "total")
    throw string("Wrong Option in section [output].\nYou should put ")
      + "\"plume\" for the liquid water content in the plume or \"total\" "
      + "for the liquid water content in the plume and in ambient air.";

  unit = lower_case(unit);
  if (unit != "a" && unit != "b")
    throw string("Wrong Unit in section [output].\nYou should put a ")
      +  "(for g/kg) or b (for g/m^3).";

  cout << " done." << endl;


  //////////
  // DATA //
  //////////


  cout << "Memory allocation for data fields...";
  cout.flush();

  /*** Input fields ***/

  // Grids.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridT_sim(t_min, Delta_t_sim, Nt_sim);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);

  // Meteorological fields.
  Data<real, 4> Temperature_in(GridT, GridZ, GridY, GridX);
  Data<real, 4> Pressure_in(GridT, GridZ, GridY, GridX);
  Data<real, 4> SpecificHumidity_in(GridT, GridZ, GridY, GridX);
  Data<real, 4> LiquidWaterContent_in(GridT, GridZ, GridY, GridX);
  // On the simulated grid.
  Data<real, 4> Temperature(GridT_sim, GridZ, GridY, GridX);
  Data<real, 4> Pressure(GridT_sim, GridZ, GridY, GridX);
  Data<real, 4> SpecificHumidity(GridT_sim, GridZ, GridY, GridX);
  Data<real, 4> LiquidWaterContent(GridT_sim, GridZ, GridY, GridX);

  // Total water content in the plume.
  Data<real, 4> PlumeWaterContent(GridT_sim, GridZ, GridY, GridX);

  // Temporary fields.
  // Liquid water potential temperature in K.
  Data<real, 4> PotentialTemperatureLW(GridT_sim, GridZ, GridY, GridX);
  // Air density in kg/m^3.
  Data<real, 4> AirDensity(GridT_sim, GridZ, GridY, GridX);
  // Liquid water potential temperature in K in the plume.
  Data<real, 4> PlumePotentialTemperatureLW(GridT_sim, GridZ, GridY, GridX);

  /*** Output fields ***/

  // Plume liquid water content (diagnosed) in the unit chosen:
  // g/kg (option a) or g/m^3 (option b).
  Data<real, 4> PlumeLiquidWaterContent(GridT_sim, GridZ, GridY, GridX);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////


  cout << "Extracting data...";
  cout.flush();

  FormatBinary<float> IO;

  // Meteorological fields.
  int step_meteo = int(date_beg.GetSecondsFrom(date_meteo) / 3600 / Delta_t
                       + 0.5);
  IO.ReadSteps(temperature_file, step_meteo, Temperature_in);
  IO.ReadSteps(pressure_file, step_meteo, Pressure_in);
  IO.ReadSteps(specific_humidity_file, step_meteo, SpecificHumidity_in);
  IO.ReadSteps(liquid_water_content_file, step_meteo, LiquidWaterContent_in);

  // Plume total water content.
  int step_sim = int(date_beg.GetSecondsFrom(date_sim) / 3600 / Delta_t_sim
                     + 0.5);
  IO.ReadSteps(plume_water_file, step_sim, PlumeWaterContent);

  cout << " done." << endl;


  ///////////////////
  // INTERPOLATION //
  //////////////////


  cout << "Performing interpolations...";
  cout.flush();

  LinearInterpolationDimension(Temperature_in, Temperature, 0);
  LinearInterpolationDimension(Pressure_in, Pressure, 0);
  LinearInterpolationDimension(SpecificHumidity_in, SpecificHumidity, 0);
  LinearInterpolationDimension(LiquidWaterContent_in, LiquidWaterContent, 0);

  cout << " done." << endl;


  ///////////////
  // DIAGNOSES //
  ///////////////


  cout << "Performing diagnosis...";
  cout.flush();

  ComputeAirDensity(Temperature, Pressure, AirDensity);

  // To g.m^{-3}.
  PlumeWaterContent.Mlt(factor);

  // From concentrations in g.m^{-3}
  // to (kg of total water)/(kg of dry air with water vapor).
  PlumeWaterContent.GetArray() =
    1.e-3 * PlumeWaterContent.GetArray() / AirDensity.GetArray();

  ComputePotentialTemperatureLW(Temperature, Pressure, LiquidWaterContent,
                                PotentialTemperatureLW);

  // Similarity assumption between water content and plume temperature.
  for (h = 0; h < Nt_sim; h++)
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          PlumePotentialTemperatureLW(h, k, j, i)
            = PlumeWaterContent(h, k, j, i) / source_water_content *
            (source_temperature - PotentialTemperatureLW(h, k, j, i))
            + PotentialTemperatureLW(h, k, j, i);

  // Computes PlumeLiquidWaterContent.
  WaterDiagnosis(Pressure, Temperature, LiquidWaterContent,
                 PlumePotentialTemperatureLW, PlumeWaterContent,
                 PlumeLiquidWaterContent, option);

  cout << " done." << endl;

  double mean;
  // Conversion of units.
  if (unit == "a")
    {
      PlumeLiquidWaterContent.Mlt(1.e3);
      cout << "Output in g/kg ";
    }
  else
    {
      PlumeLiquidWaterContent.GetArray() = PlumeLiquidWaterContent.GetArray()
        * 1000. * AirDensity.GetArray();
      cout << "Output in g/m^3 ";
    }
  if (option == "plume")
    cout << "for the plume alone: " << endl;
  else
    cout << "for the plume and the ambient air: " << endl;

  PlumeLiquidWaterContent.Mean(mean);

  cout << "   Max: " << PlumeLiquidWaterContent.GetMax() << endl;
  cout << "   Min: " << PlumeLiquidWaterContent.GetMin() << endl;
  cout << "   Mean: " << mean << endl;

  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////


  cout << "Writing data...";
  cout.flush();

  IO.Append(PlumeLiquidWaterContent, plume_liquid_file);

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
