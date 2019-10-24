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


// This program computes biogenic emissions.


//////////////
// INCLUDES //

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("bio.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Input domain.
  int Nt, Nz, Ny, Nx;
  real Delta_t, Delta_y, Delta_x;
  real t_min, y_min, x_min;

  config.SetSection("[domain]");

  config.PeekValue("Nz", "> 0", Nz);
  config.PeekValue("Ny", "> 0", Ny);
  config.PeekValue("Nx", "> 0", Nx);
  config.PeekValue("Delta_t", "> 0", Delta_t);
  config.PeekValue("Delta_y", "> 0", Delta_y);
  config.PeekValue("Delta_x", "> 0", Delta_x);
  config.PeekValue("y_min", y_min);
  config.PeekValue("x_min", x_min);

  // Dates.
  string date_meteo_str;
  config.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);
  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw string("\nThe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  int days = date_beg.GetDaysFrom(date_meteo);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  // Files.
  string dir_out, surface_temperature_file,
    PAR_file, LUC_file, land_data_file;
  real Delta_t_bio;
  int Nc, Nt_bio;

  config.SetSection("[paths]");

  config.PeekValue("SurfaceTemperature", surface_temperature_file);
  config.PeekValue("PAR", PAR_file);
  // An extra step available for PAR?
  bool extra_step = (int(file_size(PAR_file) / (Nt * Nx * Ny * sizeof(float)))
                     > days + 1);

  config.PeekValue("LUC_file", LUC_file);
  // Number of land use categories.
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny * Nx);
  config.PeekValue("Land_data", land_data_file);
  config.PeekValue("Directory_bio", dir_out);

  config.SetSection("[biogenic]");

  config.PeekValue("Delta_t", "> 0", Delta_t_bio);
  Nt_bio = int(real(Nt) * Delta_t / Delta_t_bio);

  bool save_rates;
  config.PeekValue("Rates", save_rates);

  int Nterpenes, Nratios;
  vector<string> terpenes_names;
  vector<real> terpenes_ratios;

  config.Find("Terpenes");
  split(config.GetLine(), terpenes_names);
  Nterpenes = int(terpenes_names.size());

  config.Find("Terpenes_ratios");
  split(config.GetLine(), terpenes_ratios);
  Nratios = int(terpenes_ratios.size());

  if (Nterpenes != Nratios)
    {
      throw string("Pb in configuration file: Terpenes") +
        string(" and Terpenes_ratios have not the same number of elements.");
    }


  /////////////////////////
  // METEOROLOGICAL DATA //
  /////////////////////////

  cout << "Reading and interpolating meteorological data...";
  cout.flush();

  // Input fields.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridT_cumulated(t_min - Delta_t / 2.0, Delta_t,
                                    extra_step ? Nt + 1 : Nt);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridC(Nc);

  Data<real, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<real, 3> PAR(GridT_cumulated, GridY, GridX);
  Data<real, 3> LUC(GridC, GridY, GridX);

  FormatBinary<float> Meteo;
  Meteo.ReadSteps(surface_temperature_file, step, SurfaceTemperature);
  Meteo.ReadSteps(PAR_file, step, PAR);
  Meteo.Read(LUC_file, LUC);

  // Output fields.
  RegularGrid<real> GridT_out(t_min, Delta_t_bio, Nt_bio);

  Data<real, 3> SurfaceTemperature_out(GridT_out, GridY, GridX);
  Data<real, 3> PAR_out(GridT_out, GridY, GridX);

  // Temporal linear interpolation.
  LinearInterpolationDimension(SurfaceTemperature, SurfaceTemperature_out, 0);
  LinearInterpolationDimension(PAR, PAR_out, 0);
  PAR_out.ThresholdMin(0.);

  cout << "  done." << endl;
  cout << endl;


  ///////////////////
  // EMISSION DATA //
  ///////////////////

  cout << "Reading biogenic data...";
  cout.flush();

  ifstream LandStream(land_data_file.c_str());
  string line;
  Data<float, 1> Density(Nc), EF_isoprene(Nc), EF_terpenes(Nc), EF_NO(Nc);

  FormatFormattedText LandData("<c 0 54><e><e><e><e>");

  LandData.Read(land_data_file, "1", Density);
  LandData.Read(land_data_file, "2", EF_isoprene);
  LandData.Read(land_data_file, "3", EF_terpenes);
  LandData.Read(land_data_file, "4", EF_NO);

  EF_isoprene.Mlt(1. / 3600.);
  EF_terpenes.Mlt(1. / 3600.);
  EF_NO.Mlt(1.e-3);

  cout << " done." << endl;


  ///////////////
  // EMISSIONS //
  ///////////////

  cout << "Computing biogenic emissions...";
  cout.flush();

  Data<float, 3, real> Isoprene(GridT_out, GridY, GridX);
  Data<float, 3, real> Terpenes(GridT_out, GridY, GridX);
  Data<float, 3, real> NO(GridT_out, GridY, GridX);


  // Emissions without environmental correction factors.
  Data<float, 2, real> IsopreneRates(GridY, GridX);
  Data<float, 2, real> TerpenesRates(GridY, GridX);
  Data<float, 2, real> NORates(GridY, GridX);

  Isoprene.SetZero();
  Terpenes.SetZero();
  NO.SetZero();

  IsopreneRates.SetZero();
  TerpenesRates.SetZero();
  NORates.SetZero();

  ComputeBiogenicRates(LUC, Density, EF_isoprene, EF_terpenes, EF_NO,
                       IsopreneRates, TerpenesRates, NORates);
  ComputeBiogenicEmissions(SurfaceTemperature_out, PAR_out, LUC, Density,
                           EF_isoprene, EF_terpenes, EF_NO, Isoprene,
                           Terpenes, NO);

  cout << " done." << endl;


  cout << "Writing emissions...";
  cout.flush();
  FormatBinary<float> EmissionsFormat;
  EmissionsFormat.Append(Isoprene, dir_out + "Isoprene.bin");
  EmissionsFormat.Append(NO, dir_out + "NO.bin");

  // Terpenes speciation
  for (int nt = 0; nt < Nterpenes; nt++)
    if (terpenes_ratios[nt] > 0)
      {
        Data<float, 3, real> Terpenes_tmp(Terpenes);
        Terpenes_tmp.Mlt(terpenes_ratios[nt]);
        EmissionsFormat.Append(Terpenes_tmp,
                               dir_out + terpenes_names[nt] + ".bin");
      }
  if (save_rates)
    {
      EmissionsFormat.Write(IsopreneRates, dir_out + "IsopreneRates.bin");
      EmissionsFormat.Write(TerpenesRates, dir_out + "TerpenesRates.bin");
      EmissionsFormat.Write(NORates, dir_out + "NORates.bin");
    }



  cout << " done." << endl;


  //////////////
  // ANALYSIS //
  //////////////

  cout << "\n-- Analysis --\n" << endl;

  DISP(Isoprene.GetMax());
  DISP(Isoprene.Mean());
  cout << endl;

  DISP(NO.GetMax());
  DISP(NO.Mean());
  cout << endl;

  for (int nt = 0; nt < Nterpenes; nt++)
    if (terpenes_ratios[nt] > 0)
      {
        Data<float, 3, real> Terpenes_tmp(Terpenes);
        Terpenes_tmp.Mlt(terpenes_ratios[nt]);
        cout << "Max " <<  terpenes_names[nt]
             << " = " << Terpenes_tmp.GetMax() << endl;
        cout << "Mean " <<  terpenes_names[nt]
             << " = " << Terpenes_tmp.Mean() << endl;
        cout << endl;
      }

  END;

  return 0;

}
