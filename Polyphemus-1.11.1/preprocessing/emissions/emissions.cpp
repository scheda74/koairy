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


// This program reads EMEP total annual emissions and put it on Polair3D Grid.

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <map>

using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

int main(int argc, char* argv[])
{

  TRY;

  cout << endl;

  typedef float real;

  string configuration_file, sec_config_file, default_name("emissions.cfg");
  Date date;

  parse_argument(argc, argv, configuration_file, sec_config_file, date,
                 default_name);

  cout << "Date: " << date.GetDate() << endl << endl;


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  int h, i, j, k, l, m, s;

  const real avogadro = 6.02213e23;

  // Configuration file.
  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  /*** Output ***/

  int Nt = 24;
  int Nx, Ny, Nz;
  real Delta_x, Delta_y;
  real x_min, y_min;
  string levels;
  bool divide_by_heights;
  string unit, molecular_weights_file;
  map<string, real> molecular_weights;

  /*** Options ***/

  configuration.SetSection("[options]");

  configuration.PeekValue("Divide_by_heights", divide_by_heights);
  configuration.PeekValue("Output_unit", "mass | number", unit);
  lower_case(unit);

  if (unit == "number")
    {
      // Molecular weights.
      configuration.PeekValue("Molecular_weights", molecular_weights_file);
      ExtStream molecular_weights_stream(molecular_weights_file);
      string species;
      while (!molecular_weights_stream.IsEmpty())
        {
          species = molecular_weights_stream.GetElement();
          molecular_weights_stream.GetNumber(molecular_weights[species]);
        }
    }

  /*** Domain ***/

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("x_min", x_min);
  configuration.PeekValue("y_min", y_min);

  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridT(Nt);
  RegularGrid<real> GridWeek(7);
  RegularGrid<real> GridMonth(12);

  RegularGrid<real> GridZ_interf(Nz + 1);
  if (divide_by_heights)
    {
      configuration.PeekValue("Vertical_levels", levels);
      FormatText Levels;
      Levels.Read(levels, GridZ_interf);
    }

  /*** EMEP ***/

  int Nx_emep, Ny_emep, Ncountries, Nsp_emis, Nsectors(10);
  string EMEP_directory, hourly_factors_file,
    weekdays_factors_file, monthly_factors_file,
    time_zones_file;
  real Ratio_urb, Ratio_for, Ratio_oth;

  configuration.SetSection("[EMEP]");

  configuration.PeekValue("Input_directory", EMEP_directory);

  configuration.PeekValue("Nx_emep", "> 0", Nx_emep);
  configuration.PeekValue("Ny_emep", "> 0", Ny_emep);
  configuration.PeekValue("Ncountries", "positive", Ncountries);

  configuration.PeekValue("Hourly_factors", hourly_factors_file);
  configuration.PeekValue("Weekdays_factors", weekdays_factors_file);
  configuration.PeekValue("Monthly_factors", monthly_factors_file);
  configuration.PeekValue("Time_zones", time_zones_file);

  configuration.PeekValue("Urban_ratio", "positive", Ratio_urb);
  configuration.PeekValue("Forest_ratio", "positive", Ratio_for);
  configuration.PeekValue("Other_ratio", "positive", Ratio_oth);

  string vertical_distribution_file;
  configuration.PeekValue("Polair_vertical_distribution",
                          vertical_distribution_file);

  // Species names.
  vector<string> Sp_emis_names;
  vector<string>::iterator iter_emis;

  configuration.Find("Species");
  split(configuration.GetLine(), Sp_emis_names);
  Nsp_emis = Sp_emis_names.size();

  for (iter_emis = Sp_emis_names.begin(); iter_emis != Sp_emis_names.end();
       iter_emis++)
    if (*iter_emis == "NMVOC" && iter_emis != Sp_emis_names.end() - 1)
      {
        Sp_emis_names.erase(iter_emis--);
        Sp_emis_names.push_back("NMVOC");
      }

  RegularGrid<real> GridY_emep(Ny_emep);
  RegularGrid<real> GridX_emep(Nx_emep);
  RegularGrid<real> Grid_countries(Ncountries);
  RegularGrid<real> GridSectors(Nsectors);
  RegularGrid<real> GridSp_emis(Nsp_emis);

  /*** LUC ***/

  int Nx_luc, Ny_luc;
  real x_min_luc, y_min_luc;
  real Delta_x_luc, Delta_y_luc;
  string LUC_file;
  string LUC_type;

  configuration.SetSection("[LUC]");

  configuration.PeekValue("Nx", "> 0", Nx_luc);
  configuration.PeekValue("Ny", "> 0", Ny_luc);
  configuration.PeekValue("x_min", x_min_luc);
  configuration.PeekValue("y_min", y_min_luc);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_luc);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_luc);

  configuration.PeekValue("File", LUC_file);
  configuration.PeekValue("LUC_type", LUC_type);

  RegularGrid<real> GridX_luc(x_min_luc, Delta_x_luc, Nx_luc);
  RegularGrid<real> GridY_luc(y_min_luc, Delta_y_luc, Ny_luc);

  if ( (x_min_luc > x_min) || (y_min_luc > y_min) || \
       (x_min_luc + Delta_x_luc * Nx_luc < x_min + Delta_x * Nx) || \
       (y_min_luc + Delta_y_luc * Ny_luc < y_min + Delta_y * Ny))
    throw string("[ERROR] \"") + LUC_file + "\" does not cover the domain "
      "as described by the '[domain]' section.";

  /*** Species ***/

  int Nsp_model;
  string aggregation_file, speciation_directory;
  real deposition_factor_nh3;

  configuration.SetSection("[Species]");

  configuration.PeekValue("N", "positive", Nsp_model);
  configuration.PeekValue("Aggregation", aggregation_file);
  configuration.PeekValue("Speciation_directory", speciation_directory);
  configuration.PeekValue("Deposition_factor_NH3", "positive",
                          deposition_factor_nh3);

  RegularGrid<real> GridSp_model(Nsp_model);

  /*** Files ***/

  string Dir_out_emis, Dir_out_src;

  configuration.SetSection("[paths]");

  configuration.PeekValue("Directory_surface_emissions", Dir_out_emis);
  configuration.PeekValue("Directory_volume_emissions", Dir_out_src);


  //////////
  // DATA //
  //////////

  Data<list<EmepCountryEmission<real> >, 4, real> Emis_land(GridSp_emis,
                                                            GridY_emep,
                                                            GridX_emep,
                                                            GridSectors);
  Data<list<EmepCountryEmission<real> >, 4, real> Emis_water(GridSp_emis,
                                                             GridY_emep,
                                                             GridX_emep,
                                                             GridSectors);

  Data<int, 2, real> Ntot_polair(GridY, GridX);

  Data<int, 2, real> Nurb_emep(GridY_emep, GridX_emep);
  Data<int, 2, real> Nwat_emep(GridY_emep, GridX_emep);
  Data<int, 2, real> Nfor_emep(GridY_emep, GridX_emep);
  Data<int, 2, real> Noth_emep(GridY_emep, GridX_emep);

  Data<int, 2, real> LUC(GridY_luc, GridX_luc);
  Data<list<EmepCountryEmission<real> >, 4, real> Emis_out(GridSp_emis,
                                                           GridSectors,
                                                           GridY, GridX);

  Data<real, 4> Emissions_out(GridSp_model, GridT, GridY, GridX);
  Data<real, 3> Species_factor(GridSp_model, GridSp_emis, GridSectors);
  Data<real, 3> MonthlyFactors(Grid_countries, GridSectors, GridMonth);
  Data<real, 3> DailyFactors(Grid_countries, GridSectors, GridWeek);
  Data<real, 2> HourlyFactors(GridSectors, GridT);

  Data<real, 1> Ground_part(GridSectors);
  Data<real, 2> vertical_distribution_out(GridSectors, GridZ);
  Data<real, 5> Sources_out(GridSp_model, GridT, GridZ, GridY, GridX);

  vector<string> Sp_model_names(Nsp_model);

  // LUC.
  FormatBinary<int> FormatLUC;
  FormatLUC.Read(LUC_file, LUC);


  //////////////////////
  // READS INPUT DATA //
  //////////////////////

  cout << "Reading monthly factors...";
  cout.flush();
  GetTemporalFactors(monthly_factors_file, MonthlyFactors);
  cout << " done." << endl;

  cout << "Reading daily factors...";
  cout.flush();
  GetTemporalFactors(weekdays_factors_file, DailyFactors);
  cout << " done." << endl;

  cout << "Grid correspondences EMEP / LUC / Polair3D...";
  cout.flush();

  GridCorrespondences(LUC, LUC_type, Nurb_emep, Nwat_emep,
                      Nfor_emep, Noth_emep, Ntot_polair);
  cout << " done." << endl;

  cout << "Reading aggregation/speciation...";
  cout.flush();
  SpeciationAggregation(Sp_emis_names, aggregation_file, speciation_directory,
                        Sp_model_names, Species_factor);

  while (Sp_model_names[Nsp_model - 1].empty() && Nsp_model >= 1)
    Nsp_model--;

  cout << " done." << endl;

  cout << "Reading vertical distribution...";
  cout.flush();
  FormatFormattedText distribution_format(string("<e><e ") + to_str(Nz) + ">");
  distribution_format.Read(vertical_distribution_file, "0", Ground_part);
  distribution_format.Read(vertical_distribution_file, "1",
                           vertical_distribution_out);

  // Percentages mapped to [0, 1].
  Ground_part.Mlt(0.01);
  vertical_distribution_out.Mlt(0.01);

  if (divide_by_heights)
    DividesByHeights(GridZ_interf, vertical_distribution_out);

  // maximum height.
  int k_max = Nz;
  bool zero(true);
  while (zero && k_max > 0)
    {
      k_max--;
      for (s = 0; s < Nsectors; s++)
        zero = vertical_distribution_out(s, k_max) == 0. && zero;
    }
  if (!zero)
    k_max++;
  cout << " done." << endl;

  /*** Reads hourly profiles ***/

  cout << "Reading hourly factors...";
  cout.flush();
  FormatFormattedText HourlyFactors_format(string("<e><e ")
                                           + to_str(HourlyFactors.GetLength(1))
                                           + ">");
  HourlyFactors_format.Read(hourly_factors_file, "1", HourlyFactors);
  // Local times associated with countries.
  TimeZone local_time(Ncountries);
  local_time.Init(time_zones_file);
  cout << " done." << endl;

  /*** Reads input Emissions ***/

  cout << "Reading input emissions...";
  cout.flush();
  ReadEmep(date, Sp_emis_names, EMEP_directory, time_zones_file,
           MonthlyFactors, DailyFactors, deposition_factor_nh3, Emis_land, Emis_water);
  cout << " done." << endl;


  ////////////////////////
  // COMPUTES EMISSIONS //
  ////////////////////////

  cout << "Computes emissions on Polair3D grid...";
  cout.flush();
  EmepToLatLon(LUC, LUC_type, Ratio_urb, Ratio_for, Ratio_oth, Nurb_emep, Nwat_emep,
               Nfor_emep, Noth_emep, Ntot_polair, Emis_land, Emis_water,
               Emis_out);
  cout << " done." << endl;

  cout << "Computes final emissions...";
  cout.flush();

  // Winter/summer time.
  int inc, time_index;
  // Searches for the last sunday in March.
  Date winter_limit;
  winter_limit.SetYear(date.GetYear());
  winter_limit.SetMonth(3);
  winter_limit.SetDay(31);
  while (winter_limit.GetWeekDay() != 6)
    winter_limit.AddDays(-1);
  // Searches for the last sunday in October.
  Date summer_limit;
  summer_limit.SetYear(date.GetYear());
  summer_limit.SetMonth(10);
  summer_limit.SetDay(31);
  while (summer_limit.GetWeekDay() != 6)
    summer_limit.AddDays(-1);

  if (date.GetNumberOfDays() < winter_limit.GetNumberOfDays()
      || date.GetNumberOfDays() > summer_limit.GetNumberOfDays())
    inc = 0;
  else
    inc = 1;

  real factor, emission;

  list<EmepCountryEmission<real> >::iterator iter;

  // Surface emissions.
  Emissions_out.Fill(0.0);
  for (l = 0; l < Nsp_emis; l++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        for (s = 0; s < Nsectors; s++)
          for (iter = Emis_out(l, s, j, i).begin();
               iter != Emis_out(l, s, j, i).end(); ++iter)
            {
              emission = iter->emission_;
              for (m = 0; m < Nsp_model; m++)
                {
                  factor = Species_factor(m, l, s) * emission;
                  time_index = inc + local_time(iter->country_);
                  for (h = 0; h < Nt; h++)
                    Emissions_out(m, h, j, i) += factor
                      * Ground_part(s)
                      * HourlyFactors(s, time_index++ % 24);
                }
            }
  // From Tons/m^2/day to microg/m^2/s.
  Emissions_out.Mlt(1.e12 / (real(Nt) * 3600.));

  Sources_out.Fill(0.0);
  if (k_max != 0)
    for (l = 0; l < Nsp_emis; l++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          for (s = 0; s < Nsectors; s++)
            for (iter = Emis_out(l, s, j, i).begin();
                 iter != Emis_out(l, s, j, i).end(); ++iter)
              {
                emission = iter->emission_;
                for (m = 0; m < Nsp_model; m++)
                  {
                    factor = Species_factor(m, l, s) * emission;
                    for (k = 0; k < k_max; k++)
                      {
                        time_index = inc + local_time(iter->country_);
                        for (h = 0; h < Nt; h++)
                          Sources_out(m, h, k, j, i) += factor
                            * vertical_distribution_out(s, k)
                            * HourlyFactors(s, time_index++ % 24);
                      }
                  }
              }
  // From Tons/m^2/day/height to microg/m^3/s.
  Sources_out.Mlt(1.e12 / (3600. * real(Nt)));

  cout << " done." << endl;


  /////////////////
  // WRITES DATA //
  /////////////////

  cout << "Writing output emissions...";
  cout.flush();

  FormatBinary<float> Output;

  bool zero_surface_emissions(true);
  for (s = 0; s < Nsectors; s++)
    zero_surface_emissions = Ground_part(s) == 0. && zero_surface_emissions;

  if (!zero_surface_emissions)
    for (m = 0; m < Nsp_model; m++)
      {
        Data<real, 3> Emis_extract(Nt, Ny, Nx);

        Emis_extract.SubData(Emissions_out,
                             m, Range::all(), Range::all(), Range::all());
        // From microgramm/m^2 to molecule/cm^2.
        if (unit == "number")
          Emis_extract.Mlt(1e-10 * avogadro
                           / molecular_weights[Sp_model_names[m]]);
        Output.Append(Emis_extract, Dir_out_emis + Sp_model_names[m] + ".bin");
      }

  if (k_max != 0)
    for (m = 0; m < Nsp_model; m++)
      {
        Data<real, 4> Sources_extract(Nt, k_max, Ny, Nx);
        Sources_extract.SubData(Sources_out,
                                m, Range::all(), Range(0, k_max - 1),
                                Range::all(), Range::all());
        // From microgramm/m^2 to molecule/cm^2.
        if (unit == "number")
          Sources_extract.Mlt(1e-10 * avogadro
                              / molecular_weights[Sp_model_names[m]]);
        Output.Append(Sources_extract, Dir_out_src + Sp_model_names[m] + ".bin");
      }

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
