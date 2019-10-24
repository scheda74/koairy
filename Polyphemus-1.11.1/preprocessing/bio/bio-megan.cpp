// Copyright (C) 2013, ENPC
//    Author(s): Florian Couvidat
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


#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <map>

using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_NETCDF

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;


typedef float T;


T compute_area(T y, T dx, T dy)
{
  const T pi = 3.14159265358979323846264;
  return 3600.0 * 1852.0 * 1852.0 * cos(y * pi / 180.0) * dx * dy;
}


int main(int argc, char** argv)
{

  TRY;
  string configuration_file, sec_config_file, default_name("MEGAN.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file, date_beg,
                 date_end, default_name);

  ConfigStreams config(configuration_file);
  if (exists(sec_config_file))
    config.AddFile(sec_config_file);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  // Input domain.
  int Nt, Nz, Ny, Nx;
  T Delta_t, Delta_y, Delta_x;
  T t_min, y_min, x_min;

  config.SetSection("[domain]");
  config.PeekValue("Nz", "> 0", Nz);
  config.PeekValue("Ny", "> 0", Ny);
  config.PeekValue("Nx", "> 0", Nx);
  config.PeekValue("Delta_t", "> 0", Delta_t);
  config.PeekValue("Delta_y", "> 0", Delta_y);
  config.PeekValue("Delta_x", "> 0", Delta_x);
  config.PeekValue("y_min", y_min);
  config.PeekValue("x_min", x_min);

  //! Mapping of LAI and EF data on the simulation domain.
  int nred = 5.0 / min(Delta_y, Delta_x);
  
  // Dates.
  string date_meteo_str;
  config.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw "The date provided in command line should be after the date "
      "in the main configuration file.";
  int step = int(difference / 3600 / Delta_t + 0.5);

  int days = date_beg.GetDaysFrom(date_meteo);

  Nt = compute_Nt(date_beg, date_end, Delta_t);

  int Ndays = int(Nt * Delta_t / 24.);

  if (Ndays * 24 != Nt * Delta_t)
    throw "[ERROR] The computation has to be done on a whole number "
      "of days.";
  t_min = T(date_beg.GetHour()) + T(date_beg.GetMinutes()) / 60.
    + T(date_beg.GetSeconds()) / 3600.;

  // Files.
  string dir_out, surface_temperature_file, SoilMoisture_file,
    PAR_file, LUC_file, land_data_file, dir_EF, dir_LAI,
    Aggregation_file;
  T Delta_t_bio, Delta_x_EF, Delta_y_EF, xmin_EF, ymin_EF;
  int Nc, Nt_bio, Nx_in, Ny_in, Nx_EFin, Ny_EFin;

  config.SetSection("[paths]");

  config.PeekValue("SurfaceTemperature", surface_temperature_file);
  config.PeekValue("PAR", PAR_file);
  config.PeekValue("SoilMoisture", SoilMoisture_file);
  config.PeekValue("Aggregation_file", Aggregation_file);

  // An extra step available for PAR?
  bool extra_step = (int(file_size(PAR_file) / (Nt * Nx * Ny * sizeof(float)))
                     > days + 1);

  // Number of land use categories.
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny * Nx);
  config.PeekValue("Directory_bio", dir_out);

  config.SetSection("[LAI]");
  config.PeekValue("Directory_LAI", dir_LAI);

  config.SetSection("[emissions_factors]");
  config.PeekValue("Directory_EF", dir_EF);
  config.PeekValue("Nx", Nx_EFin);
  config.PeekValue("Ny", Ny_EFin);
  config.PeekValue("Delta_x", Delta_x_EF);
  config.PeekValue("Delta_y", Delta_y_EF);
  config.PeekValue("x_min", xmin_EF);
  config.PeekValue("y_min", ymin_EF);

  config.SetSection("[biogenic]");
  config.PeekValue("Delta_t", "> 0", Delta_t_bio);
  Nt_bio = int(T(Nt) * Delta_t / Delta_t_bio);

  vector<string> biogenic_names, output_names;
  int Nspecies, Noutspecies;
  config.SetSection("[Input_MEGAN]");
  config.Find("Species");
  split(config.GetLine(), biogenic_names);
  Nspecies = int(biogenic_names.size());
  config.Find("OutputSpecies");
  split(config.GetLine(), output_names);
  Noutspecies = int(output_names.size());

  string wpfile, withwp;

  for (int s = 0; s < Nspecies; s++)
    assert_defined_species(biogenic_names[s]);


  /////////////////////////
  // METEOROLOGICAL DATA //
  /////////////////////////


  cout << "Reading and interpolating meteorological data...";
  cout.flush();

  // Input fields.
  RegularGrid<T> GridT(t_min, Delta_t, Nt);
  RegularGrid<T> GridT_cumulated(t_min - Delta_t / 2.0, Delta_t,
                                 extra_step ? Nt + 1 : Nt);
  RegularGrid<T> GridZ(Nz);
  RegularGrid<T> GridY(y_min, Delta_y, Ny);
  RegularGrid<T> GridX(x_min, Delta_x, Nx);
  RegularGrid<T> GridC(Nc);
  RegularGrid<T> GridY_EFin(ymin_EF, Delta_y_EF, Ny_EFin);
  RegularGrid<T> GridX_EFin(xmin_EF, Delta_x_EF, Nx_EFin);

  Data<T, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<T, 3> PAR(GridT_cumulated, GridY, GridX);
  Data<T, 3> SoilMoisture(GridT, GridY, GridX);

  FormatBinary<float> Meteo;
  Meteo.ReadSteps(surface_temperature_file, step, SurfaceTemperature);
  Meteo.ReadSteps(PAR_file, step, PAR);
  Meteo.ReadSteps(SoilMoisture_file, step, SoilMoisture);

  // Output fields.
  RegularGrid<T> GridT_out(t_min, Delta_t_bio, Nt_bio);

  Data<T, 3> SurfaceTemperature_out(GridT_out, GridY, GridX);
  Data<T, 3> PAR_out(GridT_out, GridY, GridX);
  Data<T, 3> SoilMoisture_out(GridT_out, GridY, GridX);

  // Temporal linear interpolation.
  LinearInterpolationDimension(SurfaceTemperature, SurfaceTemperature_out, 0);
  LinearInterpolationDimension(PAR, PAR_out, 0);
  LinearInterpolationDimension(SoilMoisture, SoilMoisture_out, 0);
  PAR_out.ThresholdMin(0.);

  cout << " done." << endl;


  ///////////////////
  // EMISSION DATA //
  ///////////////////


  string EF_file, EF_var;

  cout << "Reading LAI data...";
  cout.flush();

  //Monthly LAI.
  string LAI01_file(dir_LAI + "laiv200301_150sec.nc"),
    LAI02_file(dir_LAI + "laiv200302_150sec.nc"),
    LAI03_file(dir_LAI + "laiv200303_150sec.nc"),
    LAI04_file(dir_LAI + "laiv200304_150sec.nc"),
    LAI05_file(dir_LAI + "laiv200305_150sec.nc"),
    LAI06_file(dir_LAI + "laiv200306_150sec.nc"),
    LAI07_file(dir_LAI + "laiv200307_150sec.nc"),
    LAI08_file(dir_LAI + "laiv200308_150sec.nc"),
    LAI09_file(dir_LAI + "laiv200309_150sec.nc"),
    LAI10_file(dir_LAI + "laiv200310_150sec.nc"),
    LAI11_file(dir_LAI + "laiv200311_150sec.nc"),
    LAI12_file(dir_LAI + "laiv200312_150sec.nc");

  FormatNetCDF<float> LAIdata;
  LAIdata.ReadDimension(LAI01_file, "lon", 0, Nx_in);
  LAIdata.ReadDimension(LAI01_file, "lat", 0, Ny_in);

  RegularGrid<T> GridX_in(Nx_in),
    GridY_in(Ny_in);
  RegularGrid<T> GridMonth(12);
  RegularGrid<T> GridSpecies(Nspecies);
  RegularGrid<T> GridDays(int(Nt_bio / 24.));
  Data<T, 3> EF_in(GridSpecies, GridY_EFin, GridX_EFin);
  Data<T, 2> EF_temp(GridY_EFin, GridX_EFin);
  Data<T, 3> EF(GridSpecies, GridY, GridX);
  Data<T, 2> EF2(GridY, GridX);
  Data<T, 2> WP(GridY, GridX);
  Data<T, 2> LAI01(GridY_in, GridX_in),
    LAI02(GridY_in, GridX_in),
    LAI03(GridY_in, GridX_in),
    LAI04(GridY_in, GridX_in),
    LAI05(GridY_in, GridX_in),
    LAI06(GridY_in, GridX_in),
    LAI07(GridY_in, GridX_in),
    LAI08(GridY_in, GridX_in),
    LAI09(GridY_in, GridX_in),
    LAI10(GridY_in, GridX_in),
    LAI11(GridY_in, GridX_in),
    LAI12(GridY_in, GridX_in),
    LAIp(GridY_in, GridX_in),
    LAIc(GridY_in, GridX_in);
  Data<T, 3> LAI(GridMonth, GridY, GridX),
    LAI_in(GridMonth, GridY_in, GridX_in);
  Data<T, 3> DailyTemperature(GridDays, GridY, GridX),
    DailyPAR(GridDays, GridY, GridX);
  Data<T, 4> emis_out(GridSpecies, GridY, GridX, GridT_out);

  LAIdata.Read(LAI01_file, "lon", GridX_in);
  LAIdata.Read(LAI01_file, "lat", GridY_in);

  LAIdata.Read(LAI01_file, "LAI_for_Jan_2003_(m2_per_m2)", LAI01);
  LAIdata.Read(LAI02_file, "LAI_for_Feb_2003_(m2_per_m2)", LAI02);
  LAIdata.Read(LAI03_file, "LAI_for_Mar_2003_(m2_per_m2)", LAI03);
  LAIdata.Read(LAI04_file, "LAI_for_Apr_2003_(m2_per_m2)", LAI04);
  LAIdata.Read(LAI05_file, "LAI_for_May_2003_(m2_per_m2)", LAI05);
  LAIdata.Read(LAI06_file, "LAI_for_Jun_2003_(m2_per_m2)", LAI06);
  LAIdata.Read(LAI07_file, "LAI_for_Jul_2003_(m2_per_m2)", LAI07);
  LAIdata.Read(LAI08_file, "LAI_for_Aug_2003_(m2_per_m2)", LAI08);
  LAIdata.Read(LAI09_file, "LAI_for_Sep_2003_(m2_per_m2)", LAI09);
  LAIdata.Read(LAI10_file, "LAI_for_Oct_2003_(m2_per_m2)", LAI10);
  LAIdata.Read(LAI11_file, "LAI_for_Nov_2003_(m2_per_m2)", LAI11);
  LAIdata.Read(LAI12_file, "LAI_for_Dec_2003_(m2_per_m2)", LAI12);

  for (int i = 0; i < Nx_in; i++)
    for (int j = 0; j < Ny_in; j++)
      {
        LAI_in(0, j, i) = LAI01(j, i) / 1000.0;
        LAI_in(1, j, i) = LAI02(j, i) / 1000.0;
        LAI_in(2, j, i) = LAI03(j, i) / 1000.0;
        LAI_in(3, j, i) = LAI04(j, i) / 1000.0;
        LAI_in(4, j, i) = LAI05(j, i) / 1000.0;
        LAI_in(5, j, i) = LAI06(j, i) / 1000.0;
        LAI_in(6, j, i) = LAI07(j, i) / 1000.0;
        LAI_in(7, j, i) = LAI08(j, i) / 1000.0;
        LAI_in(8, j, i) = LAI09(j, i) / 1000.0;
        LAI_in(9, j, i) = LAI10(j, i) / 1000.0;
        LAI_in(10, j, i) = LAI11(j, i) / 1000.0;
        LAI_in(11, j, i) = LAI12(j, i) / 1000.0;
        for (int k = 0; k < 12; k++)
          if (LAI_in(k, j, i) < 0.0)
            LAI_in(k, j, i) = 0.0;
      }

  LAI_in.ReverseData(1);

  T Delta_x_LAI, Delta_y_LAI, xmin_LAI, ymin_LAI;
  T total_area;
  int imin, imax, jmin, jmax;

  Delta_x_LAI = GridX_in(1) - GridX_in(0);
  Delta_y_LAI = GridY_in(0) - GridY_in(1);
  xmin_LAI = GridX_in(0);
  ymin_LAI = GridY_in(Ny_in - 1);

  RegularGrid<T> GridY_LAI(ymin_LAI, Delta_y_LAI, Ny_in);
  RegularGrid<T> GridX_LAI(xmin_LAI, Delta_x_LAI, Nx_in);

  T Delta_y_LAI2, Delta_x_LAI2;
  int Ny_LAIin2, Nx_LAIin2;
  T xmin_LAI2, ymin_LAI2;

  imin = (GridX(0) - Delta_x / 2 - xmin_LAI) / Delta_x_LAI;
  imax = (GridX(Nx - 1) + Delta_x / 2 - xmin_LAI) / Delta_x_LAI + 1;
  jmin = (GridY(0) - Delta_y / 2 - ymin_LAI) / Delta_y_LAI;
  jmax = (GridY(Ny - 1) + Delta_y / 2 - ymin_LAI) / Delta_y_LAI + 1;

  xmin_LAI2 = GridX_LAI(imin);
  ymin_LAI2 = GridY_LAI(jmin);

  Delta_x_LAI2 = Delta_x_LAI / nred;
  Delta_y_LAI2 = Delta_y_LAI / nred;
  Nx_LAIin2 = (imax - imin + 1) * nred;
  Ny_LAIin2 = (jmax - jmin + 1) * nred;

  RegularGrid<T> GridX_LAI2(xmin_LAI2, Delta_x_LAI2, Nx_LAIin2),
    GridY_LAI2(ymin_LAI2, Delta_y_LAI2, Ny_LAIin2);
  Data<T, 3> LAI_in2(GridMonth, GridY_LAI2, GridX_LAI2);
  LAI_in2.Fill(0.0);

  for (int i = 0; i < Nx_LAIin2; i++)
    for (int j = 0; j < Ny_LAIin2; j++)
      for (int k = 0; k < 12; k++)
        LAI_in2(k, j, i) = LAI_in(k, jmin + j / nred, imin + i / nred);

  LAI.Fill(0.0);
  for (int k = 0; k < Nx; k++)
    for (int l = 0; l < Ny; l++)
      {
        imin = (GridX(k) - Delta_x / 2. - xmin_LAI2) / (Delta_x_LAI2);
        imax = (GridX(k) + Delta_x / 2. - xmin_LAI2) / (Delta_x_LAI2) + 1.;
        jmin = (GridY(l) - Delta_y / 2. - ymin_LAI2) / (Delta_y_LAI2);
        jmax = (GridY(l) + Delta_y / 2. - ymin_LAI2) / (Delta_y_LAI2) + 1.;
        for (int i = imin; i <= imax; i++)
          for (int j = jmin; j <= jmax; j++)
            if ((GridX(k) - Delta_x / 2. <= GridX_LAI2(i)) &&
                (GridX(k) + Delta_x / 2. > GridX_LAI2(i)) &&
                (GridY(l) - Delta_y / 2. <= GridY_LAI2(j)) &&
                (GridY(l) + Delta_y / 2. > GridY_LAI2(j)))
              {
                T area = compute_area(GridY_LAI2(j), Delta_x_LAI2,
                                      Delta_y_LAI2);
                for (int s = 0; s < 12; s++)
                  LAI(s, l, k) += LAI_in2(s, j, i) * area;
              }
        total_area = compute_area(GridY(l), Delta_x, Delta_y);
        for (int s = 0; s < 12; s++)
          LAI(s, l, k) /= total_area;
      }

  cout << " done." << endl;
  cout << "Reading Emissions Factors...";
  cout.flush();

  FormatBinary<float> EFdata;

  for (int s = 0; s < Nspecies; s++)
    {
      EF_file = dir_EF + biogenic_names[s] + ".bin";
      ifstream EF_stream(EF_file.c_str());
      EFdata.Read(EF_stream, EF_temp);
      for (int i = 0; i < Nx_EFin; i++)
        for (int j = 0; j < Ny_EFin; j++)
          EF_in(s, j, i) = EF_temp(j, i);
    }

  T Delta_x_EF2, Delta_y_EF2, xmin_EF2, ymin_EF2;
  int Nx_EFin2, Ny_EFin2;

  imin = (GridX(0) - Delta_x / 2 - xmin_EF) / Delta_x_EF;
  imax = (GridX(Nx - 1) + Delta_x / 2 - xmin_EF) / Delta_x_EF + 1;
  jmin = (GridY(0) - Delta_y / 2 - ymin_EF) / Delta_y_EF;
  jmax = (GridY(Ny - 1) + Delta_y / 2 - ymin_EF) / Delta_y_EF + 1;

  xmin_EF2 = GridX_EFin(imin);
  ymin_EF2 = GridY_EFin(jmin);

  Delta_x_EF2 = Delta_x_EF / nred;
  Delta_y_EF2 = Delta_y_EF / nred;
  Nx_EFin2 = (imax - imin + 1) * nred;
  Ny_EFin2 = (jmax - jmin + 1) * nred;

  RegularGrid<T> GridX_EFin2(xmin_EF2, Delta_x_EF2, Nx_EFin2),
    GridY_EFin2(ymin_EF2, Delta_y_EF2, Ny_EFin2);
  Data<T, 3> EF_in2(GridSpecies, GridY_EFin2, GridX_EFin2);

  EF_in2.Fill(0.0);
  for (int i = 0; i < Nx_EFin2; i++)
    for (int j = 0; j < Ny_EFin2; j++)
      for (int s = 0; s < Nspecies; s++)
        EF_in2(s, j, i) = EF_in(s, jmin + j / nred, imin + i / nred);

  EF.Fill(0.0);
  for (int k = 0; k < Nx; k++)
    for (int l = 0; l < Ny; l++)
      {
        imin = (GridX(k) - Delta_x / 2 - xmin_EF2) / (Delta_x_EF2);
        imax = (GridX(k) + Delta_x / 2 - xmin_EF2) / (Delta_x_EF2) + 1;
        jmin = (GridY(l) - Delta_y / 2 - ymin_EF2) / (Delta_y_EF2);
        jmax = (GridY(l) + Delta_y / 2 - ymin_EF2) / (Delta_y_EF2) + 1;
        for (int i = imin; i <= imax; i++)
          for (int j = jmin; j <= jmax; j++)
            if ((GridX(k) - Delta_x / 2 <= GridX_EFin2(i)) &&
                (GridX(k) + Delta_x / 2 > GridX_EFin2(i)) &&
                (GridY(l) - Delta_y / 2 <= GridY_EFin2(j)) &&
                (GridY(l) + Delta_y / 2 > GridY_EFin2(j)))
              {
                T area = compute_area(GridY_EFin2(j), Delta_x_EF2,
                                      Delta_y_EF2);
                for (int s = 0; s < Nspecies; s++)
                  EF(s, l, k) += EF_in2(s, j, i) * area;
              }
        total_area = compute_area(GridY(l), Delta_x, Delta_y);
        for (int s = 0; s < Nspecies; s++)
          EF(s, l, k) /= total_area;
      }

  cout << " done." << endl;


  ///////////////
  // EMISSIONS //
  ///////////////


  cout << "Computing biogenic emissions...";
  cout.flush();
  Date date_str;
  date_str = date_meteo_str;

  DailyTemperature.Fill(0.);
  DailyPAR.Fill(0.);
  for (int k = 0; k < Nx; k++)
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nt_bio; i++)
        {
          DailyTemperature(i / 24, j, k)
            += SurfaceTemperature_out(i, j, k) / 24.0;
          DailyPAR(i / 24, j, k)
            += PAR_out(i, j, k) / 24.0;
        }

  T gamma_temp, gamma_light, gamma_AGE, gamma_SM, gamma_LAI;
  T zenith_angle;
  T ut;
  T localhour, lon, lat, sinangle, LDF;
  const T rhocanopy = 0.96;

  // Length of the time step (days) between the previous LAI (LAIp) and the
  // current LAI (LAIc).
  const T tstlen = 30.0;

  sinangle = 0.0;
  LDF = 0.0;

  Data<T, 3> angle(GridT_out, GridY, GridX);

  for (int i = 0; i < Nt_bio; i++)
    {
      for (int k = 0; k < Nx; k++)
        {
          lon = GridX(k);
          localhour = date_beg.GetHour() + lon / 15;
          ut = T(date_beg.GetHour());
          for (int j = 0; j < Ny; j++)
            {
              lat = y_min + Delta_y * j;
              if (date_beg.GetMonth() == 1)
                LAIp(j, k) = LAI(11, j, k);
              else
                LAIp(j, k) = LAI(date_beg.GetMonth() - 2, j, k);
              LAIc(j, k) = LAI(date_beg.GetMonth() - 1, j, k);
              gamma_light = 1.;
              gamma_LAI = 1.;
              gamma_SM = 1.;
              gamma_LAI = GammaLAI(LAIc(j, k));

              zenith_angle = ZenithAngle(GridX(k), GridY(j),
                                         date_beg, ut) * 3.14159 / 180.0;
              angle(i, j, k) = ZenithAngle(GridX(k), GridY(j), date_beg, ut);
              sinangle = sin(1.5708 - zenith_angle);

              gamma_light = GammaLight(date_beg.GetOrdinalDay(), sinangle,
                                       PAR_out(i, j, k), DailyPAR(i / 24, j, k));
              for (int s = 0; s < Nspecies; s++)
                {
                  emis_out(s, j, k, i) = 0.0;
                  gamma_temp = 1.;
                  gamma_AGE = 1.;
                  if (biogenic_names[s] == "ISOP")
                    gamma_temp = GammaTempISOP(SurfaceTemperature_out(i, j, k),
                                               DailyTemperature(i / 24, j, k));
                  else
                    gamma_temp = GammaTemp(biogenic_names[s],
                                           SurfaceTemperature_out(i, j, k));
                  gamma_AGE = GammaAGE(biogenic_names[s], LAIp(j, k),
                                       LAIc(j, k), DailyTemperature(i / 24, j, k), tstlen);
                  LDF = light_dependent_factor<T>(biogenic_names[s]);
                  emis_out(s, j, k, i) = EF(s, j, k)
                    * gamma_AGE * gamma_temp * gamma_SM * gamma_LAI
                    * rhocanopy * (1. - LDF + gamma_light * LDF) / 3600.0;
                }
            }
        }

      date_beg.AddHours(Delta_t_bio);
    }

  cout << " done." << endl;


  /////////////////
  // WRITES DATA //
  /////////////////


  cout << "Writing output emissions...";
  cout.flush();

  RegularGrid<T> GridOutSpecies(Noutspecies);
  Data<T, 4> emissions(GridOutSpecies, GridT_out, GridY, GridX);

  MeganAggregation(emis_out, Aggregation_file, biogenic_names, output_names,
                   emissions);

  FormatBinary<float> Output;
  for (int s = 0; s < Noutspecies; s++)
    {
      Data<T, 3> Emis_extract(Nt_bio, Ny, Nx);

      Emis_extract.SubData(emissions,
                           s, Range::all(), Range::all(), Range::all());
      Output.Append(Emis_extract, dir_out + output_names[s] + ".bin");
    }

  cout << " done." << endl;

  END;

  return 0;
}
