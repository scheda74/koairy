// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).

///////////////
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


// INCLUDES //
/////////////

int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file,
    default_name("MM5-meteo-castor.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, date_end, default_name);

  Date date_prev = date_beg;
  date_prev.AddDays(-1);

  ///////////////////
  // CONFIGURATION //
  ///////////////////


  typedef float real;

  int h, i, j, k;

  // Constants.
  const real pi = 3.14159265358979323846264;

  cout << "Reading configuration...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (sec_config_file != "")
    configuration.AddFile(sec_config_file);


  /*** MM5 ***/

  // Input dimensions for meteological data.
  int Nt_in, Nz_in, Ny_in, Nx_in;
  // Input-data steps.
  real Delta_t_in, Delta_y_in, Delta_x_in;
  real y_min_in, x_min_in;

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

  int projection_type;
  configuration.PeekValue("projection_type", "= 1 2 3", projection_type);

  string horizontal_interpolation;
  configuration.PeekValue("Horizontal_interpolation", "latlon | MM5",
                          horizontal_interpolation);

  string dot_coordinate_file;
  if (horizontal_interpolation == "latlon")
    configuration.PeekValue("Dot_coordinates", dot_coordinate_file);

  /*** Simulation domain ***/

  // Output dimensions.
  int Nt_out, Nz_out, Ny_out, Nx_out;
  // Output-data steps.
  real Delta_t_out, Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx_out);
  configuration.PeekValue("Ny", "> 0", Ny_out);
  configuration.PeekValue("Nz", "> 0", Nz_out);

  configuration.PeekValue("Delta_t", "> 0", Delta_t_out);
  configuration.PeekValue("Delta_y", "> 0", Delta_y_out);
  configuration.PeekValue("Delta_x", "> 0", Delta_x_out);

  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);

  string hybrid_coefficient_file;
  configuration.PeekValue("Vertical_levels", hybrid_coefficient_file);

  Nt_out = compute_Nt(date_beg, date_end, Delta_t_out);

  date_end.AddSeconds(- Delta_t_out * 3600);

  /*** Meteorological options ***/

  configuration.SetSection("[meteo]");

  // Relative humidity threshold for cloud formation.
  real relative_humidity_threshold;
  configuration.PeekValue("Relative_humidity_threshold", ">= 0 | <= 1",
                          relative_humidity_threshold);

  // Low cloud maximum height (m).
  real low_cloud_top_max;
  configuration.PeekValue("Low_cloud_top_max", "positive", low_cloud_top_max);

  /*** Vertical diffusion ***/

  // Kz thresholds.
  real Kz_min_dry;
  real Kz_min_wet;
  real Kz_min_up;
  real Kz_max;

  configuration.SetSection("[Kz]");

  configuration.PeekValue("Min_dry", "positive", Kz_min_dry);
  configuration.PeekValue("Min_wet", "positive", Kz_min_wet);
  configuration.PeekValue("Min_above_PBLH", "positive", Kz_min_up);
  configuration.PeekValue("Max", "positive", Kz_max);

  /*** Path ***/

  string directory_out, file_in, Roughness_file;

  configuration.SetSection("[paths]");

  configuration.PeekValue("Database_MM5-meteo", file_in);

  // Processes the input file name to prepare date expansion.
  file_in = find_replace(file_in, "&D", "%y-%m-%d");
  file_in = find_replace(file_in, "&y", "%y");
  file_in = find_replace(file_in, "&m", "%m");
  file_in = find_replace(file_in, "&d", "%d");
  file_in = date_beg.GetDate(file_in);

  if (!exists(file_in))
    throw string("Unable to find MM5 file \"") + file_in + "\".";

  Date date_beg_meteo = read_date_MM5(file_in);
  Date date_end_meteo = date_beg_meteo;
  date_end_meteo.AddSeconds(Nt_in * Delta_t_in * 3600);

  if (date_beg_meteo > date_beg ||  date_end_meteo < date_end)
    throw string("The MM5 file provided ") + file_in
      + " does not contain data for all dates.";

  real t_min_in, t_min_out, differences;
  t_min_in = real(date_beg_meteo.GetHour())
    + real(date_beg_meteo.GetMinutes()) / 60.
    + real(date_beg_meteo.GetSeconds()) / 3600.;
  differences = real(date_beg.GetSecondsFrom(date_beg_meteo));
  t_min_out = t_min_in + differences / 3600.;

  configuration.PeekValue("Roughness_file", Roughness_file);

  configuration.PeekValue("Directory_meteo", directory_out);

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////


  cout << "Memory allocation for grids...";
  cout.flush();

  /*** Input ***/

  // Input grids.
  RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
  // Grid for Z in height. Z depends on Z, Y, X.
  GeneralGrid<real, 3> GridZ_in(shape(Nz_in, Ny_in - 1, Nx_in - 1),
                                1, shape(1, 2, 3));
  // Vertical levels are shared.
  GridZ_in.SetVariable(1);
  GridZ_in.SetDuplicate(false);

  // Vertical levels for dot grids are not the same.
  GeneralGrid<real, 3> GridZ_Dot_in(shape(Nz_in, Ny_in, Nx_in),
                                    1, shape(1, 2, 3));
  GridZ_Dot_in.SetVariable(1);
  GridZ_Dot_in.SetDuplicate(false);


  // Input grid in MM5 coordinates (z, x, y).
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
  RegularGrid<real> GridZ_ground_in(Nz_in + 1);

  // Interfaces.
  RegularGrid<real> GridX_interf_in(x_min_in - Delta_x_in / 2.,
                                    Delta_x_in, Nx_in);
  RegularGrid<real> GridY_interf_in(y_min_in - Delta_y_in / 2.,
                                    Delta_y_in, Ny_in);

  /*** Output ***/

  // Output grids.
  RegularGrid<real> GridT_out(t_min_out, Delta_t_out, Nt_out);
  RegularGrid<real> GridZ_out(Nz_out);
  // Data may be provided on interfaces.
  RegularGrid<real> GridZ_interf_out(Nz_out + 1);

  // Latitudes and longitudes.
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  GeneralGrid<real, 2> GridX_2D_out(shape(Ny_out, Nx_out), 2, shape(1, 2));
  GeneralGrid<real, 2> GridY_2D_out(shape(Ny_out, Nx_out), 1, shape(1, 2));
  GeneralGrid<real, 2> GridX_3D_out(shape(Ny_out, Nx_out), 3, shape(2, 3));
  GeneralGrid<real, 2> GridY_3D_out(shape(Ny_out, Nx_out), 2, shape(2, 3));

  GridX_2D_out.SetVariable(2);
  GridX_2D_out.SetDuplicate(false);
  GridY_2D_out.SetVariable(1);
  GridY_2D_out.SetDuplicate(false);
  GridX_3D_out.SetVariable(3);
  GridX_3D_out.SetDuplicate(false);
  GridY_3D_out.SetVariable(2);
  GridY_3D_out.SetDuplicate(false);

  cout << " done." << endl;


  /////////////////
  // OUTPUT DATA //
  /////////////////


  cout << "Memory allocation for output data fields...";
  cout.flush();

  Data<real, 3> RoughnessHeight(12, Ny_out, Nx_out);

  FormatBinary<float>().Read(Roughness_file, RoughnessHeight);

  // 3D Output fields.
  // Cross fields.
  Data<real, 4> Temperature_out(GridT_out, GridZ_out,
                                GridY_3D_out, GridX_3D_out);
  Data<real, 4> Pressure_out(GridT_out, GridZ_out,
                             GridY_3D_out, GridX_3D_out);
  Data<real, 4> SpecificHumidity_out(GridT_out, GridZ_out,
                                     GridY_3D_out, GridX_3D_out);
  Data<real, 4> LiquidWaterContent_out(GridT_out, GridZ_out,
                                       GridY_3D_out, GridX_3D_out);
  Data<real, 3> Attenuation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 4> AirDensity_out(GridT_out, GridZ_out,
                               GridY_3D_out, GridX_3D_out);
  // Dot fields.
  Data<real, 4> MeridionalWind_out(GridT_out, GridZ_out,
                                   GridY_3D_out, GridX_3D_out);
  Data<real, 4> ZonalWind_out(GridT_out, GridZ_out,
                              GridY_3D_out, GridX_3D_out);

  // Diagnosed fields.
  Data<real, 4> WindModule_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> PotentialTemperature_out(GridT_out, GridZ_out,
                                         GridY_out, GridX_out);
  Data<real, 4> Kz_out(GridT_out, GridZ_interf_out, GridY_out, GridX_out);

  // 2D Output fields.
  Data<real, 3> SurfacePressure_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SurfaceTemperature_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> BoundaryHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> CloudBaseHeight_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FrictionModule_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> SensibleHeat_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> Evaporation_out(GridT_out, GridY_2D_out, GridX_2D_out);
  Data<real, 3> FirstLevelWindModule_out(GridT_out, GridY_2D_out,
                                         GridX_2D_out);

  cout << " done." << endl;


  ////////////////////////////////
  // SIGMA TO HEIGHT CONVERSION //
  ////////////////////////////////


  cout << "Conversion from sigma levels to altitudes...";
  cout.flush();

  FormatMM5 InputMeteo;

  // Data to convert Z from sigma to heights.
  Data<real, 2> SigmaH(Nt_in, Nz_in);
  // For full sigma levels.
  Data<real, 1> SigmaF(Nz_in + 1);

  // Reads data to convert Z from sigma levels to heights.
  InputMeteo.ReadWholeField(file_in, "SIGMAH", SigmaH);
  SigmaH.ReverseData(1);
  // Reference pressure.
  Data<real, 3> ReferencePressure(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "PSTARCRS", ReferencePressure);
  ReferencePressure.SwitchDimensions(shape(0, 2, 1),
                                     GridT_in, GridY_in, GridX_in);

  // Big header data and general constants.
  Array<int, 2> BHI;
  Array<float, 2> BHR;
  Array<string, 2> BHIC;
  Array<string, 2> BHRC;
  real p00, ptop, Ts0, A;
  real refpress, reftemp;
  const real g(9.81);
  const real R(287.04);

  InputMeteo.ReadBigHeader(file_in, BHI, BHR, BHIC, BHRC);

  p00 = BHR(4, 1);
  ptop = BHR(1, 1);
  Ts0 = BHR(4, 2);
  A = BHR(4, 3);

  SigmaF(0) = 1.;
  for (k = 0; k < SigmaH.GetLength(1); k++)
    SigmaF(k + 1) = max(0., 2. * SigmaH(0, k) - SigmaF(k));

  // At cross points.
  real altitude, thickness;
  for (j = 0; j < Ny_in - 1; j++)
    for (i = 0; i < Nx_in - 1; i++)
      {
        altitude = 0.;
        for (k = 0; k < SigmaH.GetLength(1); k++)
          {
            refpress = ReferencePressure(0, j, i) * SigmaH(0, k) + ptop;
            reftemp = Ts0 + A * log(refpress / p00);
            thickness = R / g * reftemp
              * log((ReferencePressure(0, j, i) * SigmaF(k) + ptop)
                    / (ReferencePressure(0, j, i) * SigmaH(0, k) + ptop));
            GridZ_in.Value(0, k, j, i) = altitude + thickness;
            thickness = R / g * reftemp
              * log((ReferencePressure(0, j, i) * SigmaF(k) + ptop)
                    / (ReferencePressure(0, j, i) * SigmaF(k + 1) + ptop));
            altitude += thickness;
          }
      }

  // At dot points.
  for (k = 0; k < SigmaH.GetLength(1); k++)
    for (j = 0; j < Ny_in; j++)
      for (i = 0; i < Nx_in; i++)
        // Warning: this is only an approximation. 'GridZ_Dot_in' should not
        // be used as such.
        GridZ_Dot_in.Value(0, k, j, i)
          = GridZ_in.Value(0, k, min(j, Ny_in - 2), min(i, Nx_in - 2));

  cout << " done." << endl;

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
    }

  cout << " done." << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////


  cout << "Computing pressure...";
  cout.flush();
  cout.flush();

  // MM5 3D fields.
  Data<real, 4> PressurePerturbation(GridT_in, GridZ_MM5_in,
                                     GridX_in, GridY_in);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);

  InputMeteo.ReadWholeField(file_in, "PP", PressurePerturbation);
  PressurePerturbation.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                        GridY_in, GridX_in);
  PressurePerturbation.ReverseData(1);

  Data<real, 4> Temperature(GridT_in, GridZ_MM5_in, GridX_in, GridY_in);

  InputMeteo.ReadWholeField(file_in, "T", Temperature);
  Temperature.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                               GridY_in, GridX_in);
  Temperature.ReverseData(1);

  for (h = 0; h < Nt_in; h++)
    for (k = 0; k < Nz_in; k++)
      for (j = 0; j < Ny_in - 1; j++)
        for (i = 0; i < Nx_in - 1; i++)
          Pressure(h, k, j, i) = ReferencePressure(h, j, i) * SigmaH(0, k)
            + ptop + PressurePerturbation(h, k, j, i);

  PressurePerturbation.Resize();
  ReferencePressure.Resize();

  cout << " done." << endl;

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

  cout << " done." << endl;


  ///////////////////
  // WIND ROTATION //
  ///////////////////


  // Interpolations.
  cout << "Wind rotation...";
  cout.flush();

  // Winds.
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

  // Cosine and sine associated with rotation angles.
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
      Cosine(i, Ny_in - 1) = 2. * Cosine(i, Ny_in - 2) - Cosine(i, Ny_in - 3);
      Sine(i, Ny_in - 1) = 2. * Sine(i, Ny_in - 2) - Sine(i, Ny_in - 3);
    }
  for (j = 0; j < Ny_in; j++)
    {
      Cosine(0, j) = 2. * Cosine(1, j) - Cosine(2, j);
      Sine(0, j) = 2. * Sine(1, j) - Sine(2, j);
      Cosine(Nx_in - 1, j) = 2. * Cosine(Nx_in - 2, j) - Cosine(Nx_in - 3, j);
      Sine(Nx_in - 1, j) = 2. * Sine(Nx_in - 2, j) - Sine(Nx_in - 3, j);
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

  MeridionalWind.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_Dot_in,
                                  GridY_interf_in, GridX_interf_in);
  MeridionalWind.ReverseData(1);

  ZonalWind.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_Dot_in,
                             GridY_interf_in, GridX_interf_in);
  ZonalWind.ReverseData(1);

  cout << " done." << endl;


  //////////////////////////////
  // HORIZONTAL INTERPOLATION //
  //////////////////////////////


  // Interpolations.
  cout << "Horizontal interpolations...";
  cout.flush();

  // Fields interpolated on output horizontal grid, but not interpolated on
  // the vertical. '_h_out' is appended to their name.
  Data<real, 4> MeridionalWind_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> ZonalWind_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> Pressure_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> Temperature_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 3> Temperature_2m_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> SoilMoisture_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> U_10m_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> V_10m_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> WindModule_10m_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> PotentialTemperature_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> VirtualPotentialTemperature_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);

  Data<real, 4> SpecificHumidity_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> RelativeHumidity_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 3> RelativeHumidityMax_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> LiquidWaterContent_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> IceContent_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> RainWater_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);
  Data<real, 4> LWPath_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);

  // Heat fluxes.
  Data<real, 3> SensibleHeat_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> LatentHeat_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PotentialHeatFlux_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> PotentialThermal_h_out(GridT_out, GridY_out, GridX_out);

  Data<real, 3> FrictionModule_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> BoundaryHeight_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> LowCloudTop_h_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> LowCloudTopW_h_out(GridT_out, GridY_out, GridX_out);

  // Diagnosed air density.
  Data<real, 4> AirDensity_h_out(GridT_out, GridZ_ground_in, GridY_out, GridX_out);

  // Altitudes.
  Data<real, 3> Altitude_h_out(GridZ_ground_in, GridY_out, GridX_out);

  // Indices and weights for horizontal interpolation of cross and dot points.
  Array<int, 2> Index_i_in(Ny_out, Nx_out);
  Array<int, 2> Index_j_in(Ny_out, Nx_out);
  Array<int, 2> Index_dot_i_in(Ny_out, Nx_out);
  Array<int, 2> Index_dot_j_in(Ny_out, Nx_out);
  Array<real, 2> Weight_i_in(Ny_out, Nx_out);
  Array<real, 2> Weight_j_in(Ny_out, Nx_out);
  Array<real, 2> Weight_dot_i_in(Ny_out, Nx_out);
  Array<real, 2> Weight_dot_j_in(Ny_out, Nx_out);

  if (horizontal_interpolation == "MM5")
    {
      InterpolationRegularInput(x_min_in, Delta_x_in, y_min_in, Delta_y_in,
                                GridX_2D_out.GetArray(),
                                GridY_2D_out.GetArray(),
                                Index_i_in, Index_j_in,
                                Weight_i_in, Weight_j_in);
      InterpolationRegularInput(x_min_in - Delta_x_in / real(2.), Delta_x_in,
                                y_min_in - Delta_y_in / real(2.), Delta_y_in,
                                GridX_2D_out.GetArray(),
                                GridY_2D_out.GetArray(),
                                Index_dot_i_in, Index_dot_j_in,
                                Weight_dot_i_in, Weight_dot_j_in);
    }
  else // horizontal_interpolation == "latlon"
    {
      Array<real, 3> Longitude_tmp(Nt_in, Nx_in - 1, Ny_in - 1);
      Array<real, 3> Latitude_tmp(Nt_in, Nx_in - 1, Ny_in - 1);

      InputMeteo.ReadWholeField(file_in, "LONGICRS", Longitude_tmp);
      InputMeteo.ReadWholeField(file_in, "LATITCRS", Latitude_tmp);

      Array<real, 2> Longitude_in(Nx_in - 1, Ny_in - 1);
      Array<real, 2> Latitude_in(Nx_in - 1, Ny_in - 1);

      for (i = 0; i < Nx_in - 1; i++)
        for (j = 0; j < Ny_in - 1; j++)
          {
            Longitude_in(i, j) = Longitude_tmp(0, i, j);
            Latitude_in(i, j) = Latitude_tmp(0, i, j);
          }

      // For dot points, coordinates are temporarily read in a text file.
      Array<real, 2> Longitude_dot_in(Nx_in, Ny_in);
      Array<real, 2> Latitude_dot_in(Nx_in, Ny_in);

      ifstream coordinate_stream(dot_coordinate_file.c_str());
      for (j = 0; j < Ny_in; j++)
        for (i = 0; i < Nx_in; i++)
          {
            coordinate_stream >> Longitude_dot_in(i, j);
            coordinate_stream >> Latitude_dot_in(i, j);
          }
      if (!coordinate_stream.good())
        throw string("Reading MM5 dot coordinates in file \"")
          + dot_coordinate_file + "\" failed.";

      InterpolationChimere(Longitude_in, Latitude_in,
                           x_min_out, Delta_x_out, Nx_out,
                           y_min_out, Delta_y_out, Ny_out,
                           Index_i_in, Index_j_in,
                           Weight_i_in, Weight_j_in);

      InterpolationChimere(Longitude_dot_in, Latitude_dot_in,
                           x_min_out, Delta_x_out, Nx_out,
                           y_min_out, Delta_y_out, Ny_out,
                           Index_dot_i_in, Index_dot_j_in,
                           Weight_dot_i_in, Weight_dot_j_in);
    }

  // Horizontal interpolation of altitudes.
  Altitude_h_out() = 0.;
  Array<real, 3> Altitude_h_out_extract(Altitude_h_out(), Range(1, Nz_in));
  HorizontalInterpolation(GridZ_in.GetArray(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, Altitude_h_out_extract);

  // Horizontal interpolation of pressure and temperature.
  Data<real, 4> Pressure_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(Pressure(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 Pressure_tmp());

  LinearInterpolationDimension(Pressure_tmp, Pressure_h_out, 0);
  Pressure_tmp.Resize();

  Data<real, 4> Temperature_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(Temperature(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 Temperature_tmp());

  LinearInterpolationDimension(Temperature_tmp, Temperature_h_out, 0);
  Temperature_tmp.Resize();

  // 2m temperature.
  Data<real, 3> Temperature_2m(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "T2", Temperature_2m);
  Temperature_2m.SwitchDimensions(shape(0, 2, 1),
                                  GridT_in, GridY_in, GridX_in);

  Data<real, 3> Temperature_2m_tmp(GridT_in, GridY_out, GridX_out);
  HorizontalInterpolation(Temperature_2m(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, Temperature_2m_tmp());
  Temperature_2m.Resize();
  LinearInterpolationDimension(Temperature_2m_tmp, Temperature_2m_h_out, 0);
  Temperature_2m_tmp.Resize();

  // Soil moisture (first level).
  Data<real, 3> SoilMoisture(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "SOIL M 1", SoilMoisture);
  SoilMoisture.SwitchDimensions(shape(0, 2, 1), GridT_in, GridY_in, GridX_in);
  Data<real, 3> SoilMoisture_tmp(GridT_in, GridY_out, GridX_out);
  HorizontalInterpolation(SoilMoisture(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, SoilMoisture_tmp());
  SoilMoisture.Resize();
  LinearInterpolationDimension(SoilMoisture_tmp, SoilMoisture_h_out, 0);
  SoilMoisture_tmp.Resize();

  // 10m wind module.
  Data<real, 3> U_10m(GridT_in, GridX_in, GridY_in);
  Data<real, 3> V_10m(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "U10", U_10m);
  InputMeteo.ReadWholeField(file_in, "V10", V_10m);
  U_10m.SwitchDimensions(shape(0, 2, 1), GridT_in, GridY_in, GridX_in);
  V_10m.SwitchDimensions(shape(0, 2, 1), GridT_in, GridY_in, GridX_in);
  Data<real, 3> U_10m_tmp(GridT_in, GridY_out, GridX_out);
  Data<real, 3> V_10m_tmp(GridT_in, GridY_out, GridX_out);
  HorizontalInterpolation(U_10m(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, U_10m_tmp());
  HorizontalInterpolation(V_10m(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, V_10m_tmp());
  U_10m.Resize();
  V_10m.Resize();
  LinearInterpolationDimension(U_10m_tmp, U_10m_h_out, 0);
  LinearInterpolationDimension(V_10m_tmp, V_10m_h_out, 0);
  U_10m_tmp.Resize();
  V_10m_tmp.Resize();

  ComputeModule(U_10m_h_out, V_10m_h_out, WindModule_10m_h_out);

  // Special correction on temperature.
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        Temperature_h_out(h, 0, j, i) = Temperature_2m_h_out(h, j, i)
          - 2. * (Temperature_h_out(h, 1, j, i)
                  - Temperature_2m_h_out(h, j, i))
          / (Altitude_h_out(1, j, i) - 2.);

  // Horizontal interpolation of meridional wind.
  Data<real, 4> MeridionalWind_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  MeridionalWind_tmp.SetZero();
  Array<real, 4>
    MeridionalWind_extract(MeridionalWind_tmp(),
                           Range::all(), Range(1, Nz_in));
  HorizontalInterpolation(MeridionalWind(), Index_dot_i_in, Index_dot_j_in,
                          Weight_dot_i_in, Weight_dot_j_in,
                          MeridionalWind_extract);
  LinearInterpolationDimension(MeridionalWind_tmp, MeridionalWind_h_out, 0);
  MeridionalWind_tmp.Resize();

  // Horizontal interpolation of zonal wind.
  Data<real, 4> ZonalWind_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  ZonalWind_tmp.SetZero();
  Array<real, 4> ZonalWind_extract(ZonalWind_tmp(),
                                   Range::all(), Range(1, Nz_in));
  HorizontalInterpolation(ZonalWind(), Index_dot_i_in, Index_dot_j_in,
                          Weight_dot_i_in, Weight_dot_j_in,
                          ZonalWind_extract);
  LinearInterpolationDimension(ZonalWind_tmp, ZonalWind_h_out, 0);
  ZonalWind_tmp.Resize();

  Data<real, 4> Altitude_out(Nt_out, Nz_out + 1, Ny_out, Nx_out);
  FormatFormattedText InputLevels("<e><e>");
  Data<real, 1> alpha(Nz_out);
  Data<real, 1> beta(Nz_out);
  InputLevels.Read(hybrid_coefficient_file, "0", alpha);
  InputLevels.Read(hybrid_coefficient_file, "1", beta);

  Altitude_out.SetZero();
  int k_in;
  real level_pressure;
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        {
          k_in = 1;
          for (k = 0; k < Nz_out; k++)
            {
              level_pressure = alpha(k) * 1.e5
                + beta(k) * Pressure_h_out(h, 0, j, i);
              while (Pressure_h_out(h, k_in, j, i) >= level_pressure)
                k_in++;
              Altitude_out(h, k + 1, j, i) = Altitude_h_out(k_in - 1, j, i)
                + (Altitude_h_out(k_in, j, i)
                   - Altitude_h_out(k_in - 1, j, i))
                * (level_pressure - Pressure_h_out(h, k_in - 1, j, i))
                / (Pressure_h_out(h, k_in, j, i)
                   - Pressure_h_out(h, k_in - 1, j, i));
            }
        }

  /*** Humidity ***/
  Data<real, 4> SpecificHumidity(GridT_in, GridZ_MM5_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "Q", SpecificHumidity);
  SpecificHumidity.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                    GridY_in, GridX_in);
  SpecificHumidity.ReverseData(1);
  Data<real, 4> SpecificHumidity_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(SpecificHumidity(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 SpecificHumidity_tmp(), real(1.e-10));
  SpecificHumidity_tmp.ThresholdMin(0.);
  SpecificHumidity.Resize();

  LinearInterpolationDimension(SpecificHumidity_tmp, SpecificHumidity_h_out, 0);
  SpecificHumidity_tmp.Resize();

  // Water content.
  Data<real, 4> LiquidWaterContent(GridT_in, GridZ_MM5_in,
                                   GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "CLW", LiquidWaterContent);
  LiquidWaterContent.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                                      GridY_in, GridX_in);
  LiquidWaterContent.ReverseData(1);
  Data<real, 4> LiquidWaterContent_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(LiquidWaterContent(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 LiquidWaterContent_tmp(), real(1.e-10));
  LiquidWaterContent_tmp.ThresholdMin(0.);
  LiquidWaterContent.Resize();

  LinearInterpolationDimension(LiquidWaterContent_tmp, LiquidWaterContent_h_out, 0);
  LiquidWaterContent_tmp.Resize();

  Data<real, 4> IceContent(GridT_in, GridZ_MM5_in,
                           GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "ICE", IceContent);
  IceContent.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                              GridY_in, GridX_in);
  IceContent.ReverseData(1);
  Data<real, 4> IceContent_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(IceContent(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 IceContent_tmp(), real(1.e-10));
  IceContent_tmp.ThresholdMin(0.);
  IceContent.Resize();

  LinearInterpolationDimension(IceContent_tmp, IceContent_h_out, 0);
  IceContent_tmp.Resize();

  Data<real, 4> RainWater(GridT_in, GridZ_MM5_in,
                          GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "RNW", RainWater);
  RainWater.SwitchDimensions(shape(0, 1, 3, 2), GridT_in, GridZ_in,
                             GridY_in, GridX_in);
  RainWater.ReverseData(1);

  Data<real, 4> RainWater_tmp(GridT_in, GridZ_ground_in, GridY_out, GridX_out);
  HorizontalInterpolation_ground(RainWater(), Index_i_in, Index_j_in,
                                 Weight_i_in, Weight_j_in, Altitude_h_out(),
                                 RainWater_tmp(), real(1.e-10));

  RainWater_tmp.ThresholdMin(0.);
  RainWater.Resize();

  LinearInterpolationDimension(RainWater_tmp, RainWater_h_out, 0);
  RainWater_tmp.Resize();


  // Computes relative humidity and liquid water path.
  real P_sat, additional;
  for (h = 0; h < Nt_out; h++)
    for (k = 0; k < Nz_in + 1; k++)
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            P_sat = 611. * exp(17.27 *
                               (Temperature_h_out(h, k, j, i) - 273.15)
                               / (Temperature_h_out(h, k, j, i) - 35.86));

            RelativeHumidity_h_out(h, k, j, i)
              = SpecificHumidity_h_out(h, k, j, i) / (.622 * P_sat)
              * (Pressure_h_out(h, k, j, i) - P_sat);

            AirDensity_h_out(h, k, j, i) = 7.2868e16
              * Pressure_h_out(h, k, j, i) / Temperature_h_out(h, k, j, i);

            if (k == 0)
              LWPath_h_out(h, k, j, i) = 0.;
            else
              {
                additional =
                  .5 * ((LiquidWaterContent_h_out(h, k, j, i)
                         + RainWater_h_out(h, k, j, i))
                        * AirDensity_h_out(h, k, j, i)
                        + (LiquidWaterContent_h_out(h, k - 1, j, i)
                           + RainWater_h_out(h, k - 1, j, i))
                        * AirDensity_h_out(h, k - 1, j, i)) * 1.8e2
                  + .5 * (IceContent_h_out(h, k, j, i)
                          * AirDensity_h_out(h, k, j, i)
                          + IceContent_h_out(h, k - 1, j, i)
                          * AirDensity_h_out(h, k - 1, j, i)) * 60. / .9;
                LWPath_h_out(h, k, j, i) = LWPath_h_out(h, k - 1, j, i)
                  + additional / 7.2868e16 / 287.04
                  * (Altitude_h_out(k, j, i) - Altitude_h_out(k - 1, j, i));
              }
          }

  /*** Clouds ***/

  real relative_humidity_max;
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        {
          relative_humidity_max = relative_humidity_threshold;
          LowCloudTop_h_out(h, j, i) = -1.;
          for (k = 0; k < Nz_in; k++)
            if (Altitude_h_out(k, j, i) <= low_cloud_top_max)
              {
                relative_humidity_max
                  = min(1., max(relative_humidity_max,
                                RelativeHumidity_h_out(h, k, j, i)));
                if (RelativeHumidity_h_out(h, k, j, i)
                    >= relative_humidity_threshold)
                  {
                    LowCloudTop_h_out(h, j, i)
                      = Altitude_h_out(k + 1, j, i);
                    if (RelativeHumidity_h_out(h, k + 1, j, i)
                        < relative_humidity_threshold)
                      LowCloudTop_h_out(h, j, i)
                        = Altitude_h_out(k, j, i)
                        + (Altitude_h_out(k + 1, j, i)
                           - Altitude_h_out(k, j, i))
                        * (RelativeHumidity_h_out(h, k, j, i)
                           - relative_humidity_threshold)
                        / (RelativeHumidity_h_out(h, k, j, i)
                           - RelativeHumidity_h_out(h, k + 1, j, i));
                  }
              }
          LowCloudTop_h_out(h, j, i)
            = min(low_cloud_top_max, LowCloudTop_h_out(h, j, i));
          LowCloudTopW_h_out(h, j, i)
            = (relative_humidity_max - relative_humidity_threshold)
            / (1. - relative_humidity_threshold);

          RelativeHumidityMax_h_out(h, j, i) = relative_humidity_max;
        }

  /*** Heat fluxes ***/

  // Friction module.
  Data<real, 3> FrictionModule(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "UST", FrictionModule);
  FrictionModule.SwitchDimensions(shape(0, 2, 1),
                                  GridT_in, GridY_in, GridX_in);
  Data<real, 3> FrictionModule_tmp(GridT_in, GridY_out, GridX_out);
  HorizontalInterpolation(FrictionModule(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, FrictionModule_tmp());
  FrictionModule.Resize();
  LinearInterpolationDimension(FrictionModule_tmp, FrictionModule_h_out, 0);
  FrictionModule_tmp.Resize();

  // Virtual potential temperature.
  for (h = 0; h < Nt_out; h++)
    for (k = 0; k < Nz_in + 1; k++)
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            PotentialTemperature_h_out(h, k, j, i)
              = Temperature_h_out(h, k, j, i)
              * pow(1.e5 / Pressure_h_out(h, k, j, i), 0.2857);
            VirtualPotentialTemperature_h_out(h, k, j, i)
              = PotentialTemperature_h_out(h, k, j, i)
              * (1. + .61 * SpecificHumidity_h_out(h, k, j, i)
                 - LiquidWaterContent_h_out(h, k, j, i)
                 - IceContent_h_out(h, k, j, i)
                 - RainWater_h_out(h, k, j, i));
          }

  // PBLH.
  Data<real, 3> BoundaryHeight(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "PBL HGT", BoundaryHeight);
  BoundaryHeight.SwitchDimensions(shape(0, 2, 1),
                                  GridT_in, GridY_in, GridX_in);
  Data<real, 3> BoundaryHeight_tmp(GridT_in, GridY_2D_out, GridX_2D_out);
  HorizontalInterpolation(BoundaryHeight(),
                          Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in,
                          BoundaryHeight_tmp());
  BoundaryHeight_tmp.ThresholdMin(20.);
  LinearInterpolationDimension(BoundaryHeight_tmp, BoundaryHeight_h_out, 0);
  BoundaryHeight_tmp.Resize();

  // Clouds.
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        if (LowCloudTop_h_out(h, j, i) > BoundaryHeight_h_out(h, j, i))
          BoundaryHeight_h_out(h, j, i)
            = (1. - LowCloudTopW_h_out(h, j, i))
            * BoundaryHeight_h_out(h, j, i)
            + LowCloudTopW_h_out(h, j, i) * LowCloudTop_h_out(h, j, i);
  BoundaryHeight.Resize();

  // Thermals.
  real thermals_start_height = 25.;
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        for (k = 0; k < Nz_out; k++)
          if (Altitude_h_out(k, j, i) < thermals_start_height
              && Altitude_h_out(k + 1, j, i) >= thermals_start_height)
            {
              PotentialThermal_h_out(h, j, i)
                = VirtualPotentialTemperature_h_out(h, k, j, i)
                + (VirtualPotentialTemperature_h_out(h, k + 1, j, i)
                   - VirtualPotentialTemperature_h_out(h, k, j, i))
                * (thermals_start_height - Altitude_h_out(k, j, i))
                / (Altitude_h_out(k + 1, j, i)
                   - Altitude_h_out(k, j, i));
            }

  // Sensible heat.
  Data<real, 3> SensibleHeat(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "SHFLUX", SensibleHeat);
  SensibleHeat.SwitchDimensions(shape(0, 2, 1),
                                GridT_in, GridY_in, GridX_in);
  Data<real, 3> SensibleHeat_tmp(GridT_in, GridY_2D_out, GridX_2D_out);
  HorizontalInterpolation(SensibleHeat(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, SensibleHeat_tmp());
  SensibleHeat.Resize();
  LinearInterpolationDimension(SensibleHeat_tmp, SensibleHeat_h_out, 0);
  SensibleHeat_tmp.Resize();

  // Latent heat.
  Data<real, 3> LatentHeat(GridT_in, GridX_in, GridY_in);
  InputMeteo.ReadWholeField(file_in, "LHFLUX", LatentHeat);
  LatentHeat.SwitchDimensions(shape(0, 2, 1),
                              GridT_in, GridY_in, GridX_in);
  Data<real, 3> LatentHeat_tmp(GridT_in, GridY_out, GridX_out);
  HorizontalInterpolation(LatentHeat(), Index_i_in, Index_j_in,
                          Weight_i_in, Weight_j_in, LatentHeat_tmp());
  LatentHeat.Resize();
  LinearInterpolationDimension(LatentHeat_tmp, LatentHeat_h_out, 0);
  LatentHeat_tmp.Resize();

  const real r_P = 287.04 / 1005.;
  const real r_L = 287.04 / 2.45e6;
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        {
          SensibleHeat_h_out(h, j, i) *= r_P * Temperature_2m_h_out(h, j, i)
            / Pressure_h_out(h, 0, j, i);
          LatentHeat_h_out(h, j, i) *= r_L * Temperature_2m_h_out(h, j, i)
            / Pressure_h_out(h, 0, j, i);
          // The last term below is always 0. But it would be better to
          // activate the correction above urban areas. The best would be to
          // replace "0."  with the amount of urban area in the cell.
          PotentialHeatFlux_h_out(h, j, i) = SensibleHeat_h_out(h, j, i)
            * (1. + 0.61 * SpecificHumidity_h_out(h, 0, j, i))
            + LatentHeat_h_out(h, j, i) * 0.61
            * PotentialTemperature_h_out(h, 0, j, i)
            + 0. * r_P * Temperature_2m_h_out(h, j, i)
            / Pressure_h_out(h, 0, j, i);
        }

  cout << " done." << endl;


  ////////////////////////
  // VERTICAL DIFFUSION //
  ////////////////////////

  cout << "Vertical diffusion...";
  cout.flush();

  Data<real, 3> MoninObukov_h_out(Nt_out, Ny_out, Nx_out);
  Data<real, 3> ConvectiveVelocity_h_out(Nt_out, Ny_out, Nx_out);
  Data<real, 3> AerodynamicResistance_out(Nt_out, Ny_out, Nx_out);

  Date current_date = date_beg;
  real friction, heat_flux, eta, eta0;
  for (h = 0; h < Nt_out; h++)
    {
      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            int month = current_date.GetMonth() - 1;
            friction = FrictionModule_h_out(h, j, i);
            heat_flux = PotentialHeatFlux_h_out(h, j, i) > 0 ?
              max(PotentialHeatFlux_h_out(h, j, i), 1.e-30)
              : min(PotentialHeatFlux_h_out(h, j, i), -1.e-30);
            MoninObukov_h_out(h, j, i) = -PotentialThermal_h_out(h, j, i)
              * friction * friction * friction / (.4 * 9.81 * heat_flux);
            if (MoninObukov_h_out(h, j, i) >= 0.)
              {
                AerodynamicResistance_out(h, j, i)
                  = log(.5 * Altitude_out(h, 1, j, i)
                        / RoughnessHeight(month, j, i))
                  + 4.7 * (.5 * Altitude_out(h, 1, j, i)
                           - RoughnessHeight(month, j, i))
                  / MoninObukov_h_out(h, j, i);
                AerodynamicResistance_out(h, j, i)
                  /= 0.4 * friction;
                ConvectiveVelocity_h_out(h, j, i) = 0.;
              }
            else
              {
                eta = pow(1. - 7.5 * Altitude_out(h, 1, j, i)
                          / MoninObukov_h_out(h, j, i), 0.25);
                eta0 = pow(1. - 15. * RoughnessHeight(month, j, i)
                           / MoninObukov_h_out(h, j, i), 0.25);
                AerodynamicResistance_out(h, j, i)
                  = (log(.5 * Altitude_out(h, 1, j, i)
                         / RoughnessHeight(month, j, i))
                     + log((eta0 * eta0 + 1.) * (eta0 + 1) * (eta0 + 1)
                           / ((eta * eta + 1.) * (eta + 1) * (eta + 1)))
                     + 2. * (atan(eta) - atan(eta0)))
                  / (.4 * friction);
                ConvectiveVelocity_h_out(h, j, i)
                  = pow(9.81 * max(PotentialHeatFlux_h_out(h, j, i), 1.e-6)
                        * BoundaryHeight_h_out(h, j, i)
                        / PotentialThermal_h_out(h, j, i), 1. / 3.);
              }
          }
      current_date.AddSeconds(Delta_t_out * 3600);
    }

  /*** Troen & Mahrt ***/

  real height_ratio, velocity_scale, convective, kz_min;
  real weight, coeff;
  real dz, grad, du, dv, dk;
  real richardson;
  real chi, temperature_mean, specific_humidity_mean;
  const real Lv = 2.45e6;
  const real Cp = 1005;
  const real Rv = 461.5;

  Kz_out.SetZero();

  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        for (k = 0; k < Nz_out; k++)
          {
            friction = FrictionModule_h_out(h, j, i);
            convective = ConvectiveVelocity_h_out(h, j, i);
            height_ratio = Altitude_out(h, k + 1, j, i)
              / BoundaryHeight_h_out(h, j, i);
            if (MoninObukov_h_out(h, j, i) > 0.)
              velocity_scale = friction / (1. + 4.7 * Altitude_out(h, k + 1, j, i)
                                           / MoninObukov_h_out(h, j, i));
            else
              velocity_scale = pow(friction * friction * friction
                                   + min(height_ratio, 0.1) * 2.8 * convective
                                   * convective * convective, 1. / 3.);
            if (height_ratio <= 1.)
              // Within boundary layer.
              {
                Kz_out(h, k + 1, j, i) = 0.4 * velocity_scale
                  * BoundaryHeight_h_out(h, j, i) * height_ratio
                  * (1. - height_ratio) * (1. - height_ratio);
                kz_min = Kz_min_dry + (Kz_min_wet - Kz_min_dry)
                  * (RelativeHumidityMax_h_out(h, j, i)
                     - relative_humidity_threshold)
                  / (1. - relative_humidity_threshold);
                Kz_out(h, k + 1, j, i) = max(Kz_out(h, k + 1, j, i), kz_min);
                if (k < Nz_out - 2
                    && Altitude_out(h, k + 2, j, i)
                    >= BoundaryHeight_h_out(h, j, i))
                  {
                    weight
                      = (BoundaryHeight_h_out(h, j, i)
                         - Altitude_out(h, k + 1, j, i))
                      / (Altitude_out(h, k + 2, j, i) - Altitude_out(h, k + 1, j, i));
                    Kz_out(h, k + 1, j, i) = weight * Kz_out(h, k + 1, j, i)
                      + (1. - weight) * Kz_min_up;
                  }
              }
            else
              {
                k_in = 1;
                while (k_in != Nz_in
                       && Altitude_h_out(k_in, j, i) < Altitude_out(h, k + 1, j, i))
                  k_in++;
                dz = Altitude_h_out(k_in, j, i)
                  - Altitude_h_out(k_in - 1, j, i);
                du = ZonalWind_h_out(h, k_in, j, i)
                  - ZonalWind_h_out(h, k_in - 1, j, i);
                dv = MeridionalWind_h_out(h, k_in, j, i)
                  - MeridionalWind_h_out(h, k_in - 1, j, i);
                grad = 1.e-6 + (du * du + dv * dv) / (dz * dz);
                richardson = 9.81
                  * (VirtualPotentialTemperature_h_out(h, k_in, j, i)
                     - VirtualPotentialTemperature_h_out(h, k_in - 1, j, i))
                  / (dz * grad) / PotentialTemperature_h_out(h, 0, j, i);
                if (RelativeHumidity_h_out(h, k_in, j, i)
                    > relative_humidity_threshold
                    || RelativeHumidity_h_out(h, k_in - 1, j, i)
                    > relative_humidity_threshold)
                  {
                    specific_humidity_mean
                      = 0.5 * (SpecificHumidity_h_out(h, k_in, j, i)
                               + SpecificHumidity_h_out(h, k_in - 1, j, i));
                    temperature_mean
                      = 0.5 * (Temperature_h_out(h, k_in, j, i)
                               + Temperature_h_out(h, k_in - 1, j, i));
                    coeff = Lv * specific_humidity_mean
                      / (R * temperature_mean);
                    chi  = Lv * Lv * specific_humidity_mean
                      / (Cp * Rv * temperature_mean * temperature_mean);
                    richardson = (1. + coeff)
                      * (richardson - 9.81 * 9.81 * (chi - coeff)
                         / (1. + chi) / grad / Cp / temperature_mean);
                  }
                dk = 1. / (0.4 * Altitude_out(h, k + 1, j, i)) + 1. / 150.;
                dk = sqrt(grad) / (dk * dk);

                if (richardson < 0.)
                  dk *= 1. - 8. * richardson
                    / (1. + 1.286 * sqrt(-richardson));
                else
                  {
                    coeff = 1. + 5. * richardson;
                    dk /= coeff * coeff;
                  }
                Kz_out(h, k + 1, j, i) = max(Kz_min_up, dk);
              }
          }
  Kz_out.ThresholdMax(Kz_max);

  cout << " done." << endl;


  /////////////////
  // ATTENUATION //
  /////////////////


  cout << "Computing attenuation...";
  cout.flush();

  int level_low, level_medium, level_high;
  real top_low(2500.), top_medium(6000.), top_high(20000.);
  real path_low, path_medium, path_high;
  real depth_low, depth_medium, depth_high;
  Data<real, 1> tmp(Nz_in);
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        {
          level_low = Nz_in;
          while (level_low != 0
                 && Altitude_h_out(level_low, j, i) >= top_low)
            level_low--;

          level_medium = Nz_in;
          while (level_medium != 0
                 && Altitude_h_out(level_medium, j, i) >= top_medium)
            level_medium--;

          level_high = Nz_in;
          while (level_high != 0
                 && Altitude_h_out(level_high, j, i) >= top_high)
            level_high--;

          if (level_low < Nz_in)
            path_low = LWPath_h_out(h, level_low, j, i)
              + (LWPath_h_out(h, level_low + 1, j, i)
                 - LWPath_h_out(h, level_low, j, i))
              * (top_low - Altitude_h_out(level_low, j, i))
              / (Altitude_h_out(level_low + 1, j, i)
                 - Altitude_h_out(level_low, j, i));
          else
            path_low = LWPath_h_out(h, level_low, j, i);

          if (level_medium < Nz_in)
            path_medium = LWPath_h_out(h, level_medium, j, i)
              + (LWPath_h_out(h, level_medium + 1, j, i)
                 - LWPath_h_out(h, level_medium, j, i))
              * (top_medium - Altitude_h_out(level_medium, j, i))
              / (Altitude_h_out(level_medium + 1, j, i)
                 - Altitude_h_out(level_medium, j, i));
          else
            path_medium = LWPath_h_out(h, level_medium, j, i);

          if (level_high < Nz_in)
            path_high = LWPath_h_out(h, level_high, j, i)
              + (LWPath_h_out(h, level_high + 1, j, i)
                 - LWPath_h_out(h, level_high, j, i))
              * (top_high - Altitude_h_out(level_high, j, i))
              / (Altitude_h_out(level_high + 1, j, i)
                 - Altitude_h_out(level_high, j, i));
          else
            path_high = LWPath_h_out(h, level_high, j, i);

          depth_medium = path_medium - path_low;
          depth_high = path_high - path_medium;

          tmp() = 0.;
          for (k = 1; k < Nz_in; k++)
            tmp(k) = tmp(k - 1)
              + max(0., 0.5 * (RelativeHumidity_h_out(h, k, j, i)
                               + RelativeHumidity_h_out(h, k - 1, j, i))
                    - 0.85)
              * 0.025 * (Altitude_h_out(k, j, i) - Altitude_h_out(k - 1, j, i))
              / (1. - .85);
          if (level_low < Nz_in)
            depth_low = tmp(level_low)
              + (tmp(level_low + 1) - tmp(level_low))
              * (top_low - Altitude_h_out(level_low, j, i))
              / (Altitude_h_out(level_low + 1, j, i)
                 - Altitude_h_out(level_low, j, i));
          else
            depth_low = tmp(level_low);

          // TODO: pow base can be negative whereas the exponent is real => invalid operation.
          Attenuation_out(h, j, i)
            = exp(-.11 * pow(depth_low + depth_medium
                             + depth_high, real(.67)));

        }

  cout << " done." << endl;


  ///////////////////////
  // VERTICAL AVERAGES //
  ///////////////////////


  cout << "Vertical averages...";
  cout.flush();

  // Wind averages along z.
  VerticalAverage(Altitude_h_out.GetArray(), MeridionalWind_h_out(),
                  Altitude_out(), MeridionalWind_out());
  VerticalAverage(Altitude_h_out.GetArray(), ZonalWind_h_out(),
                  Altitude_out(), ZonalWind_out());

  // Temperature and pressure average along z.
  VerticalAverage(Altitude_h_out.GetArray(), Pressure_h_out(),
                  Altitude_out(), Pressure_out());
  VerticalAverage(Altitude_h_out.GetArray(), Temperature_h_out(),
                  Altitude_out(), Temperature_out());

  // Other variables.
  VerticalAverage(Altitude_h_out.GetArray(), SpecificHumidity_h_out(),
                  Altitude_out(), SpecificHumidity_out());
  VerticalAverage(Altitude_h_out.GetArray(), LiquidWaterContent_h_out(),
                  Altitude_out(), LiquidWaterContent_out());
  VerticalAverage(Altitude_h_out.GetArray(), AirDensity_h_out(),
                  Altitude_out(), AirDensity_out());

  cout << " done." << endl;


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////


  FormatBinary<float> OutputMeteo;

  cout << "Writing data...";
  cout.flush();

  OutputMeteo.Append(Altitude_out, directory_out + "Altitude.bin");
  OutputMeteo.Append(Pressure_out, directory_out + "Pressure.bin");
  OutputMeteo.Append(Temperature_out, directory_out + "Temperature.bin");
  OutputMeteo.Append(SpecificHumidity_out,
                     directory_out + "SpecificHumidity.bin");
  OutputMeteo.Append(LiquidWaterContent_out,
                     directory_out + "LiquidWaterContent.bin");
  OutputMeteo.Append(ZonalWind_out, directory_out + "ZonalWind.bin");
  OutputMeteo.Append(MeridionalWind_out, directory_out
                     + "MeridionalWind.bin");
  OutputMeteo.Append(AirDensity_out, directory_out + "AirDensity.bin");
  OutputMeteo.Append(Kz_out, directory_out + "Kz.bin");
  OutputMeteo.Append(BoundaryHeight_h_out, directory_out + "PBLH.bin");
  Data<real, 3> RelativeHumidity_tmp(Nt_out, Ny_out, Nx_out);
  for (h = 0; h < Nt_out; h++)
    for (j = 0; j < Ny_out; j++)
      for (i = 0; i < Nx_out; i++)
        RelativeHumidity_tmp(h, j, i) = RelativeHumidity_h_out(h, 0, j, i);
  RelativeHumidity_tmp.Threshold(0., 1.);
  OutputMeteo.Append(RelativeHumidity_tmp,
                     directory_out + "SurfaceRelativeHumidity.bin");
  OutputMeteo.Append(Temperature_2m_h_out,
                     directory_out + "Temperature_2m.bin");
  OutputMeteo.Append(FrictionModule_h_out,
                     directory_out + "FrictionModule.bin");
  OutputMeteo.Append(AerodynamicResistance_out,
                     directory_out + "AerodynamicResistance.bin");
  OutputMeteo.Append(Attenuation_out, directory_out + "Attenuation.bin");
  OutputMeteo.Append(WindModule_10m_h_out,
                     directory_out + "WindModule_10m.bin");
  OutputMeteo.Append(ConvectiveVelocity_h_out,
                     directory_out + "ConvectiveVelocity.bin");
  OutputMeteo.Append(SoilMoisture_h_out, directory_out + "SoilMoisture.bin");

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
