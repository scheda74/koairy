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

// INCLUDES //
//////////////


template <class T>
void exp_(T& x)
{
  x = exp(x);
}


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("meteo.cfg");
  Date date_beg;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, default_name);

  Date date_end(date_beg);
  date_end.AddDays(1);

  Date date_prev(date_beg);
  date_prev.AddDays(-1);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  int h, i, j, k;

  // Constants.
  const real P0 = 101325.;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

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
  string directory_in, directory_out;

  configuration.SetSection("[paths]");

  configuration.PeekValue("Database_meteo", directory_in);
  configuration.PeekValue("Directory_attenuation", directory_out);

  string file_in = directory_in + date_beg.GetDate("ECMWF-%y%m%d.grb");
  string file_in_prev = directory_in + date_prev.GetDate("ECMWF-%y%m%d.grb");

  // Accumulated data
  int accumulated_time, accumulated_index;

  configuration.SetSection("[accumulated_data]");

  configuration.PeekValue("Accumulated_time", "positive", accumulated_time);
  configuration.PeekValue("Accumulated_index", "positive", accumulated_index);

  // Attenuation parameterization.
  int atte_type, option_critical_relative_humidity;

  configuration.SetSection("[photolysis_rates]");

  configuration.PeekValue("Attenuation_Type", "= 1 2",  atte_type);

  // Clouds.
  real min_height;

  configuration.SetSection("[clouds]");

  configuration.PeekValue("Min_height", "> 0", min_height);

  configuration.PeekValue("Critical_relative_humidity", "= 1 2",
                          option_critical_relative_humidity);

  cout << " done." << endl;


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

  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels, GridZ_interf_out);
  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_interf_out(k) + GridZ_interf_out(k + 1)) / 2.0;


  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 4> Temperature(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> SurfaceTemperature(GridT_in, GridY_in, GridX_in);
  Data<real, 4> Pressure(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Pressure_interf(GridT_in, GridZ_interf_in,
                                GridY_in, GridX_in);
  Data<real, 3> SurfacePressure(GridT_in, GridY_in, GridX_in);
  Data<real, 4> SpecificHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> RelativeHumidity(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> LiquidWaterContent(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> LowCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> LowCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> MediumCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness(GridT_in, GridY_in, GridX_in);
  Data<real, 3> HighCloudiness_diag(GridT_in, GridY_in, GridX_in);
  Data<real, 3> BoundaryLayerHeight(GridT_in, GridY_in, GridX_in);
  Data<real, 4> CloudFraction(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> CRH(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> Attenuation(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 3> CloudBaseHeight(GridT_in, GridY_in, GridX_in);
  Data<real, 3> CloudTopHeight(GridT_in, GridY_in, GridX_in);

  Data<real, 3> ConvectiveRain(GridT_in, GridY_in, GridX_in);
  Data<real, 3> ConvectiveRain_prev(GridT_in, GridY_in, GridX_in);
  Data<real, 3> LargeScaleRain(GridT_in, GridY_in, GridX_in);
  Data<real, 3> LargeScaleRain_prev(GridT_in, GridY_in, GridX_in);

  // Output fields.
  Data<real, 4> Attenuation_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 3> CloudBaseHeight_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> CloudTopHeight_out(GridT_out, GridY_out, GridX_out);
  Data<real, 4> CloudFraction_out(GridT_out, GridZ_out, GridY_out, GridX_out);

  Data<real, 3> ConvectiveRain_out(GridT_out, GridY_out, GridX_out);
  Data<real, 3> LargeScaleRain_out(GridT_out, GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  FormatGrib InputMeteo;

  cout << "Extracting data...";
  cout.flush();

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

  InputMeteo.Read(file_in, 246, LiquidWaterContent);
  LiquidWaterContent.ReverseData(1);
  LiquidWaterContent.ReverseData(2);

  InputMeteo.Read(file_in, 186, LowCloudiness);
  LowCloudiness.ReverseData(1);

  InputMeteo.Read(file_in, 187, MediumCloudiness);
  MediumCloudiness.ReverseData(1);

  InputMeteo.Read(file_in, 188, HighCloudiness);
  HighCloudiness.ReverseData(1);

  InputMeteo.Read(file_in, 143, ConvectiveRain);
  ConvectiveRain.ReverseData(1);

  InputMeteo.Read(file_in_prev, 143, ConvectiveRain_prev);
  ConvectiveRain_prev.ReverseData(1);

  InputMeteo.Read(file_in, 142, LargeScaleRain);
  LargeScaleRain.ReverseData(1);

  InputMeteo.Read(file_in_prev, 142, LargeScaleRain_prev);
  LargeScaleRain_prev.ReverseData(1);

  InputMeteo.Read(file_in, 159, BoundaryLayerHeight);
  BoundaryLayerHeight.ReverseData(1);

  // Transformations.

  SurfacePressure.Apply(exp_);

  SpecificHumidity.Threshold(0., 1.);
  ConvectiveRain.ThresholdMin(0.);
  ConvectiveRain_prev.ThresholdMin(0.);
  LargeScaleRain.ThresholdMin(0.);
  LargeScaleRain_prev.ThresholdMin(0.);

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

  cout << " done." << endl;


  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////

  // Computes level heights with pressure levels.
  cout << "Computing level heights in meter...";
  cout.flush();

  ExtStream file_hybrid_coefficients(directory_in
                                     + "hybrid_coefficients.dat");
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
  ComputeMiddleHeight(Pressure_interf, Temperature,
                      GridZ_interf_in, GridZ_in);

  cout << " done." << endl;


  /////////////////
  // ATTENUATION //
  /////////////////

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

  ComputeCloudFraction(BoundaryLayerHeight, RelativeHumidity, CRH,
                       CloudFraction);

  Data<int, 4> LowIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> MediumIndices(Nt_in, Ny_in, Nx_in, 2);
  Data<int, 4> HighIndices(Nt_in, Ny_in, Nx_in, 2);

  ComputeCloudiness(CloudFraction, Pressure, GridZ_interf_in, LowIndices,
                    MediumIndices, HighIndices, LowCloudiness_diag,
                    MediumCloudiness_diag, HighCloudiness_diag);
  ComputeCloudBaseHeight(LowIndices, MediumIndices, HighIndices,
                         GridZ_interf_in, CloudBaseHeight);
  ComputeCloudTopHeight(LowIndices, MediumIndices, HighIndices,
                        GridZ_interf_in, CloudTopHeight);

  cout << " done." << endl;

  cout << "Computing attenuation...";
  cout.flush();

  if (atte_type == 1)
    ComputeAttenuation_LWC(LiquidWaterContent,
                           LowIndices, MediumIndices, HighIndices,
                           MediumCloudiness, HighCloudiness,
                           date_beg, Delta_t_out, Attenuation);
  else if (atte_type == 2)
    ComputeAttenuation_ESQUIF(MediumCloudiness, HighCloudiness,
                              RelativeHumidity, Attenuation);

  cout << " done." << endl;
  cout << endl;


  ///////////////////////////
  // LINEAR INTERPOLATIONS //
  ///////////////////////////

  cout << "Linear interpolations...";
  cout.flush();
  LinearInterpolationOneGeneral(Attenuation, Attenuation_out, 1);
  Attenuation_out.Threshold(0., 2.);
  LinearInterpolationRegular(CloudBaseHeight, CloudBaseHeight_out);
  CloudBaseHeight_out.ThresholdMin(min_height);
  LinearInterpolationRegular(CloudTopHeight, CloudTopHeight_out);
  CloudTopHeight_out.ThresholdMin(min_height);

  LinearInterpolationRegular(ConvectiveRain, ConvectiveRain_out);
  LinearInterpolationRegular(LargeScaleRain, LargeScaleRain_out);

  ConvectiveRain_out.ThresholdMin(0.);
  LargeScaleRain_out.ThresholdMin(0.);
  cout << " done." << endl;

  Data<real, 4> RelativeHumidity_out(GridT_out, GridZ_out, GridY_out, GridX_out);
  LinearInterpolationOneGeneral(RelativeHumidity, RelativeHumidity_out, 1);
  LinearInterpolationOneGeneral(CloudFraction, CloudFraction_out, 1);


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  FormatBinary<float> OutputMeteo;

  cout << "Writing data...";
  cout.flush();
  OutputMeteo.Append(Attenuation_out, directory_out + "Attenuation_ref.bin");
  OutputMeteo.Append(RelativeHumidity_out, directory_out + "RelativeHumidity_ref.bin");
  OutputMeteo.Append(CloudFraction_out, directory_out + "CloudFraction_ref.bin");
  OutputMeteo.Append(CloudBaseHeight_out, directory_out + "CloudBaseHeight.bin");
  OutputMeteo.Append(CloudTopHeight_out,
                     directory_out + "CloudTopHeight.bin");
  OutputMeteo.Append(ConvectiveRain_out,
                     directory_out + "ConvectiveRain.bin");
  LargeScaleRain_out.GetArray() = LargeScaleRain_out.GetArray()
    + ConvectiveRain_out.GetArray();
  OutputMeteo.Append(LargeScaleRain_out, directory_out + "Rain.bin");
  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
