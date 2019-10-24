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
#include <cstring>
#include<string>
#include<sstream>
#include <iostream>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////

#define _cls cls_

extern "C"
{
  void _cls(float*, float*, float*, float*, float*,
            float*, int*, float*, float*, float*, float*,
            float*, float*, float*, float*, float*,
            float*);
}


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("meteo.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, date_end, default_name);

  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  // Warning: 'real' must be 'float' because of Fortran functions
  // that are called.
  typedef float real;

  int h, i, j, k;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);

  // Output domain.
  int Nt, Nz, Ny, Nx;
  real Delta_t, Delta_y, Delta_x;
  real t_min, y_min, x_min;
  string vertical_levels;

  configuration.SetSection("[domain]");

  configuration.PeekValue("Nx", "> 0", Nx);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Delta_t", "> 0", Delta_t);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);
  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);
  configuration.PeekValue("Vertical_levels", vertical_levels);

  string date_str;
  configuration.PeekValue("Date", date_str);
  Date date_pre(date_str);

  double difference = date_beg.GetSecondsFrom(date_pre);
  if (difference < 0)
    throw string("\nThe date you provide in command line")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  t_min = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;


  // Paths.
  string LUC_file, roughness_file;
  int Nc, sea_index, urban_index;
  string dir_in, dir_out, file_in;

  configuration.SetSection("[paths]");

  configuration.PeekValue("LUC_file", LUC_file);
  if (!exists(LUC_file))
    throw "Unable to open land use cover file \"" + LUC_file + "\".";
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny * Nx);
  configuration.PeekValue("Sea_index", ">= 0 | < " + to_str(Nc), sea_index);
  configuration.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                          urban_index);

  configuration.PeekValue("Roughness_file", roughness_file);

  configuration.PeekValue("Directory_meteo", dir_in);
  configuration.PeekValue("File_Kz", file_in);
  configuration.PeekValue("Directory_Kz_TM", dir_out);

  // Vertical diffusion.
  real Kz_min, Kz_min_urban, Kz_max, perturbed_BL;
  bool apply_vert;

  configuration.SetSection("[Kz]");

  configuration.PeekValue("Min", "positive", Kz_min);
  configuration.PeekValue("Min_urban", "positive", Kz_min_urban);
  configuration.PeekValue("Max", "positive", Kz_max);
  configuration.PeekValue("Apply_vert", apply_vert);
  configuration.PeekValue("Perturbed_BL", "positive", perturbed_BL);

  // Troen & Mahrt coefficients.
  real p;
  real C;
  // Ratio between the SBL and the PBL.
  real SBL;
  // Critical Richardson number.
  real Ric;

  configuration.PeekValue("p", "positive", p);
  configuration.PeekValue("C", "positive", C);
  configuration.PeekValue("SBL", ">= 0 | <= 1", SBL);
  configuration.PeekValue("Ric", "positive", Ric);

  bool fluxes_diagnosed;
  int BL_diag;
  configuration.PeekValue("Fluxes_diagnosed", fluxes_diagnosed);
  configuration.PeekValue("BL_diag", "= 1 2 3", BL_diag);

  bool TM_stable;
  configuration.PeekValue("TM_stable", TM_stable);

  cout << " done." << endl;

  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  // Vertical levels depend on t, z, y and x.
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  // On in_interfaces.
  RegularGrid<real> GridZ_interf(Nz + 1);
  RegularGrid<real> GridY_interf(y_min - Delta_y / 2., Delta_y, Ny + 1);
  RegularGrid<real> GridX_interf(x_min - Delta_x / 2., Delta_x, Nx + 1);
  RegularGrid<real> GridC(Nc);

  // Reads vertical-heights.
  FormatText Heights;
  Heights.Read(vertical_levels, GridZ_interf);
  // Sets values at nodes.
  for (k = 0; k < Nz; k++)
    GridZ(k) = (GridZ_interf(k) + GridZ_interf(k + 1)) / 2.0;


  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 3> LUC(GridC, GridY, GridX);
  Data<real, 4> Temperature(GridT, GridZ, GridY, GridX);
  Data<real, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<real, 3> SkinTemperature(GridT, GridY, GridX);
  Data<real, 4> Pressure(GridT, GridZ, GridY, GridX);
  Data<real, 3> SurfacePressure(GridT, GridY, GridX);
  Data<real, 4> SpecificHumidity(GridT, GridZ, GridY, GridX);
  Data<real, 4> MeridionalWind(GridT, GridZ, GridY_interf, GridX);
  Data<real, 4> ZonalWind(GridT, GridZ, GridY, GridX_interf);
  Data<real, 4> WindModule(GridT, GridZ, GridY, GridX);
  Data<real, 3> FrictionModule(GridT, GridY, GridX);
  Data<real, 3> SoilWater(GridT, GridY, GridX);
  Data<real, 3> BoundaryHeight(GridT, GridY, GridX);
  Data<real, 2> Roughness(GridY, GridX);

  // Output fields.
  Data<real, 4> Richardson(GridT, GridZ, GridY, GridX);
  Data<real, 3> SaturationHumidity(GridT, GridY, GridX);
  Data<real, 3> SurfaceSpecificHumidity(GridT, GridY, GridX);
  Data<real, 3> SurfaceRelativeHumidity(GridT, GridY, GridX);
  Data<real, 3> LMO(GridT, GridY, GridX);
  Data<real, 3> Evaporation(GridT, GridY, GridX);
  Data<real, 3> SensibleHeat(GridT, GridY, GridX);
  Data<real, 3> Fm(GridT, GridY, GridX);
  Data<real, 3> Ft(GridT, GridY, GridX);
  Data<real, 4> Kz(GridT, GridZ_interf, GridY, GridX);
  Data<real, 4> PotentialTemperature(GridT, GridZ, GridY, GridX);
  Data<real, 3> SurfacePotentialTemperature(GridT, GridY, GridX);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  // Polair format.
  FormatBinary<float> InputMeteo;

  InputMeteo.Read(LUC_file, LUC);
  if (fluxes_diagnosed)
    InputMeteo.Read(roughness_file, Roughness);

  cout << "Extracting fields...";
  cout.flush();
  InputMeteo.ReadSteps(dir_in + "SurfaceTemperature.bin",
                       step, SurfaceTemperature);
  InputMeteo.ReadSteps(dir_in + "Temperature.bin", step, Temperature);
  InputMeteo.ReadSteps(dir_in + "SkinTemperature.bin",
                       step, SkinTemperature);
  InputMeteo.ReadSteps(dir_in + "Richardson.bin", step, Richardson);
  InputMeteo.ReadSteps(dir_in + "SoilWater.bin", step, SoilWater);
  InputMeteo.ReadSteps(dir_in + "SurfacePressure.bin",
                       step, SurfacePressure);
  InputMeteo.ReadSteps(dir_in + "Pressure.bin", step, Pressure);
  InputMeteo.ReadSteps(dir_in + "SpecificHumidity.bin", step,
                       SpecificHumidity);
  InputMeteo.ReadSteps(dir_in + "WindModule.bin", step, WindModule);
  InputMeteo.ReadSteps(dir_in + "MeridionalWind.bin",
                       step, MeridionalWind);
  InputMeteo.ReadSteps(dir_in + "ZonalWind.bin", step, ZonalWind);
  InputMeteo.ReadSteps(dir_in + "BoundaryHeight.bin",
                       step, BoundaryHeight);

  if (!fluxes_diagnosed)
    {
      InputMeteo.ReadSteps(dir_in + "FrictionModule.bin",
                           step, FrictionModule);
      InputMeteo.ReadSteps(dir_in + "Evaporation.bin", step, Evaporation);
      InputMeteo.ReadSteps(dir_in + "SensibleHeat.bin",
                           step, SensibleHeat);
    }

  InputMeteo.ReadSteps(file_in, step, Kz);

  cout << " done." << endl;
  cout << endl;


  ////////////////////////////
  // OUTPUT DATA PROCESSING //
  ////////////////////////////

  cout << "Computing Kz...";
  cout.flush();

  /*** First computes meteorological fields ***/

  ComputePotentialTemperature(Temperature, Pressure, PotentialTemperature);
  ComputePotentialTemperature(SkinTemperature, SurfacePressure,
                              SurfacePotentialTemperature);

  ComputeSaturationHumidity(SkinTemperature, SurfacePressure,
                            SaturationHumidity);

  ComputeSurfaceHumidity_diag(SpecificHumidity, SaturationHumidity, SoilWater,
                              LUC, sea_index, SurfaceSpecificHumidity);

  SurfaceRelativeHumidity.GetArray() = SurfaceSpecificHumidity.GetArray()
    / SaturationHumidity.GetArray();

  /*** Fluxes ***/

  // Computes the Monin-Obukhov length and the friction velocity.
  if (fluxes_diagnosed)
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            // Roughness.
            real zo = min(Roughness(j, i), 2.);
            real zot = zo / 10.;
            int is_sea;
            if (LUC(sea_index, j, i) > 0.5)
              is_sea = 1;
            else
              is_sea = 0;
            _cls(&ZonalWind(h, 0, j, i), &MeridionalWind(h, 0, j, i),
                 &PotentialTemperature(h, 0, j, i),
                 &SurfacePotentialTemperature(h, j, i),
                 &SpecificHumidity(h, 0, j, i),
                 &SurfaceSpecificHumidity(h, j, i),
                 &is_sea, &GridZ(0), &zo, &zot, &GridZ(0),
                 &LMO(h, j, i), &FrictionModule(h, j, i),
                 &SensibleHeat(h, j, i), &Evaporation(h, j, i),
                 &Fm(h, j, i), &Ft(h, j, i));
          }
  else
    ComputeLMO(FrictionModule, SurfacePotentialTemperature,
               PotentialTemperature, SensibleHeat, Evaporation,
               LMO);

  SensibleHeat.GetArray() += 0.608 * Evaporation.GetArray()
    * SurfaceTemperature.GetArray();

  /*** Boundary layer height ***/

  if (BL_diag == 1)
    /*** T&M diagnostic ***/
    ComputePBLH_TM(SurfaceTemperature, SurfacePotentialTemperature,
                   PotentialTemperature, FrictionModule, WindModule,
                   SensibleHeat, LMO, GridZ_interf, BoundaryHeight,
                   SBL, Ric, C);
  else if (BL_diag == 2)
    /*** Richardson diagnostic ***/
    ComputePBLH_Richardson(Richardson, GridZ_interf, BoundaryHeight, Ric);
  else
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            BoundaryHeight(h, j, i)
              = max(GridZ_interf(0), BoundaryHeight(h, j, i));
            BoundaryHeight(h, j, i)
              = min(GridZ_interf(Nz), BoundaryHeight(h, j, i));
          }

  BoundaryHeight() *= perturbed_BL;

  /*** Kz ***/

  ComputeTM_Kz(SurfaceTemperature, FrictionModule, SensibleHeat, LMO,
               BoundaryHeight, Kz, TM_stable, SBL, p);

  // Thresholds.
  real local_min;
  for (j = 0; j < Ny; j++)
    for (i = 0; i < Nx; i++)
      {
        local_min = Kz_min * (1. - LUC(urban_index, j, i))
          + Kz_min_urban * LUC(urban_index, j, i);
        if (apply_vert)
          {
            for (h = 0; h < Nt; h++)
              for (k = 0; k < Nz + 1; k++)
                if (Kz(h, k, j, i) < local_min)
                  Kz(h, k, j, i) = local_min;
          }
        else
          for (h = 0; h < Nt; h++)
            if (Kz(h, 1, j, i) < local_min)
              Kz(h, 1, j, i) = local_min;
      }

  Kz.ThresholdMax(Kz_max);

  cout << " done." << endl;


  ////////////////////////
  // WRITES OUTPUT DATA //
  ////////////////////////

  cout << "Writing output files...";
  cout.flush();

  FormatBinary<float> OutputMeteo;

  OutputMeteo.Append(Kz, dir_out + "Kz_TM.bin");

  OutputMeteo.Append(SurfaceRelativeHumidity,
                     dir_out + "SurfaceRelativeHumidity.bin");

  if (!fluxes_diagnosed)
    OutputMeteo.Append(LMO, dir_out + "LMO.bin");

  if (fluxes_diagnosed)
    {
      OutputMeteo.Append(Evaporation, dir_out + "Evaporation_diag.bin");
      OutputMeteo.Append(SensibleHeat, dir_out + "SensibleHeat_diag.bin");
      OutputMeteo.Append(FrictionModule, dir_out + "FrictionModule_diag.bin");
      OutputMeteo.Append(Fm, dir_out + "Fm.bin");
      OutputMeteo.Append(Ft, dir_out + "Ft.bin");
      OutputMeteo.Append(LMO, dir_out + "LMO_diag.bin");
    }

  if (BL_diag == 1 || BL_diag == 2)
    OutputMeteo.Append(BoundaryHeight, dir_out + "BoundaryHeight_diag.bin");

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
