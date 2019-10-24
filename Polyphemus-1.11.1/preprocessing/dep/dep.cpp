// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Chi-Sian Soh and Vivien Mallet
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
#include <limits>
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "SeldonData.hxx"
using namespace SeldonData;

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////


template<class T>
void Cut(T& x)
{
  x = max(T(0), x);
  x = min(T(1e+30), x);
}

// Sets infinite resistance values to 1e+30.
template<class T>
void Infinity(T& x)
{
  if (x == T(9999))
    x = T(1e+30);
}


int main(int argc, char** argv)
{
  TRY;

  cout << endl;

  string configuration_file, sec_config_file, default_name("dep.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, date_end, default_name);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;

  // Units: SI.
  const real pi = 3.14159265358979323846264;
  const real Karman = 0.4;
  const real nu = 0.15;
  const real D_H2O = 0.25;
  const real Pr = 0.74;
  const real M_vap = 18.015e-3;
  const real M_air = 28.97e-3;

  int h, i, j, k, l, s;

  // Input and output directories.
  string dir_in, dir_out;

  // Polair Inputs
  string SurfaceTemperature_file, SurfaceRichardson_file, SolarRadiation_file,
    WindModule_file, U_file, V_file, PARdiff_file, PARdir_file,
    SpecificHumidity_file, SurfacePressure_file, FrictionVelocity_file,
    CanopyWetness_file, Rain_file, RoughnessHeight_file;

  // Output domain.
  int Nt, Nz, Ny, Nx;
  // Output steps.
  real Delta_t, Delta_y, Delta_x;
  real y_min, x_min;
  // Output altitudes.
  string levels;

  // Type of landuse classification.
  string landuse_type;
  // Landuse file.
  string LUC_file;
  // Number of LUC categories.
  int Nc;
  // Land Data Files
  string midsummer, autumn, late_autumn, snow, spring;

  // Number of species.
  int Ns_out;
  // Species data.
  string species_data;

  // Options.
  string option_Ra, option_Rb, option_Rc;
  bool option_save_resistances, option_Roughness;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files...";
  cout.flush();

  ConfigStreams configuration(configuration_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);


  configuration.SetSection("[domain]");

  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nx", "> 0", Nx);

  configuration.PeekValue("Delta_t", "> 0", Delta_t);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);

  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);
  configuration.PeekValue("Vertical_levels", levels);


  // Dates.
  string date_meteo_str;
  configuration.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw string("\nthe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  int days = date_beg.GetDaysFrom(date_meteo);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  real t_min = date_beg.GetHour() + real(date_beg.GetMinutes()) / 60.
    + real(date_beg.GetSeconds()) / 3600.;

  // Input/output files.
  configuration.SetSection("[paths]");

  configuration.PeekValue("Directory_dep", dir_out);
  configuration.PeekValue("SurfaceTemperature", SurfaceTemperature_file);
  configuration.PeekValue("SurfaceRichardson", SurfaceRichardson_file);
  configuration.PeekValue("SolarRadiation", SolarRadiation_file);

  // An extra step available for accumulated data?
  bool extra_step = (int(file_size(Rain_file) / (Nt * Nx * Ny
                                                 * sizeof(float)))
                     > days + 1);

  configuration.PeekValue("WindModule", WindModule_file);
  configuration.PeekValue("PARdiff", PARdiff_file);
  configuration.PeekValue("PARdir", PARdir_file);
  configuration.PeekValue("SpecificHumidity", SpecificHumidity_file);
  configuration.PeekValue("SurfacePressure", SurfacePressure_file);
  configuration.PeekValue("FrictionVelocity", FrictionVelocity_file);
  configuration.PeekValue("CanopyWetness", CanopyWetness_file);
  configuration.PeekValue("Rain", Rain_file);
  configuration.PeekValue("RoughnessHeight", RoughnessHeight_file);

  // LUC.
  configuration.PeekValue("Type", landuse_type);

  // Species.
  configuration.PeekValue("Data", species_data);
  configuration.SetSection("[Species]");
  configuration.PeekValue("Ns", "positive", Ns_out);

  // Options.
  configuration.SetSection("[Options]");

  configuration.PeekValue("CellRoughness", option_Roughness);
  configuration.PeekValue("Ra", "fh | fm | diag", option_Ra);
  configuration.PeekValue("Rb", "friction | diag", option_Rb);
  configuration.PeekValue("Rc", "zhang | wesely", option_Rc);
  configuration.PeekValue("Save_resistances", option_save_resistances);

  // Reads the configuration file associated with the landuse type.
  ConfigStreams classification(landuse_type);
  classification.AddFile(configuration_file);
  classification.AddFile(sec_config_file);

  classification.PeekValue("File", LUC_file);
  Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny * Nx);

  classification.PeekValue("Midsummer", midsummer);
  classification.PeekValue("Autumn", autumn);
  classification.PeekValue("Late_autumn", late_autumn);
  classification.PeekValue("Snow", snow);
  classification.PeekValue("Spring", spring);

  cout << " done." << endl;
  cout << endl;


  ///////////
  // GRIDS //
  ///////////

  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input grids.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridT_cumulated(t_min - Delta_t / 2.0,
                                    Delta_t,
                                    extra_step ? Nt + 1 : Nt);
  RegularGrid<real> GridZ(Nz);
  RegularGrid<real> GridZ_interf(Nz + 1);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridY_interf(y_min - Delta_y / 2.,
                                 Delta_y, Ny + 1);
  RegularGrid<real> GridX_interf(x_min - Delta_x / 2.,
                                 Delta_x, Nx + 1);
  RegularGrid<real> GridC(Nc);


  // Species.
  RegularGrid<real> GridS(Ns_out);

  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(levels.c_str(), GridZ_interf);
  // Sets values at nodes.
  for (k = 0; k < Nz; k++)
    GridZ(k) = (GridZ_interf(k) + GridZ_interf(k + 1)) / 2.0;


  //////////
  // DATA //
  //////////

  // Input fields.
  Data<real, 3> LandUse(GridC, GridY, GridX);
  Data<real, 2> RoughnessHeight_cell(GridY, GridX);

  Data<real, 3> SurfaceRichardson(GridT, GridY, GridX);
  Data<real, 4> WindModule(GridT, GridZ, GridY, GridX);
  Data<real, 3> SolarRadiation(GridT_cumulated, GridY, GridX);
  Data<real, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<real, 3> PARdiff(GridT_cumulated, GridY, GridX);
  Data<real, 3> PARdir(GridT_cumulated, GridY, GridX);
  Data<real, 4> SpecificHumidity(GridT, GridZ, GridY, GridX);
  Data<real, 3> SurfacePressure(GridT, GridY, GridX);
  Data<real, 3> FrictionVelocity(GridT, GridY, GridX);
  Data<real, 3> CanopyWetness(GridT, GridY, GridX);
  Data<real, 3> Rain(GridT_cumulated, GridY, GridX);

  Data<string, 1, real> SpeciesNames(GridS);
  Data<real, 1> MolarMass(GridS);
  Data<real, 1> Henry(GridS);
  Data<real, 1> Diffusivity(GridS);
  Data<real, 1> Reactivity(GridS);
  Data<real, 1> Alpha(GridS);
  Data<real, 1> Beta(GridS);
  Data<real, 1> Rm(GridS);

  Data<real, 1> LeafLength(GridC);
  Data<real, 1> rst_min(GridC);
  Data<real, 1> brs(GridC);
  Data<real, 1> T_min(GridC);
  Data<real, 1> T_max(GridC);
  Data<real, 1> T_opt(GridC);
  Data<real, 1> b_vpd(GridC);
  Data<real, 1> psi_c1(GridC);
  Data<real, 1> psi_c2(GridC);
  Data<real, 1> H(GridC);
  Data<real, 1> LAI(GridC);
  Data<real, 1> RoughnessHeight(GridC);
  Data<real, 1> ZeroDisplacementHeight(GridC);
  Data<real, 1> Rg_SO2(GridC);
  Data<real, 1> Rg_O3(GridC);
  Data<real, 1> Rc_srf_SO2(GridC);
  Data<real, 1> Rc_srf_O3(GridC);
  Data<real, 1> R_ac0(GridC);
  Data<real, 1> Ri(GridC);
  Data<real, 1> Rlu(GridC);
  Data<real, 1> WithVegetation(GridC);
  Data<real, 1> TallVegetation(GridC);

  // Output fields.
  Data<real, 3> Velocities(GridT, GridY, GridX);

  Data<real, 4> DRa(GridT, GridY, GridX, GridC);
  Data<real, 4> DRb(GridT, GridY, GridX, GridC);
  Data<real, 4> DRc(GridT, GridY, GridX, GridC);

  Data<real, 3> SurfaceRelativeHumidity(GridT, GridY, GridX);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////

  cout << "Extracting input data...";
  cout.flush();

  FormatText InputParam;

  // Species data.
  ifstream SpeciesData(species_data.c_str());
  for (i = 0; i < Ns_out; i++)
    SpeciesData >> SpeciesNames(i);
  InputParam.Read(SpeciesData, MolarMass);
  InputParam.Read(SpeciesData, Henry);
  InputParam.Read(SpeciesData, Diffusivity);
  InputParam.Read(SpeciesData, Reactivity);
  InputParam.Read(SpeciesData, Alpha);
  InputParam.Read(SpeciesData, Beta);
  InputParam.Read(SpeciesData, Rm);
  SpeciesData.close();

  // Land use resistances and roughness heights.
  ifstream LandData;
  int MM = date_beg.GetMonth();
  if (MM == 5 || MM == 6 || MM == 7 || MM == 8)
    LandData.open(midsummer.c_str());
  else if (MM == 9 || MM == 10)
    LandData.open(autumn.c_str());
  else if (MM == 11 || MM == 12 || MM == 1 || MM == 2)
    LandData.open(late_autumn.c_str());
  else if (MM == 3 || MM == 4)
    LandData.open(spring.c_str());
  else
    LandData.open(snow.c_str());

  InputParam.Read(LandData, LeafLength);
  InputParam.Read(LandData, rst_min);
  InputParam.Read(LandData, brs);
  InputParam.Read(LandData, T_min);
  InputParam.Read(LandData, T_max);
  InputParam.Read(LandData, T_opt);
  InputParam.Read(LandData, b_vpd);
  InputParam.Read(LandData, psi_c1);
  InputParam.Read(LandData, psi_c2);
  InputParam.Read(LandData, H);
  InputParam.Read(LandData, LAI);
  InputParam.Read(LandData, RoughnessHeight);
  InputParam.Read(LandData, ZeroDisplacementHeight);
  InputParam.Read(LandData, Rg_SO2);
  InputParam.Read(LandData, Rg_O3);
  InputParam.Read(LandData, Rc_srf_SO2);
  InputParam.Read(LandData, Rc_srf_O3);
  InputParam.Read(LandData, R_ac0);
  InputParam.Read(LandData, Ri);
  InputParam.Read(LandData, Rlu);
  InputParam.Read(LandData, WithVegetation);
  InputParam.Read(LandData, TallVegetation);
  LandData.close();

  LeafLength.Apply(Infinity<real>);
  rst_min.Apply(Infinity<real>);
  brs.Apply(Infinity<real>);
  T_min.Apply(Infinity<real>);
  T_max.Apply(Infinity<real>);
  T_opt.Apply(Infinity<real>);
  b_vpd.Apply(Infinity<real>);
  psi_c1.Apply(Infinity<real>);
  psi_c2.Apply(Infinity<real>);
  H.Apply(Infinity<real>);
  LAI.Apply(Infinity<real>);
  Rg_SO2.Apply(Infinity<real>);
  Rg_O3.Apply(Infinity<real>);
  Rc_srf_SO2.Apply(Infinity<real>);
  Rc_srf_O3.Apply(Infinity<real>);
  R_ac0.Apply(Infinity<real>);
  Ri.Apply(Infinity<real>);
  Rlu.Apply(Infinity<real>);

  // Polair inputs.
  FormatBinary<float> InputPolair;

  InputPolair.Read(LUC_file, LandUse);
  InputPolair.Read(RoughnessHeight_file, RoughnessHeight_cell);

  InputPolair.ReadSteps(SurfaceTemperature_file, step, SurfaceTemperature);
  InputPolair.ReadSteps(SurfaceRichardson_file, step, SurfaceRichardson);
  InputPolair.ReadSteps(SolarRadiation_file, step, SolarRadiation);
  InputPolair.ReadSteps(WindModule_file, step, WindModule);
  InputPolair.ReadSteps(PARdiff_file, step, PARdiff);
  InputPolair.ReadSteps(PARdir_file, step, PARdir);
  InputPolair.ReadSteps(SpecificHumidity_file, step, SpecificHumidity);
  InputPolair.ReadSteps(SurfacePressure_file, step, SurfacePressure);
  InputPolair.ReadSteps(FrictionVelocity_file, step, FrictionVelocity);
  InputPolair.ReadSteps(CanopyWetness_file, step, CanopyWetness);
  InputPolair.ReadSteps(Rain_file, step, Rain);

  cout << " done." << endl;
  cout << endl;


  /////////////////////////
  // METEOROLOGICAL DATA //
  /////////////////////////

  cout << "Computing meteorological fields...";
  cout.flush();

  SurfaceTemperature.Add(-273.15);

  // Calculation of relative humidity.
  real P_sat_vap, P_vap;
  for (h = 0; h < Nt; h++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        {
          P_sat_vap = 611.2 * exp(17.67 * SurfaceTemperature(h, j, i)
                                  / (SurfaceTemperature(h, j, i) + 243.50));
          P_vap = SpecificHumidity(h, 0, j, i) * SurfacePressure(h, j, i)
            / (SpecificHumidity(h, 0, j, i) +
               (1 - SpecificHumidity(h, 0, j, i)) * (M_vap / M_air));
          SurfaceRelativeHumidity(h, j, i) = P_vap / P_sat_vap;
        }
  SurfaceRelativeHumidity.Threshold(0.0, 1.0);

  CanopyWetness.SetZero();

  // PAR to W.m^-2.
  PARdiff.Mlt(1. / 4.6);
  PARdir.Mlt(1. / 4.6);

  cout << " done." << endl;
  cout << endl;


  ///////////////////////////
  // DEPOSITION VELOCITIES //
  ///////////////////////////

  cout << "Computing deposition velocities...";
  cout.flush();
  cout << endl;

  real Rsx, Rmx, Rlux, Rdc, Rclx, Rgx;
  real r0, r1, wind(0.), friction_wind, tmp;
  real Gs, f_T, f_D, f_psi, Wst(0.), bt, D, e_sat, e_amb, psi;
  real lon, lat, ut, alpha, theta;
  real PAR_shade, PAR_sun, F_sun, F_shade, R_dir, R_diff;
  real Rcut, Rg, Rst, Rns, Rac;
  real stab, rich, Cd(0.), Ra, Rb, Rc(0.);
  string surface_condition;
  real Sc;
  real zszo, zszot, zot, al_zszo, al_zszot, a2, alt, alu, c, drib;
  real fm(0.), fh(0.);
  real b1 = 5., c1 = 5., d1 = 5.;
  real roughness_height;


  for (s = 0; s < Ns_out; s++)
    {
      Date current_date = date_beg;
      Velocities.SetZero();
      Sc = nu / Diffusivity(s);
      for (h = 0; h < Nt; h++)
        {
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {

                // Surface condition.
                if (CanopyWetness(h, j, i) <= 0.1 &&
                    SurfaceRelativeHumidity(h, j, i) < 0.8)
                  surface_condition = "dry";
                else if (CanopyWetness(h, j, i) <= 0.1 &&
                         SurfaceRelativeHumidity(h, j, i) >= 0.9)
                  surface_condition = "high humidity";
                else if (CanopyWetness(h, j, i) > 0.8)
                  surface_condition = "dew";
                else if (CanopyWetness(h, j, i) > 0.8 && Rain(h, j, i) > 0.01)
                  surface_condition = "rain";
                else
                  surface_condition = "dry";

                // Loop for the different LUCs.
                for (l = 0; l < Nc; l++)
                  {
                    if (option_Roughness)
                      roughness_height = RoughnessHeight_cell(j, i);
                    else
                      roughness_height = RoughnessHeight(l);


                    /////////////
                    // Ra & Rb //
                    /////////////

                    if (option_Ra == "diag" || option_Rb == "diag")
                      {
                        // Stability function.
                        r0 = log(25. / roughness_height);
                        r0 = Karman * Karman / (r0 * r0);
                        r1 = 9.4 * 7.4 * r0 * sqrt(25. / roughness_height);
                        rich = SurfaceRichardson(h, j, i);
                        if (rich < 0.0)
                          stab = 1. - 9.4 * rich / (1. + r1 * sqrt(-rich));
                        else
                          stab = 1. / (pow(1. + 4.7 * rich, 2.0));
                        // Drag coefficient Cd.
                        Cd = r0 * stab;
                        // Wind module.
                        wind = max(WindModule(h, 0, j, i), 0.1);
                      }


                    ////////
                    // Ra //
                    ////////

                    if (option_Ra != "diag")
                      {
                        zszo = (GridZ(0) + roughness_height)
                          / roughness_height;
                        zot = roughness_height / 10;
                        zszot = (GridZ(0) + zot) / zot;
                        al_zszo = log(zszo);
                        al_zszot = log(zszot);
                        alu = Karman / al_zszo;
                        alt = Karman / al_zszot;

                        if (SurfaceRichardson(h, j, i) > 0)
                          {
                            drib = sqrt(1. + d1 * SurfaceRichardson(h, j, i));
                            if (option_Ra == "fh")
                              fh = 1. / (1. + 3. * b1
                                         * SurfaceRichardson(h, j, i) * drib);
                            else
                              fm = 1. / (1. + 2. * b1
                                         * SurfaceRichardson(h, j, i) / drib);
                          }
                        else
                          {
                            c = alu * alt * b1 * c1 * sqrt(1. - 1. / zszot)
                              * pow(double(pow(double(zszot), double(1. / 3.))
                                           - 1.), double(3. / 2.));
                            drib = 1. + 3. * c
                              * sqrt(-SurfaceRichardson(h, j, i));
                            if (option_Ra == "fh")
                              fh = 1. - 3. * b1 * SurfaceRichardson(h, j, i)
                                / drib;
                            else
                              fm = 1. - 2. * b1 * SurfaceRichardson(h, j, i)
                                / drib;
                          }

                        if (option_Ra == "fh")
                          Ra = 1. / (alu * alt * WindModule(h, 0, j, i) * fh);
                        else
                          {
                            a2 = (Karman * Karman) / (al_zszo * al_zszo);
                            Ra = 1. / (a2 * WindModule(h, 0, j, i) * fm);
                          }
                      }
                    else // old parameterization.
                      // Aerodynamic resistance.
                      Ra = 1. / (wind * Cd);


                    ////////
                    // Rb //
                    ////////

                    if (option_Rb == "friction")
                      {
                        if (WithVegetation(l) == 1)
                          Rb = 2. / (Karman * FrictionVelocity(h, j, i))
                            * pow(double(Sc / Pr), double(2. / 3.));
                        else
                          Rb = 1. / (Karman * FrictionVelocity(h, j, i))
                            * pow(double(Sc / Pr), double(2. / 3.));
                      }
                    else
                      {
                        // Friction wind.
                        friction_wind = wind * sqrt(Cd);
                        //
                        tmp = sqrt(MolarMass(s) / 18.);
                        tmp = 2. * float(pow(double((nu / (D_H2O * Pr)) * tmp),
                                             double(2. / 3.))) / Karman;
                        // Quasilaminary sublayer resistance.
                        Rb = tmp / friction_wind;
                      }


                    ////////
                    // Rc //
                    ////////

                    if (option_Rc == "zhang")
                      {
                        /*** Fraction of stomatal blocking, Wst ***/
                        if (surface_condition == "dry" ||
                            surface_condition == "high humidity")
                          Wst = 0;
                        else if (SolarRadiation(h, j, i) <= 200)
                          Wst = 0.;
                        else if (SolarRadiation(h, j, i) > 200 &&
                                 SolarRadiation(h, j, i) <= 600)
                          Wst = (SolarRadiation(h, j, i) - 200.) / 800.;
                        else if (SolarRadiation(h, j, i) > 600)
                          Wst = 0.5;


                        /*** Stomatal Resistance, Rst ***/

                        if (rst_min(l) == real(1e+30) || brs(l) == real(1e+30)
                            || T_min(l) == real(1e+30) || T_max(l) == real(1e+30)
                            || T_opt(l) == real(1e+30) || b_vpd(l) == real(1e+30)
                            || psi_c1(l) == real(1e+30) || psi_c2(l) == real(1e+30)
                            || LAI(l) == real(1e+30))
                          Rst = 1e+30;
                        else if (SolarRadiation(h, j, i) == 0)
                          Rst = 1e+30;
                        else if (SurfaceTemperature(h, j, i) < T_min(l)
                                 || SurfaceTemperature(h, j, i) > T_max(l))
                          Rst = 1e+30;
                        else
                          {
                            // Unstressed canopy stomatal conductance.
                            lon = x_min + i * Delta_x;
                            lat = y_min + j * Delta_y;
                            ut = real(current_date.GetHour() +
                                      current_date.GetMinutes() / 60. +
                                      current_date.GetSeconds() / 3600.);
                            theta = ZenithAngle(lon, lat,
                                                current_date.GetDate(), ut)
                              / 180. * pi;
                            alpha = pi / 3.;

                            if (theta > pi / 2.)
                              {
                                F_sun = 0.;
                                F_shade = 0.;
                              }
                            else
                              {
                                F_sun = 2. * cos(theta)
                                  * (1. -  exp(-0.5 * LAI(l) / cos(theta)));
                                F_shade = LAI(l) - F_sun;
                              }

                            // R_dir and R_diff.
                            R_diff = PARdiff(h, j, i);
                            R_dir = PARdir(h, j, i);

                            if (LAI(l) < 2.5 || SolarRadiation(h, j, i) < 200.)
                              {
                                PAR_shade = R_diff
                                  * exp(-0.5 * pow(double(LAI(l)), 0.7))
                                  + 0.07 * R_dir * (1.1 - 0.1 * LAI(l))
                                  * exp(-cos(theta));
                                PAR_sun = R_dir * cos(alpha) / cos(theta)
                                  + PAR_shade;
                              }
                            else
                              {
                                PAR_shade = R_diff
                                  * exp(-0.5 * pow(double(LAI(l)), 0.8))
                                  + 0.07 * R_dir * (1.1 - 0.1 * LAI(l))
                                  * exp(-cos(theta));
                                PAR_sun = pow(double(R_dir), 0.8) * cos(alpha)
                                  / cos(theta) + PAR_shade;
                              }

                            // Computing Gs.
                            {
                              real Gs_sun = (F_sun * PAR_sun) / (rst_min(l) * (PAR_sun + brs(l)));
                              real Gs_shade = (F_shade * PAR_shade) / (rst_min(l) * (PAR_shade + brs(l)));
                              Gs = Gs_sun + Gs_shade;
                              if (Gs == 0.)
                                Gs = numeric_limits<real>::epsilon();
                            }

                            // Conductance-reducing effects of air temperature.
                            bt = (T_max(l) - T_opt(l)) / (T_opt(l) - T_min(l));
                            f_T = (SurfaceTemperature(h, j, i) - T_min(l)) /
                                         (T_opt(l) - T_min(l))
                                         * pow((T_max(l) - SurfaceTemperature(h, j, i)) /
                                               (T_max(l) - T_opt(l)), bt);

                            // Conductance-reducing effects
                            // of water vapour deficit.
                            e_amb = SpecificHumidity(h, 0, j, i)
                                         * SurfacePressure(h, j, i)
                                         / (0.622 * (1. - SpecificHumidity(h, 0, j, i))
                                            + SpecificHumidity(h, 0, j, i));
                            e_sat = 610.78 * exp(17.2694
                                                 * SurfaceTemperature(h, j, i)
                                                 / (SurfaceTemperature(h, j, i)
                                                    + 237.29));
                            D = (e_sat - e_amb) / 1000.;
                            f_D = 1. - b_vpd(l) * D;

                            // Conductance-reducing effects of water stress.
                            if (TallVegetation(l) == 1)
                              psi = (-0.72 - 0.0013 * SolarRadiation(h, j, i))
                                         * 102.;
                            else
                              psi = (-0.395 - 0.043
                                     * SurfaceTemperature(h, j, i)) * 102.;

                            if (psi > psi_c1(l))
                              f_psi = 1.;
                            else
                              f_psi = (psi - psi_c2(l))
                                / (psi_c1(l) - psi_c2(l));

                            Rst = 1. / (Gs * f_T * f_D * f_psi
                                        * Diffusivity(s) / D_H2O);
                          }
                        Cut(Rst);


                        /*** Ground resistance Rg ***/

                        if (SpeciesNames(s) == "O3")
                          Rg = Rg_O3(l);
                        else if (SpeciesNames(s) == "SO2")
                          Rg = Rg_SO2(l);
                        else
                          {
                            if (Alpha(s) == 0. && Beta(s) == 0.)
                              Rg = numeric_limits<real>::infinity();
                            else
                              Rg = 1. / (Alpha(s) / Rg_SO2(l) + Beta(s) / Rg_O3(l));
                          }
                        if (SurfaceTemperature(h, j, i) < -1.)
                          Rg = Rg * exp(0.2
                                        * (-1. - SurfaceTemperature(h, j, i)));
                        Cut(Rg);


                        /*** Local leaf cuticular resistance Rcut ***/

                        if (SpeciesNames(s) == "O3")
                          Rcut = Rc_srf_O3(l);
                        else if (SpeciesNames(s) == "SO2")
                          Rcut = Rc_srf_SO2(l);
                        else
                          {
                            if (Alpha(s) == 0. && Beta(s) == 0.)
                              Rcut = numeric_limits<real>::infinity();
                            else
                              Rcut = 1. / (Alpha(s) / Rc_srf_SO2(l) + Beta(s)
                                           / Rc_srf_O3(l));
                          }
                        if (SurfaceTemperature(h, j, i) < -1)
                          Rcut = Rcut
                            * exp(0.2 * (-1 - SurfaceTemperature(h, j, i)));
                        Cut(Rcut);


                        /*** In canopy aerodynamic resistance Rac ***/

                        Rac = R_ac0(l) * pow(double(LAI(l)), 0.25)
                          / pow(double(FrictionVelocity(h, j, i)), 2.);
                        Cut(Rac);


                        /*** Rc ***/

                        Rns = 1. / (1. / (Rac + Rg) + 1. / Rcut);

                        Rc = (Rst + Rm(s)) / ((1. - Wst) + (Rst + Rm(s)) / Rns);
                      }
                    else if (option_Rc == "wesely")
                      {
                        /*** Stomatal resistance ***/
                        if ((SurfaceTemperature(h, j, i) < 0.01) ||
                            (SurfaceTemperature(h, j, i) > 39.99))
                          Rsx = 1e+30;
                        else
                          Rsx = Ri(l) * (1. + 40000.
                                         / float(pow(double(SolarRadiation(h, j, i)
                                                            + 0.1), 2.0)))
                            * (400. / (SurfaceTemperature(h, j, i)
                                       * (40. - SurfaceTemperature(h, j, i))));
                        Rsx = Rsx * D_H2O / Diffusivity(s);

                        /*** Mesophyll resistance ***/
                        Rmx = 1. / (Henry(s) / 3000. + 100.*Reactivity(s));

                        /*** Cuticle resistance ***/
                        Rlux = Rlu(l) / (1.e-5 * Henry(s) + Reactivity(s));

                        /*** Buoyant-convection resistance ***/
                        Rdc = 100. * (1. + 1000. /
                                      (SolarRadiation(h, j, i) + 10.));

                        /*** Rclx ***/
                        if (SpeciesNames(s) == "O3")
                          Rclx = Rc_srf_O3(l);
                        else if (SpeciesNames(s) == "SO2")
                          Rclx = Rc_srf_SO2(l);
                        else
                          Rclx = 1. / (Alpha(s) / Rc_srf_SO2(l)
                                       + Beta(s) / Rc_srf_O3(l));

                        /*** Rgx ***/
                        if (SpeciesNames(s) == "O3")
                          Rgx = Rg_O3(l);
                        else if (SpeciesNames(s) == "SO2")
                          Rgx = Rg_SO2(l);
                        else
                          Rgx = 1. / (Alpha(s) / Rg_SO2(l)
                                      + Beta(s) / Rg_O3(l));

                        /*** Rc ***/
                        Cut(Rsx);
                        Cut(Rmx);
                        Cut(Rlux);
                        Cut(Rdc);
                        Cut(Rclx);
                        Cut(Rgx);

                        Rc = 1. / (1. / (Rsx + Rmx) + 1. / Rlux
                                   + 1. / (Rdc + Rclx) + 1. / (R_ac0(l) + Rgx));
                      }


                    if (Rc > 9999)
                      Rc = 9999.;
                    else if (Rc < 10)
                      Rc = 10.;

                    /////////////////////////
                    // DEPOSITION VELOCITY //
                    /////////////////////////

                    Velocities(h, j, i) += LandUse(l, j, i) / (Rc + Rb + Ra);

                    DRa(h, j, i, l) = Ra;
                    DRb(h, j, i, l) = Rb;
                    DRc(h, j, i, l) = Rc;
                  }
              }
          current_date.AddHours(int(Delta_t));
        }

      // Writes output files.
      FormatBinary<float> PolairOut;
      PolairOut.Append(Velocities, dir_out + SpeciesNames(s) + ".bin");
      cout << "Done for " << SpeciesNames(s) << endl;

      if (option_save_resistances)
        {
          PolairOut.Append(DRa, dir_out + "Ra.bin");
          PolairOut.Append(DRb, dir_out + "Rb.bin");
          PolairOut.Append(DRc, dir_out + "Rc.bin");
        }

    }

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
