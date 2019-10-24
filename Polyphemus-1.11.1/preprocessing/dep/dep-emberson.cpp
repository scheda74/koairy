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


// This program is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).
// Base program: "depvel.f90".


// This program computes deposition velocities based on Emberson et
// al. (2000).


// Warning: it is assumed that 't_min' and 'Delta_t' are integers.


#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string configuration_file, sec_config_file,
    default_name("dep-emberson.cfg");
  Date date_beg, date_end;

  parse_argument(argc, argv, configuration_file, sec_config_file,
                 date_beg, date_end, default_name);


  // Configuration files.
  ConfigStreams configuration(configuration_file);
  if (sec_config_file != "")
    configuration.AddFile(sec_config_file);


  ///////////////////
  // CONFIGURATION //
  ///////////////////

  typedef float real;

  const real pi = 3.14159265358979323846264;

  int h, i, j, s, v, c;

  cout << "Reading configuration files...";
  cout.flush();

  /*** Domain ***/

  configuration.SetSection("[domain]");

  // Main dimensions.
  int Nt, Nz, Ny, Nx;
  // Associated steps.
  real Delta_t, Delta_y, Delta_x;
  real y_min, x_min;

  configuration.PeekValue("Nz", "> 0", Nz);
  configuration.PeekValue("Ny", "> 0", Ny);
  configuration.PeekValue("Nx", "> 0", Nx);

  configuration.PeekValue("Delta_t", "> 0", Delta_t);
  configuration.PeekValue("Delta_y", "> 0", Delta_y);
  configuration.PeekValue("Delta_x", "> 0", Delta_x);

  configuration.PeekValue("y_min", y_min);
  configuration.PeekValue("x_min", x_min);

  // Dates.
  string date_meteo_str;
  configuration.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  double difference = date_beg.GetSecondsFrom(date_meteo);
  if (difference < 0)
    throw string("\nthe date you provide in command line ")
      + "should be after the date in the main configuration file.";
  int step = int(difference  / 3600 / Delta_t + 0.5);

  Nt = compute_Nt(date_beg, date_end, Delta_t);
  int t_min = date_beg.GetHour();

  /*** Paths ***/

  // Input/output files.
  configuration.SetSection("[paths]");

  // Output directory.
  string dir_out;
  configuration.PeekValue("Directory_dep", dir_out);

  // Paths to input meteorological files.
  string SurfaceTemperature_file, FrictionVelocity_file,
    CanopyWetness_file, SurfaceRelativeHumidity_file,
    AerodynamicResistance_file, Attenuation_file, Altitude_file;

  configuration.PeekValue("SurfaceTemperature", SurfaceTemperature_file);
  configuration.PeekValue("SurfaceRelativeHumidity",
                          SurfaceRelativeHumidity_file);
  configuration.PeekValue("Altitude", Altitude_file);
  configuration.PeekValue("AerodynamicResistance",
                          AerodynamicResistance_file);
  configuration.PeekValue("Attenuation", Attenuation_file);
  configuration.PeekValue("FrictionVelocity", FrictionVelocity_file);

  int Nveg, Nc;
  string LUC_file, land_data_file;

  // Number of vegetation classes.
  configuration.PeekValue("Nveg", "positive", Nveg);
  // Number of land use categories.
  configuration.PeekValue("Nc", "positive", Nc);
  configuration.PeekValue("LUC_file", LUC_file);
  configuration.PeekValue("Land_data", land_data_file);

  // Species data.
  string species_data;
  configuration.PeekValue("Species_data", species_data);

  /*** Species ***/

  configuration.SetSection("[species]");

  // Number of species.
  int Ns;
  configuration.PeekValue("Ns", "positive", Ns);

  cout << " done." << endl;


  ///////////
  // GRIDS //
  ///////////


  cout << "Memory allocation for data fields...";
  cout.flush();

  // Input grids.
  RegularGrid<real> GridT(t_min, Delta_t, Nt);
  RegularGrid<real> GridZ(Nz);
  // Includes ground level.
  RegularGrid<real> GridZ_ext(Nz + 1);
  RegularGrid<real> GridY(y_min, Delta_y, Ny);
  RegularGrid<real> GridX(x_min, Delta_x, Nx);
  RegularGrid<real> GridY_interf(y_min - Delta_y / 2., Delta_y, Ny + 1);
  RegularGrid<real> GridX_interf(x_min - Delta_x / 2., Delta_x, Nx + 1);
  RegularGrid<real> GridC(Nc);
  RegularGrid<real> GridV(Nveg);

  // Species.
  RegularGrid<real> GridS(Ns);


  //////////
  // DATA //
  //////////


  // Input fields.
  Data<real, 3> LandUse(GridC, GridY, GridX);
  Data<real, 2> RoughnessHeight_cell(GridY, GridX);

  Data<real, 4> Altitude(GridT, GridZ_ext, GridY, GridX);
  Data<real, 3> SurfaceTemperature(GridT, GridY, GridX);
  Data<real, 3> SurfaceRelativeHumidity(GridT, GridY, GridX);
  Data<real, 3> AerodynamicResistance(GridT, GridY, GridX);
  Data<real, 3> Attenuation(GridT, GridY, GridX);
  Data<real, 3> FrictionVelocity(GridT, GridY, GridX);

  Data<string, 1, real> SpeciesNames(GridS);
  Data<real, 1> MolarMass(GridS);
  Data<real, 1> Henry(GridS);
  Data<real, 1> Diffusivity(GridS);
  Data<real, 1> Reactivity(GridS);
  Data<real, 1> Alpha(GridS);
  Data<real, 1> Beta(GridS);
  Data<real, 1> Rm(GridS);

  Data<real, 1> gmax(GridV);
  Data<real, 1> fmin(GridV);
  Data<real, 1> deptmin(GridV);
  Data<real, 1> deptopt(GridV);
  Data<real, 1> deptmax(GridV);
  Data<real, 1> depalph(GridV);
  Data<real, 1> depvpd1(GridV);
  Data<real, 1> depvpd2(GridV);
  Data<real, 1> depsgs(GridV);
  Data<real, 1> depegs(GridV);
  Data<real, 1> depsgl(GridV);
  Data<real, 1> depegl(GridV);
  Data<real, 1> deplai1(GridV);
  Data<real, 1> deplai2(GridV);
  Data<real, 1> depphe0(GridV);
  Data<real, 1> depphe1(GridV);
  Data<real, 1> depphe2(GridV);
  Data<real, 1> zcanopy(GridV);
  Data<real, 1> RGSO3(GridV);
  Data<real, 1> RGSSO2(GridV);
  Data<real, 1> so2rh(GridV);

  // Output fields.
  Data<real, 4> DepositionVelocity(GridS, GridT, GridY, GridX);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////


  cout << "Extracting input data...";
  cout.flush();

  FormatFormattedText InputParam("<e><e><e><e><e>");

  // Species data.
  InputParam.Read(species_data, "0", SpeciesNames);
  InputParam.Read(species_data, "1", MolarMass);
  InputParam.Read(species_data, "2", Henry);
  InputParam.Read(species_data, "3", Reactivity);

  int NH3_index = -1;
  for (s = 0; s < Ns; s++)
    if (SpeciesNames(s) == "NH3")
      NH3_index = s;
  int HNO3_index = -1;
  for (s = 0; s < Ns; s++)
    if (SpeciesNames(s) == "HNO3")
      HNO3_index = s;

  Data<real, 1> factRb(Ns);
  for (s = 0; s < Ns; s++)
    factRb(s) = pow(.15 * sqrt(MolarMass(s) / 18.) / (.25 * .72), 2. / 3.)
      * 2. / .4;
  for (s = 0; s < Ns; s++)
    Rm(s) = .01 / (Henry(s) / 3000. + 100. * Reactivity(s));

  // Land data.
  ifstream InputData(land_data_file.c_str());
  string tmp;
  getline(InputData, tmp);
  FormatText InputDataFormat;

  InputDataFormat.Read(InputData, gmax);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, fmin);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, deptmin);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, deptopt);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, deptmax);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depalph);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depvpd1);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depvpd2);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depsgs);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depegs);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depsgl);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depegl);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, deplai1);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, deplai2);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depphe0);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depphe1);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, depphe2);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, zcanopy);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, RGSO3);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, RGSSO2);
  getline(InputData, tmp);
  InputDataFormat.Read(InputData, so2rh);

  // Polair inputs.
  FormatBinary<float> InputPolair;

  InputPolair.Read(LUC_file, LandUse);

  InputPolair.ReadSteps(SurfaceTemperature_file, step, SurfaceTemperature);
  InputPolair.ReadSteps(SurfaceRelativeHumidity_file, step,
                        SurfaceRelativeHumidity);
  InputPolair.ReadSteps(Altitude_file, step, Altitude);
  InputPolair.ReadSteps(FrictionVelocity_file, step, FrictionVelocity);
  InputPolair.ReadSteps(AerodynamicResistance_file, step,
                        AerodynamicResistance);
  InputPolair.ReadSteps(Attenuation_file, step, Attenuation);

  // Warning: land use and vegetation descriptions cannot be changed until the
  // following lines are read in the configuration.
  Data<real, 4> fveg(Nc, Nveg, Ny, Nx);
  fveg.Fill(0.);
  for (j = 0; j < Ny; j++)
    for (i = 0; i < Nx; i++)
      {
        fveg(0, 4, j, i) = 0.7;
        fveg(0, 5, j, i) = 0.3;
        fveg(1, 8, j, i) = 1.0;
        fveg(2, 12, j, i) = 1.0;
        fveg(3, 13, j, i) = 1.0;
        fveg(4, 15, j, i) = 1.0;
        fveg(5, 9, j, i) = 1.0;
        fveg(6, 0, j, i) = 1.0;
        fveg(7, 1, j, i) = 1.0;
        fveg(8, 13, j, i) = 1.0;
      }

  cout << " done." << endl;


  ///////////////////////////
  // DEPOSOTION VELOCITIES //
  ///////////////////////////

  cout << "Computing deposition velocities...";
  cout.flush();

  SurfaceTemperature.Add(-273.15);

  real P_sat_vap;
  real Rlow;
  real dsnow = 0.;

  real zenith_angle;
  real cloud_factor;

  real surf_temp;

  real sol_time, hour, loc_time, dec;
  Date current_date(date_beg);
  Date sol_date;
  sol_date.SetYear(date_beg.GetYear());
  sol_date.SetMonth(12);
  sol_date.SetDay(21);
  sol_date.SetHour(12);

  real pardif, pardir, facQ, beamrad, dirrad, difrad;
  real sunlai, shalai, shapar, sunpar, Glisun, Glisha,
    O3Gext, Rin;
  real fsn, frh, SO2Rdry, SO2Rwet, Gsto, Gnst;
  real Rb, Rc;

  Data<real, 1> ftem(Nveg), fvpd(Nveg), fswp(Nveg),
    dlai(Nveg), dsai(Nveg), fphen(Nveg), flig(Nveg),
    O3Gsto(Nveg), O3Gnst(Nveg), SO2Gnst(Nveg), NH3Gnst(Nveg);
  real djul;

  NH3Gnst() = 0.;

  DepositionVelocity() = 0.;

  for (h = 0; h < Nt; h++)
    {
      sol_time = current_date.GetSecondsFrom(sol_date) / 3600.;
      hour = real(current_date.GetHour());
      djul = real(current_date.GetNumberOfDays()) + hour / 24.;
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            // Zenith angle.
            dec = -0.409105176 * cos(7.168657935e-04 * sol_time);
            loc_time = hour - 12.
              + 0.066666666666667 * (x_min + Delta_x * real(i));
            loc_time *= 0.261799387;
            zenith_angle =
              sin((y_min + Delta_y * real(j)) / 180. * pi) * sin(dec)
              + cos((y_min + Delta_y * real(j)) / 180. * pi)
              * cos(dec) * cos(loc_time);

            cloud_factor = 6. * Attenuation(h, j, i);

            if (zenith_angle > 1.e-3)
              {
                pardir = exp(-0.16 / zenith_angle) * zenith_angle;
                pardif = 0.4 * (1. - pardir) * zenith_angle;
                facQ = 1.
                  + 1. / ((cloud_factor * (pardir + pardif) + 5.e-4)
                          * (cloud_factor * (pardir + pardif) + 5.e-4));
                beamrad = 4.57 * Attenuation(h, j, i) * 1069.
                  * exp(-0.21 / zenith_angle);
                dirrad = beamrad * zenith_angle;
                difrad = 0.13 * beamrad;
              }
            else
              {
                zenith_angle = 1.e-3;
                dirrad = 0.;
                difrad = 0.;
                facQ = 1. + 2000. * 2000.;
              }

            P_sat_vap = 611. * exp(17.27 * SurfaceTemperature(h, j, i)
                                   / (SurfaceTemperature(h, j, i) + 237.29));
            P_sat_vap *= 1.e-3 * (1. - SurfaceRelativeHumidity(h, j, i));

            Rlow = 10. * exp(-SurfaceTemperature(h, j, i) - 4.);

            for (v = 0; v < Nveg; v++)
              {
                surf_temp = max(deptmin(v) + 0.1,
                                SurfaceTemperature(h, j, i));
                surf_temp = min(surf_temp,
                                2. * deptopt(v) - deptmin(v) - 0.1);

                ftem(v) = (surf_temp - deptopt(v)) / (deptopt(v) - deptmin(v));
                ftem(v) = 1. - ftem(v) * ftem(v);

                if (P_sat_vap <= depvpd1(v))
                  fvpd(v) = 1.;
                else if (P_sat_vap >= depvpd2(v))
                  fvpd(v) = .1;
                else
                  fvpd(v) = 1. - (depvpd1(v) - P_sat_vap) * .9
                    / (depvpd1(v) - depvpd2(v));

                fswp(v) = 1.;

                dlai(v) = 0.;
                dsai(v) = 0.;
                fphen(v) = 0.;

                if (djul >= depsgs(v) && djul <= depegs(v))
                  {
                    if (djul < depsgs(v) + depsgl(v))
                      {
                        dlai(v) = deplai1(v) + (djul - depsgs(v))
                          * (deplai2(v) - deplai1(v)) / depsgl(v);
                        dsai(v) = 5. * dlai(v) / 3.5;
                      }
                    else if (djul > depegs(v) - depegl(v))
                      {
                        dlai(v) = deplai1(v) + (depegs(v) - djul)
                          * (deplai2(v) - deplai1(v)) / depegl(v);
                        dsai(v) = dlai(v) + 1.5;
                      }
                    else
                      {
                        dlai(v) = deplai2(v);
                        dsai(v) = dlai(v) + 1.5;
                      }

                    if (djul < depsgs(v) + depphe1(v))
                      fphen(v) = depphe0(v) + (djul - depsgs(v))
                        * (1. - depphe0(v)) / depphe1(v);
                    else if (djul > depegs(v) - depphe2(v))
                      fphen(v) = depphe0(v) + (depegs(v) - djul)
                        * (1. - depphe0(v)) / depphe2(v);
                    else
                      fphen(v) = 1.;
                  }

                if (v <= 3)
                  dsai(v) = dlai(v) + 1.;
                else if (v >= 7)
                  dsai(v) = dlai(v);

                sunlai = (1. - exp(-.5 * dlai(v) / zenith_angle))
                  * 2. * zenith_angle;
                shalai = dlai(v) - sunlai;
                shapar = difrad * exp(-.5 * pow(dlai(v), real(0.7)))
                  + 0.07 * dirrad * (1.1 - 0.1 * dlai(v))
                  * exp(-zenith_angle);
                sunpar = .5 * dirrad / zenith_angle + shapar;
                Glisun = 1. - exp(-depalph(v) * sunpar);
                Glisha = 1. - exp(-depalph(v) * shapar);
                flig(v) = (Glisun * sunlai + Glisha * shalai)
                  / max(dlai(v), 1.e-20);

                O3Gsto(v) = dlai(v) * gmax(v) * fphen(v) * flig(v)
                  * max(fmin(v), ftem(v) * fvpd(v) * fswp(v)) / 410.;

                O3Gext = dsai(v) / 25.;

                Rin = 14. * dsai(v) * zcanopy(v)
                  / max(FrictionVelocity(h, j, i), 1.e-10);

                O3Gnst(v) = O3Gext
                  + 1. / (Rin + .01 * (RGSO3(v) + Rlow + 20. * dsnow));

                fsn = exp(-2.);

                if (SurfaceRelativeHumidity(h, j, i) > so2rh(v))
                  frh = (SurfaceRelativeHumidity(h, j, i) - so2rh(v))
                    / (1. - so2rh(v));
                else
                  frh = 0.;

                SO2Rdry = 1.8 * fsn + Rlow + 20. * dsnow;
                SO2Rwet = 1.0 * fsn + Rlow + 20. * dsnow;

                if (v <= 9)
                  SO2Gnst(v) = (1. - frh) / SO2Rdry + frh / SO2Rwet;
                else
                  SO2Gnst(v) = 100. / RGSSO2(v);

                real ts = SurfaceTemperature(h, j, i);
                if (v <= 9)
                  if (ts <= -5.)
                    NH3Gnst(v) = .1;
                  else if (ts <= 0.)
                    NH3Gnst(v) = .5;
                  else
                    {
                      NH3Gnst(v) = 100.
                        / (0.0455 * (ts + 2.)
                           * exp((1. - SurfaceRelativeHumidity(h, j, i))
                                 / 0.07)
                           * pow(10., -1.1099 * 0. + 1.6769));
                      NH3Gnst(v) = max(min(10., NH3Gnst(v)), .5);
                    }
              }

            for (c = 0; c < Nc; c++)
              for (s = 0; s < Ns; s++)
                {
                  Rb = factRb(s) / FrictionVelocity(h, j, i);

                  for (v = 0; v < Nveg; v++)
                    {

                      Gsto = O3Gsto(v)
                        / (sqrt(MolarMass(s) / 48.) + Rm(s) * O3Gsto(v));
                      Gnst = 1.e-5 * Henry(s) * SO2Gnst(v)
                        + Reactivity(s) * O3Gnst(v);
                      if (s == NH3_index)
                        Gnst = NH3Gnst(v);

                      Rc = 1. / (Gsto + Gnst);

                      if (s == HNO3_index)
                        Rc = max(1.e-2, Rlow);

                      DepositionVelocity(s, h, j, i)
                        += fveg(c, v, j, i) * LandUse(c, j, i)
                        / (AerodynamicResistance(h, j, i)
                           + Rb + 100. * Rc) / Altitude(h, 1, j, i);
                    }
                }
          }

      current_date.AddHours(int(Delta_t));
    }

  cout << " done." << endl;


  //////////////////
  // WRITING DATA //
  //////////////////


  cout << "Writing output data...";
  cout.flush();

  for (s = 0; s < Ns; s++)
    {
      Array<real, 3> Dep
        = DepositionVelocity()(s, Range::all(), Range::all(), Range::all());
      FormatBinary<float>().Append(Dep, dir_out + SpeciesNames(s) + ".bin");
    }

  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
