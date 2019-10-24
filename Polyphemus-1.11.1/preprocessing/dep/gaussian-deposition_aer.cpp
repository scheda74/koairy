// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "SeldonData.hxx"
using namespace SeldonData ;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////


//////////////////////
// FORTRAN FUNCTION //
//////////////////////


#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _compute_gravitational_settling compute_gravitational_settling__
#else
#define _compute_gravitational_settling compute_gravitational_settling_
#endif

extern "C"
{
  void _compute_gravitational_settling(double*, double*,
                                       double*, double*, double*, double*);
}


///////////////////
// MAIN FUNCTION //
///////////////////


int main(int argc, char **argv)

{
  TRY;

  cout << endl;

  string configuration_file;

  parse_argument(argc, argv, configuration_file);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef double real;

  // Meteorological data.
  vector<map<string, real> > meteo_data;
  // Number of meteorological situations.
  int Nmeteo;
  // Stability conditions.
  vector<string> stability;

  // List of gaseous species.
  vector<string> species_list;
  // Number of gas species.
  int Ns;

  // List of aerosol species.
  vector<string> species_list_aer;
  // Number of aerosol species.
  int Ns_aer;

  // Diameter list.
  vector<real> diameter_list;
  // Number of diameters.
  int Ndiam;

  // Files.
  string meteo_file, diameter_file, species_file, output_file;

  // Option: are comments to be written?
  bool with_comment;

  // Parameterization types to compute scavenging coefficients and deposition
  // velocities.
  string scavenging_type, velocity_type;
  // Parameterization types for aerosol species.
  string scavenging_type_aer, velocity_type_aer, velocity_part;
  // Files where the values are given for each situation (in case the type is
  // "file").
  string scavenging_file, deposition_file;
  string scavenging_file_aer, deposition_file_aer;

  // Scavenging coefficients (s^(-1)).
  Array<real, 2> scavenging_coefficient;
  // Deposition velocities (m/s).
  Array<real, 2> deposition_velocity;

  // Scavenging coefficients for aerosol species (s^(-1)).
  Array<real, 3> scavenging_coefficient_aer;
  // Deposition velocities for aerosol species (m/s).
  Array<real, 3> deposition_velocity_aer;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////


  cout << "Reading configuration file...";
  cout.flush();

  ConfigStream config(configuration_file);
  config.SetSection("[data]");
  config.PeekValue("Species", species_file);
  config.PeekValue("Meteo", meteo_file);
  config.PeekValue("Diameter", diameter_file);

  config.SetSection("[scavenging]");
  config.PeekValue("Type", "none | constant | belot | file", scavenging_type);
  config.PeekValue("Type_aer", "none | constant | slinn | file",
                   scavenging_type_aer);
  if (scavenging_type == "file")
    config.PeekValue("File", scavenging_file);
  if (scavenging_type_aer == "file")
    config.PeekValue("File_aer", scavenging_file_aer);

  config.SetSection("[deposition]");
  config.PeekValue("Type", "none | constant | file", velocity_type);
  config.PeekValue("Type_aer", "none | constant | file", velocity_type_aer);
  config.PeekValue("Velocity_part", "diffusive | total", velocity_part);
  if (velocity_type == "file")
    config.PeekValue("File", deposition_file);
  if (velocity_type_aer == "file")
    config.PeekValue("File_aer", deposition_file_aer);

  config.SetSection("[output]");
  config.PeekValue("With_comment", with_comment);
  config.PeekValue("Output_file", output_file);

  cout << " done." << endl;


  ////////////////
  // DATA FILES //
  ////////////////


  cout << "Reading meteorological data...";
  cout.flush();

  // Browses all sections "[situation]".
  ConfigStream meteo(meteo_file);
  Nmeteo = 0;
  while (!meteo.IsEmpty())
    if (split(meteo.GetLine())[0] == "[situation]")
      {
        map<string, real> data;
        string stability_;

        meteo.PeekValue("Temperature", data["temperature"]);
        meteo.PeekValue("Wind_angle", data["wind_angle"]);
        meteo.PeekValue("Wind", data["wind"]);
        meteo.PeekValue("Boundary_height", data["boundary_height"]);
        meteo_data.push_back(data);
        meteo.PeekValue("Stability", stability_);
        stability.push_back(stability_);
        Nmeteo++;
      }
  meteo.Rewind();
  cout << " done." << endl;

  cout << "Reading species...";
  cout.flush();
  ConfigStream species(species_file);

  // Section "[species]" contains all species names.
  species.SetSection("[species]");
  while (!species.IsEmpty())
    species_list.push_back(species.GetElement());
  Ns = int(species_list.size());

  // Section "[species_aer]" contains all aerosol species names.
  species.SetSection("[aerosol_species]");
  while (!species.IsEmpty())
    species_list_aer.push_back(species.GetElement());
  Ns_aer = int(species_list_aer.size());
  cout << " done." << endl;

  cout << "Reading diameter...";
  cout.flush();
  ConfigStream diameter_stream(diameter_file);
  diameter_stream.SetSection("[diameter]");
  while (!diameter_stream.IsEmpty())
    {
      real diameter = to_num<real>(diameter_stream.GetElement());
      diameter_list.push_back(diameter);
    }
  Ndiam = int(diameter_list.size());
  cout << " done." << endl;


  ////////////////
  // SCAVENGING //
  ////////////////


  cout << "Computation of the scavenging coefficients...";
  cout.flush();

  /*** Reads species for which scavenging occurs ***/

  // List of the species for which scavenging occurs (gaseous species).
  vector<string> scavenging_list;
  // Number of species for which scavenging occurs.
  int Nscav;

  species.SetSection("[scavenging]");
  while (!species.IsEmpty())
    scavenging_list.push_back(species.GetElement());
  Nscav = int(scavenging_list.size());

  vector<bool> species_scavenging(Ns);
  for (int i = 0; i < Ns; i++)
    {
      bool in_list = false;
      for (int j = 0; j < Nscav; j++)
        if (species_list[i] == scavenging_list[j])
          in_list = true;
      species_scavenging[i] = in_list;
    }

  /*** Scavenging coefficients for gaseous species ***/

  scavenging_coefficient.resize(Nmeteo, Ns);

  // Scavenging type: none. All coefficients for gaseous species are set to 0.
  if (scavenging_type == "none")
    scavenging_coefficient = 0;

  // Scavenging type: constant. All coefficients for gaseous species are
  // read in the section [scavenging_constant] of the species file.
  else  if (scavenging_type == "constant")
    {
      Array<real, 1> scav_coeff(Ns);
      scav_coeff = 0.;

      species.SetSection("[scavenging_constant]");
      for (int i = 0; i < Ns; i++)
        if (species_scavenging[i])
          species.PeekValue(species_list[i], scav_coeff(i));

      for (int i = 0; i < Nmeteo; i++)
        for (int j = 0; j < Ns; j++)
          scavenging_coefficient(i, j) = scav_coeff(j);
    }

  // Scavenging type: belot. All coefficients for gaseous species are computed
  // with Belot parameterization. The Belot coefficients are read in the
  // section [scavenging_belot] of the species_file.
  else if (scavenging_type == "belot")
    {
      Array<real, 2> belot_coefficient(Ns, 2);
      for (int i = 0; i < Ns; i++)
        if (species_scavenging[i])
          {
            species.SetSection("[scavenging_belot]");
            species.Find(species_list[i]);
            species.GetNumber(belot_coefficient(i, 0));
            species.GetNumber(belot_coefficient(i, 1));
          }
        else
          {
            belot_coefficient(i, 0) = 0.;
            belot_coefficient(i, 1) = 0.;
          }

      for (int i = 0; i < Nmeteo; i++)
        {
          real rainfall_rate;
          meteo.Find("[situation]");
          meteo.PeekValue("Rainfall_rate", rainfall_rate);

          for (int j = 0; j < Ns; j++)
            scavenging_coefficient(i, j) =  belot_coefficient(j, 0) *
              pow(rainfall_rate, belot_coefficient(j, 1));
        }
      meteo.Rewind();
    }

  // Scavenging type: file. The file containing the list of scavenging
  // coefficients for all species and all meteorological situations is read.
  else if (scavenging_type == "file")
    {
      ConfigStream scav(scavenging_file);
      int i = 0;
      while (!scav.IsEmpty())
        if (split(scav.GetLine())[0] == "[situation]")
          {
            for (int j = 0; j < Ns; j++)
              if (species_scavenging[j])
                scav.PeekValue(species_list[j], scavenging_coefficient(i, j));
            i++;
          }
    }

  /*** Scavenging coefficient for aerosol species ***/

  scavenging_coefficient_aer.resize(Nmeteo, Ns_aer, Ndiam);

  // Scavenging type_aer: none. All coefficients for aerosol species are set
  // to 0.
  if (scavenging_type_aer == "none")
    scavenging_coefficient_aer = 0;

  // Scavenging type_aer: constant. All coefficients for aerosol species are
  // read in the section [scavenging_constant_aer] of the species file.
  else  if (scavenging_type_aer == "constant")
    {
      Array<real, 1> scav_coeff(Ndiam);

      species.SetSection("[scavenging_constant_aer]");
      for (int i = 0; i < Ndiam; i++)
        species.PeekValue(to_str(i), scav_coeff(i));
      for (int i = 0; i < Nmeteo; i++)
        for (int j = 0; j < Ns_aer; j++)
          for (int k = 0; k < Ndiam; k++)
            scavenging_coefficient_aer(i, j, k) = scav_coeff(k);
    }

  // Scavenging type_aer: slinn. Coefficients for aerosol species are computed
  // with the Slinn parameterization.
  else  if (scavenging_type_aer == "slinn")
    {
      Array<real, 2> scav_coeff(Nmeteo, Ndiam);
      scav_coeff = 0.;
      string value;
      config.SetSection("[scavenging]");
      config.PeekValue("Value", "best_estimate | conservative", value);
      for (int i = 0; i < Nmeteo; i++)
        {
          real rainfall_rate;
          meteo.Find("[situation]");
          meteo.PeekValue("Rainfall_rate", rainfall_rate);
          if (rainfall_rate != 0.)
            {
              // Mean rain droplet diameter (mm)
              real droplet_diameter = 0.7 * pow(rainfall_rate, 0.25) ;
              for (int j = 0; j < Ndiam; j++)
                {
                  real efficiency = 1.0;
                  if (value == "conservative")
                    efficiency = 1.0;
                  else if (value == "best_estimate")
                    {
                      if (diameter_list[j] <= 1.0)
                        efficiency = 0.1;
                      else if (diameter_list[j] > 1.0
                               && diameter_list[j] < 10.)
                        efficiency = 0.1 + 0.9 * (diameter_list[j] - 1.0)
                          / (10. - diameter_list[j]);
                      else
                        efficiency = 1.0;
                    }
                  scav_coeff(i, j) = 3. * efficiency * rainfall_rate
                    / (3600. * 2. * droplet_diameter);
                }
            }
        }
      meteo.Rewind();

      for (int i = 0; i < Nmeteo; i++)
        for (int j = 0; j < Ns_aer; j++)
          for (int k = 0; k < Ndiam; k++)
            scavenging_coefficient_aer(i, j, k) = scav_coeff(i, k);
    }
  else if (scavenging_type_aer == "file")
    {
      ConfigStream scav(scavenging_file_aer);
      int i = 0;
      while (!scav.IsEmpty())
        if (split(scav.GetLine())[0] == "[situation]")
          {
            for (int j = 0; j < Ns_aer; j++)
              for (int k = 0; k < Ndiam; k++)
                {
                  string species_bin = species_list_aer[j] + '_' + to_str(k);
                  scav.PeekValue(species_bin,
                                 scavenging_coefficient_aer(i, j, k));
                }
            i++;
          }
    }

  cout << " done." << endl;


  ////////////////
  // DEPOSITION //
  ////////////////


  cout << "Computation of the deposition velocities..";
  cout.flush();

  /*** Reads species for which deposition occurs ***/

  // List of the species for which deposition occurs (gas species).
  vector<string> deposition_list;
  // Number of species for which deposition occurs.
  int Ndep;

  species.SetSection("[deposition]");
  while (!species.IsEmpty())
    deposition_list.push_back(species.GetElement());
  Ndep = int(deposition_list.size());

  vector<bool> species_deposition(Ns);
  for (int i = 0; i < Ns; i++)
    {
      bool in_list = false;
      for (int j = 0; j < Ndep; j++)
        if (species_list[i] == deposition_list[j])
          in_list = true;
      if (in_list)
        species_deposition[i] = true;
      else species_deposition[i] = false;
    }

  /*** Deposition velocities for gaseous species ***/

  deposition_velocity.resize(Nmeteo, Ns);

  // Deposition type: none. All velocities for gaseous species are set to 0.
  if (velocity_type == "none")
    deposition_velocity = 0.;

  // Deposition type: constant. All velocities for gaseous species are
  // read in the section [deposition_constant] of the species file.
  else if (velocity_type == "constant")
    {
      Array<real, 1> dep_velocity(Ns);
      dep_velocity = 0.;

      species.SetSection("[deposition_constant]");
      for (int i = 0; i < Ns; i++)
        {
          if (species_deposition[i])
            species.PeekValue(species_list[i], dep_velocity(i));
        }

      for (int i = 0; i < Nmeteo; i++)
        for (int j = 0; j < Ns; j++)
          deposition_velocity(i, j) = dep_velocity(j);
    }
  else if (velocity_type == "file")
    {
      ConfigStream dep(deposition_file);
      int i = 0;
      while (!dep.IsEmpty())
        if (split(dep.GetLine())[0] == "[situation]")
          {
            for (int j = 0; j < Ns; j++)
              if (species_deposition[j])
                dep.PeekValue(species_list[j], deposition_velocity(i, j));
            i++;
          }
    }

  /*** Deposition velocities for aerosol species ***/

  deposition_velocity_aer.resize(Nmeteo, Ns_aer, Ndiam);

  // Deposition type_aer: none. All velocities for aerosol species are set to
  // 0.
  if (velocity_type_aer == "none")
    deposition_velocity_aer = 0.;

  // Deposition type_aer: constant. The deposition velocities of aerosol
  // species are read in the section [deposition_constant_aer].
  else  if (velocity_type_aer == "constant")
    {
      // Gets the deposition velocity (diffusive part) (m/s)
      Array<real, 1> dep_velocity(Ndiam);
      species.SetSection("[deposition_constant_aer]");
      for (int i = 0; i < Ndiam; i++)
        species.PeekValue(to_str(i), dep_velocity(i));

      for (int i = 0; i < Nmeteo; i++)
        for (int j = 0; j < Ns_aer; j++)
          for (int k = 0; k < Ndiam; k++)
            deposition_velocity_aer(i, j, k) = dep_velocity(k);
    }
  else if (velocity_type_aer == "file")
    {
      ConfigStream dep(deposition_file_aer);
      int i = 0;
      while (!dep.IsEmpty())
        if (split(dep.GetLine())[0] == "[situation]")
          {
            for (int j = 0; j < Ns_aer; j++)
              for (int k = 0; k < Ndiam; k++)
                {
                  string species_bin = species_list_aer[j] + '_' + to_str(k);
                  dep.PeekValue(species_bin,
                                deposition_velocity_aer(i, j, k));
                }
            i++;
          }
    }

  // If velocity_part is "all", the given deposition velocities are supposed
  // to be the total deposition velocity; if it is "diffusive", the
  // gravitational settling velocity is computed and added.
  if (velocity_part == "diffusive")
    {
      // Gets the aerosol density (kg/m^3)
      vector<real> density_list(Ns_aer);
      species.SetSection("[density_aer]");
      for (int i = 0; i < Ns_aer; i++)
        species.PeekValue(species_list_aer[i], density_list[i]);

      // Computes the gravitational settling and the deposition velocity (m/s)
      for (int i = 0; i < Nmeteo; i++)
        {
          meteo.Find("[situation]");
          meteo.PeekValue("Pressure", meteo_data[i]["pressure"]);

          for (int j = 0; j < Ns_aer; j++)
            for (int k = 0; k < Ndiam; k++)
              {
                real vd_tot = 0.;
                real air_path = 0.;
                real settling_velocity = 0.;
                real temperature_kelvin = meteo_data[i]["temperature"]
                  + 273.15;
                real diameter_m = diameter_list[k] * 1.e-6;
                _compute_gravitational_settling(&temperature_kelvin,
                                                &meteo_data[i]["pressure"],
                                                &density_list[j],
                                                &diameter_m,
                                                &air_path,
                                                &settling_velocity);
                vd_tot = settling_velocity
                  / (1. - exp(-settling_velocity
                              / deposition_velocity_aer(i, j, k)));
                deposition_velocity_aer(i, j, k) = vd_tot;
              }
        }
      meteo.Rewind();
    }

  cout << "done." << endl;


  ////////////
  // OUTPUT //
  ////////////


  cout << "Writing data...";
  cout.flush();

  ofstream meteo_stream(output_file.c_str());
  if (!meteo_stream.is_open())
    {
      cout << "unable to open file \"" << output_file << "\"." << endl;
      return 1;
    }

  for (int i = 0; i < Nmeteo; i++)
    {
      meteo_stream << "[situation]" << endl << endl;

      if (with_comment)
        meteo_stream << "# Temperature (Celsius degrees)\n";
      meteo_stream << "Temperature = " << meteo_data[i]["temperature"]
                   << endl;

      if (with_comment)
        meteo_stream << endl << "# Wind angle (degrees)\n";
      meteo_stream << "Wind_angle = " <<  meteo_data[i]["wind_angle"]
                   << endl;

      if (with_comment)
        meteo_stream << endl << "# Wind speed (m/s)\n";
      meteo_stream << "Wind = " << meteo_data[i]["wind"] << endl;

      if (with_comment)
        meteo_stream << endl << "# Boundary height (m)\n";
      meteo_stream << "Boundary_height = "
                   << meteo_data[i]["boundary_height"] << endl;

      if (with_comment)
        meteo_stream << endl << "# Stability class \n";
      meteo_stream << "Stability = " << stability[i] << endl;

      if (with_comment)
        meteo_stream << endl
                     << "# Scavenging coefficients of"
                     << " gaseous species (s^-1)\n";
      meteo_stream << "Scavenging_coefficient =\n";
      for (int j = 0; j < Ns; j++)
        meteo_stream << species_list[j] << "  "
                     << scavenging_coefficient(i, j) << "\t";
      meteo_stream << endl;

      if (with_comment)
        meteo_stream << endl
                     << "# Deposition velocities "
                     << "of gaseous species (m/s)\n";
      meteo_stream << "Deposition_velocity =\n";
      for (int j = 0; j < Ns; j++)
        meteo_stream << species_list[j] << "  "
                     << deposition_velocity(i, j) << "\t";
      meteo_stream << endl;

      if (with_comment)
        meteo_stream << endl
                     << "# Scavenging coefficients "
                     << "of aerosol species (s^-1)\n";
      meteo_stream << "Scavenging_coefficient_aer =\n";
      for (int j = 0; j < Ns_aer; j++)
        for (int k = 0; k < Ndiam; k++)
          meteo_stream << species_list_aer[j] << "_" << k << "  "
                       << scavenging_coefficient_aer(i, j, k) << "\t";
      meteo_stream << endl;

      if (with_comment)
        meteo_stream << endl
                     << "# Deposition velocities "
                     << "of the aerosol species (m/s)\n";
      meteo_stream << "Deposition_velocity_aer =\n";
      for (int j = 0; j < Ns_aer; j++)
        for (int k = 0; k < Ndiam; k++)
          meteo_stream << species_list_aer[j] << "_" << k << "  "
                       << deposition_velocity_aer(i, j, k) << "\t";

      meteo_stream << endl << endl << endl;
    }

  meteo_stream.close();

  cout << " done." << endl;

  END;

  return 0;
}
