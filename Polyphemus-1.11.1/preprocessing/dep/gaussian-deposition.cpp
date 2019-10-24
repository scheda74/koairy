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

  // Files.
  string meteo_file, species_file, output_file;

  // Option: are comments to be written?
  bool with_comment;

  // Parameterization types to compute scavenging coefficients and deposition
  // velocities.
  string scavenging_type, velocity_type;
  // Files where the values are given for each situation (in case the type is
  // "file").
  string scavenging_file, deposition_file;

  // Scavenging coefficients (s^(-1)).
  Array<real, 2> scavenging_coefficient;
  // Deposition velocities (m/s).
  Array<real, 2> deposition_velocity;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////


  cout << "Reading configuration file...";
  cout.flush();

  ConfigStream config(configuration_file);
  config.SetSection("[data]");
  config.PeekValue("Species", species_file);
  config.PeekValue("Meteo", meteo_file);

  config.SetSection("[scavenging]");
  config.PeekValue("Type", "none | constant | belot | file", scavenging_type);

  if (scavenging_type == "file")
    config.PeekValue("File", scavenging_file);

  config.SetSection("[deposition]");
  config.PeekValue("Type", "none | constant | file", velocity_type);
  if (velocity_type == "file")
    config.PeekValue("File", deposition_file);

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
  cout << " done." << endl;


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

      meteo_stream << endl << endl << endl;
    }

  meteo_stream.close();

  cout << " done." << endl;

  END;

  return 0;
}
