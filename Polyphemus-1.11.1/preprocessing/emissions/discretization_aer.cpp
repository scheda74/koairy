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

#include "Talos.hxx"
using namespace Talos ;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////


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

  // Source characterization.
  real rate, quantity, source_velocity, velocity, temperature, diameter;
  string species_name, source_type;

  // Files.
  string trajectory_file, output_file;

  // Option: are comments written?
  bool with_comment;

  // Number of points on the trajectory.
  int Np;

  // Time step.
  real delta_t;

  // Emission dates.
  Date date_beg;
  Date date_end;

  vector<real> puff_quantity;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////


  cout << "Reading configuration file...";
  cout.flush();
  ConfigStreams config(configuration_file);

  config.SetSection("[trajectory]");

  // File containing the trajectory data.
  config.PeekValue("Trajectory_file", trajectory_file);

  // Number of points on the trajectory.
  config.PeekValue("Np", "> 0", Np);

  // Time step  between two points on the trajectory (for mobile source).
  config.PeekValue("Delta_t", "> 0", delta_t);

  // Source data.
  config.SetSection("[aerosol_source]");

  config.PeekValue("Source_type", "puff | continuous", source_type);
  config.PeekValue("Species_name", species_name);
  config.PeekValue("Velocity", velocity);
  config.PeekValue("Temperature", temperature);
  config.PeekValue("Diameter", "> 0", diameter);
  date_beg = config.PeekValue("Date_beg");

  if (source_type == "puff")
    {
      config.SetSection("[puff-source]");
      config.PeekValue("Quantity", "positive", quantity);
      config.PeekValue("Source_velocity", "positive", source_velocity);
    }
  else if (source_type == "continuous")
    {
      config.SetSection("[plume-source]");
      config.PeekValue("Rate", "positive", rate);
      date_end = config.PeekValue("Date_end");
    }

  // Output file.
  config.SetSection("[output]");
  config.PeekValue("With_comment", with_comment);
  config.PeekValue("Source_file", output_file);

  cout << " done." << endl;


  ////////////////
  // TRAJECTORY //
  ////////////////


  map<string, real> coord;
  vector<map<string, real> >  node;

  int Nnode;

  cout << "Reading trajectory data...";
  cout.flush();

  ExtStream trajectory_stream(trajectory_file);
  if (!trajectory_stream.is_open())
    throw string("File \"") + trajectory_file + "\" does not exist.";

  real X, Y, Z;

  while (!trajectory_stream.IsEmpty())
    {
      trajectory_stream.GetNumber(X);
      trajectory_stream.GetNumber(Y);
      trajectory_stream.GetNumber(Z);

      coord["X"] = X;
      coord["Y"] = Y;
      coord["Z"] = Z;

      node.push_back(coord);
    }
  Nnode = int(node.size());

  cout << "  done." << endl;


  ////////////////////
  // DISCRETIZATION //
  ////////////////////


  map<string, real> data;
  vector<map<string, real> > source_data;
  vector<string> species;
  vector<Date> date_list;
  int Nsource;

  vector<real> segment_length(Nnode - 1);
  real length = 0.;
  real v = 0.;
  int Ntot = 0;
  vector<real> time(Nnode);
  Date date;

  if (source_type == "puff")
    {
      v = source_velocity / 3.6;
      time[0] = 0.;
    }

  // Calculation of each segment length and of the total trajectory length.
  for (int i = 0; i < Nnode - 1; i++)
    {
      segment_length[i] = sqrt(pow((node[i]["X"] - node[i + 1]["X"]), 2) +
                               pow((node[i]["Y"] - node[i + 1]["Y"]), 2) +
                               pow((node[i]["Z"] - node[i + 1]["Z"]), 2));
      length += segment_length[i];
      if (source_type == "puff" && v != 0.)
        time[i + 1] = time[i] + segment_length[i] / v;
      else if (source_type == "puff" && v == 0.)
        time[i + 1] = 0.;
    }
  cout << "Length of the trajectory: " << length << endl;

  // Calculation of the points coordinates.
  for (int i = 0; i < Nnode - 1; i++)
    {
      real  delta_X, delta_Y, delta_Z;
      real X1, X2, Y1, Y2, Z1, Z2;
      int Npi;

      X1 = node[i]["X"];
      Y1 = node[i]["Y"];
      Z1 = node[i]["Z"];
      X2 = node[i + 1]["X"];
      Y2 = node[i + 1]["Y"];
      Z2 = node[i + 1]["Z"];

      date = date_beg;
      if (source_type == "puff")
        date.AddSeconds(int(time[i]));

      // Number of points to calculate on the current segment.
      if (source_type == "continous" || v == 0.)
        {
          if (i != Nnode - 2)
            {
              Npi =  int(segment_length[i] * Np / length);
              Ntot += Npi;
            }
          else
            Npi = Np - Ntot;

          delta_X = (X2 - X1) / (Npi + 1);
          delta_Y = (Y2 - Y1) / (Npi + 1);
          delta_Z = (Z2 - Z1) / (Npi + 1);
        }
      else
        {
          Npi = max(int(segment_length[i] / (v * delta_t)), 0);
          // To avoid taking a node twice into account.
          if (Npi * v * delta_t == segment_length[i])
            Npi -= 1;
          delta_X = (X2 - X1) * v * delta_t / segment_length[i];
          delta_Y = (Y2 - Y1) * v * delta_t / segment_length[i];
          delta_Z = (Z2 - Z1) * v * delta_t / segment_length[i];
        }

      // Data of the calculated points on the current segment.
      for (int j = 0; j < Npi + 1; j++)
        {
          // Coordinates.
          data["X"] = X1 + real(j) * delta_X;
          data["Y"] = Y1 + real(j) * delta_Y;
          data["Z"] = Z1 + real(j) * delta_Z;

          // Rate.
          if (source_type == "continuous")
            data["rate"] = rate;

          // Plume rise parameters.
          data["velocity"] = velocity;
          data["temperature"] = temperature;
          data["diameter"] = diameter;

          source_data.push_back(data);
          species.push_back(species_name);

          // Release time (seconds from beginning date).
          date_list.push_back(date);
          if (source_type == "puff" && v != 0.)
            date.AddSeconds(int(delta_t));
        }
    }

  data["X"] = node[Nnode - 1]["X"];
  data["Y"] = node[Nnode - 1]["Y"];
  data["Z"] = node[Nnode - 1]["Z"];
  date = date_beg;
  if (source_type == "puff")
    date.AddSeconds(int(time[Nnode - 1]));
  if (source_type == "continuous")
    {
      data["rate"] = rate;
    }
  data["velocity"] = velocity;
  data["temperature"] = temperature;
  data["diameter"] = diameter;

  source_data.push_back(data);
  species.push_back(species_name);
  date_list.push_back(date);

  Nsource = int(source_data.size());
  cout << "Number of points on the trajectory: " << Nsource << endl;

  // Calculation of the mass released for each puff.
  if (source_type == "puff")
    for (int i = 0; i < Nsource; i++)
      puff_quantity.push_back(quantity / Nsource);


  ////////////
  // OUTPUT //
  ////////////


  cout << "Writing source data...";
  cout.flush();

  ofstream sources(output_file.c_str());

  if (!sources.is_open())
    throw string("Unable to open file \"") + output_file + "\".";

  for (int i = 0; i < Nsource; i++)
    {
      sources << "[aerosol_source]" << endl << endl;

      // Source coordinates.
      if (with_comment)
        sources << "# Source coordinates (meters)" << endl;
      sources << "Abscissa: " << source_data[i]["X"] << endl;
      sources << "Ordinate: " << source_data[i]["Y"] << endl;
      sources << "Altitude: " << source_data[i]["Z"] << endl;
      sources << endl;

      // Source species and type.
      if (with_comment)
        sources << "# Species name" << endl;
      sources << "Species_name: " << species[i] << endl;
      if (with_comment)
        sources << "# Source type " << endl;
      sources << "Type: " << source_type << endl;
      sources << endl;

      // Plume rise parameters.
      if (with_comment)
        sources << "# Source velocity (m/s)" << endl;
      sources << "Velocity: " << source_data[i]["velocity"] << endl;
      if (with_comment)
        sources << "# Source temperature (Celsius degree)" << endl;
      sources << "Temperature: " << source_data[i]["temperature"] << endl;
      if (with_comment)
        sources << "# Source diameter (m)" << endl;
      sources << "Diameter: " << source_data[i]["diameter"] << endl;
      sources << endl;

      // Beginning date.
      sources << "Date_beg: " << date_list[i].GetDate("%y-%m-%d_%h-%i-%s") << endl;

      if (source_type == "puff")
        {
          if (with_comment)
            sources << "# Total mass released (mass)" << endl;
          sources << "Quantity: " << puff_quantity[i] << endl;
        }
      else if (source_type == "continuous")
        {
          sources << "Date_end: " << date_end.GetDate("%y-%m-%d_%h-%i-%s") << endl;
          if (with_comment)
            sources << endl << "# Source rate (mass/s)" << endl;
          sources << "Rate: " << source_data[i]["rate"] << endl;
        }
      sources << endl << endl;
    }
  sources.close();

  cout << " done." << endl << endl;

  END;

  return 0;
}
