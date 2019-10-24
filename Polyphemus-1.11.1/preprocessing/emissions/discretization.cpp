// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Régis Briant
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
  string source_basename;
  real source_velocity, velocity, temperature, diameter;
  string source_type;
  vector<string> species_list;
  vector<real> rate;
  vector<real> quantity;

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

  // (Broken_line or independant_segment)
  bool with_broken_line;

  vector<vector<real> > puff_quantity;


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
  config.SetSection("[source]");

  config.PeekValue("Source", source_basename);

  config.PeekValue("Source_type", "puff | continuous", source_type);
  config.Find("Species");
  species_list =  split(config.GetLine());
  int Ns = species_list.size();
  config.PeekValue("Velocity", velocity);
  config.PeekValue("Temperature", temperature);
  config.PeekValue("Diameter", ">= 0", diameter);
  date_beg = config.PeekValue("Date_beg");

  if (source_type == "puff")
    {
      config.SetSection("[puff-source]");
      config.Find("Quantity");
      vector<string> tmp = split(config.GetLine());
      if (int(tmp.size()) != Ns)
        throw string("Error in \"PuffEmission<T>::Init()\": ")
          + " there must be one quantity per emitted species.";
      for (int i = 0; i < Ns; i++)
        quantity.push_back(to_num<real>(tmp[i]));
      config.PeekValue("Source_velocity", "positive", source_velocity);
    }
  else if (source_type == "continuous")
    {
      config.SetSection("[plume-source]");
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
  vector<vector<real> > rate_continuous;
  vector<map<string, real> >  node;
  vector<real> line_source_width;
  vector<string> coordinate_line;

  int Nnode;

  cout << "Reading trajectory data...";
  cout.flush();

  ExtStream trajectory_stream(trajectory_file);
  if (!trajectory_stream.is_open())
    throw string("File \"") + trajectory_file + "\" does not exist.";

  real X, Y, Z, Id;

  // Read the first line of the trajectory file
  coordinate_line = split(trajectory_stream.GetLine());

  with_broken_line = (int(coordinate_line.size()) == 4 + Ns);
  if (!with_broken_line)
    {
      if (int(coordinate_line.size()) != 9 + Ns)
        throw string("There must be ") + to_str(9 + Ns)
          + " coordinates per line, but "
          + to_str(coordinate_line.size()) + " were found.";

      // Adds the first point of the segment.
      coord["Id"] = to_num<int>(coordinate_line[0]);
      coord["X"] = to_num<real>(coordinate_line[1]);
      coord["Y"] = to_num<real>(coordinate_line[2]);
      coord["Z"] = to_num<real>(coordinate_line[3]);
      node.push_back(coord);
      // Adds the second point of the segment.
      coord["Id"] = to_num<int>(coordinate_line[4]);
      coord["X"] = to_num<real>(coordinate_line[5]);
      coord["Y"] = to_num<real>(coordinate_line[6]);
      coord["Z"] = to_num<real>(coordinate_line[7]);
      node.push_back(coord);

      line_source_width.push_back(to_num<real>(coordinate_line[8]));
      if (source_type == "continuous")
        {
          vector<real> tmp;
          for (int i = 0; i < Ns; ++i)
            tmp.push_back(to_num<real>(coordinate_line[9 + i]));
          rate_continuous.push_back(tmp);
        }

      while (!trajectory_stream.IsEmpty())
        {
          coordinate_line = split(trajectory_stream.GetLine());

          if (int(coordinate_line.size()) != 9 + Ns)
            throw string("There must be ") + to_str(9 + Ns)
              + string(" coordinates per line, but ")
              + to_str(coordinate_line.size()) + " were found.";

          // Adds the first point of the segment.
          coord["Id"] = to_num<int>(coordinate_line[0]);
          coord["X"] = to_num<real>(coordinate_line[1]);
          coord["Y"] = to_num<real>(coordinate_line[2]);
          coord["Z"] = to_num<real>(coordinate_line[3]);
          node.push_back(coord);
          // Adds the second point of the segment.
          coord["Id"] = to_num<int>(coordinate_line[4]);
          coord["X"] = to_num<real>(coordinate_line[5]);
          coord["Y"] = to_num<real>(coordinate_line[6]);
          coord["Z"] = to_num<real>(coordinate_line[7]);
          node.push_back(coord);

          line_source_width.push_back(to_num<real>(coordinate_line[8]));
          if (source_type == "continuous")
            {
              vector<real> tmp;
              tmp.reserve(Ns);
              for (int i = 0; i < Ns; ++i)
                tmp.push_back(to_num<real>(coordinate_line[9 + i]));
              rate_continuous.push_back(tmp);
            }
        }
    }
  else
    {
      if (int(coordinate_line.size()) != 4 + Ns)
        throw string("There must be ") + to_str(4 + Ns)
          + string(" coordinates per line, but ")
          + to_str(coordinate_line.size()) + " were found.";

      line_source_width.push_back(to_num<real>(coordinate_line[3]));
      coord["X"] = to_num<real>(coordinate_line[0]);
      coord["Y"] = to_num<real>(coordinate_line[1]);
      coord["Z"] = to_num<real>(coordinate_line[2]);
      node.push_back(coord);

      if (int(coordinate_line.size()) == int(4 + Ns))
        {
          line_source_width.push_back(to_num<real>(coordinate_line[3]));
          if (source_type == "continuous")
            {
              vector<real> tmp;
              tmp.reserve(Ns);
              for (int i = 0; i < Ns; ++i)
                tmp.push_back(to_num<real>(coordinate_line[4 + i]));
              rate_continuous.push_back(tmp);
            }
        }

      while (!trajectory_stream.IsEmpty())
        {
          coordinate_line = split(trajectory_stream.GetLine());
          if (int(coordinate_line.size()) != 3)
            {
              line_source_width.push_back(to_num<real>(coordinate_line[3]));

              if (source_type == "continuous")
                {
                  vector<real> tmp;
                  for (int i = 0; i < Ns; i++)
                    tmp.push_back(to_num<real>(coordinate_line[4 + i]));
                  rate_continuous.push_back(tmp);
                }
              if (int(coordinate_line.size()) != 4 + Ns
                  && !trajectory_stream.IsEmpty())
                throw string("There must be ") + to_str(4 + Ns)
                  + string(" coordinates per line, but ")
                  + to_str(coordinate_line.size()) + " were found.";
            }

          coord["X"] = to_num<real>(coordinate_line[0]);
          coord["Y"] = to_num<real>(coordinate_line[1]);
          coord["Z"] = to_num<real>(coordinate_line[2]);
          node.push_back(coord);

          if (int(coordinate_line.size()) == int(4 + Ns))
            {
              line_source_width.push_back(to_num<real>(coordinate_line[3]));
              if (source_type == "continuous")
                {
                  vector<real> tmp;
                  for (int i = 0; i < Ns; i++)
                    tmp.push_back(to_num<real>(coordinate_line[4 + i]));
                  rate_continuous.push_back(tmp);
                }
            }
        }
    }

  Nnode = int(node.size());
  cout << "  done." << endl;


  ////////////////////
  // DISCRETIZATION //
  ////////////////////


  map<string, real> data;
  vector<string> name_list;
  vector<map<string, real> > source_data;
  vector<Date> date_list;
  int Nsource;

  vector<real> segment_length(Nnode - 1);
  vector<real> time(Nnode);
  Date date;
  int Ntot = 0;

  real length = 0., area = 0., v = 0.;
  if (source_type == "puff")
    {
      v = source_velocity / 3.6;
      time[0] = 0.;
    }

  // Calculation of each segment length and of the total trajectory length.
  for (int id_width = 0, i = 0; i < Nnode - 1; i++)
    {
      segment_length[i] = sqrt(pow((node[i]["X"] - node[i + 1]["X"]), 2) +
                               pow((node[i]["Y"] - node[i + 1]["Y"]), 2) +
                               pow((node[i]["Z"] - node[i + 1]["Z"]), 2));
      length += segment_length[i];
      area += segment_length[i] * line_source_width[id_width];
      if (source_type == "puff" && v != 0.)
        time[i + 1] = time[i] + segment_length[i] / v;
      else if (source_type == "puff" && v == 0.)
        time[i + 1] = 0.;

      id_width++;
      // Skip an iteration (necessary for independant segment).
      if (!with_broken_line)
        i++;
    }
  cout << "Length of the trajectory: " << length << endl;
  cout << "Total area of the trajectory: " << area << endl;

  // Calculation of the points coordinates.
  vector<int> id_section;
  for (int id_width = 0, i = 0; i < Nnode - 1; i++)
    {
      real  delta_X, delta_Y, delta_Z;
      real X1, X2, Y1, Y2, Z1, Z2, cos_theta, sin_theta;
      int Npi = Np;

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
      if (source_type == "continuous" || v == 0.)
        {
          delta_X = (X2 - X1) / real(Npi);
          delta_Y = (Y2 - Y1) / real(Npi);
          delta_Z = (Z2 - Z1) / real(Npi);
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

      for (int k = 0; k <= line_source_width[id_width]; k++)
        {
          cos_theta = -(Y2 - Y1) / segment_length[i];
          sin_theta = (X2 - X1) / segment_length[i];

          for (int j = 0; j < Npi + 1; j++)
            {
              name_list.push_back(source_basename
                                  + "-" + to_str(j)
                                  + "-" + to_str(k));

              // Coordinates.
              data["X"] = X1 + real(j) * delta_X + cos_theta
                * (k - line_source_width[id_width] / 2.);
              data["Y"] = Y1 + real(j) * delta_Y + sin_theta
                * (k - line_source_width[id_width] / 2.);
              data["Z"] = Z1 + real(j) * delta_Z;

              // Plume rise parameters.
              data["velocity"] = velocity;
              data["temperature"] = temperature;
              data["diameter"] = diameter;
              if (line_source_width[id_width] != 0.)
                {
                  if (j == 0 || j == Npi)
                    data["coef_rate"] = segment_length[i]
                      * line_source_width[id_width]
                      / (2 * Npi * (line_source_width[id_width] + 1.));
                  else
                    data["coef_rate"] = segment_length[i]
                      * line_source_width[id_width]
                      / (Npi * (line_source_width[id_width] + 1.));
                }
              else
                {
                  if (j == 0 || j == Npi)
                    data["coef_rate"] = segment_length[i] / (2 * Npi);
                  else
                    data["coef_rate"] = segment_length[i] / Npi;
                }
              source_data.push_back(data);
              if (with_broken_line)
                id_section.push_back(i);
              else
                id_section.push_back(int(i / 2));

              // Release time (seconds from beginning date).
              date_list.push_back(date);
              if (source_type == "puff" && v != 0.)
                date.AddSeconds(int(delta_t));
            }
        }
      id_width ++;

      // Skip an iteration
      if (!with_broken_line)
        i++;
    }

  Nsource = int(source_data.size());
  cout << "Number of points on the trajectory: " << Nsource << endl;

  // Calculation of the mass released for each puff.
  vector<real> tmp;
  if (source_type == "puff")
    for (int s = 0; s < Ns; s++)
      {
        for (int i = 0; i < Nsource; i++)
          tmp.push_back(quantity[s] / Nsource);
        puff_quantity.push_back(tmp);
      }


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
      sources << "[source]" << endl << endl;
      sources << "Source: " << name_list[i] << endl;

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
      sources << "Species: ";
      for (int s = 0; s < Ns; s++)
        sources << species_list[s] << " ";
      sources << endl;
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
      sources << "Date_beg: "
              << date_list[i].GetDate("%y-%m-%d_%h-%i-%s") << endl;

      if (source_type == "puff")
        {
          if (with_comment)
            sources << "# Total mass released (mass)" << endl;
          sources << "Quantity: ";
          for (int s = 0; s < Ns; s++)
            sources << puff_quantity[s][i] << " ";
          sources << endl;
        }
      else if (source_type == "continuous")
        {
          sources << "Date_end: "
                  << date_end.GetDate("%y-%m-%d_%h-%i-%s") << endl;
          if (with_comment)
            sources << endl << "# Source rate (mass/s)" << endl;
          sources << "Rate: ";
          for (int s = 0; s < Ns; s++)
            sources << rate_continuous[id_section[i]][s]
              * source_data[i]["coef_rate"] << " ";
          sources << endl;
        }
      sources << endl << endl;
    }
  sources.close();

  cout << " done." << endl << endl;

  END;

  return 0;
}
