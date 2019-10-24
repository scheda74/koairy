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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_AER_CXX


#include "SaverUnitPoint_aer.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  ////////////////////
  // SAVERUNITPOINT //
  ////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitPoint_aer<T, ClassModel>
  ::SaverUnitPoint_aer():
    BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitPoint_aer<T, ClassModel>::~SaverUnitPoint_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "indices_list" or "coordinates_list".
  */
  template<class T, class ClassModel>
  string SaverUnitPoint_aer<T, ClassModel>::GetType()  const
  {
    return type;
  }


  //! First initialization.
  /*! Reads the configuration.
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex(string)
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNz()
    <li> GetConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitPoint_aer<T, ClassModel>::Init(ConfigStream& config_stream,
                                               ClassModel& Model)
  {
    config_stream.PeekValue("Type", type);

    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Vertical levels to be saved.
    if (type == "indices_list_aer")
      {
        config_stream.Find("Levels");
        split(config_stream.GetLine(), levels);
        Nlevels = int(levels.size());
        with_indices_list = 1;
      }
    else if (type == "coordinates_list_aer")
      {
        config_stream.Find("Levels_coordinates");
        split(config_stream.GetLine(), levels_coord);
        Nlevels = int(levels_coord.size());
        with_indices_list = 0;
      }
    else
      throw string("Saver unit \"") + type +  "\" is not a point list.";

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species and bins.
    string field, species, file;
    vector<string> vsplit, bounds;
    int first, last, j, k;
    if (this->species_list[0] == "all")
      {
        int i;
        int Ns_aer = Model.GetNs_aer();
        vector<string> list_aer = Model.GetSpeciesList_aer();
        vector<int> bins;
        pair<string, vector<int> > tmp;
        for (i = 0; i < Model.GetNbin_aer(); i++)
          bins.push_back(i);
        tmp.second = bins;
        for (i = 0; i < Ns_aer; i++)
          {
            tmp.first = list_aer[i];
            species_list_aer.push_back(tmp);
            // Adds a buffer to compute averaged concentrations.
            Concentration_.push_back(vector<Data<T, 1> >());
            output_file.push_back(vector<string>());
            for (k = 0; k < int(bins.size()); k++)
              {
                file = find_replace(filename, "&f", tmp.first);
                file = find_replace(file, "&n", to_str(bins[k]));
                output_file[i].push_back(file);
                Concentration_[i].push_back(Data<T, 1>());
              }
          }
      }
    else
      for (unsigned int i = 0; i < this->species_list.size(); i++)
        {
          field = this->species_list[i];
          vsplit = split(field, "{}");
          if (field[0] == '{')
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: it cannot be parsed by the output saver.";
          if (vsplit.size() == 1)
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: bins are needed by the output saver.";
          if (vsplit.size() > 2)
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: it cannot be parsed by the output saver.";

          // Searches for species index in 'species_list_aer' and
          // 'output_file'.  If the species is not found, it is added in
          // 'species_list_aer' and 'output_file'.
          species = split(field, "_")[0];
          j = 0;
          while (j < int(species_list_aer.size())
                 && species_list_aer[j].first != species)
            j++;
          if (j == int(species_list_aer.size()))
            {
              species_list_aer.push_back(pair<string, vector<int> >
                                         (species, vector<int>()));
              // Adds a buffer to compute averaged concentrations.
              Concentration_.push_back(vector<Data<T, 1> >());
              output_file.push_back(vector<string>());
            }
          bounds = split(vsplit[1], "-");
          // First bound.
          first = convert<int>(bounds[0]);
          // Last bound.
          if (bounds.size() != 1)
            last = convert<int>(bounds[1]);
          else
            last = first;
          for (k = first; k < last + 1; k++)
            {
              // Adds the bin (associated to species #j).
              species_list_aer[j].second.push_back(k);
              Concentration_[j].push_back(Data<T, 1>());

              // Field base name is replaced in the generic file name.
              file = find_replace(filename, "&f",
                                  species);

              // Field number is also replaced.
              file = find_replace(file, "&n", to_str(k));
              // Adds the corresponding output file.
              output_file[j].push_back(file);
            }
        }

    // Indices of points to be saved.
    if (with_indices_list)
      {
        config_stream.Find("Indices");
        string line;
        while (config_stream.PeekElement() != "Coordinates"
               && config_stream.PeekElement() != "Point_file")
          {
            line = config_stream.GetLine();
            if (split(line).size() == 3)
              {
                vector<int> temp;
                split(line, temp);
                point_list.push_back(temp);
              }
            else if (split(line).size() == 2)
              for (int i = 0; i < Nlevels; i++)
                {
                  vector<int> temp;
                  temp.push_back(levels[i]);
                  temp.push_back(to_num<int>(split(line)[0]));
                  temp.push_back(to_num<int>(split(line)[1]));
                  point_list.push_back(temp);
                }
            else
              throw string("Two or three indices must be written per line, ")
                + "but "  + to_str(split(line).size()) + " numbers found.";
          }
        // Writing the list of points to be saved.
        Npoint_aer = int(point_list.size());
        string point_file = config_stream.GetValue("Point_file");
        ofstream points(point_file.c_str());
        if (!points.is_open())
          throw string("Unable to open file \"") + point_file + "\".";
        for (int i = 0; i < Npoint_aer; i++)
          points << point_list[i][0] << "\t" << point_list[i][1] << "\t"
                 << point_list[i][2] << endl;
        points.close();

        // Initializing the buffer used to compute averaged concentrations.
        if (this->averaged)
          {
            int s, b, point;
            int base_s, base_b;
            int Nspecies = int(this->species_list_aer.size());
            for (s = 0; s < Nspecies; s++)
              {
                int Nbin = int(this->species_list_aer[s].second.size());
                for (b = 0; b < Nbin; b++)
                  {
                    base_s = Model.
                      GetSpeciesIndex_aer(species_list_aer[s].first);
                    base_b = species_list_aer[s].second[b];
                    Concentration_[s][b].Resize(Npoint_aer);
                    for (point = 0; point < Npoint_aer; point++)
                      Concentration_[s][b](point) = 0.5
                        * Model.GetConcentration_aer()(base_s, base_b,
                                                       point_list[point][0],
                                                       point_list[point][1],
                                                       point_list[point][2]);
                  }
              }
          }
      }
    // Coordinates of points to be saved.
    else
      {
        config_stream.Find("Coordinates");
        string line;
        while (config_stream.PeekElement() != "Point_file")
          {
            line = config_stream.GetLine();
            if (split(line).size() == 3)
              {
                vector<T> temp;
                split(line, temp);
                coord_list.push_back(temp);
              }
            else if (split(line).size() == 2)
              for (int i = 0; i < Nlevels; i++)
                {
                  vector<T> temp;
                  temp.push_back(levels_coord[i]);
                  temp.push_back(to_num<T>(split(line)[0]));
                  temp.push_back(to_num<T>(split(line)[1]));
                  coord_list.push_back(temp);
                }
            else
              throw string("Two or three coordinates must be written")
                + "per line, but " + to_str(split(line).size()) +
                " numbers found.";
          }
        // Writing the list of points to be saved.
        Npoint_aer = int(coord_list.size());
        string point_file = config_stream.GetValue("Point_file");
        ofstream points(point_file.c_str());
        if (!points.is_open())
          throw string("Unable to open file \"") + point_file + "\".";
        for (int i = 0; i < Npoint_aer; i++)
          points << coord_list[i][0] << "\t" << coord_list[i][1] << "\t"
                 << coord_list[i][2] << endl;
        points.close();

        // Initializing the buffer used to compute averaged concentrations.
        if (this->averaged)
          {
            int s, b, point;
            int base_s, base_b;
            int Nspecies = int(this->species_list_aer.size());
            for (s = 0; s < Nspecies; s++)
              {
                int Nbin = int(this->species_list_aer[s].second.size());
                for (b = 0; b < Nbin; b++)
                  {
                    base_s = Model.
                      GetSpeciesIndex_aer(species_list_aer[s].first);
                    base_b = species_list_aer[s].second[b];
                    Concentration_[s][b].Resize(Npoint_aer);
                    for (point = 0; point < Npoint_aer; point++)
                      Concentration_[s][b](point) = 0.5
                        * Model.GetConcentration_aer(base_s, base_b,
                                                     coord_list[point][0],
                                                     coord_list[point][1],
                                                     coord_list[point][2]);
                  }
              }
          }
      }


    int s, b;
    // Empties output files.
    for (s = 0; s < int(this->species_list_aer.size()); s++)
      for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
        ofstream tmp_stream(output_file[s][b].c_str());

    if (this->initial_concentration)
      this->Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitPoint_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetCurrentDate()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitPoint_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    // Indices in the model.
    int base_s, base_b;
    int s, b, point;
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
        for (s = 0; s < int(this->species_list_aer.size()); s++)
          for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
            {
              base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
              base_b = species_list_aer[s].second[b];
              for (point = 0; point < Npoint_aer; point++)
                {
                  if (with_indices_list)
                    Concentration_[s][b](point) += 0.5
                      * Model.GetConcentration_aer()(base_s, base_b,
                                                     point_list[point][0],
                                                     point_list[point][1],
                                                     point_list[point][2]);
                  else
                    Concentration_[s][b](point) += 0.5
                      * Model.GetConcentration_aer(base_s, base_b,
                                                   coord_list[point][0],
                                                   coord_list[point][1],
                                                   coord_list[point][2]);
                }

              Concentration_[s][b].GetArray() /= T(this->interval_length);

              if (Model.GetCurrentDate() >= this->date_beg
                  && Model.GetCurrentDate() <= this->date_end)
                FormatBinary<float>().Append(Concentration_[s][b],
                                             output_file[s][b]);

              for (point = 0; point < Npoint_aer; point++)
                {
                  if (with_indices_list)
                    Concentration_[s][b](point) = 0.5
                      * Model.GetConcentration_aer()(base_s, base_b,
                                                     point_list[point][0],
                                                     point_list[point][1],
                                                     point_list[point][2]);
                  else
                    Concentration_[s][b](point) = 0.5
                      * Model.GetConcentration_aer(base_s, base_b,
                                                   coord_list[point][0],
                                                   coord_list[point][1],
                                                   coord_list[point][2]);
                }
              this->counter = 0;
            }
      else
        for (s = 0; s < int(this->species_list_aer.size()); s++)
          for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
            {
              base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
              base_b = species_list_aer[s].second[b];
              for (point = 0; point < Npoint_aer; point++)
                {
                  if (with_indices_list)
                    Concentration_[s][b](point) +=
                      Model.GetConcentration_aer()(base_s, base_b,
                                                   point_list[point][0],
                                                   point_list[point][1],
                                                   point_list[point][2]);
                  else
                    Concentration_[s][b](point) +=
                      Model.GetConcentration_aer(base_s, base_b,
                                                 coord_list[point][0],
                                                 coord_list[point][1],
                                                 coord_list[point][2]);
                }
            }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (s = 0; s < int(this->species_list_aer.size()); s++)
        for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
          {
            base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
            base_b = species_list_aer[s].second[b];
            Data<T, 1> Concentration_point(Npoint_aer);
            for (int point = 0; point < Npoint_aer; point++)
              if (with_indices_list)
                Concentration_point(point) = Model.GetConcentration_aer()
                  (base_s, base_b, point_list[point][0],
                   point_list[point][1], point_list[point][2]);
              else
                Concentration_point(point) =
                  Model.GetConcentration_aer(base_s, base_b,
                                             coord_list[point][0],
                                             coord_list[point][1],
                                             coord_list[point][2]);
            FormatBinary<float>().Append(Concentration_point,
                                         output_file[s][b]);
          }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_AER_CXX
#endif
