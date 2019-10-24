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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_CXX


#include "SaverUnitPoint.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  ////////////////////
  // SAVERUNITPOINT //
  ////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitPoint<T, ClassModel>
  ::SaverUnitPoint():
    BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitPoint<T, ClassModel>::~SaverUnitPoint()
  {
  }


  //! Type of saver.
  /*!
    \return The string "indices_list" or "coordinates_list".
  */
  template<class T, class ClassModel>
  string SaverUnitPoint<T, ClassModel>::GetType()  const
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
  void SaverUnitPoint<T, ClassModel>::Init(ConfigStream& config_stream,
                                           ClassModel& Model)
  {
    config_stream.PeekValue("Type", type);
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Vertical levels to be saved.
    if (type == "coordinates_list")
      {
        with_indices_list = 0;
        // This option is temporally removed.
        // We need to check if it is useful for Gaussian models (YK/YR).
        // config_stream.PeekValue("Spatial_average", spatial_average);
        // if (spatial_average)
        //   {
        //     config_stream.PeekValue("Averaging_length", averaging_length);
        //     config_stream.PeekValue("Averaging_height", averaging_height);
        //   } 
        config_stream.Find("Levels_coordinates");
        split(config_stream.GetLine(), levels_coord);
        Nlevels = int(levels_coord.size());
      }
    else if (type == "indices_list")
      {
        config_stream.Find("Levels");
        split(config_stream.GetLine(), levels);
        Nlevels = int(levels.size());
        with_indices_list = 1;
      }
    else
      throw string("Saver unit \"") + type +  "\" is not a point list.";

    // Species to be saved.
    this->Ns = int(this->species_list.size());
    for (int s = 0; s < this->Ns; s++)
      this->species_index.push_back(Model.GetSpeciesIndex
                                    (this->species_list[s]));

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(this->Ns);
    for (unsigned int i = 0; i < this->species_list.size(); i++)
      output_file[i] = find_replace(filename, "&f", this->species_list[i]);

    // Empties output files.
    for (unsigned int s = 0; s < this->species_list.size(); s++)
      ofstream tmp_stream(output_file[s].c_str());

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
        Npoint = int(point_list.size());
        string point_file = config_stream.GetValue("Point_file");
        ofstream points(point_file.c_str());
        if (!points.is_open())
          throw string("Unable to open file \"") + point_file + "\".";
        for (int i = 0; i < Npoint; i++)
          points << point_list[i][0] << "\t" << point_list[i][1] << "\t"
                 << point_list[i][2] << endl;
        points.close();

        // Initializing the buffer used to compute averaged concentrations.
        if (this->averaged)
          {
            int s, point;
            Concentration_.Resize(this->Ns, Npoint);
            for (s = 0; s < this->Ns; s++)
              for (point = 0; point < Npoint; point++)
                Concentration_(s, point) = 0.5
                  * Model.GetConcentration()(this->species_index[s],
                                             point_list[point][0],
                                             point_list[point][1],
                                             point_list[point][2]);
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
                + " per line, but " + to_str(split(line).size()) +
                " numbers found.";
          }

        // Writing the list of points to be saved.
        Npoint = int(coord_list.size());
        string point_file = config_stream.GetValue("Point_file");
        ofstream points(point_file.c_str());
        if (!points.is_open())
          throw string("Unable to open file \"") + point_file + "\".";
        for (int i = 0; i < Npoint; i++)
          points << coord_list[i][0] << "\t" << coord_list[i][1] << "\t"
                 << coord_list[i][2] << endl;
        points.close();

        // Initializing the buffer used to compute averaged concentrations.
        if (this->averaged)
          {
            int s, point;
            Concentration_.Resize(this->Ns, Npoint);
            for (s = 0; s < this->Ns; s++)
              for (point = 0; point < Npoint; point++)
                Concentration_(s, point) = 0.5
                  * Model.GetConcentration(this->species_index[s],
                                           coord_list[point][0],
                                           coord_list[point][1],
                                           coord_list[point][2]);
          }
      }

    if (this->initial_concentration && !this->averaged)
      this->Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitPoint<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetCurrentDate()
    <li> GetConcentration(int, T, T, T)
    <li> GetIntegratedConcentration(int, T, T, T, T, T, T)
    <li> ComputeConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitPoint<T, ClassModel>::Save(ClassModel& Model)
  {
    if (this->averaged)
      {
        // Useful for models where concentrations are not automatically
        // computed on the grid. Otherwise, nothing is performed.
        if (with_indices_list)
          Model.ComputeConcentration();

        if (this->counter % this->interval_length == 0)
          {

            int s, point;
            for (s = 0; s < this->Ns; s++)
              for (point = 0; point < Npoint; point++)
                {
                  if (with_indices_list)
                    Concentration_(s, point) += 0.5
                      * Model.GetConcentration()(this->species_index[s],
                                                 point_list[point][0],
                                                 point_list[point][1],
                                                 point_list[point][2]);
                  else
                    Concentration_(s, point) += 0.5
                      * Model.GetConcentration(this->species_index[s],
                                               coord_list[point][0],
                                               coord_list[point][1],
                                               coord_list[point][2]);
                }

            Concentration_.GetArray() /= T(this->interval_length);

            if (Model.GetCurrentDate() >= this->date_beg
                && Model.GetCurrentDate() <= this->date_end)
              for (s = 0; s < this->Ns; s++)
                {
                  Data<T, 1>
                    Concentration_tmp(&Concentration_(s, 0), shape(Npoint));
                  FormatBinary<float>().Append(Concentration_tmp,
                                               output_file[s]);
                }
            for (s = 0; s < this->Ns; s++)
              for (point = 0; point < Npoint; point++)
                {
                  if (with_indices_list)
                    Concentration_(s, point) = 0.5
                      * Model.GetConcentration()(this->species_index[s],
                                                 point_list[point][0],
                                                 point_list[point][1],
                                                 point_list[point][2]);
                  else
                    Concentration_(s, point) = 0.5
                      * Model.GetConcentration(this->species_index[s],
                                               coord_list[point][0],
                                               coord_list[point][1],
                                               coord_list[point][2]);
                }
            this->counter = 0;
          }
        else
          {
            int s, point;
            for (s = 0; s < this->Ns; s++)
              for (point = 0; point < Npoint; point++)
                {
                  if (with_indices_list)
                    Concentration_(s, point) +=
                      Model.GetConcentration()(this->species_index[s],
                                               point_list[point][0],
                                               point_list[point][1],
                                               point_list[point][2]);
                  else
                    Concentration_(s, point) +=
                      Model.GetConcentration(this->species_index[s],
                                             coord_list[point][0],
                                             coord_list[point][1],
                                             coord_list[point][2]);
                }
          }
      }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      {
        // Useful for models where concentrations are not automatically
        // computed on the grid. Otherwise, nothing is performed.
        if (with_indices_list)
          Model.ComputeConcentration();

        // Instantaneous concentrations.
        for (int s = 0; s < this->Ns; s++)
          {
            Data<T, 1> Concentration_point(Npoint);
            for (int point = 0; point < Npoint; point++)
              if (with_indices_list)
                Concentration_point(point) = Model.GetConcentration()
                  (this->species_index[s], point_list[point][0],
                   point_list[point][1], point_list[point][2]);
              else
                // Computing concentration a given coordinates.
                // if (!spatial_average)
                  Concentration_point(point) =
                    Model.GetConcentration(this->species_index[s],
                                           coord_list[point][0],
                                           coord_list[point][1],
                                           coord_list[point][2]);
                // else
                //   Concentration_point(point) =
                //     Model.GetIntegratedConcentration(this->species_index[s],
                //                                      coord_list[point][0],
                //                                      coord_list[point][1],
                //                                      coord_list[point][2],
                //                                      averaging_height,
                //                                      averaging_length,
                //                                      averaging_length);
            FormatBinary<float>().Append(Concentration_point, output_file[s]);
          }
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_CXX
#endif
