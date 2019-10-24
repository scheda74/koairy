// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_CXX


#include "ContinuousLineEmission.hxx"
#include "PointEmissionUnit.cxx"


namespace Polyphemus
{


  ////////////////////
  // LINE EMISSIONS //
  ////////////////////


  //! Main constructor.
  template<class T>
  ContinuousLineEmission<T>::ContinuousLineEmission()
  {
  }


  //! Destructor.
  template<class T>
  ContinuousLineEmission<T>::~ContinuousLineEmission()
  {
    // Nothing.
  }


  //! Model initialization.
  /*! It reads the configuration. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Date beg: the release date,
    <li> Rate: the list of release rates (mass unit per seconds),
    <li> Species: list of species names,
    <li> Coordinate_file: name of the file containing coordinates
    of line sources. In this file, each line contains 8 elements: ID number of
    node 1, coordinates (x, y, z) of node 1, ID number of node 2 and
    coordinates of node 2.
    </ul>
    \param config The ConfigStream instance containing the sources parameters.
    \note The cross-sectional surface of the source is assumed to be a disc.
  */
  template<class T>
  void ContinuousLineEmission<T>::Init(ConfigStream& config,
                                       vector<string> species_list)
  {
    date_end = config.PeekValue("Date_end");
    PointEmissionUnit<T>::Init(config, species_list);
    config.PeekValue("Diameter", this->diameter_);
    config.PeekValue("Velocity", this->velocity_);
    config.PeekValue("Temperature", this->temperature_);

    string coordinate_file;
    T Id, x1, y1, z1, x2, y2, z2;
    config.PeekValue("Coordinate_file", coordinate_file);

    // Reading the coordinate file
    int id_section;
    ConfigStreams coordinate_stream(coordinate_file);
    while (!coordinate_stream.IsEmpty())
      {
        Array<T, 1> section_coordinate(6);
        vector<T> section_rate;
        vector<string> param = split(coordinate_stream.GetLine());

        if (int(param.size()) != 12 + this->Ns_emis)
          throw "Error in 'ContinuousLineEmission<T>::Init()': "
            "for each source section, there must be one rate per emitted "
            "species.";

        section_coordinate = to_num<T>(param[1]), to_num<T>(param[2]),
          to_num<T>(param[3]), to_num<T>(param[5]), to_num<T>(param[6]),
          to_num<T>(param[7]);

        Width.push_back(to_num<T>(param[8]));
        VehicleVelocity.push_back(to_num<T>(param[9]));
        Area.push_back(to_num<T>(param[10]));
        Density.push_back(to_num<T>(param[11]));

        for (int i = 0; i < this->Ns_emis; ++i)
          section_rate.push_back(to_num<T>(param[12 + i]));

        rate.push_back(section_rate);
        coordinate_list_.push_back(section_coordinate);
        ++id_section;
      }
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool ContinuousLineEmission<T>::HasEnded(Date date)
  {
    return date >= date_end;
  }


  //! Returns the ending date of emission.
  /*!
    \return The ending date of emission.
  */
  template<class T>
  Date ContinuousLineEmission<T>::GetEndDate() const
  {
    return date_end;
  }


  //! Gets the emission coordinates.
  /*!
    Gets the emission coordinates.
    \param coordinate_list: list of coordinate of line sources.
  */
  template<class T>
  void ContinuousLineEmission<T>::GetCoordinates(list<Array<T, 1> >&
                                                 coordinate_list) const
  {
    coordinate_list = coordinate_list_;
  }


  //! Gets the species rate for a given source section.
  /*!
    \return The species rate.
  */
  template<class T>
  T ContinuousLineEmission<T>::GetRate(int species, int id_section) const
  {
    return rate[id_section][species];
  }


  //! Gets the source width.
  /*!
    \return The source width.
  */
  template<class T>
  T ContinuousLineEmission<T>::GetWidth(int id_section) const
  {
    return Width[id_section];
  }


  //! Sets the emission coordinates.
  /*!
    \param[in] coordinate_list.
  */
  template<class T>
  void ContinuousLineEmission<T>::SetCoordinates(const list<Array<T, 1> >&
                                                 coordinate_list)
  {
    coordinate_list_ = coordinate_list;
  }


  //! Gets the source parameters needed for plume rise computation.
  /*!
    \param[out] velocity the source rejection rate (m/s).
    \param[out] temperature the source temperature (K).
    \param[out] diameter the source diameter (m).
  */
  template<class T>
  void ContinuousLineEmission<T>::GetPlumeRiseParam(T& velocity,
                                                    T& temperature,
                                                    T& diameter)
  {
    velocity = this->velocity_;
    temperature = this->temperature_;
    diameter = this->diameter_;
  }


  //! Gets the source VehicleVelocity.
  /*!
    \return The source VehicleVelocity.
  */
  template<class T>
  T ContinuousLineEmission<T>::GetVehicleVelocity(int id_section) const
  {
    return VehicleVelocity[id_section];
  }

  //! Gets the source Area.
  /*!
    \return The source Area.
  */
  template<class T>
  T ContinuousLineEmission<T>::GetArea(int id_section) const
  {
    return Area[id_section];
  }

  //! Gets the source Density.
  /*!
    \return The source Density.
  */
  template<class T>
  T ContinuousLineEmission<T>::GetDensity(int id_section) const
  {
    return Density[id_section];
  }
} // namespace Polyphemus.

#define POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_CXX
#endif
