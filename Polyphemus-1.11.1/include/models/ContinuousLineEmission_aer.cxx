// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_AER_CXX


#include "ContinuousLineEmission_aer.hxx"


namespace Polyphemus
{


  ////////////////////
  // LINE EMISSIONS //
  ////////////////////


  //! Main constructor.
  template<class T>
  ContinuousLineEmission_aer<T>::ContinuousLineEmission_aer():
    ContinuousLineEmission<T>()
  {
  }


  // //! Model initialization.
  // /*! It reads the configuration. Each source is described in a
  //   dedicated section "[source]" in which one finds the following entries:
  //   <ul>
  //   <li> Date beg: the release date,
  //   <li> Rate: the list of release rates (mass unit per seconds),
  //   <li> Species: list of species names,
  //   <li> Coordinate_file: name of the file containing coordinates
  //   of line sources. In this file, each line contains 8 elements: ID number of
  //   node 1, coordinates (x, y, z) of node 1, ID number of node 2 and
  //   coordinates of node 2.
  //   </ul>
  //   \param config The ConfigStream instance containing the sources parameters.
  //   \note The cross-sectional surface of the source is assumed to be a disc.
  // */
  template<class T>
  void ContinuousLineEmission_aer<T>::Init(ConfigStream& config,
                                           vector<string> species_list,
                                           vector<string> species_list_aer,
                                           int Nbin_aer)
  {
    //   // Reading the emssion file
    //   config.PeekValue("Diameter", this->diameter_);

    ContinuousLineEmission<T>::Init(config, species_list);
    //   date_end = config.PeekValue("Date_end");
    PointEmissionUnit<T>::Init(config, species_list, species_list_aer, Nbin_aer);
    config.Find("Rate_aer");
    vector<string> tmp = split(config.GetLine());
    int Nrate = tmp.size();
    if (Nrate != this->Ns_emis_aer)
      throw string("Error in \"ContinuousLineEmission<T>::Init()\": ")
        + " there must be one rate per emitted species.";

    for (int i = 0; i < Nrate; i++)
      rate_aer.push_back(to_num<T>(tmp[i]));

    //   string coordinate_file;
    //   T Id, x1, y1, z1, x2, y2, z2;
    //   config.PeekValue("Coordinate_file", coordinate_file);

    //   // Reading the coordinate file
    //   ConfigStreams coordinate_stream(coordinate_file);
    //   while(!coordinate_stream.IsEmpty())
    //     {
    //       Array<T,1> coordinate(6);
    //       coordinate_stream.GetNumber(Id);
    //       coordinate_stream.GetNumber(x1);
    //       coordinate_stream.GetNumber(y1);
    //       coordinate_stream.GetNumber(z1);
    //       coordinate_stream.GetNumber(Id);
    //       coordinate_stream.GetNumber(x2);
    //       coordinate_stream.GetNumber(y2);
    //       coordinate_stream.GetNumber(z2);
    //      coordinate = x1, y1, z1, x2, y2, z2;
    //        coordinate_list_.push_back(coordinate);
    //     }
  }


  // //! Checks if the emission has begun at the given date.
  // /*!
  //   \return True if the emission has begun at date date.
  // */
  // template<class T>
  // bool ContinuousLineEmission<T>::HasEnded(Date date)
  // {
  //   return date >= date_end;
  // }


  // //! Gets the emission coordinates.
  // /*!
  //   Gets the emission coordinates.
  //   \param coordinate_list: list of coordinate of line sources.
  // */
  // template<class T>
  // void ContinuousLineEmission<T>::GetCoordinates(list<Array<T, 1> >&
  //                                                coordinate_list) const
  // {
  //   coordinate_list = coordinate_list_;
  // }


  // //! Gets the species rate.
  // /*!
  //   \return The species rate.
  // */
  // template<class T>
  // T ContinuousLineEmission<T>::GetRate(int species) const
  // {
  //   return rate[species];
  // }


  // //! Gets the species diameter.
  // /*!
  //   \return The species diameter.
  // */
  // template<class T>
  // T ContinuousLineEmission<T>::GetDiameter() const
  // {
  //   return this->diameter_;
  // }


} // namespace Polyphemus.

#define POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_AER_CXX
#endif
