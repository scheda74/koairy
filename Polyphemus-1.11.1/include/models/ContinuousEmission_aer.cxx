// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_AER_CXX


#include "ContinuousEmission_aer.hxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  ContinuousEmission_aer<T>::ContinuousEmission_aer():
    ContinuousEmission<T>()
  {
  }


  //! Model initialization.
  /*! It reads the configuration. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Rate_aer: the list of release rates for aerosol species(mass unit per seconds),
    </ul>
    \param config The ConfigStream instance containing the sources parameters.
  */
  template<class T>
  void ContinuousEmission_aer<T>::Init(ConfigStream& config,
                                       vector<string> species_list,
                                       vector<string> species_list_aer,
                                       int Nbin_aer)
  {
    ContinuousEmission<T>::Init(config, species_list);
    PointEmissionUnit<T>::Init(config, species_list, species_list_aer, Nbin_aer);

    config.Find("Rate_aer");
    vector<string> tmp = split(config.GetLine());
    int Nrate = tmp.size();
    if (Nrate != this->Ns_emis_aer)
      throw string("Error in \"ContinuousEmission_aer<T>::Init()\": ")
        + " there must be one rate per emitted species.";
    for (int i = 0; i < Nrate; i++)
      rate_aer.push_back(to_num<T>(tmp[i]));
  }

  //! Gets emitted quantity for one species during a given time interval.
  /*!
    \param date_beg beginning date of the current time interval
    \param date_end ending date of the current time interval
    \param s species index
    \param point_emission a 2D array of size (1, 4) containing the source
    coordinates (x, y, z) and the species quantity.
  */
  template<class T>
  void ContinuousEmission_aer<T>::GetEmission_aer(Date current_date,
                                                  Date next_date, int s,
                                                  T& quantity)
  {
    T delta_t = next_date.GetSecondsFrom(current_date);
    quantity = delta_t
      - max(0., this->date_beg.GetSecondsFrom(current_date))
      - max(0., next_date.GetSecondsFrom(this->date_end));
    quantity *= rate_aer[s];
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_AER_CXX
#endif
