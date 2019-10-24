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


#ifndef POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_AER_CXX


#include "TemporalEmission_aer.hxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Model initialization.
  /*! It reads the configuration and allocates memory. In addition to the
    parameters needed for a continuous source, the following parameters must
    be provided:
    <ul>
    <li> TemporalFactor: binary file or value of the temporal factors.
    <li> Delta_t: time step of the temporal factors in the binary file.
    <li> Date_min_file: starting date of the temporal factors in the binary
    file.
    </ul>
  */
  template<class T>
  void TemporalEmission_aer<T>::Init(ConfigStream& config,
                                     vector<string> species_list,
                                     vector<string> species_list_aer,
                                     int Nbin_aer)
  {

    TemporalEmission<T>::Init(config, species_list);
    PointEmissionUnit<T>::Init(config, species_list, species_list_aer, Nbin_aer);

    config.Find("Rate_aer");
    vector<string> tmp = split(config.GetLine());
    int Nrate = tmp.size();
    if (Nrate != this->Ns_emis_aer)
      throw string("Error in \"TemporalEmission_aer<T>::Init()\": ")
        + " there must be one rate per emitted species.";
    for (int i = 0; i < Nrate; i++)
      rate_aer.push_back(to_num<T>(tmp[i]));

  }


  //! Gets emitted quantities for one species.
  /*!
    \return The emitted quantities for one species.
  */
  template<class T>
  void TemporalEmission_aer<T>::GetEmission_aer(Date current_date,
                                                Date next_date, int s,
                                                T& quantity)
  {

    T delta_t = next_date.GetSecondsFrom(current_date);
    quantity = delta_t
      - max(0., this->date_beg.GetSecondsFrom(current_date))
      - max(0., next_date.GetSecondsFrom(this->date_end));
    quantity *= rate_aer[s];

    // Reading temporal factors.
    T time_distance = current_date.GetSecondsFrom(this->date_min_file);
    int record = int(floor(time_distance / this->Delta_t_file));
    if (record != this->current_record)
      {
        if (is_num(this->input_file))
          {
            T value = to_num<T>(this->input_file);
            this->TemporalFactor(0) = value;
          }
        else
          FormatBinary<float>().ReadRecord(this->input_file, record,
                                           this->TemporalFactor);
        this->current_record = record;
      }
    // Computing quantity.
    quantity *= this->TemporalFactor(0);
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_AER_CXX
#endif
