// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_CXX


#include "TemporalEmission.hxx"
#include "PointEmissionUnit.cxx"
#include "ContinuousEmission.cxx"


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
  void TemporalEmission<T>::Init(ConfigStream& config,
                                 vector<string> species_list)
  {

    config.PeekValue("TemporalFactor", input_file);
    config.PeekValue("Delta_t", "positive", Delta_t_file);
    date_min_file = config.PeekValue("Date_min_file");
    current_record = -1;
    TemporalFactor.resize(1);

    ContinuousEmission<T>::Init(config, species_list);

  }


  //! Gets emitted quantities for one species.
  /*!
    \return The emitted quantities for one species.
  */
  template<class T>
  void TemporalEmission<T>::GetEmission(Date current_date,
                                        Date next_date, int s,
                                        Array<T, 2>& point_emission)
  {
    ContinuousEmission<T>::GetEmission(current_date, next_date,
                                       s, point_emission);

    // Reading temporal factors.
    T time_distance = current_date.GetSecondsFrom(date_min_file);
    int record = int(floor(time_distance / Delta_t_file));
    if (record != current_record)
      {
        if (is_num(input_file))
          {
            T value = to_num<T>(input_file);
            TemporalFactor(0) = value;
          }
        else
          FormatBinary<float>().ReadRecord(input_file, record,
                                           TemporalFactor);
        current_record = record;
      }
    // Computing quantity.
    point_emission(0, 3) *= TemporalFactor(0);
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_CXX
#endif
