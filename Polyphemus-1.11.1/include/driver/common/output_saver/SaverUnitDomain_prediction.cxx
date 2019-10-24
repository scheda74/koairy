// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_PREDICTION_CXX


#include "SaverUnitDomain_prediction.hxx"

#include "SaverUnitDomain.cxx"


namespace Polyphemus
{


  ////////////////////////////////
  // SAVERUNITDOMAIN_PREDICTION //
  ////////////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitDomain_prediction<T, ClassModel>
  ::SaverUnitDomain_prediction(): SaverUnitDomain<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitDomain_prediction<T, ClassModel>::~SaverUnitDomain_prediction()
  {
  }

  //! First initialization.
  /*! Reads the configuration.
    \param config_stream configuration stream.
    \param Model model.
  */
  template<class T, class ClassModel>
  void SaverUnitDomain_prediction<T, ClassModel>
  ::Init(ConfigStream& config_stream, ClassModel& Model)
  {
    SaverUnitDomain<T, ClassModel>::Init(config_stream, Model);

    // Date file name.
    date_file = config_stream.GetValue("Date_file");
    // Empties file.
    ofstream(date_file.c_str()).close();
  }

  //! Type of saver.
  /*!
    \return The string "domain_prediction".
  */
  template<class T, class ClassModel>
  string SaverUnitDomain_prediction<T, ClassModel>::GetType()  const
  {
    return "domain_prediction";
  }


  //! Group of the saver unit.
  /*!
    \return The group of the saver unit, that is, "prediction".
  */
  template<class T, class ClassModel>
  string SaverUnitDomain_prediction<T, ClassModel>::GetGroup()  const
  {
    return "prediction";
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitDomain_prediction<T, ClassModel>::Save(ClassModel& Model)
  {
    SaverUnitDomain<T, ClassModel>::Save(Model);

    if (this->counter % this->interval_length == 0)
      {
        ofstream out_file(date_file.c_str(), ios::app);
        out_file << Model.GetCurrentDate().GetDate("%y-%m-%d_%hh%i %s\n");
        out_file.close();
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_PREDICTION_CXX
#endif
