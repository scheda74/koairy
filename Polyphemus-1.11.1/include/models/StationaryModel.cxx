// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_MODELS_STATIONARYMODEL_CXX


#include "StationaryModel.hxx"
#include "Norms.hxx"

namespace Polyphemus
{


  //! Main constructor.
  /*! Builds the model and the underlying model.
    \param config_file configuration file.
  */
  template<class T, class ClassModel>
  StationaryModel<T, ClassModel>
  ::StationaryModel(string config_file):
    BaseModel<T>(config_file), Model(config_file)
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  StationaryModel<T, ClassModel>::~StationaryModel()
  {
  }


  //! Model initialization.
  /*! It reads the configuration and allocates memory.
   */
  template<class T, class ClassModel>
  void StationaryModel<T, ClassModel>::Init()
  {
    Model.Init();
    BaseModel<T>::Init();

    /*** Stationary conditions ***/

    this->config.SetSection("[stationary]");
    this->config.GetValue("Delta_t", "> 0", this->Delta_t);
    this->config.GetValue("Nt", "> 0", this->Nt);

    this->config.SetSection("[convergence]");
    this->config.GetValue("Norm", "one | two | infinity", norm);
    this->config.GetValue("Epsilon", epsilon);
    this->config.GetValue("Method", "mean | max", method);
  }


  //! Model initialization for each step.
  /*! It sets the underlying model at the right date and reads on file the
    data that is needed for the current step.
  */
  template<class T, class ClassModel>
  void StationaryModel<T, ClassModel>::InitStep()
  {
    Model.SetDate(this->GetCurrentDate());
    Model.InitStep();
  }


  //! Performs one step forward.
  /*!
   */
  template<class T, class ClassModel>
  void StationaryModel<T, ClassModel>::Forward()
  {
    int n = 0;

    convergence = false;
    /*** Time inner-loop ***/
    for (int i = 0; i < Model.GetNt(); i++)
      {
        if (!convergence)
          {
            Conc_prev.Copy(Model.GetConcentration());
            Model.Forward();
            Conc_current.Copy(Model.GetConcentration());

            // Goes back.
            Model.StepBack();

            // Initializes the model.
            this->InitStep();

            // Computation of convergence criterion.
            if (norm == "one")
              convergence = CheckConvergenceOne(Conc_prev, Conc_current,
                                                epsilon, method);
            else if (norm == "two")
              convergence = CheckConvergenceTwo(Conc_prev, Conc_current,
                                                epsilon, method);
            else
              convergence = CheckConvergenceInfinity(Conc_prev, Conc_current,
                                                     epsilon, method);
            n++;
          }
      }
    if (n < Model.GetNt())
      cout << "Convergence reached after " << n << " sub-iterations." << endl;
    else
      cout << "The inner-loop did not converge." << endl;

    this->Concentration.Copy(Model.GetConcentration());
    this->AddTime(this->Delta_t);
    this->step++;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STATIONARYMODEL_CXX
#endif
