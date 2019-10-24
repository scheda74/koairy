// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_DRIVER_RRSQRTDRIVER_HXX


#include <map>
#include "SequentialDriver.cxx"
#include "AlgebraFunctions.cxx"


/////////////////////////////
// BLAS & LAPACK INTERFACE //

extern "C"
{
#include "cblas.h"
#include "clapack.h"
}

// BLAS & LAPACK INTERFACE //
/////////////////////////////


namespace Polyphemus
{


  //////////////////
  // RRSQRTDRIVER //
  //////////////////


  /*! \brief This driver performs assimilation with the reduced rank square
    root (RRSQRT) Kalman filter.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  class RRSQRTDriver:
    public SequentialDriver < T, ClassModel, ClassOutputSaver,
                              ClassObsManager >
  {

  protected:

    //! The number of analysis modes.
    int Nmode_max;
    //! The number of model modes.
    int Nmode_model;
    //! The number of observation modes.
    int Nmode_observation;
    //! Option for mode matrix propagation (forecast step).
    string propagation_option;
    //! Perturbation value for finite differences.
    T finite_difference_perturbation;

  public:

    /*** Constructor and destructor ***/

    RRSQRTDriver(string config_file);
    virtual ~RRSQRTDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    virtual void Run();

    void Forecast(Array<T, 1>& state_vector, Array<T, 2>& ModeMatrix);
    void Analyze(Array<T, 1>& state_vector, Array<T, 2>& ModeMatrix);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_RRSQRTDRIVER_HXX
#endif
