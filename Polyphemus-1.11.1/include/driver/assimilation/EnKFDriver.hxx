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


#ifndef POLYPHEMUS_FILE_DRIVER_ENKFDRIVER_HXX


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


  ////////////////
  // ENKFDRIVER //
  ////////////////


  /*! \brief This driver performs assimilation with the ensemble Kalman
    filter.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  class EnKFDriver:
    public SequentialDriver < T, ClassModel, ClassOutputSaver,
                              ClassObsManager >
  {

  protected:

    //! The number of samples in the ensemble.
    int Nensemble;
    //! Flag that indicates whether observations are perturbed.
    bool with_observation_perturbation;
    //! Flag that indicates whether square root scheme is employed.
    bool with_square_root;
    //! Flag that indicates whether advance sampling is employed.
    bool with_advance_sampling;
    //! Flag indicates whether ensemble prediction is supported.
    bool with_ensemble_prediction;

    /*! Name list of the files that saves coresponding entire concentration
      data for each member in the ensemble.
    */
    vector<string> conc_file_list;

  public:

    /*** Constructor and destructor ***/

    EnKFDriver(string config_file);
    virtual ~EnKFDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    virtual void Run();

    void InitEnsemble(Array<T, 1>& state_vector, Array<T, 2>& Ensemble);
    void Forecast(Array<T, 2>& Ensemble);
    void Analyze(Array<T, 2>& ObsPerturbation, Array<T, 2>& Ensemble);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ENKFDRIVER_HXX
#endif
