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


#ifndef POLYPHEMUS_FILE_DRIVER_FOURDIMVARDRIVER_HXX


#include "Trajectory.cxx"
#include "VariationalDriver.cxx"
#include "Bfgs.cxx"
#include <map>


namespace Polyphemus
{


  //////////////////////
  // FOURDIMVARDRIVER //
  //////////////////////


  //! This driver performs 4D-Var.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  class FourDimVarDriver:
    public VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  {

  protected:

    //! Display precision for sensitivity checking results.
    int display_precision;
    //! With trajectory managements?
    bool with_trajectory_management;

    //! Trajectory manager.
    Trajectory<double, 4> TrajManager;
    //! Trajectory time step in seconds.
    T Trajectory_delta_t;
    //! Trajectory number of time steps.
    string Trajectory_file_name;
    //! Data array that stores trajectories.
    Data<T, 5> TrajData;

    //! Pointer to optimization solver.
    BaseOptimizer<T>* Optimizer;
    //! Type of optimization solver.
    string type;
    //! Maximal iteration for optimization solver.
    int max_iter;
    //! Display iterations during optimizing?
    bool disp_iter;

    //! Is background term taken into account in the cost function?
    bool with_background_term;

    //! Name of the file that saves Jb values during optimization.
    string Jb_file;
    //! Name of the file that saves Jo values during optimization.
    string Jo_file;
    //! Name of the file that saves gradient norms during optimization.
    string grad_norm_file;

    /*! \brief Should the inverse of the background error covariance matrix be
      read from disk? */
    bool read_inverse_background_matrix;
    /*! \brief Name of the file that stores the inverse of the background
      error covariance matrix. */
    string file_background_inverse_matrix;
    /*! \brief Name of the file that stores the background error covariance
      matrix. */
    string file_background_matrix;

  public:

    /*** Constructor and destructor ***/

    FourDimVarDriver(string config_file);
    virtual ~FourDimVarDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    void ComputeInverseBgCovMatrix(Array<T, 2>& matrix);

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_FOURDIMVARDRIVER_HXX
#endif
