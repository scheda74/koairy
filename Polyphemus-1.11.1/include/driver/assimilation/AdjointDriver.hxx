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


#ifndef POLYPHEMUS_FILE_DRIVER_ADJOINTDRIVER_HXX


#include "Trajectory.cxx"
#include "AssimilationDriver.cxx"
#include <map>


namespace Polyphemus
{


  ///////////////////
  // ADJOINTDRIVER //
  ///////////////////


  /*! This driver checks sensitivity of model output with respect to model
    initial conditions using adjoint model.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  class AdjointDriver:
    public AssimilationDriver < T, ClassModel, ClassOutputSaver,
                                ClassObsManager >
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

    /*! Species name of the selected point in model domain for sensitivity
      calculation. */
    string point_species_name;
    //! x-index of the selected point for sensitivity calculation.
    int point_nx;
    //! y-index of the selected point for sensitivity calculation.
    int point_ny;
    //! z-index of the selected point for sensitivity calculation.
    int point_nz;

    /*** Initial perturbation settings for finite difference calculation ***/

    //! With random directions for the perturbation?
    bool with_random_perturbation;
    //! Norm of the initial perturbation vector.
    T norm_scale;

    //! Decreasing ratio of the sequence of perturbation vectors.
    T decreasing_root;
    //! Index for the calculation of the first decreasing ratio.
    int start_index;
    //! Index for the calculation of the last decreasing ratio.
    int end_index;

    /*** Checking options ***/

    //! Checking left-side finite difference results?
    bool with_left_finite_difference_checking;
    //! Display sensitivity results for the decreasing perturbation sequences?
    bool display_sensitivity;

  public:

    /*** Constructor and destructor ***/

    AdjointDriver(string config_file);
    virtual ~AdjointDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ADJOINTDRIVER_HXX
#endif
