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


#ifndef POLYPHEMUS_FILE_DRIVER_GRADIENTDRIVER_HXX


#include "Trajectory.cxx"
#include "AssimilationDriver.cxx"
#include <map>


namespace Polyphemus
{


  ////////////////////
  // GRADIENTDRIVER //
  ////////////////////


  /*! \brief This driver checks gradient calculation using adjoint model. The
    cost function is chosen to be the sum of least squares of the difference
    between model simulations and synthetic observations.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  class GradientDriver:
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
    /*! \brief Display cost function values for the decreasing perturbation
      sequences? */
    bool display_cost;

  public:

    /*** Constructor and destructor ***/

    GradientDriver(string config_file);
    virtual ~GradientDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_GRADIENTDRIVER_HXX
#endif
