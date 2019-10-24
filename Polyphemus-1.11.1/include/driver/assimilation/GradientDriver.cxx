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


#ifndef POLYPHEMUS_FILE_DRIVER_GRADIENTDRIVER_CXX


#include "GradientDriver.hxx"
#include "newran.h"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  GradientDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::GradientDriver(string config_file):
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
    (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  GradientDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~GradientDriver()
  {
  }


  //! Reads configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void GradientDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    // Reads "gradient" section.
    this->config.SetSection("[adjoint]");

    this->config.PeekValue("Display_precision", display_precision);
    this->config.PeekValue("With_trajectory_management",
                           with_trajectory_management);
    if (with_trajectory_management)
      {
        this->config.PeekValue("Trajectory_delta_t", Trajectory_delta_t);
        this->config.PeekValue("Trajectory_file", Trajectory_file_name);
      }

    this->config.PeekValue("Norm_perturbation_vector", norm_scale);
    this->config.PeekValue("With_random_perturbation",
                           with_random_perturbation);

    this->config.PeekValue("Decreasing_root", decreasing_root);
    this->config.PeekValue("Start_index", start_index);
    if (start_index < 0) start_index = 0;
    this->config.PeekValue("End_index", end_index);

    this->config.PeekValue("With_left_finite_difference_checking",
                           with_left_finite_difference_checking);
    this->config.PeekValue("Display_cost", display_cost);
  }


  //! Driver initialization.
  /*! It reads configurations */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void GradientDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  //! Checks the gradient calculation using adjoint model.
  /*! The cost function is chosen to be the sum of least squares of the
    difference between model simulations and synthetic observations, therefore
    no specifications of error covariance matrices, i.e; B or R, are needed.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void GradientDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int t, s, z, y, x;


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    this->Model.Init();
    this->ObsManager.Init(this->Model);
    Init();

    if (with_trajectory_management)
      TrajManager.Init(this->Model.GetConcentration().GetArray().shape(),
                       this->Model.GetDate_min(),
                       double(Trajectory_delta_t),
                       Trajectory_file_name,
                       true);
    else
      TrajData.Resize(this->Model.GetNt() + 1,
                      this->Model.GetNs(), this->Model.GetNz(),
                      this->Model.GetNy(), this->Model.GetNx());

    // Data array that stores synthetic observations.
    Data<T, 5> Observation(this->Model.GetNt() + 1,
                           this->Model.GetNs(), this->Model.GetNz(),
                           this->Model.GetNy(), this->Model.GetNx());

    // Sets output precision.
    cout.precision(streamsize(display_precision));


    ////////////////////////////////
    // GENERATION OF OBSERVATIONS //
    ////////////////////////////////


    cout << "  Generation of observations ..." << endl;
    for (t = 0; t < this->Model.GetNt(); t++)
      {
        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
               << endl;

        Array<T, 4> Obs(&Observation.GetArray()(t, 0, 0, 0, 0),
                        shape(this->Model.GetNs(), this->Model.GetNz(),
                              this->Model.GetNy(), this->Model.GetNx()),
                        neverDeleteData);

        Obs = this->Model.GetConcentration().GetArray();

        this->Model.InitStep();
        this->Model.Forward();
      }
    if (this->option_display["show_date"])
      cout << "Current date: " <<
        this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

    Array<T, 4> LastObs(&Observation.GetArray()(this->Model.GetNt(),
                                                0, 0, 0, 0),
                        shape(this->Model.GetNs(), this->Model.GetNz(),
                              this->Model.GetNy(), this->Model.GetNx()),
                        neverDeleteData);
    LastObs = this->Model.GetConcentration().GetArray();


    //////////////////
    // FORWARD LOOP //
    //////////////////


    // In the forward loop, the simulation starts from first guess; the state
    // trajectories are stored for backward loop; and the cost function is
    // calculated.
    //
    // The first guess is set to be the average of all synthetic observations,
    // and the cost function is simply the sum of squares of differences
    // between synthetic observations and model simulations starting from the
    // first guess.


    // Sets observation availabilities for synthetic observations.
    Array<int, 1> obs_availability(this->Model.GetNt() + 1);
    obs_availability = 1;

    // Sets first guess.
    Data<T, 4> conc_guess(this->Model.GetNs(), this->Model.GetNz(),
                          this->Model.GetNy(), this->Model.GetNx());

    for (s = 0; s < this->Model.GetNs(); s++)
      for (z = 0; z < this->Model.GetNz(); z++)
        for (y = 0; y < this->Model.GetNy(); y++)
          for (x = 0; x < this->Model.GetNx(); x++)
            {
              T tmp = 0;
              for (t = 0; t < this->Model.GetNt() + 1; t++)
                tmp += Observation(t, s, z, y, x);
              tmp /= this->Model.GetNt() + 1;

              conc_guess(s, z, y, x) = tmp;
            }
    this->Model.GetConcentration().GetArray() = conc_guess.GetArray();

    /*** Integrations ***/

    this->Model.SetDate(this->Model.GetDate_min());

    cout << "  Forward loop ..." << endl;
    T cost = 0.;
    for (t = 0; t < this->Model.GetNt(); t++)
      {
        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        if (with_trajectory_management)
          TrajManager.Append(this->Model.GetConcentration(),
                             this->Model.GetCurrentDate());
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                             shape(this->Model.GetNs(), this->Model.GetNz(),
                                   this->Model.GetNy(), this->Model.GetNx()),
                             neverDeleteData);

            Traj = this->Model.GetConcentration().GetArray();
          }

        // Calculates cost functions.
        for (s = 0; s < this->Model.GetNs(); s++)
          for (z = 0; z < this->Model.GetNz(); z++)
            for (y = 0; y < this->Model.GetNy(); y++)
              for (x = 0; x < this->Model.GetNx(); x++)
                cost += obs_availability(t) *
                  pow(double(Observation(t, s, z, y, x) -
                             this->Model.GetConcentration()(s, z, y, x)), 2.0)
                  / 2.0;

        this->Model.InitStep();
        this->Model.Forward();
      }

    if (this->option_display["show_date"])
      cout << "Current date: " <<
        this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

    if (with_trajectory_management)
      TrajManager.Append(this->Model.GetConcentration(),
                         this->Model.GetCurrentDate());
    else
      {
        Array<T, 4> LastTraj(&TrajData.
                             GetArray()(this->Model.GetNt(), 0, 0, 0, 0),
                             shape(this->Model.GetNs(),
                                   this->Model.GetNz(),
                                   this->Model.GetNy(),
                                   this->Model.GetNx()),
                             neverDeleteData);

        LastTraj = this->Model.GetConcentration().GetArray();
      }

    // Calculates cost functions.
    for (s = 0; s < this->Model.GetNs(); s++)
      for (z = 0; z < this->Model.GetNz(); z++)
        for (y = 0; y < this->Model.GetNy(); y++)
          for (x = 0; x < this->Model.GetNx(); x++)
            cost += obs_availability(this->Model.GetNt()) *
              pow(double(Observation(this->Model.GetNt(), s, z, y, x)
                         -  this->Model
                         .GetConcentration()(s, z, y, x)), 2.0)
              / 2.0;

    // Displays reference cost function values.
    cout << "  Cost : " << cost << endl;


    ///////////////////
    // Backward loop //
    ///////////////////


    if (with_trajectory_management)
      TrajManager.InitLoad();

    this->Model.SetBackward(true);
    this->Model.GetConcentration_ccl().GetArray() = T(0.);

    Date cur_date = this->Model.GetDate_min();
    cur_date.AddSeconds(this->Model.GetNt() * this->Model.GetDelta_t());

    cout << "  Backward loop ..." << endl;
    for (t = this->Model.GetNt(); t > 0; t--)
      {

        // Loads state values of the last model intergration results.
        if (with_trajectory_management)
          TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                             shape(this->Model.GetNs(), this->Model.GetNz(),
                                   this->Model.GetNy(), this->Model.GetNx()),
                             neverDeleteData);

            this->Model.GetConcentration().GetArray() = Traj;
          }

        // Adds forcing terms.
        for (s = 0; s < this->Model.GetNs(); s++)
          for (z = 0; z < this->Model.GetNz(); z++)
            for (y = 0; y < this->Model.GetNy(); y++)
              for (x = 0; x < this->Model.GetNx(); x++)
                this->Model.GetConcentration_ccl()(s, z, y, x) +=
                  obs_availability(t) *
                  (this->Model.GetConcentration()(s, z, y, x)
                   - Observation(t, s, z, y, x));

        // Loads state values before the last model intergration.
        cur_date.AddSeconds(-this->Model.GetDelta_t());
        if (with_trajectory_management)
          TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t - 1, 0, 0, 0, 0),
                             shape(this->Model.GetNs(), this->Model.GetNz(),
                                   this->Model.GetNy(), this->Model.GetNx()),
                             neverDeleteData);

            this->Model.GetConcentration().GetArray() = Traj;
          }

        this->Model.SetDate(cur_date);
        this->Model.InitStep();

        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        this->Model.Backward();

        this->Model.SubtractTime(this->Model.GetDelta_t());
        cur_date = this->Model.GetCurrentDate();
      }

    // Treatments of observations at intial time.
    cout << "Forcing item for initial condition. " << endl;

    if (with_trajectory_management)
      TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
    else
      {
        Array<T, 4> Traj(&TrajData.GetArray()(0, 0, 0, 0, 0),
                         shape(this->Model.GetNs(), this->Model.GetNz(),
                               this->Model.GetNy(), this->Model.GetNx()),
                         neverDeleteData);

        this->Model.GetConcentration().GetArray() = Traj;
      }

    for (s = 0; s < this->Model.GetNs(); s++)
      for (z = 0; z < this->Model.GetNz(); z++)
        for (y = 0; y < this->Model.GetNy(); y++)
          for (x = 0; x < this->Model.GetNx(); x++)
            this->Model.GetConcentration_ccl()(s, z, y, x) +=
              obs_availability(0) *
              (this->Model.GetConcentration()(s, z, y, x)
               - Observation(0, s, z, y, x));

    // Saves gradients.
    Data<T, 4> grad(this->Model.GetNs(), this->Model.GetNz(),
                    this->Model.GetNy(), this->Model.GetNx());
    grad.GetArray() = this->Model.GetConcentration_ccl().GetArray();


    //////////////////////
    // CHECKS GRADIENTS //
    //////////////////////


    /*** Generates initial perturbation vector ***/

    Data<T, 4> delta_h, pert_h, conc_pert;
    delta_h.Resize(this->Model.GetNs(), this->Model.GetNz(),
                   this->Model.GetNy(), this->Model.GetNx());
    pert_h.Resize(this->Model.GetNs(), this->Model.GetNz(),
                  this->Model.GetNy(), this->Model.GetNx());
    conc_pert.Resize(this->Model.GetNs(), this->Model.GetNz(),
                     this->Model.GetNy(), this->Model.GetNx());

    srand((unsigned long)time(0));
    double seed = rand() / double(RAND_MAX);
    NEWRAN::LGM_mixed* urng = new NEWRAN::LGM_mixed(seed);
    NEWRAN::Random::Set(*urng);
    // Gaussian generator with 0 mean and unitary variance.
    NEWRAN::Normal N1;

    if (with_random_perturbation)
      for (s = 0; s < this->Model.GetNs(); s++)
        for (z = 0; z < this->Model.GetNz(); z++)
          for (y = 0; y < this->Model.GetNy(); y++)
            for (x = 0; x < this->Model.GetNx(); x++)
              delta_h(s, z, y, x) = N1.Next()
                * conc_guess.GetArray()(s, z, y, x);
    else
      delta_h.GetArray() = conc_guess.GetArray();

    T delta_norm = delta_h.Norm2();
    delta_h.GetArray() = norm_scale / delta_norm * delta_h.GetArray();

    // Temporary variables.
    Array<T, 1> ratio_right(end_index - start_index);
    Array<T, 1> cost_right(end_index - start_index);
    Array<T, 1> ratio_left(end_index - start_index);
    Array<T, 1> cost_left(end_index - start_index);

    // Forward mode.
    this->Model.SetBackward(false);

    cout << "  Checks gradients ..." << endl;
    for (int i = start_index; i < end_index; i++)
      {
        // Decreases perturbations.
        T alpha = pow(double(decreasing_root), double(-i));
        pert_h.GetArray() = alpha * delta_h.GetArray();

        /*** right-side finite difference ***/

        Date cur_date = this->Model.GetDate_min();
        this->Model.SetDate(cur_date);

        // Calculates sum of squares of differences between observations and
        // simulations with respect to perturbed initial concentrations.
        conc_pert.GetArray() = pert_h.GetArray() + conc_guess.GetArray();
        this->Model.GetConcentration().GetArray() = conc_pert.GetArray();

        cost_right(i - start_index) = 0.;
        for (t = 0; t < this->Model.GetNt(); t++)
          {
            for (s = 0; s < this->Model.GetNs(); s++)
              for (z = 0; z < this->Model.GetNz(); z++)
                for (y = 0; y < this->Model.GetNy(); y++)
                  for (x = 0; x < this->Model.GetNx(); x++)
                    cost_right(i - start_index) += obs_availability(t) *
                      (pow(double(Observation(t, s, z, y, x) -
                                  this->Model.
                                  GetConcentration()(s, z, y, x)), 2.0))
                      / 2.0;

            this->Model.InitStep();
            this->Model.Forward();
          }

        for (s = 0; s < this->Model.GetNs(); s++)
          for (z = 0; z < this->Model.GetNz(); z++)
            for (y = 0; y < this->Model.GetNy(); y++)
              for (x = 0; x < this->Model.GetNx(); x++)
                cost_right(i - start_index)
                  += obs_availability(this->Model.GetNt()) *
                  (pow(double(Observation(this->Model.GetNt(), s, z, y, x) -
                              this->Model.
                              GetConcentration()(s, z, y, x)), 2.0))
                  / 2.0;

        // Calculates inner product of gradient and initial perturbation
        // vector.
        T tmp = T(0.);
        for (s = 0; s < this->Model.GetNs(); s++)
          for (z = 0; z < this->Model.GetNz(); z++)
            for (y = 0; y < this->Model.GetNy(); y++)
              for (x = 0; x < this->Model.GetNx(); x++)
                tmp += grad(s, z, y, x) * delta_h(s, z, y, x);

        // Ratio of the results obtained by finite difference and by adjoint
        // model.
        ratio_right(i - start_index)
          = (cost_right(i - start_index) - cost) / (alpha * tmp);

        /*** left-side finite difference ***/

        if (with_left_finite_difference_checking)
          {

            cur_date = this->Model.GetDate_min();
            this->Model.SetDate(cur_date);
            conc_pert.GetArray() = conc_guess.GetArray() - pert_h.GetArray();
            this->Model.GetConcentration().GetArray() = conc_pert.GetArray();

            cost_left(i - start_index) = 0.;
            for (t = 0; t < this->Model.GetNt(); t++)
              {
                for (s = 0; s < this->Model.GetNs(); s++)
                  for (z = 0; z < this->Model.GetNz(); z++)
                    for (y = 0; y < this->Model.GetNy(); y++)
                      for (x = 0; x < this->Model.GetNx(); x++)
                        cost_left(i - start_index) += obs_availability(t) *
                          (pow(double(Observation(t, s, z, y, x) -
                                      this->Model.
                                      GetConcentration()(s, z, y, x)), 2.0))
                          / 2.0;

                this->Model.InitStep();
                this->Model.Forward();
              }

            for (s = 0; s < this->Model.GetNs(); s++)
              for (z = 0; z < this->Model.GetNz(); z++)
                for (y = 0; y < this->Model.GetNy(); y++)
                  for (x = 0; x < this->Model.GetNx(); x++)
                    cost_left(i - start_index)
                      += obs_availability(this->Model.GetNt()) *
                      (pow(double(Observation(this->Model.GetNt(), s, z, y, x)
                                  - this->Model.
                                  GetConcentration()(s, z, y, x)), 2.0))
                      / 2.0;

            ratio_left(i - start_index)
              = (cost - cost_left(i - start_index)) / (alpha * tmp);
          }
      }

    /*** Display checking results ***/

    cout << "  Gradient check results : " << endl;

    if (display_cost)
      {
        cout << "  Right-side cost function values :" << endl;
        for (int i = 0; i < end_index - start_index; i++)
          cout << cost_right(i) << endl;
        if (with_left_finite_difference_checking)
          {
            cout << "  Left-side cost function values :" << endl;
            for (int i = 0; i < end_index - start_index; i++)
              cout << cost_left(i) << endl;
          }
      }

    cout << "  Right-side ratio values :" << endl;
    for (int i = 0; i < end_index - start_index; i++)
      cout << ratio_right(i) << endl;

    if (with_left_finite_difference_checking)
      {
        cout << "  Left-side ratio values :" << endl;
        for (int i = 0; i < end_index - start_index; i++)
          cout << ratio_left(i) << endl;
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_GRADIENTDRIVER_CXX
#endif
