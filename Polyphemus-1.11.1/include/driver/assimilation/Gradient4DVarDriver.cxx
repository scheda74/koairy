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


#ifndef POLYPHEMUS_FILE_DRIVER_GRADIENT4DVARDRIVER_CXX


#include "Gradient4DVarDriver.hxx"
#include "newran.h"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  Gradient4DVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Gradient4DVarDriver(string config_file):
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
    (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  Gradient4DVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~Gradient4DVarDriver()
  {
  }


  //! Reads configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void Gradient4DVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
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
  void Gradient4DVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  //! Checks the gradient calculation using adjoint model.
  /*! The cost function is chosen to be the 4D-Var one.
   */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void Gradient4DVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int t;


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

    // Sets output precision.
    cout.precision(display_precision);


    //////////////////
    // FORWARD LOOP //
    //////////////////


    // In the forward loop, the simulation starts from first guess; the state
    // trajectories are stored for backward loop, and the cost function is
    // calculated.
    //
    // The first guess is set to initial conditions.

    // Sets first guess.
    Data<T, 4> conc_guess(this->Model.GetNs(), this->Model.GetNz(),
                          this->Model.GetNy(), this->Model.GetNx());
    Data<T, 1> conc_guess_vec;

    // "With_initial_condition" in main configuration file should be
    // "yes".
    conc_guess_vec.Resize(this->Model.GetNstate());
    this->Model.GetState(conc_guess_vec.GetArray());

    // Saves all concentration array.
    conc_guess.GetArray() = this->Model.GetConcentration().GetArray();

    /*** Integrations ***/

    T cost = 0.;
    cout << "  Forward loop ..." << endl;
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

        // Calculates cost function.
        this->ObsManager.SetDate(this->Model.GetCurrentDate());
        if (this->ObsManager.IsAvailable())
          {
            int Nobs = this->ObsManager.GetNobs();
            int Nstate = this->Model.GetNstate();
            Array<T, 1> departure(Nobs);
            Array<T, 1> state_vector(Nstate);

            this->Model.GetState(state_vector);
            for (int r = 0; r < Nobs; r++)
              {
                departure(r) = this->ObsManager.GetObservation()(r) -
                  this->ObsManager.MltObsOperator(r, state_vector);
                cost += 1 / this->ObsManager.GetCovariance()(r, r) *
                  departure(r) * departure(r) / 2.0;
              }
          }

        this->Model.InitStep();
        this->Model.Forward();
      }

    // Saves state vector after the last forward intergration.
    if (with_trajectory_management)
      TrajManager.Append(this->Model.GetConcentration(),
                         this->Model.GetCurrentDate());
    else
      {
        Array<T, 4> LastTraj(&TrajData.
                             GetArray()(this->Model.GetNt(), 0, 0, 0, 0),
                             shape(this->Model.GetNs(), this->Model.GetNz(),
                                   this->Model.GetNy(), this->Model.GetNx()),
                             neverDeleteData);

        LastTraj = this->Model.GetConcentration().GetArray();
      }

    // Calculates cost functions.
    this->ObsManager.SetDate(this->Model.GetCurrentDate());
    if (this->ObsManager.IsAvailable())
      {
        int Nobs = this->ObsManager.GetNobs();
        int Nstate = this->Model.GetNstate();
        Array<T, 1> departure(Nobs);
        Array<T, 1> state_vector(Nstate);

        this->Model.GetState(state_vector);
        for (int r = 0; r < Nobs; r++)
          {
            departure(r) = this->ObsManager.GetObservation()(r) -
              this->ObsManager.MltObsOperator(r, state_vector);
            cost += 1 / this->ObsManager.GetCovariance()(r, r) *
              departure(r) * departure(r) / 2.0;
          }
      }

    // Displays reference cost function values.
    cout << "  Cost : " << cost << endl;


    ///////////////////
    // BACKWARD LOOP //
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
        // Loads state vector of the last model integration results.
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
        this->ObsManager.SetDate(cur_date);
        if (this->ObsManager.IsAvailable())
          {
            int Nobs = this->ObsManager.GetNobs();
            int Nstate = this->Model.GetNstate();
            Array<T, 1> normailized_departure(Nobs);
            Array<T, 1> state_vector(Nstate);
            Array<T, 1> state_ccl(Nstate);

            this->Model.GetState(state_vector);
            this->Model.GetState_ccl(state_ccl);
            for (int r = 0; r < Nobs; r++)
              normailized_departure(r) =
                (this->ObsManager.GetObservation()(r) - this->ObsManager.
                 MltObsOperator(r, state_vector)) /
                this->ObsManager.GetCovariance()(r, r);
            for (int i = 0; i < Nstate; i++)
              state_ccl(i) += - this->ObsManager.
                ADJMltObsOperator(i, normailized_departure);

            this->Model.SetState_ccl(state_ccl);
          }

        // Loads state vector before the last model intergration.
        cur_date.AddSeconds(-this->Model.GetDelta_t());
        if (with_trajectory_management)
          TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t - 1, 0, 0, 0, 0),
                             shape(this->Model.GetNs(),
                                   this->Model.GetNz(),
                                   this->Model.GetNy(),
                                   this->Model.GetNx()),
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

    this->ObsManager.SetDate(cur_date);
    if (this->ObsManager.IsAvailable())
      {
        int Nobs = this->ObsManager.GetNobs();
        int Nstate = this->Model.GetNstate();
        Array<T, 1> normailized_departure(Nobs);
        Array<T, 1> state_vector(Nstate);
        Array<T, 1> state_ccl(Nstate);

        this->Model.GetState(state_vector);
        this->Model.GetState_ccl(state_ccl);
        for (int r = 0; r < Nobs; r++)
          normailized_departure(r) =
            (this->ObsManager.GetObservation()(r) - this->ObsManager.
             MltObsOperator(r, state_vector)) /
            this->ObsManager.GetCovariance()(r, r);
        for (int i = 0; i < Nstate; i++)
          state_ccl(i) += - this->ObsManager.
            ADJMltObsOperator(i, normailized_departure);

        this->Model.SetState_ccl(state_ccl);
      }

    // Saves gradients.
    Data<T, 1> grad(this->Model.GetNstate());
    this->Model.GetState_ccl(grad.GetArray());


    //////////////////////
    // CHECKS GRADIENTS //
    //////////////////////


    /*** Generates initial perturbation vector ***/

    Data<T, 1> delta_h_vec, pert_h_vec, conc_pert_vec;
    delta_h_vec.Resize(this->Model.GetNstate());
    pert_h_vec.Resize(this->Model.GetNstate());
    conc_pert_vec.Resize(this->Model.GetNstate());

    srand((unsigned long)time(0));
    double seed = rand() / double(RAND_MAX);
    NEWRAN::LGM_mixed* urng = new NEWRAN::LGM_mixed(seed);
    NEWRAN::Random::Set(*urng);
    // Gaussian generator with 0 mean and unitary variance.
    NEWRAN::Normal N1;

    if (with_random_perturbation)
      for (int i = 0; i < this->Model.GetNstate(); i++)
        delta_h_vec(i) = N1.Next() * conc_guess_vec(i);
    else
      delta_h_vec.GetArray() = conc_guess_vec.GetArray();

    T delta_norm = delta_h_vec.Norm2();
    delta_h_vec.GetArray() = norm_scale / delta_norm * delta_h_vec.GetArray();

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
        pert_h_vec.GetArray() = alpha * delta_h_vec.GetArray();

        /*** right-side finite difference ***/

        Date cur_date = this->Model.GetDate_min();
        this->Model.SetDate(cur_date);
        this->Model.GetConcentration().GetArray() = conc_guess.GetArray();

        conc_pert_vec.GetArray()
          = pert_h_vec.GetArray() + conc_guess_vec.GetArray();
        this->Model.SetState(conc_pert_vec.GetArray());

        // Calculates cost functions.
        cost_right(i - start_index) = 0.;
        for (t = 0; t < this->Model.GetNt(); t++)
          {
            this->ObsManager.SetDate(this->Model.GetCurrentDate());
            if (this->ObsManager.IsAvailable())
              {
                int Nobs = this->ObsManager.GetNobs();
                int Nstate = this->Model.GetNstate();
                Array<T, 1> departure(Nobs);
                Array<T, 1> state_vector(Nstate);

                this->Model.GetState(state_vector);
                for (int r = 0; r < Nobs; r++)
                  {
                    departure(r) = this->ObsManager.GetObservation()(r) -
                      this->ObsManager.MltObsOperator(r, state_vector);
                    cost_right(i - start_index)
                      += 1 / this->ObsManager.GetCovariance()(r, r) *
                      departure(r) * departure(r) / 2.0;
                  }
              }

            this->Model.InitStep();
            this->Model.Forward();
          }

        // Calculates the last term of the cost function.
        this->ObsManager.SetDate(this->Model.GetCurrentDate());
        if (this->ObsManager.IsAvailable())
          {
            int Nobs = this->ObsManager.GetNobs();
            int Nstate = this->Model.GetNstate();
            Array<T, 1> departure(Nobs);
            Array<T, 1> state_vector(Nstate);

            this->Model.GetState(state_vector);
            for (int r = 0; r < Nobs; r++)
              {
                departure(r) = this->ObsManager.GetObservation()(r) -
                  this->ObsManager.MltObsOperator(r, state_vector);
                cost_right(i - start_index)
                  += 1 / this->ObsManager.GetCovariance()(r, r) *
                  departure(r) * departure(r) / 2.0;
              }
          }

        // Calculates inner product of gradient and initial perturbation
        // vector.
        T tmp = T(0.);
        for (int j = 0; j < this->Model.GetNstate(); j++)
          tmp += grad(j) * delta_h_vec(j);

        // Ratio of the results obtained by finite difference and by adjoint
        // model.
        ratio_right(i - start_index)
          = (cost_right(i - start_index) - cost) / (alpha * tmp);

        /*** left-side finite difference ***/

        if (with_left_finite_difference_checking)
          {

            cur_date = this->Model.GetDate_min();
            this->Model.SetDate(cur_date);
            this->Model.GetConcentration().GetArray() = conc_guess.GetArray();

            conc_pert_vec.GetArray()
              = conc_guess_vec.GetArray() - pert_h_vec.GetArray();
            this->Model.SetState(conc_pert_vec.GetArray());

            // Calculates cost function.
            cost_left(i - start_index) = 0.;
            for (t = 0; t < this->Model.GetNt(); t++)
              {
                this->ObsManager.SetDate(this->Model.GetCurrentDate());
                if (this->ObsManager.IsAvailable())
                  {
                    int Nobs = this->ObsManager.GetNobs();
                    int Nstate = this->Model.GetNstate();
                    Array<T, 1> departure(Nobs);
                    Array<T, 1> state_vector(Nstate);

                    this->Model.GetState(state_vector);
                    for (int r = 0; r < Nobs; r++)
                      {
                        departure(r) = this->ObsManager.GetObservation()(r) -
                          this->ObsManager.MltObsOperator(r, state_vector);
                        cost_left(i - start_index)
                          += 1 / this->ObsManager.GetCovariance()(r, r) *
                          departure(r) * departure(r) / 2.0;
                      }
                  }

                this->Model.InitStep();
                this->Model.Forward();
              }

            // Calculates the last term of the cost function.
            this->ObsManager.SetDate(this->Model.GetCurrentDate());
            if (this->ObsManager.IsAvailable())
              {
                int Nobs = this->ObsManager.GetNobs();
                int Nstate = this->Model.GetNstate();
                Array<T, 1> departure(Nobs);
                Array<T, 1> state_vector(Nstate);

                this->Model.GetState(state_vector);
                for (int r = 0; r < Nobs; r++)
                  {
                    departure(r) = this->ObsManager.GetObservation()(r) -
                      this->ObsManager.MltObsOperator(r, state_vector);
                    cost_left(i - start_index)
                      += 1 / this->ObsManager.GetCovariance()(r, r) *
                      departure(r) * departure(r) / 2.0;
                  }
              }

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


#define POLYPHEMUS_FILE_DRIVER_GRADIENT4DVARDRIVER_CXX
#endif
