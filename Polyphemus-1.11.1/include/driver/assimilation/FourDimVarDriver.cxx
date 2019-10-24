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


#ifndef POLYPHEMUS_FILE_DRIVER_FOURDIMVARDRIVER_CXX


#include "FourDimVarDriver.hxx"


//////////////////////
// LAPACK INTERFACE //
//////////////////////


extern "C"
{
#include "clapack.h"
}

namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::FourDimVarDriver(string config_file):
    VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
    (config_file),
    Optimizer(NULL)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~FourDimVarDriver()
  {
    if (Optimizer != NULL)
      delete Optimizer;
  }


  //! Reads configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    // Reads "4DVar" section.
    this->config.SetSection("[4DVar]");

    this->config.PeekValue("Display_precision", display_precision);
    this->config.PeekValue("With_trajectory_management",
                           with_trajectory_management);
    if (with_trajectory_management)
      {
        this->config.PeekValue("Trajectory_delta_t", Trajectory_delta_t);
        this->config.PeekValue("Trajectory_file", Trajectory_file_name);
      }

    this->config.PeekValue("With_background_term", with_background_term);

    if (with_background_term)
      {
        this->config.PeekValue("Jb_file", Jb_file);
        this->config.PeekValue("Read_inverse_background_matrix",
                               read_inverse_background_matrix);
        this->config.PeekValue("File_background_inverse_matrix",
                               file_background_inverse_matrix);
        this->config.PeekValue("File_background_matrix",
                               file_background_matrix);
      }

    this->config.PeekValue("Jo_file", Jo_file);
    this->config.PeekValue("Gradient_norm_file", grad_norm_file);

    // Reads "optimizer" section.
    this->config.SetSection("[optimizer]");
    this->config.PeekValue("Type", type);
    this->config.PeekValue("Maximal_iteration", max_iter);
    this->config.PeekValue("Display_iterations", disp_iter);
  }


  //! Driver initialization.
  /*! It reads configurations */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  //! Computes the inversion of the background error covariance matrix.
  /*!
    \param matrix (input/output) on entry, the matrix with size equal to the
    background error covariance matrix; on exit, the inversion of the
    background error covariance matrix.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ComputeInverseBgCovMatrix(Array<T, 2>& matrix)
  {
    // Computes B.
    for (int r = 0; r < this->Nstate; r++)
      matrix(r, Range::all()) = this->Model.BackgroundErrorCovariance(r);

    if (file_background_matrix != "empty")
      FormatBinary<double>().Write(matrix, file_background_matrix);

    // Computes B^{-1}.
    Array<int, 1> permutations(this->Nstate);
    int info;
    dgetrf_(&this->Nstate, &this->Nstate, matrix.data(), &this->Nstate,
            permutations.data(), &info);
    if (info != 0)
      throw Error("FourDimVarDriver::ComputeInverseBgCovMatrix",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGETRF.");

    int lwork = this->Nstate * 64;
    Array<T, 1> work(lwork);
    dgetri_(&this->Nstate, matrix.data(), &this->Nstate,
            permutations.data(), work.data(), &lwork, &info);
    if (info != 0)
      throw Error("FourDimVarDriver::ComputeInverseBgCovMatrix",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGETRI.");
  }


  //! Performs four-dimensional variational assimilations.
  /*! It minimizes the 4D-Var cost function, the gradient with respect to
    initial conditions is calculated using adjoint model. */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void FourDimVarDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int t;


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    this->Model.Init();
    this->ObsManager.Init(this->Model);
    Init();

    // Gets inverse of background error covariance matrix that is not
    // diagonal.
    string background_error_model = this->Model.GetErrorCovarianceModel();
    bool is_background_diagonal
      = background_error_model == "diagonal_constant";

    Array<T, 2> B_inv;
    if (with_background_term)
      {
        // Gets the inversion of background error covariance matrix.
        if (!is_background_diagonal)
          {
            B_inv.resize(this->Nstate, this->Nstate);

            if (read_inverse_background_matrix)
              FormatBinary<double>().
                Read(file_background_inverse_matrix, B_inv);

            // Balgovind parameterization.
            if (background_error_model == "Balgovind" &&
                !read_inverse_background_matrix)
              {
                ComputeInverseBgCovMatrix(B_inv);
                if (file_background_inverse_matrix != "empty")
                  FormatBinary<double>().
                    Write(B_inv, file_background_inverse_matrix);
              }
          }
      }


    // Sets output precision.
    cout.precision(display_precision);


    ///////////////////////
    // 4DVAR APPLICATION //
    ///////////////////////


    // The first guess is set to initial condition (background), and the cost
    // function is chosen to be the classical 4D-Var one with B set to
    // diagonal matrix.

    // Saves first guess.
    Data<T, 4> conc_guess(this->Model.GetNs(), this->Model.GetNz(),
                          this->Model.GetNy(), this->Model.GetNx());
    conc_guess.GetArray() = this->Model.GetConcentration().GetArray();
    Array<T, 1> state_guess(this->Nstate);
    this->Model.GetState(state_guess);

    // Optimizer initializations.
    if (type == "BFGS")
      Optimizer = new Bfgs<T>(this->config.GetFileName(),
                              this->Nstate, disp_iter);

    // Sets initial parameter values, and boundary conditions for optimization
    // solver.
    this->Model.GetState(Optimizer->param);

    // For calling Bfgs member functions.
    Bfgs<T>* pOpt = (Bfgs<T> *)Optimizer;
    pOpt->SetBoundary(0., 10000., 2);

    Optimizer->Init();

    // Arrays for optimization results.
    Array<T, 1> Jb(max_iter), Jo(max_iter);
    Array<T, 1> grad_norm2(max_iter);


    //////////////////
    // ASSIMILATION //
    //////////////////


    int opt_iter = 0;
    cout << "  Start optimizing..." << endl;
    while (!Optimizer->IsStop())
      {
        // Trajectory managements during iterations.
        if (with_trajectory_management)
          TrajManager.Init(this->Model.GetConcentration().GetArray().shape(),
                           this->Model.GetDate_min(),
                           double(Trajectory_delta_t),
                           Trajectory_file_name,
                           true);
        else
          TrajData.Resize(this->Nt_assim + 1,
                          this->Model.GetNs(), this->Model.GetNz(),
                          this->Model.GetNy(), this->Model.GetNx());


        //////////////////
        // FORWARD LOOP //
        //////////////////


        // Reintializes date and intial concentration array.
        this->Model.SetDate(this->Model.GetDate_min());
        this->Model.GetConcentration().GetArray() = conc_guess.GetArray();
        // Sets state variables to iterative values.
        this->Model.SetState(Optimizer->param);

        /*** Integrations ***/

        T cost = 0.;
        cout << "  Forward loop ..." << endl;
        for (t = 0; t < this->Nt_assim; t++)
          {
            if (this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            if (with_trajectory_management)
              TrajManager.Append(this->Model.GetConcentration(),
                                 this->Model.GetCurrentDate());
            else
              {
                Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                                 shape(this->Model.GetNs(),
                                       this->Model.GetNz(),
                                       this->Model.GetNy(),
                                       this->Model.GetNx()),
                                 neverDeleteData);

                Traj = this->Model.GetConcentration().GetArray();
              }

            // Calculates cost function.
            this->ObsManager.SetDate(this->Model.GetCurrentDate());
            if (this->ObsManager.IsAvailable())
              {
                int Nobs = this->ObsManager.GetNobs();
                Array<T, 1> departure(Nobs);
                Array<T, 1> state_vector(this->Nstate);

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
                                 GetArray()(this->Nt_assim, 0, 0, 0, 0),
                                 shape(this->Model.GetNs(),
                                       this->Model.GetNz(),
                                       this->Model.GetNy(),
                                       this->Model.GetNx()),
                                 neverDeleteData);

            LastTraj = this->Model.GetConcentration().GetArray();
          }

        // Calculates cost functions.
        this->ObsManager.SetDate(this->Model.GetCurrentDate());
        if (this->ObsManager.IsAvailable())
          {
            int Nobs = this->ObsManager.GetNobs();
            Array<T, 1> departure(Nobs);
            Array<T, 1> state_vector(this->Nstate);

            this->Model.GetState(state_vector);
            for (int r = 0; r < Nobs; r++)
              {
                departure(r) = this->ObsManager.GetObservation()(r) -
                  this->ObsManager.MltObsOperator(r, state_vector);
                cost += 1 / this->ObsManager.GetCovariance()(r, r) *
                  departure(r) * departure(r) / 2.0;
              }
          }

        Jo(opt_iter) = cost;

        // Calculates J_b item.
        if (with_background_term)
          {
            int i, j;
            Array<T, 1> state_increment(this->Nstate);
            for (i = 0; i < this->Nstate; i++)
              state_increment(i) = Optimizer->param(i) - state_guess(i);

            if (is_background_diagonal)
              {
                int Nstate_s = this->Model.GetStateNz() * this->Model.GetNy()
                  * this->Model.GetNx();

                for (i = 0; i < this->Nstate; i++)
                  cost += state_increment(i) * state_increment(i) /
                    this->Model.GetBackgroundErrVar()(i / Nstate_s) / 2.0;
              }
            else
              {
                T sum = T(0.);
                for (j = 0; j < this->Nstate; j++)
                  {
                    T tmp = T(0.);
                    for (i = 0; i < this->Nstate; i++)
                      tmp += B_inv(i, j) * state_increment(i);

                    sum += tmp * state_increment(j);
                  }
                cost += sum;
              }

            Jb(opt_iter) = cost - Jo(opt_iter);
          }

        // Sets cost function values.
        Optimizer->cost = double(cost);


        ///////////////////
        // BACKWARD LOOP //
        ///////////////////


        if (with_trajectory_management)
          TrajManager.InitLoad();

        this->Model.SetBackward(true);
        this->Model.GetConcentration_ccl().GetArray() = T(0.);

        Date cur_date = this->Model.GetDate_min();
        cur_date.AddSeconds(this->Nt_assim * this->Model.GetDelta_t());

        cout << "  Backward loop ..." << endl;
        for (t = this->Nt_assim; t > 0; t--)
          {

            // Loads state vector of the last model integration results.
            if (with_trajectory_management)
              TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
            else
              {
                Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                                 shape(this->Model.GetNs(),
                                       this->Model.GetNz(),
                                       this->Model.GetNy(),
                                       this->Model.GetNx()),
                                 neverDeleteData);

                this->Model.GetConcentration().GetArray() = Traj;
              }


            // Adds forcing terms.
            this->ObsManager.SetDate(cur_date);
            if (this->ObsManager.IsAvailable())
              {
                int Nobs = this->ObsManager.GetNobs();
                Array<T, 1> normalized_departure(Nobs);
                Array<T, 1> state_vector(this->Nstate);
                Array<T, 1> state_ccl(this->Nstate);

                this->Model.GetState(state_vector);
                this->Model.GetState_ccl(state_ccl);
                for (int r = 0; r < Nobs; r++)
                  normalized_departure(r) =
                    (this->ObsManager.GetObservation()(r) - this->ObsManager.
                     MltObsOperator(r, state_vector)) /
                    this->ObsManager.GetCovariance()(r, r);
                for (int i = 0; i < this->Nstate; i++)
                  state_ccl(i) += - this->ObsManager.
                    ADJMltObsOperator(i, normalized_departure);

                this->Model.SetState_ccl(state_ccl);
              }

            // Loads state vector before the last model integration.
            cur_date.AddSeconds(-this->Model.GetDelta_t());
            if (with_trajectory_management)
              TrajManager.Load(cur_date, this->Model.GetConcentration(),
                               true);
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
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            this->Model.Backward();

            this->Model.SubtractTime(this->Model.GetDelta_t());
            cur_date = this->Model.GetCurrentDate();
          }

        // Treatments of observations at initial time.
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
            Array<T, 1> normalized_departure(Nobs);
            Array<T, 1> state_vector(this->Nstate);
            Array<T, 1> state_ccl(this->Nstate);

            this->Model.GetState(state_vector);
            this->Model.GetState_ccl(state_ccl);
            for (int r = 0; r < Nobs; r++)
              normalized_departure(r) =
                (this->ObsManager.GetObservation()(r) - this->ObsManager.
                 MltObsOperator(r, state_vector)) /
                this->ObsManager.GetCovariance()(r, r);
            for (int i = 0; i < this->Nstate; i++)
              state_ccl(i) += - this->ObsManager.
                ADJMltObsOperator(i, normalized_departure);

            this->Model.SetState_ccl(state_ccl);
          }

        // Sets gradients for Jo item.
        this->Model.GetState_ccl(Optimizer->grad);

        // Adds gradients for Jb item.
        if (with_background_term)
          {
            int i, j;
            Array<T, 1> state_increment(this->Nstate);
            for (i = 0; i < this->Nstate; i++)
              state_increment(i) = Optimizer->param(i) - state_guess(i);

            if (is_background_diagonal)
              {
                int Nstate_s = this->Model.GetStateNz() * this->Model.GetNy()
                  * this->Model.GetNx();

                for (i = 0; i < this->Nstate; i++)
                  Optimizer->grad(i) += state_increment(i) /
                    this->Model.GetBackgroundErrVar()(i / Nstate_s);
              }
            else
              {
                for (i = 0; i < this->Nstate; i++)
                  {
                    T sum = T(0.);
                    for (j = 0; j < this->Nstate; j++)
                      sum += B_inv(i, j) * state_increment(j);

                    Optimizer->grad(i) += sum;
                  }
              }
          }
        // Calculates the norm of gradient.
        grad_norm2(opt_iter) = T(0.);
        for (int i = 0; i < this->Nstate; i++)
          grad_norm2(opt_iter) += Optimizer->grad(i) * Optimizer->grad(i);
        grad_norm2(opt_iter) = sqrt(grad_norm2(opt_iter));

        if (!Optimizer->Optimize())
          break;

        if (opt_iter >= max_iter - 1)
          {
            cout << "  Maximal iteration number reached, stops optimization."
                 << endl;
            break;
          }

        opt_iter++;
      }

    // Optimization ending information.
    cout << "  Return information of optimization solver : " << pOpt->task
         << endl;

    // Records number of iterations for the optimization process.
    max_iter = opt_iter + 1;

    // Records optimization results.
    if (with_background_term)
      FormatBinary<T>().Write(Jb, Jb_file);
    FormatBinary<T>().Write(Jo, Jo_file);
    FormatBinary<T>().Write(grad_norm2, grad_norm_file);

    // Saves assimilation results.
    cout << "  Saving assimilation results ..." << endl;

    this->Model.SetDate(this->Model.GetDate_min());
    this->Model.GetConcentration().GetArray() = conc_guess.GetArray();
    // Set state variables to optimized values.
    this->Model.SetState(Optimizer->param);

    this->OutputSaver.Init(this->Model);

    for (t = 0; t < this->Nt_assim; t++)
      {
        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        this->Model.InitStep();
        this->OutputSaver.InitStep(this->Model);

        this->Model.Forward();

        this->OutputSaver.SetGroup("analysis");
        this->OutputSaver.Save(this->Model);

        this->OutputSaver.SetGroup("forecast");
        this->OutputSaver.Save(this->Model);
      }


    ////////////////
    // PREDICTION //
    ////////////////


    if (this->Nt_predict > 0)
      cout << "  Predicting ..." << endl;

    for (t = 0; t < this->Nt_predict; t++)
      {
        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        this->Model.InitStep();
        this->OutputSaver.InitStep(this->Model);
        this->Model.Forward();

        this->OutputSaver.SetGroup("prediction");
        this->OutputSaver.Save(this->Model);

        this->OutputSaver.SetGroup("forecast");
        this->OutputSaver.Save(this->Model);
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_FOURDIMVARDRIVER_CXX
#endif
