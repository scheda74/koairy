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


#ifndef POLYPHEMUS_FILE_DRIVER_OPTIMALINTERPOLATIONDRIVER_CXX


#include "OptimalInterpolationDriver.hxx"


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
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  OptimalInterpolationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::OptimalInterpolationDriver(string config_file):
    SequentialDriver < T, ClassModel, ClassOutputSaver,
    ClassObsManager > (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  OptimalInterpolationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~OptimalInterpolationDriver()
  {
  }


  /////////////
  // METHODS //
  /////////////


  //! Reads the configuration.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void
  OptimalInterpolationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    SequentialDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    this->config.SetSection("[data_assimilation]");

    this->config.PeekValue("Analyze_first_step", analyze_first_step);
  }


  //! Performs a simulation with optimal interpolation.
  /*! Initializes the model, the output saver, and the observation manager,
    then performs the time loop. At each time step, whenever observations are
    avaible, it assimilates them using optimal interpolation.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void
  OptimalInterpolationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int rank;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Init();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    cout.precision(20);

    Array<T, 1> state_vector;

    /*** Initializations ***/

    this->Model.Init();
    this->ObsManager.Init(this->Model);
    this->Init();

    if (rank == 0)
      cout << "\nASSIMILATION\n" << endl;

    if (analyze_first_step)
      {
        // Retrieves observations.
        this->ObsManager.SetDate(this->Model.GetCurrentDate());

        if (this->ObsManager.IsAvailable())
          {
            if (rank == 0)
              {
                cout << "Performing optimal interpolation at date ["
                     << this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                     << "]...";
                cout.flush();
              }

            this->Model.GetState(state_vector);
            Analyze(state_vector);
            this->Model.SetState(state_vector);

            if (rank == 0)
              {
                this->OutputSaver.SetGroup("analysis");
                this->OutputSaver.Save(this->Model);
                cout << " done." << endl;
              }
          }
      }

    if (rank == 0)
      this->OutputSaver.Init(this->Model);

    /*** Time loop ***/

    string flag_control;
    for (int i = 0; i < this->Model.GetNt(); i++)
      {
        if (i < this->Nt_assim)
          flag_control = "assimilation";
        else
          flag_control = "prediction";

        if (flag_control == "assimilation")
          {

            if (rank == 0)
              {
                if (this->option_display["show_iterations"])
                  cout << "Performing iteration #" << i << endl;
                if (this->option_display["show_date"])
                  cout << "Current date: " <<
                    this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                       << endl;
              }

            this->Model.InitStep();
            if (rank == 0)
              this->OutputSaver.InitStep(this->Model);

            this->Model.Forward();

            if (rank == 0)
              {
                this->OutputSaver.SetGroup("forecast");
                this->OutputSaver.Save(this->Model);
              }

            // Retrieves observations.
            this->ObsManager.SetDate(this->Model.GetCurrentDate());

            if (this->ObsManager.IsAvailable())
              {
                if (rank == 0)
                  {
                    cout << "Performing optimal interpolation at date ["
                         << this->
                      Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                         << "]...";
                    cout.flush();
                  }

                this->Model.GetState(state_vector);
                Analyze(state_vector);
                this->Model.SetState(state_vector);

                if (rank == 0)
                  {
                    this->OutputSaver.SetGroup("analysis");
                    this->OutputSaver.Save(this->Model);
                    cout << " done." << endl;
                  }
              }
          }
        else
          {
            if (rank == 0 && i == this->Nt_assim)
              // First prediction step.
              cout << "\nPREDICTION\n" << endl;

            if (rank == 0 && this->option_display["show_iterations"])
              cout << "Performing iteration #" << i << endl;
            if (rank == 0 && this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            this->Model.InitStep();
            if (rank == 0)
              this->OutputSaver.InitStep(this->Model);

            this->Model.Forward();

            if (rank == 0)
              {
                this->OutputSaver.SetGroup("forecast");
                this->OutputSaver.Save(this->Model);
                this->OutputSaver.SetGroup("prediction");
                this->OutputSaver.Save(this->Model);
              }
          }
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Finalize();
#endif
  }


  //! Performs analyze.
  /*! The state is updated by the combination of background state and
    innovations. The weight of the combination is computed according to
    optimal interpolation algorithm.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void
  OptimalInterpolationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Analyze(Array<T, 1>& state_vector)
  {
    int r, c, j;

    // Number of observations at current date.
    this->Nobs = this->ObsManager.GetNobs();

    // One row of background matrix B.
    Array<T, 1> error_covariance_row(this->Nstate);

    /*** Analyze of optimal interpolation ***/

    // Temporary arrays.
    Array<T, 2> working_matrix;
    working_matrix.resize(this->Nobs, this->Nobs);
    working_matrix = 0.;

    Array<T, 1> row(this->Nobs);

    // Computes HBH'.
    T H_entry;
    for (j = 0; j < this->Nstate; j++)
      {
        error_covariance_row = this->Model.BackgroundErrorCovariance(j);
        // Computes the j-th row of BH'.
        for (r = 0; r < this->Nobs; r++)
          row(r)
            = this->ObsManager.TLMMltObsOperator(r, error_covariance_row);
        // Keeps on building HBH'.
        for (r = 0; r < this->Nobs; r++)
          {
            H_entry = this->ObsManager.Linearized(r, j);
            for (c = 0; c < this->Nobs; c++)
              working_matrix(r, c) += H_entry * row(c);
          }
      }

    // Computes (HBH' + R).
    for (r = 0; r < this->Nobs; r++)
      for (c = 0; c < this->Nobs; c++)
        working_matrix(r, c) += this->ObsManager.GetCovariance()(r, c);

    // Computes (HBH' + R)^{-1}.
    Array<int, 1> permutations(this->Nobs);
    int info;
    dgetrf_(&this->Nobs, &this->Nobs, working_matrix.data(), &this->Nobs,
            permutations.data(), &info);
    if (info != 0)
      throw Error("OptimalInterpolationDriver::Analyze",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGETRF.");

    int lwork = this->Nobs * 64;
    Array<T, 1> work(lwork);
    dgetri_(&this->Nobs, working_matrix.data(), &this->Nobs,
            permutations.data(), work.data(), &lwork, &info);
    if (info != 0)
      throw Error("OptimalInterpolationDriver::Analyze",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGETRI.");

    // Computes innovation.
    Array<T, 1> innovation;
    innovation.resize(this->Nobs);
    innovation = 0.;
    for (r = 0; r < this->Nobs; r++)
      {
        innovation(r) = this->ObsManager.GetObservation()(r);
        innovation(r) -= this->ObsManager.MltObsOperator(r, state_vector);
      }

    // Computes working_matrix * innovation.
    Array<T, 1> working_vector;
    working_vector.resize(this->Nobs);
    working_vector = 0.;
    for (r = 0; r < this->Nobs; r++)
      for (j = 0; j < this->Nobs; j++)
        working_vector(r) += working_matrix(r, j) * innovation(j);

    // Computes new state.
    Array<T, 1> working_vector0;
    working_vector0.resize(this->Nobs);
    working_vector0 = 0;
    for (r = 0; r < this->Nstate; r++)
      {
        for (c = 0; c < this->Nobs; c++)
          {
            error_covariance_row = this->Model.BackgroundErrorCovariance(r);
            working_vector0(c)
              = this->ObsManager.TLMMltObsOperator(c, error_covariance_row);
          }
        for (j = 0; j < this->Nobs; j++)
          state_vector(r) += working_vector0(j) * working_vector(j);
      }

    // Positivity requirement.
    if (this->with_positivity_requirement)
      for (r = 0; r < this->Nstate; r++)
        if (state_vector(r) < 0.)
          state_vector(r) = 0.;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_OPTIMALINTERPOLATIONDRIVER_CXX
#endif
