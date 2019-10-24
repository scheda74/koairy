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


#ifndef POLYPHEMUS_FILE_DRIVER_RRSQRTDRIVER_CXX


#include "RRSQRTDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::RRSQRTDriver(string config_file):
    SequentialDriver < T, ClassModel, ClassOutputSaver,
    ClassObsManager > (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~RRSQRTDriver()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    SequentialDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    this->config.SetSection("[RRSQRT]");
    this->config.PeekValue("Number_analysis_mode", Nmode_max);
    this->config.PeekValue("Number_model_mode", Nmode_model);
    this->config.PeekValue("Number_observation_mode", Nmode_observation);
    this->config.PeekValue("Propagation_option", propagation_option);
    if (propagation_option == "finite_difference")
      this->config.PeekValue("Finite_difference_perturbation",
                             finite_difference_perturbation);
    else
      throw string("Option \"") + propagation_option + "\" not supported.";
  }


  ////////////////////
  // INITIALIZATION //
  ////////////////////


  //! Driver initialization.
  /*! It reads the configuration.
   */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    SequentialDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  /////////////
  // METHODS //
  /////////////


  //! Performs a simulation with reduced-rank square root Kalman filtering.
  /*! Initializes the model, the output saver, and the observation manager,
    then performs the time loop. At each time step, whenever observations are
    available, it assimilates them using RRSQRT algorithm.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {

    /*** Initializations ***/

    this->Model.Init();
    this->OutputSaver.Init(this->Model);
    this->ObsManager.Init(this->Model);
    Init();
    this->PerturbManager.Init(this->Model);

    /*** Allocates field arrays and set perturbation random numbers ****/

    this->AllocateFieldArray();
    this->SetPerturbNumEns(Nmode_model);

    // Mode matrix initialized to zero.
    Array<T, 2> ModeMatrix(this->Nstate, Nmode_max, ColumnMajorArray<2>());
    ModeMatrix = 0.;

    Array<T, 1> state_vector;

    /*** Time loop ***/

    cout << "\nASSIMILATION\n" << endl;

    string flag_control;
    for (int i = 0; i < this->Model.GetNt(); i++)
      {
        if (i < this->Nt_assim)
          flag_control = "assimilation";
        else
          flag_control = "prediction";

        if (flag_control == "assimilation")
          {
            if (this->option_display["show_iterations"])
              cout << "Performing iteration #" << i << endl;
            if (this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            this->Model.InitStep();
            this->OutputSaver.InitStep(this->Model);

            // Model integration over one time step
            // and propagation of the mode matrix.
            Forecast(state_vector, ModeMatrix);

            this->OutputSaver.SetGroup("forecast");
            this->OutputSaver.Save(this->Model);

            // Retrieves observations.
            this->ObsManager.SetDate(this->Model.GetCurrentDate());

            if (this->ObsManager.IsAvailable())
              {
                cout << "Performing analysis at date ["
                     << this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                     << "]...";
                cout.flush();

                Analyze(state_vector, ModeMatrix);
                this->Model.SetState(state_vector);

                this->OutputSaver.SetGroup("analysis");
                this->OutputSaver.Save(this->Model);

                cout << " done." << endl;
              }
          }
        else
          {
            if (i == this->Nt_assim)
              // First prediction step.
              cout << "\nPREDICTION\n" << endl;

            if (this->option_display["show_iterations"])
              cout << "Performing iteration #" << i << endl;
            if (this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            this->Model.InitStep();
            this->OutputSaver.InitStep(this->Model);

            this->Model.Forward();

            this->OutputSaver.SetGroup("forecast");
            this->OutputSaver.Save(this->Model);
            this->OutputSaver.SetGroup("prediction");
            this->OutputSaver.Save(this->Model);
          }
      } // end time loop.

    /*** Frees field arrays ***/

    this->DeallocateFieldArray();
  }


  //! RRSQRT forecast step.
  /*! Propagates the state and the modes.
    \param state_vector (output) the forecast state.
    \param ModeMatrix (input/output) on entry, the mode matrix at current
    date; on exit, the forecast mode matrix.
    \warning works for a single observed species and for ground observations.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Forecast(Array<T, 1>& state_vector, Array<T, 2>& ModeMatrix)
  {
    int i, j;

    // Stores model current state.
    Array<T, 1> old_state;
    this->Model.GetState(old_state);

    // Stores all model concentrations so as to be able to get back in time.
    Array<T, 4> concentration(this->Model.GetNs(), this->Model.GetNz(),
                              this->Model.GetNy(), this->Model.GetNx());
    concentration = this->Model.GetConcentration().GetArray();

    // Saves field data into field array before perturbations.
    this->SetFieldArray();

    /*** Forecasts from previous analysis results ***/

    this->Model.Forward();
    // Forecast state computed without perturbation.
    this->Model.GetState(state_vector);
    // Goes back.
    this->Model.StepBack(concentration);

    // Number of modes in the current analyzed mode matrix.
    int Nmode_analyzed = ModeMatrix.columns();

    Array<T, 1> working_vector(this->Nstate);

    if (propagation_option == "finite_difference")
      {
        for (j = 0; j < Nmode_analyzed; j++)
          {
            for (i = 0; i < this->Nstate; i++)
              working_vector(i) = old_state(i) +
                finite_difference_perturbation * ModeMatrix(i, j);

            this->Model.SetDate(this->Model.GetCurrentDate());
            this->Model.InitStep();

            // Perturbed state to be propagated.
            this->Model.SetState(working_vector);

            this->Model.Forward();
            this->Model.GetState(working_vector);
            this->Model.StepBack(concentration);

            // New mode computed on the basis of the perturbation.
            for (i = 0; i < this->Nstate; i++)
              ModeMatrix(i, j) = (working_vector(i) - state_vector(i))
                / finite_difference_perturbation;
          }
      }
    else // only finite difference is supported.
      throw Error("RRSQRTDriver::Forecast()",
                  "Tangent linear model is not implemented.");

    /*** Adds columns from square root of Q to the mode matrix ***/

    Array<T, 2> ModelErrorCovariance(this->Nstate, Nmode_model,
                                     ColumnMajorArray<2>());
    Array<T, 1> error_covariance_column(this->Nstate);

    for (j = 0; j < Nmode_model; j++)
      {
        this->Model.SetDate(this->Model.GetCurrentDate());
        this->Model.InitStep();

        // Disturbs model simulation by perturbing field data.
        this->PerturbManager.GenerateField(this->Model,
                                           this->PerturbNum_ens[j]);
        this->Model.Forward();
        this->Model.GetState(working_vector);

        // Form column of square root of Q
        for (i = 0; i < this->Nstate; i++)
          error_covariance_column(i) = (working_vector(i) - state_vector(i));
        for (i = 0; i < this->Nstate; i++)
          ModelErrorCovariance(i, j) = error_covariance_column(i);

        this->Model.StepBack(concentration);
        this->SetModelFieldArray();
      }

    error_covariance_column.free();
    concentration.free();
    old_state.free();
    working_vector.free();

    // Extends the mode matrix with columns from square root of Q.
    int Nmode_total = Nmode_analyzed + Nmode_model;
    ModeMatrix.resizeAndPreserve(this->Nstate, Nmode_total);
    for (j = Nmode_analyzed; j < Nmode_total; j++)
      for (i = 0; i < this->Nstate; i++)
        ModeMatrix(i, j) = ModelErrorCovariance(i, j - Nmode_analyzed);

    ModelErrorCovariance.free();

    /*** Truncate to Nmode_max columns ***/

    if (Nmode_total > Nmode_max)
      {
        Array<T, 2> TruncateModelCovMatrix(this->Nstate, Nmode_max,
                                           ColumnMajorArray<2>());
        int Nmode_forecast;
        Nmode_forecast = TruncateMatrixByColumn(ModeMatrix, Nmode_max,
                                                TruncateModelCovMatrix);
        if (Nmode_forecast < Nmode_max)
          cout << "Warning: there less modes in the forecast mode matrix"
               << " than the expected rank (" << Nmode_max << ")." << endl;
        ModeMatrix.resize(TruncateModelCovMatrix.shape());
        ModeMatrix = TruncateModelCovMatrix;

        TruncateModelCovMatrix.free();
      }

    /*** Forward for the next step ***/

    this->Model.Forward();
  }


  //! RRSQRT analyze step.
  /*! The state is updated by the combination of background state and
    innovations. The weight of the combination is computed according to RRSQRT
    algorithm.
    \param state_vector (input/output) on entry, the forecast state; on exit,
    the analyzed state.
    \param ModeMatrix (input/output) on entry, the forecast mode matrix; on
    exit, the analyzed mode matrix.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void RRSQRTDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Analyze(Array<T, 1>& state_vector, Array<T, 2>& ModeMatrix)
  {
    int r, l, k;
    double zero = 0.;
    double one = 1.;
    char Transpose = 'T';
    char NoTranspose = 'N';
    int info;

    int Nmode_forecast = ModeMatrix.columns();
    this->Nobs = this->ObsManager.GetNobs();

    // Computes H times L where H is the observation operator and L is the
    // mode matrix.
    Array<T, 2> HL(this->Nobs, Nmode_forecast, ColumnMajorArray<2>());
    Array<T, 1> working_vector(this->Nstate);
    for (r = 0; r < this->Nobs; r++)
      for (l = 0; l < Nmode_forecast; l++)
        {
          for (k = 0; k < this->Nstate; k++)
            working_vector(k) = ModeMatrix(k, l);
          HL(r, l) = this->ObsManager.TLMMltObsOperator(r, working_vector);
        }
    working_vector.free();

    // Read R.
    Array<T, 2> working_matrix(this->Nobs, this->Nobs, ColumnMajorArray<2>());
    working_matrix = this->ObsManager.GetCovariance();

    // Computes R^{1/2}.
    Array<T, 2> obs_sqrt_matrix(this->Nobs, Nmode_observation,
                                ColumnMajorArray<2>());
    ComputeSqrtMatrix(working_matrix, Nmode_observation, obs_sqrt_matrix);

    // working_matrix stores HLL'H' + R.
    dgemm_(&NoTranspose, &Transpose, &this->Nobs, &this->Nobs,
           &Nmode_forecast, &one, HL.data(), &this->Nobs, HL.data(),
           &this->Nobs, &one, working_matrix.data(), &this->Nobs);

    // Computes innovation.
    Array<T, 1> innovation(this->Nobs);
    for (r = 0; r < this->Nobs; r++)
      innovation(r) = this->ObsManager.GetObservation()(r)
        - this->ObsManager.MltObsOperator(r, state_vector);

    // Auxiliary matrix that stores [innovation, HL, R^{1/2}].
    int num_aux = 1 + Nmode_forecast + Nmode_observation;
    Array<T, 2> aux_matrix(this->Nobs, num_aux, ColumnMajorArray<2>());
    for (r = 0; r < this->Nobs; r++)
      {
        aux_matrix(r, 0) = innovation(r);
        for (l = 1; l < 1 + Nmode_forecast; l++)
          aux_matrix(r, l) = HL(r, l - 1);
        for (l = 1 + Nmode_forecast; l < num_aux; l++)
          aux_matrix(r, l) = obs_sqrt_matrix(r, l - 1 - Nmode_forecast);
      }

    obs_sqrt_matrix.free();
    innovation.free();

    // Solves (HLL'H' + R) * aux_matrix = [innovation, HL, R^{1/2}], that
    // is, aux_matrix = (HLL'H' + R)^{-1} * [innovation, HL, R^{1/2}].
    Array<int, 1> permutation(this->Nobs);
    dgesv_(&this->Nobs, &num_aux, working_matrix.data(), &this->Nobs,
           permutation.data(), aux_matrix.data(), &this->Nobs, &info);
    if (info != 0)
      throw Error("RRSQRTDriver::Analyze",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESV.");
    working_matrix.free();
    permutation.free();

    // working_matrix0 stores LL'H'.
    Array<T, 2> working_matrix0(this->Nstate, this->Nobs,
                                ColumnMajorArray<2>());
    dgemm_(&NoTranspose, &Transpose, &this->Nstate, &this->Nobs,
           &Nmode_forecast, &one, ModeMatrix.data(), &this->Nstate,
           HL.data(), &this->Nobs, &zero, working_matrix0.data(),
           &this->Nstate);
    HL.free();

    // aux_matrix0 stores [K * innovation, KHL, K * R^{1/2}].
    // Maximal storage is this->Nstate * (this->Nobs + num_aux).
    Array<T, 2> aux_matrix0(this->Nstate, num_aux, ColumnMajorArray<2>());
    dgemm_(&NoTranspose, &NoTranspose, &this->Nstate, &num_aux, &this->Nobs,
           &one, working_matrix0.data(), &this->Nstate, aux_matrix.data(),
           &this->Nobs, &zero, aux_matrix0.data(), &this->Nstate);
    working_matrix0.free();

    // Updates the state, state_vector += K * innovation.
    for (k = 0; k < this->Nstate; k++)
      state_vector(k) += aux_matrix0(k, 0);

    // Computes L = (I - KH) * L.
    for (k = 0; k < this->Nstate; k++)
      for (l = 0; l < Nmode_forecast; l++)
        ModeMatrix(k, l) -= aux_matrix0(k, l + 1);

    // Augments L columns with the observation modes, and truncates L if
    // necessary.
    int Nmode_total = Nmode_forecast + Nmode_observation;
    ModeMatrix.resizeAndPreserve(this->Nstate, Nmode_total);

    for (k = 0; k < this->Nstate; k++)
      for (l = 0; l < Nmode_observation; l++)
        ModeMatrix(k, l + Nmode_forecast) =
          aux_matrix0(k, l + 1 + Nmode_forecast);

    aux_matrix0.free();

    if (Nmode_total > Nmode_max)
      {
        Array<T, 2> TruncateModelCovMatrix(this->Nstate, Nmode_max,
                                           ColumnMajorArray<2>());
        int Nmode_analyzed;

        Nmode_analyzed = TruncateMatrixByColumn(ModeMatrix, Nmode_max,
                                                TruncateModelCovMatrix);
        if (Nmode_analyzed < Nmode_max)
          cout << "Warning: there less modes in the forecast mode matrix"
               << " than the expected rank (" << Nmode_max << ")." << endl;

        ModeMatrix.resize(TruncateModelCovMatrix.shape());
        ModeMatrix = TruncateModelCovMatrix;

        TruncateModelCovMatrix.free();
      }

    // Positivity requirement.
    if (this->with_positivity_requirement)
      for (k = 0; k < this->Nstate; k++)
        if (state_vector(k) < 0.)
          state_vector(k) = 0.;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_RRSQRTDRIVER_CXX
#endif
