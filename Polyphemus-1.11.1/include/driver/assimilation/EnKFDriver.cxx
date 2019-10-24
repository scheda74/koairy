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


#ifndef POLYPHEMUS_FILE_DRIVER_ENKFDRIVER_CXX


#include "EnKFDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::EnKFDriver(string config_file):
    SequentialDriver < T, ClassModel, ClassOutputSaver,
    ClassObsManager > (config_file),
    Nensemble(0)
  {
  }


  //! Destructor.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~EnKFDriver()
  {
    // Empties temporal concentration files.
    for (int i = 0; i < int(conc_file_list.size()); i++)
      {
        string filename = conc_file_list[i];
        ofstream(filename.c_str()).close();
      }
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    SequentialDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    // Reads EnKF section.
    this->config.SetSection("[EnKF]");

    this->config.PeekValue("Number_ensemble", Nensemble);
    this->config.PeekValue("With_observation_perturbation",
                           with_observation_perturbation);
    this->config.PeekValue("With_square_root", with_square_root);
    this->config.PeekValue("With_advance_sampling", with_advance_sampling);
    this->config.PeekValue("With_ensemble_prediction",
                           with_ensemble_prediction);

    // Names for temporal files saving the ensemble of concentration data.
#ifdef WIN32
    string filename = "concentration_file_";
#else
    string filename = "/tmp/concentration_file_";
#endif
    string tmp;
    conc_file_list.clear();
    int index = -1;
    for (int i = 0; i < Nensemble; i++)
      {
        while (exists(filename + to_str(++index)));
        tmp = filename + to_str(index);
        ofstream(tmp.c_str()).close();
        conc_file_list.push_back(tmp);
      }
  }


  ////////////////////
  // INITIALIZATION //
  ////////////////////


  //! Driver initialization.
  /*! It reads configurations and gets the dimension of model state vector.
   */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    SequentialDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  /////////////
  // METHODS //
  /////////////


  //! Performs a simulation with ensemble Kalman filtering.
  /*! Initializes the model, the output saver, and the observation manager,
    then performs the time loop. At each time step, whenever observations are
    available, it assimilates them using EnKF algorithm.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int i, k, l;

    /*** Initializations ***/

    this->Model.Init();
    this->OutputSaver.Init(this->Model);
    this->ObsManager.Init(this->Model);
    Init();
    this->PerturbManager.Init(this->Model);

    /*** Allocates field arrays and set perturbation random numbers ****/

    this->AllocateFieldArray();
    this->SetPerturbNumEns(Nensemble);

    /*** Ensemble initialization ***/

    Array<T, 1> state_vector(this->Nstate);
    this->Model.GetState(state_vector);
    Array<T, 2> Ensemble(this->Nstate, Nensemble, ColumnMajorArray<2>());

    InitEnsemble(state_vector, Ensemble);

    /*** Time loop ***/

    cout << "\nASSIMILATION\n" << endl;

    string flag_control;
    for (i = 0; i < this->Model.GetNt(); i++)
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

            // Model integration over one time step and propagation of the
            // ensemble.
            Forecast(Ensemble);

            // Saves ensemble forecasts.
            this->OutputSaver.SetGroup("ensemble_forecast");
            for (l = 0; l < Nensemble; l++)
              {
                this->Model.SetState(Ensemble(Range::all(), l));
                this->OutputSaver.Save(this->Model);
              }

            // Sets state to ensemble mean.
            state_vector = 0.;
            for (k = 0; k < this->Nstate; k++)
              {
                for (l = 0; l < Nensemble; l++)
                  state_vector(k) += Ensemble(k, l);
                state_vector(k) = state_vector(k) / T(Nensemble);
              }
            this->Model.SetState(state_vector);

            // Saves ensemble forecast mean of concentration data.
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

                Array<T, 2> ObsPerturbation(this->ObsManager.GetNobs(),
                                            Nensemble, ColumnMajorArray<2>());
                if (with_observation_perturbation)
                  this->PerturbManager.
                    GenerateRandomColumnVector(this->ObsManager.
                                               GetCovariance(),
                                               ObsPerturbation,
                                               this->PerturbManager.
                                               GetObsMaximumSpread());
                else
                  ObsPerturbation = 0.;

                Analyze(ObsPerturbation, Ensemble);

                // Saves analyzed ensemble.
                this->OutputSaver.SetGroup("ensemble_analysis");
                for (l = 0; l < Nensemble; l++)
                  {
                    this->Model.SetState(Ensemble(Range::all(), l));
                    this->OutputSaver.Save(this->Model);
                  }

                // Sets state to ensemble mean.
                state_vector = 0.;
                for (k = 0; k < this->Nstate; k++)
                  {
                    for (l = 0; l < Nensemble; l++)
                      state_vector(k) += Ensemble(k, l);
                    state_vector(k) = state_vector(k) / T(Nensemble);
                  }
                this->Model.SetState(state_vector);

                // Saves ensemble analysis mean of concentration data.
                this->OutputSaver.SetGroup("analysis");
                this->OutputSaver.Save(this->Model);

                cout << " done." << endl;
              }
          }
        else
          {
            if (i == this->Nt_assim)
              {
                // First prediction step.
                cout << "\nPREDICTION\n" << endl;

                // Calculates average concentrations for predictions.
                if (!with_ensemble_prediction)
                  {
                    /*** Average concentration array for predictions ***/
                    Array<T, 4> average_conc(this->Model.GetNs(),
                                             this->Model.GetNz(),
                                             this->Model.GetNy(),
                                             this->Model.GetNx());
                    Array<T, 4> conc(this->Model.GetNs(),
                                     this->Model.GetNz(),
                                     this->Model.GetNy(),
                                     this->Model.GetNx());
                    average_conc = T(0.);
                    for (l = 0; l < Nensemble; l++)
                      {
                        FormatBinary<float>()
                          .Read(conc_file_list[l], conc);
                        average_conc += conc;
                      }
                    average_conc /= Nensemble;

                    this->Model.GetState(state_vector);
                    this->Model.GetConcentration().GetArray() = average_conc;
                    this->Model.SetState(state_vector);
                  }
              }

            if (this->option_display["show_iterations"])
              cout << "Performing iteration #" << i << endl;
            if (this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i")
                   << endl;

            this->Model.InitStep();
            this->OutputSaver.InitStep(this->Model);

            if (with_ensemble_prediction)
              {
                Forecast(Ensemble);

                // Joints the ensemble predictions to ensemble forecast in the
                // analysis step.
                this->OutputSaver.SetGroup("ensemble_forecast");
                for (l = 0; l < Nensemble; l++)
                  {
                    this->Model.SetState(Ensemble(Range::all(), l));
                    this->OutputSaver.Save(this->Model);
                  }

                // Sets prediction to ensemble mean.
                state_vector = 0.;
                for (k = 0; k < this->Nstate; k++)
                  {
                    for (l = 0; l < Nensemble; l++)
                      state_vector(k) += Ensemble(k, l);
                    state_vector(k) = state_vector(k) / T(Nensemble);
                  }
                this->Model.SetState(state_vector);
              }
            else
              this->Model.Forward();

            // Savings.
            this->OutputSaver.SetGroup("forecast");
            this->OutputSaver.Save(this->Model);
            this->OutputSaver.SetGroup("prediction");
            this->OutputSaver.Save(this->Model);
          }
      } // end time loop.

    /*** Frees field arrays ***/

    this->DeallocateFieldArray();
  }


  //! Ensemble initialization.
  /*! Initializes the ensemble. In the present implementation, all the members
    in the ensemble are the same as the the initial condition. The entire
    concentration data is saved in the files for each member in the ensemble.
    \param state_vector the initial condition.
    \param Ensemble (output) the initialized ensemble.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::InitEnsemble(Array<T, 1>& state_vector, Array<T, 2>& Ensemble)
  {
    // Sets ensemble members identical to initial state vector.
    for (int k = 0; k < this->Nstate; k++)
      for (int l = 0; l < Nensemble; l++)
        Ensemble(k, l) = state_vector(k);

    // Saves entire concentration data to files.
    Array<T, 4> conc(this->Model.GetNs(), this->Model.GetNz(),
                     this->Model.GetNy(), this->Model.GetNx());
    conc = this->Model.GetConcentration().GetArray();
    for (int l = 0; l < Nensemble; l++)
      FormatBinary<float>().Write(conc, conc_file_list[l]);
  }


  //! EnKF forecast step.
  /*! Propagates the ensemble.
    \param Ensemble (input/output) on entry, the ensemble matrix at current
    date; on exit, the forecast ensemble matrix.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Forecast(Array<T, 2>& Ensemble)
  {
    int l, k;
    Array<T, 1> working_vector(this->Nstate);
    Array<T, 4> conc(this->Model.GetNs(), this->Model.GetNz(),
                     this->Model.GetNy(), this->Model.GetNx());

    // Saves field data into field array before perturbations.
    this->SetFieldArray();

    // Reads entire initial concentration data for the forecast of the first
    // sample.
    FormatBinary<float>().Read(conc_file_list[0], conc);
    this->Model.GetConcentration().GetArray() = conc;

    /*** Forecasts of the first (Nensemble-1) samples with perturbations ***/

    for (l = 0; l < Nensemble; l++)
      {
        this->Model.SetDate(this->Model.GetCurrentDate());
        this->Model.InitStep();

        // Disturbs model simulation by perturbing field data.
        this->PerturbManager.GenerateField(this->Model,
                                           this->PerturbNum_ens[l]);

        // Model forecast with the perturbed field from the l-th sample.
        for (k = 0; k < this->Nstate; k++)
          working_vector(k) = Ensemble(k, l);
        this->Model.SetState(working_vector);

        this->Model.Forward();

        conc = this->Model.GetConcentration().GetArray();
        FormatBinary<float>().Write(conc, conc_file_list[l]);

        this->Model.GetState(working_vector);
        for (k = 0; k < this->Nstate; k++)
          Ensemble(k, l) = working_vector(k);

        // Moves one step back with the corresponding data for the first
        // (Nensemble - 1) members in the ensemble.
        if (l < Nensemble - 1)
          {
            // Sets back concentration data, and steps back one time step.
            FormatBinary<float>().Read(conc_file_list[l + 1], conc);
            this->Model.StepBack(conc);
          }

        // Sets back the field data before perturbation.
        this->SetModelFieldArray();
      }

    /*** Free of temporal arrays ***/

    conc.free();
    working_vector.free();
  }


  //! EnKF analyze step.
  /*! The ensemble is updated by the combination of background state and
    innovations. The weight of the combination is computed according to the
    classical EnKF algorithm with observation perturbations.
    \param ObsPerturbation the perturbation matrix for observations.
    \param Ensemble (input/output) on entry, the forecast ensemble matrix; on
    exit, the analyzed ensemble matrix.
  */
  template < class T, class ClassModel,
             class ClassOutputSaver, class ClassObsManager >
  void EnKFDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Analyze(Array<T, 2>& ObsPerturbation, Array<T, 2>& Ensemble)
  {
    int r, l, k;
    double zero = 0.;
    double one = 1.;
    char Transpose = 'T';
    char NoTranspose = 'N';
    int info;

    this->Nobs = this->ObsManager.GetNobs();

    // Constructs observation ensemble D.
    Array<T, 1> obs_vector(this->Nobs);
    for (r = 0; r < this->Nobs; r++)
      obs_vector(r) = this->ObsManager.GetObservation()(r);
    Array<T, 2> obs_ensemble(this->Nobs, Nensemble, ColumnMajorArray<2>());
    for (l = 0; l < Nensemble; l++)
      for (r = 0; r < this->Nobs; r++)
        obs_ensemble(r, l) = obs_vector(r) + ObsPerturbation(r, l);
    obs_vector.free();

    // Computes innovation matrix D' = D - HA, where H is the observation
    // operator.
    Array<T, 2> innovation_matrix(this->Nobs, Nensemble,
                                  ColumnMajorArray<2>());
    Array<T, 1> working_vector(this->Nstate);
    for (l = 0; l < Nensemble; l++)
      {
        for (k = 0; k < this->Nstate; k++)
          working_vector(k) = Ensemble(k, l);
        for (r = 0; r < this->Nobs; r++)
          innovation_matrix(r, l) = obs_ensemble(r, l)
            - this->ObsManager.MltObsOperator(r, working_vector);
      }
    working_vector.free();

    // Constructs state ensemble perturbation matrix L.
    Array<T, 1> mean_state(this->Nstate);
    Array<T, 2> ensemble_perturbation(this->Nstate, Nensemble,
                                      ColumnMajorArray<2>());
    mean_state = 0.;
    for (k = 0; k < this->Nstate; k++)
      {
        for (l = 0; l < Nensemble; l++)
          mean_state(k) += Ensemble(k, l);
        mean_state(k) = mean_state(k) / T(Nensemble);
        for (l = 0; l < Nensemble; l++)
          ensemble_perturbation(k, l) = Ensemble(k, l) - mean_state(k);
      }
    mean_state.free();

    // Computes H times L where H is the tangent linear observation operator.
    Array<T, 2> HL(this->Nobs, Nensemble, ColumnMajorArray<2>());
    working_vector.resize(this->Nstate);
    for (r = 0; r < this->Nobs; r++)
      for (l = 0; l < Nensemble; l++)
        {
          for (k = 0; k < this->Nstate; k++)
            working_vector(k) = ensemble_perturbation(k, l);
          HL(r, l) = this->ObsManager.TLMMltObsOperator(r, working_vector);
        }
    working_vector.free();

    // Read R.
    Array<T, 2> working_matrix(this->Nobs, this->Nobs, ColumnMajorArray<2>());
    working_matrix = this->ObsManager.GetCovariance();

    // working_matrix stores HLL'H' + R.
    dgemm_(&NoTranspose, &Transpose, &this->Nobs, &this->Nobs,
           &Nensemble, &one, HL.data(), &this->Nobs, HL.data(),
           &this->Nobs, &one, working_matrix.data(), &this->Nobs);

    // Computes (HLL'H' + R)^{-1} * D'.
    Array<int, 1> permutation(this->Nobs);
    dgesv_(&this->Nobs, &Nensemble, working_matrix.data(), &this->Nobs,
           permutation.data(), innovation_matrix.data(), &this->Nobs, &info);
    if (info != 0)
      throw Error("EnKFDriver::Analyze",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESV.");
    working_matrix.free();
    permutation.free();

    // working_matrix0 stores LL'H'.
    Array<T, 2> working_matrix0(this->Nstate, this->Nobs,
                                ColumnMajorArray<2>());
    dgemm_(&NoTranspose, &Transpose, &this->Nstate, &this->Nobs,
           &Nensemble, &one, ensemble_perturbation.data(), &this->Nstate,
           HL.data(), &this->Nobs, &zero, working_matrix0.data(),
           &this->Nstate);
    HL.free();

    // Computes LL'H' * (HLL'H' + R)^{-1} * D'
    // Maximal storage is this->Nstate * (this->Nobs + num_aux).
    Array<T, 2> aux_matrix0(this->Nstate, Nensemble, ColumnMajorArray<2>());
    dgemm_(&NoTranspose, &NoTranspose, &this->Nstate, &Nensemble, &this->Nobs,
           &one, working_matrix0.data(), &this->Nstate,
           innovation_matrix.data(), &this->Nobs, &zero, aux_matrix0.data(),
           &this->Nstate);
    working_matrix0.free();

    // Updates the ensemble A += K * D'.
    for (k = 0; k < this->Nstate; k++)
      for (l = 0; l < Nensemble; l++)
        Ensemble(k, l) += aux_matrix0(k, l);

    aux_matrix0.free();

    // Positivity requirement.
    if (this->with_positivity_requirement)
      for (k = 0; k < this->Nstate; k++)
        for (l = 0; l < Nensemble; l++)
          if (Ensemble(k, l) < 0.)
            Ensemble(k, l) = 0.;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ENKFDRIVER_CXX
#endif
