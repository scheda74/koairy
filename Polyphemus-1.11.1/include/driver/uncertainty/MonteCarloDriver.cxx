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


#ifndef POLYPHEMUS_FILE_DRIVER_MONTECARLODRIVER_CXX


#include "MonteCarloDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::MonteCarloDriver(string config_file):
    BaseDriver<T, ClassModel, ClassOutputSaver>(config_file),
    Nensemble(0)
  {
  }


  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::~MonteCarloDriver()
  {
    // Empties temporal concentration files.
    for (int i = 0; i < (int) conc_file_list.size(); i++)
      {
        string filename = conc_file_list[i];
        ofstream(filename.c_str()).close();
      }
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>::ReadConfiguration()
  {
    // Reads Monte Carlo section.
    this->config.SetSection("[MonteCarlo]");

    this->config.PeekValue("Number_ensemble", Nensemble);
    this->config.PeekValue("File_random_number", file_random_number);

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
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>::Init()
  {
    this->ReadConfiguration();

    Nstate = this->Model.GetNstate();
  }


  /////////////
  // METHODS //
  /////////////


  //! Performs Monte Carlo simulations.
  /*! Initializes the model, the output saver, and then performs the time
    loop.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>::Run()
  {
    int i, l;

    /*** Initializations ***/

    this->Model.Init();
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Assuming the MPI initialization was called in the model.
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
#else
    this->rank = 0;
#endif

    if (this->rank == 0)
      this->OutputSaver.Init(this->Model);
    Init();
    if (this->rank == 0)
      this->PerturbManager.Init(this->Model);

    /*** Allocates field arrays and set perturbation random numbers ****/

    if (this->rank == 0)
      {
        this->AllocateFieldArray();
        this->SetPerturbNumEns(Nensemble);
      }

    /*** Ensemble initialization ***/

    Array<T, 1> state_vector;
    if (this->rank == 0)
      {
        state_vector.resize(Nstate);
        this->Model.GetState(state_vector);
      }
    Array<T, 2> Ensemble(Nstate, Nensemble, ColumnMajorArray<2>());

    if (this->rank == 0)
      InitEnsemble(state_vector, Ensemble);

    /*** Time loop ***/

    string flag_control;
    for (i = 0; i < this->Model.GetNt(); i++)
      {
        if (this->rank == 0)
          {
            if (this->option_display["show_iterations"])
              cout << "Performing iteration #" << i << endl;
            if (this->option_display["show_date"])
              cout << "Current date: " <<
                this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;
          }

        this->Model.InitStep();
        if (this->rank == 0)
          this->OutputSaver.InitStep(this->Model);

        Forecast(Ensemble);

        if (this->rank == 0)
          {
            // Saves ensemble forecasts.
            this->OutputSaver.SetGroup("ensemble_forecast");
            for (l = 0; l < Nensemble; l++)
              {
                this->Model.SetState(Ensemble(Range::all(), l));
                this->OutputSaver.Save(this->Model);
              }
          }
      } // end time loop.

    /*** Frees field arrays ***/

    if (this->rank == 0)
      this->DeallocateFieldArray();
  }


  //! Ensemble initialization.
  /*! Initializes the ensemble. In the present implementation, all the members
    in the ensemble are the same as the the initial condition. The entire
    concentration data is saved in the files for each member in the ensemble.
    \param state_vector the initial condition.
    \param Ensemble (output) the initialized ensemble.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::InitEnsemble(Array<T, 1>& state_vector, Array<T, 2>& Ensemble)
  {
    // Sets ensemble members identical to initial state vector.
    for (int k = 0; k < Nstate; k++)
      for (int l = 0; l < Nensemble; l++)
        Ensemble(k, l) = state_vector(k);

    // Saves entire concentration data to files.
    Array<T, 4> conc(this->Model.GetNs(), this->Model.GetNz(),
                     this->Model.GetNy(), this->Model.GetNx());
    conc = this->Model.GetConcentration().GetArray();
    for (int l = 0; l < Nensemble; l++)
      FormatBinary<float>().Write(conc, conc_file_list[l]);
  }


  //! Monte Carlo forecast step.
  /*! Propagates the ensemble.
    \param Ensemble (input/output) on entry, the ensemble matrix at current
    date; on exit, the forecast ensemble matrix.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::Forecast(Array<T, 2>& Ensemble)
  {
    int l, k;
    Array<T, 1> working_vector(Nstate);
    Array<T, 4> conc(this->Model.GetNs(), this->Model.GetNz(),
                     this->Model.GetNy(), this->Model.GetNx());

    // Saves field data into field array before perturbations.
    this->SetFieldArray();

    // Reads entire initial concentration data for the forecast of the first
    // sample.
    FormatBinary<float>().Read(conc_file_list[0], conc);
    this->Model.GetConcentration().GetArray() = conc;

    /*** Forecasts of the samples with perturbations ***/

    for (l = 0; l < Nensemble; l++)
      {
        if (l != 0)
          {
            this->Model.SetDate(this->Model.GetCurrentDate());
            this->Model.InitStep();
          }

        if (this->rank == 0)
          {
            // Disturbs model simulation by perturbing field data.
            this->PerturbManager.GenerateField(this->Model,
                                               this->PerturbNum_ens[l]);

            // Model forecast with the perturbed field from the l-th sample.
            for (k = 0; k < Nstate; k++)
              working_vector(k) = Ensemble(k, l);
            this->Model.SetState(working_vector);
          }

        this->Model.Forward();

        conc = this->Model.GetConcentration().GetArray();
        FormatBinary<float>().Write(conc, conc_file_list[l]);

        this->Model.GetState(working_vector);
        for (k = 0; k < Nstate; k++)
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


  //////////////////////////////////////////
  // STORAGE MANAGEMENTS FOR PERTURBATION //
  //////////////////////////////////////////


  //! Allocates field arrays.
  /*! It allocates the arrays for the safeguard of the unperturbed field data.
   */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::AllocateFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      if (field_list[i] != "WindAngle")
        {
          if (this->Model.GetFieldDimension(field_list[i]) == 2)
            A2_perturb_map[field_list[i]] =
              new Array<T, 2>(this->Model.A2(field_list[i]).shape());
          if (this->Model.GetFieldDimension(field_list[i]) == 3)
            A3_perturb_map[field_list[i]] =
              new Array<T, 3>(this->Model.A3(field_list[i]).shape());
          if (this->Model.GetFieldDimension(field_list[i]) == 4)
            A4_perturb_map[field_list[i]] =
              new Array<T, 4>(this->Model.A4(field_list[i]).shape());
          if (this->Model.GetFieldDimension(field_list[i]) == 5)
            A5_perturb_map[field_list[i]] =
              new Array<T, 5>(this->Model.A5(field_list[i]).shape());
        }
  }


  //! Deallocates field arrays.
  /*! It frees the arrays allocated for the safeguard of the unperturbed field
    data.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::DeallocateFieldArray()
  {
    typename map<string, Array<T, 2>* >::iterator iter2;
    if (!A2_perturb_map.empty())
      for (iter2 = A2_perturb_map.begin();
           iter2 != A2_perturb_map.end(); iter2++)
        delete iter2->second;
    typename map<string, Array<T, 3>* >::iterator iter3;
    if (!A3_perturb_map.empty())
      for (iter3 = A3_perturb_map.begin();
           iter3 != A3_perturb_map.end(); iter3++)
        delete iter3->second;
    typename map<string, Array<T, 4>* >::iterator iter4;
    if (!A4_perturb_map.empty())
      for (iter4 = A4_perturb_map.begin();
           iter4 != A4_perturb_map.end(); iter4++)
        delete iter4->second;
    typename map<string, Array<T, 5>* >::iterator iter5;
    if (!A5_perturb_map.empty())
      for (iter5 = A5_perturb_map.begin();
           iter5 != A5_perturb_map.end(); iter5++)
        delete iter5->second;
  }


  //! Sets field arrays to model field data.
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>::SetFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      if (field_list[i] != "WindAngle")
        {
          if (this->Model.GetFieldDimension(field_list[i]) == 2)
            *A2_perturb_map[field_list[i]]
              = this->Model.A2(field_list[i]);
          if (this->Model.GetFieldDimension(field_list[i]) == 3)
            *A3_perturb_map[field_list[i]]
              = this->Model.A3(field_list[i]);
          if (this->Model.GetFieldDimension(field_list[i]) == 4)
            *A4_perturb_map[field_list[i]]
              = this->Model.A4(field_list[i]);
          if (this->Model.GetFieldDimension(field_list[i]) == 5)
            *A5_perturb_map[field_list[i]]
              = this->Model.A5(field_list[i]);
        }
  }


  //! Sets model field data from data stored in field arrays.
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>::SetModelFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      if (field_list[i] != "WindAngle")
        {
          if (this->Model.GetFieldDimension(field_list[i]) == 2)
            this->Model.A2(field_list[i])
              = *A2_perturb_map[field_list[i]];
          if (this->Model.GetFieldDimension(field_list[i]) == 3)
            this->Model.A3(field_list[i])
              = *A3_perturb_map[field_list[i]];
          if (this->Model.GetFieldDimension(field_list[i]) == 4)
            this->Model.A4(field_list[i])
              = *A4_perturb_map[field_list[i]];
          if (this->Model.GetFieldDimension(field_list[i]) == 5)
            this->Model.A5(field_list[i])
              = *A5_perturb_map[field_list[i]];
        }
  }


  //! Generates random numbers for the field data perturbation.
  /*! It generates random numbers for the field data perturbation. For each
    member in the ensemble, the generated randoms numbers are stored in a map
    that maps parameter names in their field to the corresponding
    random numbers.
    \param Nens ensemble number.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void MonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::SetPerturbNumEns(int Nens)
  {
    map<string, vector<T> > PerturbNum_map;
    vector<string> perturb_field_list = PerturbManager.GetPerturbFieldList();
    map<string, vector<PolairParam<float> > >& Param_map
      = PerturbManager.GetPerturbationMap();

    PerturbNum_ens.clear();

    for (int s = 0; s < Nens; s++)
      {
        PerturbNum_map = PerturbManager.GeneratePerturbNumMap();
        PerturbNum_ens.push_back(PerturbNum_map);
      }

    // Saves the random numbers.
    ofstream output_stream(file_random_number.c_str());
    output_stream.precision(16);
    for (int s = 0; s < Nens; s++)
      {
        output_stream << "Member " << s << endl;

        for (int l = 0; l < int(perturb_field_list.size()); l++)
          {
            bool additional_field
              = perturb_field_list[l] == "AdditionalField";

            for (int p = 0; p < int(Param_map[perturb_field_list[l]].size());
                 p++)
              {
                PolairParam<float>& param
                  = Param_map[perturb_field_list[l]][p];

                string fieldname;
                if (additional_field)
                  fieldname = param.GetName();
                else
                  fieldname = param.GetFieldName();

                T pert = PerturbNum_ens[s][perturb_field_list[l]][p];

                output_stream << '\t' << fieldname << '\t' << pert << endl;
              }
          }
      }

    output_stream.close();
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_MONTECARLODRIVER_CXX
#endif
