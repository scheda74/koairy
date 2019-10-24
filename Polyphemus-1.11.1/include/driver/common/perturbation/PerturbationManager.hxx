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


#ifndef POLYPHEMUS_FILE_PERTURBATION_PERTURBATIONMANAGER_HXX


#include "PolairParam.cxx"
#include "RandGenerator.cxx"


namespace Polyphemus
{


  /////////////////////////
  // PERTURBATIONMANAGER //
  /////////////////////////


  /*! \brief This class is for the perturbation managements. It reads the
    perturbation configurations, and performs perturbations. The concerning
    fields are then updated according to perturbation results for new model
    simulations in diverse applications such as data assimilation and ensemble
    predictions.
  */
  template<class T>
  class PerturbationManager
  {

  protected:

    //! Seed number.
    NEWRAN::LGM_mixed urng_;

    //! Random generator.
    RandGenerator<NEWRAN::LGM_mixed, float>* randomization;

    //! The total number of fields to be perturbed.
    int Nfield;
    //! Name list of the fields to be perturbed.
    vector<string> perturb_field_list;
    //! Name list of the model field data to be perturbed.
    vector<string> perturb_data_list;

    //! Map from field names to their parameter lists.
    map<string, vector<PolairParam<float> > > Param_map;

    //! Name of the directory that contains newran seed files.
    string rand_seed_directory;
    /*! Every field random number cannot exceed the mean plus or minus
      'field_maximum_spread' times the standard deviation. */
    T field_maximum_spread;
    /*! Every observation random number cannot exceed the mean plus or minus
      'obs_maximum_spread' times the standard deviation. */
    T obs_maximum_spread;

  public:

    /*** Constructor ***/

    PerturbationManager();
    virtual ~PerturbationManager();

    /*** Configuration ***/

    template<class ClassModel>
    void ReadConfiguration(ClassModel& Model);

    /*** Initialization ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    /*** Access methods ***/

    map<string, vector<PolairParam<float> > >& GetPerturbationMap();
    vector<string>& GetPerturbDataList();
    vector<string>& GetPerturbFieldList();
    T GetFieldMaximumSpread();
    T GetObsMaximumSpread();

    /*** Perturbation methods ***/

    map<string, vector<T> > GeneratePerturbNumMap();

    void GenerateRandomColumnVector(Array<T, 2>& Matrix,
                                    Array<T, 2>& Covariance,
                                    T maximum_spread = 2.);
    void GenerateRandomRowVector(Array<T, 2>& Matrix,
                                 Array<T, 2>& Covariance,
                                 T maximum_spread = 2.);

    template<int N>
    void PerturbArray(Array<T, N>& array, double random_number, string pdf);

    template<class ClassModel>
    void GenerateField(ClassModel& Model,
                       map<string, vector<T> > PerturbNum_map);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_PERTURBATIONMANAGER_HXX
#endif
