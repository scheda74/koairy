// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu
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


#ifndef POLYPHEMUS_FILE_DRIVER_ASSIMILATIONDRIVER_HXX


#include <map>
#include "BaseDriver.cxx"
#include "PerturbationManager.cxx"


namespace Polyphemus
{


  ////////////////////////
  // ASSIMILATIONDRIVER //
  ////////////////////////


  //! This class is the base class for most data-assimilation drivers.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  class AssimilationDriver: public BaseDriver<T, ClassModel, ClassOutputSaver>
  {

  protected:

    //! Observation manager.
    ClassObsManager ObsManager;

    //! Dimension of the state.
    int Nstate;
    //! Number of observations.
    int Nobs;

    //! Perturbation manager.
    PerturbationManager<T> PerturbManager;
    //! Perturbation random numbers of the fields for all ensemble members.
    vector<map<string, vector<T> > > PerturbNum_ens;

    //! Map from 2D perturbation field names to their arrays.
    map<string, Array<T, 2>* > A2_perturb_map;
    //! Map from 3D perturbation field names to their arrays.
    map<string, Array<T, 3>* > A3_perturb_map;
    //! Map from 4D perturbation field names to their arrays.
    map<string, Array<T, 4>* > A4_perturb_map;
    //! Map from 5D perturbation field names to their arrays.
    map<string, Array<T, 5>* > A5_perturb_map;

    /*** Experiment settings ***/

    //! Assimilation interval in number of forward integrations.
    int Nt_assim;
    //! Prediction interval in number of forward integrations.
    int Nt_predict;

  public:

    AssimilationDriver(string config_file);
    virtual ~AssimilationDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Storage managements for perturbations ***/

    void AllocateFieldArray();
    void DeallocateFieldArray();
    void SetFieldArray();
    void SetModelFieldArray();
    void SetPerturbNumEns(int Nens);

    /*** Methods ***/

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ASSIMILATIONDRIVER_HXX
#endif
