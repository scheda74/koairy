// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_DRIVER_MONTECARLODRIVER_HXX


#include "BaseDriver.cxx"
#include "PerturbationManager.cxx"
#include "GroundObservationManager.cxx"
#include <map>


namespace Polyphemus
{


  //////////////////////
  // MONTECARLODRIVER //
  //////////////////////


  /*! \brief This driver performs Monte Carlo simulations.
   */
  template<class T, class ClassModel, class ClassOutputSaver>
  class MonteCarloDriver:
    public BaseDriver<T, ClassModel, ClassOutputSaver>
  {

  protected:

    //! Dimension of the state.
    int Nstate;

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

    //! The number of samples in the ensemble.
    int Nensemble;

    //! File into which the actual random numbers should be stored.
    string file_random_number;

    /*! Name list of the files that saves coresponding entire concentration
      data for each member in the ensemble.
    */
    vector<string> conc_file_list;

  public:

    /*** Constructor and destructor ***/

    MonteCarloDriver(string config_file);
    virtual ~MonteCarloDriver();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    virtual void Init();

    /*** Methods ***/

    virtual void Run();

    void Forecast(Array<T, 2>& Ensemble);
    void InitEnsemble(Array<T, 1>& state_vector, Array<T, 2>& Ensemble);

    /*** Storage managements for perturbations ***/

    void AllocateFieldArray();
    void DeallocateFieldArray();
    void SetFieldArray();
    void SetModelFieldArray();
    void SetPerturbNumEns(int Nens);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_MONTECARLODRIVER_HXX
#endif
