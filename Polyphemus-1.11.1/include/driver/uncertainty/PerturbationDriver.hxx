// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Damien Garaud
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


#ifndef POLYPHEMUS_FILE_DRIVER_PERTURBATIONDRIVER_HXX


#include <vector>
#include <string>
#include <map>
#include <utility>
#include "AtmoData.hxx"

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
#include <mpi.h>
#endif

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////
  // BASEDRIVER //
  ////////////////


  /*! \brief Driver that performs a forward time-integration with simple
    perturbations on the input data. */
  /*! The driver is responsible for the model initialization, for the time
    loop and for the calls to the output saver. Its reference floating-point
    precision is 'T'. The model is an instance of 'ClassModel' and the output
    saver is an instance of 'ClassOutputSaver'.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  class PerturbationDriver
  {

  protected:

    /*** Main components ***/

    //! Underlying model.
    ClassModel Model;
    //! Output saver.
    ClassOutputSaver OutputSaver;

    /*** Configuration ***/

    //! Configuration stream.
    ConfigStream config;
    //! Display options.
    map<string, bool> option_display;

    /*** Perturbations ***/

    //! List of the fields (independent of the species) to be perturbed.
    vector<string> additional_field;
    //! List of perturbation types for the fields independent of the species.
    vector<string> additional_field_type;
    //! List of perturbation values for the fields independent of the species.
    vector<T> additional_field_perturbation;

    //! List of the fields (dependent of the species) to be perturbed.
    vector<string> field_list;
    /*! List of the fields (dependent of the species) and associated species
      to be perturbed. */
    map<string, vector<string> > field;
    /*! List of pertubation types per field (dependent of the species) and per
      species. */
    map<string, vector<string> > field_type;
    /*! List of pertubation values per field (dependent of the species) and
      per species. */
    map<string, vector<T> > field_perturbation;

  public:

    /*** Constructors and destructor ***/

    PerturbationDriver(string config_file);
    virtual ~PerturbationDriver();
    virtual void Run();
    template<int N>
    void Apply(Array<T, N>& A, string type, T value);
    void WindAngle(string type, T value);
    void WindModule(string type, T value);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PERTURBATIONDRIVER_HXX
#endif
