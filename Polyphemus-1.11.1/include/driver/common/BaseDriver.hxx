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


#ifndef POLYPHEMUS_FILE_DRIVER_BASEDRIVER_HXX


#include <vector>
#include <string>
#include <map>
#include "AtmoDataHeader.hxx"

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


  //! Simple driver that performs a forward time-integration.
  /*! The driver is responsible for the model initialization, for the time
    loop and for the calls to the output saver. Its reference floating-point
    precision is 'T'. The model is an instance of 'ClassModel' and the output
    saver is an instance of 'ClassOutputSaver'.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  class BaseDriver
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

    //! MPI rank.
    int rank;

  public:

    /*** Constructors and destructor ***/

    BaseDriver(string config_file);
    virtual ~BaseDriver();
    virtual void Run();

    /*** Access method ***/

    ClassModel& GetModel();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_BASEDRIVER_HXX
#endif
