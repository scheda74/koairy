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


#ifndef POLYPHEMUS_FILE_DRIVER_BASEDRIVER_CXX


#include "BaseDriver.hxx"

#include "AtmoData.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*! Builds the driver and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  BaseDriver<T, ClassModel, ClassOutputSaver>::BaseDriver(string config_file):
    Model(config_file), config(config_file)
  {

    /*** Display options ***/

    config.SetSection("[display]");
    // Should iterations be displayed on screen?
    config.PeekValue("Show_iterations", option_display["show_iterations"]);
    // Should current date be displayed on screen?
    config.PeekValue("Show_date", option_display["show_date"]);
  }


  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  BaseDriver<T, ClassModel, ClassOutputSaver>::~BaseDriver()
  {
  }


  //! Performs the simulation.
  /*! Initializes the model and the output saver, and then performs the time
    loop with calls to the model and to the output saver.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void BaseDriver<T, ClassModel, ClassOutputSaver>::Run()
  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Init();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    /*** Initializations ***/

    Model.Init();
    if (rank == 0)
      OutputSaver.Init(Model);

    /*** Time loop ***/

    for (int i = 0; i < Model.GetNt(); i++)
      {
        if (rank == 0 && option_display["show_iterations"])
          cout << "Performing iteration #" << i << endl;

        if (rank == 0 && option_display["show_date"])
          cout << "Current date: "
               << Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        Model.InitStep();
        if (rank == 0)
          OutputSaver.InitStep(Model);

        Model.Forward();
        if (rank == 0)
          OutputSaver.Save(Model);
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Finalize();
#endif
  }


  //! Returns the model.
  /*!
    \return The model.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  ClassModel& BaseDriver<T, ClassModel, ClassOutputSaver>::GetModel()
  {
    return Model;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_BASEDRIVER_CXX
#endif
