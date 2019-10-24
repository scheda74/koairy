// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Irène Korsakissok
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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEDRIVER_CXX


#include "PlumeDriver.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*! Builds the driver and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  PlumeDriver<T, ClassModel, ClassOutputSaver>
  ::PlumeDriver(string config_file):
    Model(config_file), config(config_file)
  {

    /*** Display options ***/

    config.SetSection("[display]");
    // Should meteorological data be displayed on screen?
    config.PeekValue("Show_meteorological_data",
                     option_display["show_meteo"]);
    // Should iterations be displayed on screen?
    config.PeekValue("Show_iterations", option_display["show_iterations"]);

    /*** Reads meteorological conditions ***/

    config.SetSection("[gaussian]");
    config.GetValue("File_meteo", file_meteo);
  }


  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  PlumeDriver<T, ClassModel, ClassOutputSaver>::~PlumeDriver()
  {
  }


  //! Performs the simulation.
  /*! Initializes the model and the output saver, and then performs the loop
    over all meteorological conditions with calls to the model and to the
    output saver.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void PlumeDriver<T, ClassModel, ClassOutputSaver>::Run()
  {
    string line;

    /*** Initializations ***/

    Model.Init();
    OutputSaver.Init(Model);
    Nmeteo = 0;

    /*** Meteorological conditions ***/

    if (option_display["show_meteo"])
      cout << "\tTemperature\tWind angle\tWind velocity\tStability" << endl;

    ConfigStream meteo(file_meteo);
    while (!meteo.IsEmpty())
      {
        line = meteo.GetLine();
        if (split(line)[0] == "[situation]")
          {
            if (option_display["show_iterations"])
              cout << "Case #" << Nmeteo << endl;

            Model.InitMeteo(meteo, option_display["show_meteo"]);
            OutputSaver.InitStep(Model);

            Model.Compute();

            OutputSaver.Save(Model);

            Nmeteo++;
          }
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEDRIVER_CXX
#endif
