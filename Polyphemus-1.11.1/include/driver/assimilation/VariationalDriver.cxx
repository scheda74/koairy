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


#ifndef POLYPHEMUS_FILE_DRIVER_VARIATIONALDRIVER_CXX


#include "VariationalDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::VariationalDriver(string config_file):
    AssimilationDriver < T, ClassModel, ClassOutputSaver,
    ClassObsManager > (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~VariationalDriver()
  {
  }


  //! Read configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>::
      ReadConfiguration();
  }


  //! Driver initialization.
  /*! It reads configurations */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>::
      Init();
  }


  //! Empty method.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void VariationalDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    throw string("No assimilation method is implemented in class")
      + string(" VariationalDriver.");
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_VARIATIONALDRIVER_CXX
#endif
