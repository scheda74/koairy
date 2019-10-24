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


#ifndef POLYPHEMUS_FILE_DRIVER_SEQUENTIALDRIVER_HXX


#include "AssimilationDriver.cxx"
#include <map>


namespace Polyphemus
{


  //////////////////////
  // SEQUENTIALDRIVER //
  //////////////////////


  /*! \brief This class is the base class for most sequential
    data-assimilation drivers. */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  class SequentialDriver:
    public AssimilationDriver < T, ClassModel, ClassOutputSaver,
                                ClassObsManager >
  {

  protected:

    /*! \brief Flag that indicates whether the positivity of the analyzed data
      is required.
    */
    bool with_positivity_requirement;

  public:

    /*** Constructor and destructor ***/

    SequentialDriver(string config_file);
    virtual ~SequentialDriver();

    /*** Methods ***/

    virtual void ReadConfiguration();

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_SEQUENTIALDRIVER_HXX
#endif
