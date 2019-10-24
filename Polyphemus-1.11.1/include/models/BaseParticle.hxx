// Copyright (C) 2009, ENPC - INRIA - EDF R&D
// Author(s): Pierre Tran
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

// This file is part of a Lagrangian model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_BASEPARTICLE_HXX


namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


#include <vector>
#include <cmath>
#include "AtmoData.hxx"

  using namespace std;
  using namespace AtmoData;


  ///////////////////
  // BASE PARTICLE //
  ///////////////////


  /*! Base class for particles to be used by the class LagrangianTransport
    and its derived classes.
  */
  template<class T>
  class BaseParticle
  {

  public:

    /*** Destructor ***/

    virtual ~BaseParticle() {};

    /*** Methods ***/

    virtual bool IsOutside(T x_min, T x_max, T y_min, T y_max) const = 0;
    virtual T GetConcentrationContributionOnPoint(int s, T z, T lat, T lon)
      const = 0;

    template<class ClassModel>
    void Transport(ClassModel& Model)
    {
      throw string("\"BaseParticle<T>::Transport(ClassModel&)\" is not"
                   " defined.");
    }

    virtual void Update(T time_elapsed) = 0;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_BASEPARTICLE_HXX
#endif
