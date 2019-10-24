// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Mohamed Aissaoui, Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_PERTURBATION_RANDGENERATOR_HXX


#include "newran.h"


namespace Polyphemus
{


  ///////////////////
  // RANDGENERATOR //
  ///////////////////


  //! Applies a random perturbation to a parameter.
  template<class T, class Y>
  class RandGenerator
  {

  protected:

    //! Directory where the seed is stored.
    string seed_directory;
    //! Seed number.
    Y seed_number;
    //! Seed number or directory?
    bool with_seed_number;

    //! Seed number.
    T* urng_;

    /*! Every random number cannot exceed the mean plus or minus
      'maximum_spread' times the standard deviation. */
    Y maximum_spread;

  public:

    /*** Constructor ***/

    RandGenerator(string seed_directory_ = "");
    virtual ~RandGenerator();


    /*** Access ***/

    Y GetMaximumSpread();
    void SetMaximumSpread(Y spread);

    /***  Methods ***/

    void Init(string seed_directory_ = "");
    Y GenRandomNumber(PolairParam<Y>& param);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_RANDGENERATOR_HXX
#endif
