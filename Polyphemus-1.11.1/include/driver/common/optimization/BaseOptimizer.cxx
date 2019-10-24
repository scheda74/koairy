// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_OPTIMIZATION_BASEOPTIMIZER_CXX


#include "BaseOptimizer.hxx"


namespace Polyphemus
{


  ///////////////////
  // BASEOPTIMIZER //
  ///////////////////


  //! Main constructor.
  template<class T>
  BaseOptimizer<T>::BaseOptimizer(string config_file, int Npar):
    config(config_file), Nparam(Npar), param(Nparam), cost(0.), grad(Nparam)
  {
  }


  //! Destructor.
  template<class T>
  BaseOptimizer<T>::~BaseOptimizer()
  {
  }


  //! Initializations.
  template<class T>
  void BaseOptimizer<T>::Init()
  {
  }

  //! Array allocations.
  template<class T>
  void BaseOptimizer<T>::AllocateArray(int Npar)
  {
    Nparam = Npar;
    param.resize(Nparam);
    grad.resize(Nparam);
    cost = 0.;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OPTIMIZATION_BASEOPTIMIZER_CXX
#endif
