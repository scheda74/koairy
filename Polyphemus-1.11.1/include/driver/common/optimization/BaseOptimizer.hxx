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


#ifndef POLYPHEMUS_FILE_OPTIMIZATION_BASEOPTIMIZER_HXX


namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


#include <iostream>
  using namespace std;


  ///////////////////
  // BASEOPTIMIZER //
  ///////////////////


  /*! \brief This class is the base class for optimization solvers. It defines
    the interface of all optimization solvers.
  */
  template<class T>
  class BaseOptimizer
  {

  protected:

    //! Configuration stream.
    ConfigStream config;

  public:

    //! The number of parameters to be optimized.
    int Nparam;
    /*! \brief The array that stores parameter values. On entry, the initial
      parameter array; on exit, it returns the optimized parameters. */
    Array<T, 1> param;
    /*! \brief The value of cost function for given parameter values. On first
      entry, set to zero; on final exit, cost function value for optimized
      parameters */
    double cost;
    /*! \brief The array that stores gradient values. On first entry,
      unspecified; on final exit, it returns the gradient array for optimized
      parameters. */
    Array<double, 1> grad;

  public:

    /*** Constructor ***/

    BaseOptimizer(string config_file, int Npar);
    virtual ~BaseOptimizer();

    /*** Initialization ***/

    virtual void AllocateArray(int Npar);
    virtual void Init();

    /*** Methods ***/

    virtual bool Optimize() = 0;
    virtual bool IsStop() = 0;

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OPTIMIZATION_BASEOPTIMIZER_HXX
#endif
