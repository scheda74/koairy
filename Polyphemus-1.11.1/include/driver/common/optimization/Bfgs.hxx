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


#ifndef POLYPHEMUS_FILE_OPTIMIZATION_BFGS_HXX


#include "BaseOptimizer.cxx"


namespace Polyphemus
{


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////


#define _setulb setulb_

  extern "C"
  {
    void _setulb(const int*, const int*, double*, double*, double*,
                 int*, double*, double*, double*, double*,
                 double*, int*, char*, int*, char*, int*, int*, double*,
                 int*, int*);
  }


  //////////
  // BFGS //
  //////////


  //! Optimization solver using L_BFGS_B algorithm.
  template<class T>
  class Bfgs: public BaseOptimizer<T>
  {

  public:

    //! Maximum number of variable metric corrections used to define
    //! the limited memory matrix.
    const int Nmetric;
    //! Lower bound.
    Array<double, 1> lbound;
    //! Upper bound.
    Array<double, 1> ubound;
    //! Type of bounds imposed.
    Array<int, 1> bound_type;
    //! Solution accurary.
    double accuracy;
    //! Projected gradient tolerance.
    double pgtol;
    //! Tasks.
    string task;
    //! Output? (always set to -1)
    int iprint;
    //! Display iterations during optimizing?
    bool display_iter;

  private:

    /*** Working variables ***/

    int fint;
    int work_arr_length;
    Array<double, 1> work_arr;
    Array<int, 1> work_iarr;
    Array<char, 1> work_carr;
    Array<int, 1> work_larr;
    Array<int, 1> work_iarr0;
    Array<double, 1> work_arr0;

  public:

    /*** Constructor ***/

    Bfgs(string config_file, int Npar, bool disp_iter);
    virtual ~Bfgs();

    /*** Initialization ***/

    virtual void AllocateArray(int Npar);
    virtual void Init();

    void SetBoundary(Array<double, 1>& lb, Array<double, 1>& ub,
                     Array<int, 1>& bt);
    void SetBoundary(double lb, double ub, int bt);

    /*** Methods ***/

    virtual bool Optimize();
    virtual bool IsStop();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OPTIMIZATION_BFGS_HXX
#endif
