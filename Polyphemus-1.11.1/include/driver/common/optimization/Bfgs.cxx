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


#ifndef POLYPHEMUS_FILE_OPTIMIZATION_BFGS_CXX


#include "Bfgs.hxx"


namespace Polyphemus
{


  //////////
  // BFGS //
  //////////


  //! Main constructor.
  template<class T>
  Bfgs<T>::Bfgs(string config_file, int Npar, bool disp_iter) :
    BaseOptimizer<T>(config_file, Npar),
    Nmetric(5), lbound(this->Nparam),
    ubound(this->Nparam), bound_type(this->Nparam), accuracy(1e+1), pgtol(0.),
    task(' ', 61), iprint(-1), display_iter(disp_iter), fint(1),
    work_arr_length((2 * Nmetric + 4)*this->Nparam + (11 * Nmetric + 8)*Nmetric),
    work_arr(work_arr_length),
    work_iarr(3 * this->Nparam), work_carr(60), work_larr(4), work_iarr0(44),
    work_arr0(29)
  {
    this->param = 0.;
  }


  //! Destructor.
  template<class T>
  Bfgs<T>::~Bfgs()
  {
  }


  //! Array allocations.
  template<class T>
  void Bfgs<T>::AllocateArray(int Npar)
  {
    BaseOptimizer<T>::AllocateArray(Npar);

    lbound.resize(this->Nparam);
    ubound.resize(this->Nparam);
    bound_type.resize(this->Nparam);

    work_arr_length = (2 * Nmetric + 4) * this->Nparam + (11 * Nmetric + 8) * Nmetric;
    work_arr.resize(work_arr_length);
    work_iarr.resize(3 * this->Nparam);

  }


  //! First initialization.
  template<class T>
  void Bfgs<T>::Init()
  {
    task = fill("START", 60, ' ', ios_base::left);

    _setulb(&this->Nparam, &Nmetric, this->param.data(), lbound.data(),
            ubound.data(), bound_type.data(), &this->cost, this->grad.data(),
            &accuracy, &pgtol, work_arr.data(), work_iarr.data(),
            (char*)task.c_str(), &iprint, work_carr.data(), work_larr.data(),
            work_iarr0.data(), work_arr0.data(), &fint, &fint);
  }


  //! Boundary condition settings.
  /*!
    \param lb the lower bound array for parameter values.
    \param ub the upper bound array for parameter values.
    \param bt the bound type array for parameters.
    bt(i) = 0 if param(i) is unbounded,
    bt(i) = 1 if param(i) has only a lower bound,
    bt(i) = 2 if param(i) has both lower and upper bounds, and
    bt(i) = 3 if param(i) has only an upper bound.
  */
  template<class T>
  void Bfgs<T>::SetBoundary(Array<double, 1>& lb, Array<double, 1>& ub,
                            Array<int, 1>& bt)
  {
    lbound = lb;
    ubound = ub;
    bound_type = bt;
  }


  //! Boundary condition settings.
  /*!
    \param lb the lower bound for all parameter values.
    \param ub the upper bound for all parameter values.
    \param bt the bound type for all parameters.
    bt = 0 if all parameters are unbounded,
    bt = 1 if all parameters have only lower bounds,
    bt = 2 if all parameters have both lower and upper bounds, and
    bt = 3 if all parameters have only upper bounds.
  */
  template<class T>
  void Bfgs<T>::SetBoundary(double lb, double ub, int bt)
  {
    lbound = lb;
    ubound = ub;
    bound_type = bt;
  }


  //! Is optimization stopped?
  /*!
    \return True for optimization stopped; otherwise false.
  */
  template<class T>
  bool Bfgs<T>::IsStop()
  {
    return task.substr(0, 4) == "STOP";
  }


  //! Performs optimization using L_BFGS_B algorithm.
  template<class T>
  bool Bfgs<T>::Optimize()
  {
    _setulb(&this->Nparam, &Nmetric, this->param.data(), lbound.data(),
            ubound.data(), bound_type.data(), &this->cost, this->grad.data(),
            &accuracy, &pgtol, work_arr.data(), work_iarr.data(),
            (char*)task.c_str(), &iprint, work_carr.data(), work_larr.data(),
            work_iarr0.data(), work_arr0.data(), &fint, &fint);

    if (display_iter)
      {
        static int iter = 0;
        iter++;

        int ndisp = this->Nparam > 10 ? 10 : this->Nparam;

        cout << "    Iteration : " << iter << endl;
        cout << "    Cost : " << this->cost << endl;
        cout << "    Nparam : " << this->Nparam << endl;

        cout << "    Param :\t";
        cout.flush();
        for (int i = 0; i < ndisp; i++)
          {
            cout << this->param(i) << "\t";
            cout.flush();
          }
        if (this->Nparam > ndisp)
          cout << "...";
        cout << endl;

        cout << "    Grad :\t";
        cout.flush();
        for (int i = 0; i < ndisp; i++)
          {
            cout << this->grad(i) << "\t";
            cout.flush();
          }
        if (this->Nparam > ndisp)
          cout << "...";
        cout << endl;
      }

    if (task.substr(0, 5) == "NEW_X")
      _setulb(&this->Nparam, &Nmetric, this->param.data(), lbound.data(),
              ubound.data(), bound_type.data(), &this->cost,
              this->grad.data(), &accuracy, &pgtol, work_arr.data(),
              work_iarr.data(), (char*)task.c_str(), &iprint,
              work_carr.data(), work_larr.data(), work_iarr0.data(),
              work_arr0.data(), &fint, &fint);

    return task.substr(0, 2) == "FG";
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OPTIMIZATION_BFGS_CXX
#endif
