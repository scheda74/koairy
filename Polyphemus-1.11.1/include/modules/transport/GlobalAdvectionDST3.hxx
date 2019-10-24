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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_GLOBALADVECTIONDST3_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModule.cxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////


#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _global_advection global_advection__
#else
#define _global_advection global_advection_
#endif

  extern "C"
  {
    void _global_advection(int*, int*, int*, double*, double*,
                           double*, double*, double*, double*, int*,
                           double*, double*, double*, double*, double*);
  }


  ///////////////////
  // ADVECTIONDST3 //
  ///////////////////


  //! This class is a numerical solver for advection.
  /*! It uses a Direct-Space-Time (DST) algorithm of third order with a
    Koren-like flux limiter.
  */
  template<class T>
  class GlobalAdvectionDST3: public BaseModule
  {

  public:

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);
    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    void Forward(Data<T, 3>& ZonalWind, Data<T, 3>& MeridionalWind,
                 Data<T, 3>& VerticalWind, int with_bc,
                 Array<T, 1>& CellWidth_x, Array<T, 1>& CellWidth_y,
                 Array<T, 1>& CellWidth_z, int Nx, int Ny, int Nz,
                 Data<T, 3>& BoundaryCondition_x,
                 Data<T, 3>& BoundaryCondition_y,
                 Data<T, 2>& BoundaryCondition_z, T Delta_t,
                 Data<T, 3>& Concentration);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_GLOBALADVECTIONDST3_HXX
#endif
