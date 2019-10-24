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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_DIFFUSIONROS2_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;



  ///////////////////////
  // FORTRAN FUNCTIONS //
  ///////////////////////


#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _jacdiffX jacddiffdc_x__
#define _diffX diff_x__
#define _rosdiffX rosdiff_x__
#define _diffXcl diff_xcl__
#define _jacdiffY jacddiffdc_y__
#define _diffY diff_y__
#define _rosdiffY rosdiff_y__
#define _diffYcl diff_ycl__
#define _diffZ diff_z__
#define _rosdiffZ rosdiff_z__
#define _diffZcl diff_zcl__
#else
#define _jacdiffX jacddiffdc_x_
#define _diffX diff_x_
#define _rosdiffX rosdiff_x_
#define _diffXcl diff_xcl_
#define _jacdiffY jacddiffdc_y_
#define _diffY diff_y_
#define _rosdiffY rosdiff_y_
#define _diffYcl diff_ycl_
#define _diffZ diff_z_
#define _rosdiffZ rosdiff_z_
#define _diffZcl diff_zcl_
#endif
  
  extern "C"
  {
    void _jacdiffX(int*, double*, double*, double*, double*, double*,
		   double*, double*);
    void _diffX(int*, int*, int*, double*, double*, double*,
		double*, double*, double*, double*);
    void _rosdiffX(int*, double*, double*, double*, double*, double*,
		   double*, double*, double*, double*, double*, int*);
    void _diffXcl(int*, int*, int*, double*, double*, double*,
		  double*, double*, double*, double*, double*);

    void _jacdiffY(int*, double*, double*, double*, double*, double*,
		   double*, double*);
    void _diffY(int*, int*, int*, double*, double*, double*,
		double*, double*, double*, double*);
    void _rosdiffY(int*, double*, double*, double*, double*, double*,
		   double*, double*, double*, double*, double*, int*);
    void _diffYcl(int*, int*, int*, double*, double*, double*,
		  double*, double*, double*, double*, double*);
    void _diffZ(int*, int*, int*, double*, double*,
		double*, double*, double*,
		double*, double*, double*,
		double*, double*, double*, double*);
    void _rosdiffZ(int*, double*, double*, double*, double*, double*, double*,
		   double*, double*, double*, double*, double*, double*);
    void _diffZcl(int*, int*, int*, double*, double*,
		  double*, double*, double*,
		  double*, double*, double*,
		  double*, double*, double*, double*, double*);
  }

  
  ///////////////////
  // DIFFUSIONROS2 //
  ///////////////////


  //! This class is a numerical solver for diffusion.
  /*! It uses a second-order Rosenbrock method.
   */
  template<class T>
  class DiffusionROS2: public BaseModuleParallel
  {

  protected:

    //! Is zonal diffusion taken into account?
    bool zonal_diffusion_;
    //! Is meridional diffusion taken into account?
    bool meridional_diffusion_;
    //! Is vertical diffusion taken into account?
    bool vertical_diffusion_;

  public:

    /*** Constructor ***/

    DiffusionROS2();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);
    template<class ClassModel>
    void Backward(ClassModel& Model);
    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    template<class ClassModel>
    void ForwardXY(ClassModel& Model);

    template<class ClassModel>
    void ForwardZ(ClassModel& Model);

    void Forward(T current_time, T Delta_t,
		 Data<T, 3>& ZonalDiffusionCoefficient,
		 Data<T, 3>& MeridionalDiffusionCoefficient,
		 Data<T, 3>& VerticalDiffusionCoefficient_i,
		 Data<T, 3>& VerticalDiffusionCoefficient_f,
		 Array<T, 2>& DepositionVelocity_i,
		 Array<T, 2>& DepositionVelocity_f,
		 Array<T, 1>& CellWidth_x,
		 Array<T, 1>& CellWidth_y,
		 Array<T, 1>& CellWidth_z,
		 Array<T, 1>& ModifiedCellCenterDistance_x,
		 Array<T, 1>& ModifiedCellCenterDistance_y,
		 Array<T, 1>& ModifiedCellCenterDistance_z,
		 Array<T, 2>& SurfaceEmission_i,
		 Array<T, 2>& SurfaceEmission_f,
		 int Nx, int Ny, int Nz,
		 Data<T, 3>& AirDensity, Array<T, 3>& Concentration);

    void ForwardParallel(T current_time, T Delta_t,
			 Data<T, 3>& ZonalDiffusionCoefficient,
			 Data<T, 3>& MeridionalDiffusionCoefficient,
			 Data<T, 3>& VerticalDiffusionCoefficient_i,
			 Data<T, 3>& VerticalDiffusionCoefficient_f,
			 Array<T, 2>& DepositionVelocity_i,
			 Array<T, 2>& DepositionVelocity_f,
			 Array<T, 1>& CellWidth_x,
			 Array<T, 1>& CellWidth_y,
			 Array<T, 1>& CellWidth_z,
			 Array<T, 1>& ModifiedCellCenterDistance_x,
			 Array<T, 1>& ModifiedCellCenterDistance_y,
			 Array<T, 1>& ModifiedCellCenterDistance_z,
			 Array<T, 2>& SurfaceEmission_i,
			 Array<T, 2>& SurfaceEmission_f,
			 int Nx, int Ny, int Nz,
			 Data<T, 3>& AirDensity, Array<T, 3>& Concentration);

    void Backward(T current_time, T Delta_t,
		  Data<T, 3>& ZonalDiffusionCoefficient,
		  Data<T, 3>& MeridionalDiffusionCoefficient,
		  Data<T, 3>& VerticalDiffusionCoefficient_i,
		  Data<T, 3>& VerticalDiffusionCoefficient_f,
		  Data<T, 2>& DepositionVelocity_i,
		  Data<T, 2>& DepositionVelocity_f,
		  Array<T, 1>& CellWidth_x,
		  Array<T, 1>& CellWidth_y,
		  Array<T, 1>& CellWidth_z,
		  Array<T, 1>& ModifiedCellCenterDistance_x,
		  Array<T, 1>& ModifiedCellCenterDistance_y,
		  Array<T, 1>& ModifiedCellCenterDistance_z,
		  Data<T, 2>& SurfaceEmission_i,
		  Data<T, 2>& SurfaceEmission_f,
		  int Nx, int Ny, int Nz,
		  Data<T, 3>& AirDensity,
		  Data<T, 3>& Concentration,
		  Data<T, 3>& Concentration_ccl);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_DIFFUSIONROS2_HXX
#endif
