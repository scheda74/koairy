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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_DIFFUSIONROS2_CXX


#include "DiffusionROS2Split.hxx"


namespace Polyphemus
{


  //! Default constructor.
  /*!
    \note Zonal, meridional and vertical diffusion are taken into account in
    the default configuration.
  */
  template<class T>
  DiffusionROS2<T>::DiffusionROS2():
    zonal_diffusion_(true), meridional_diffusion_(true),
    vertical_diffusion_(true)
  {
  }


  //! Initialization of the scheme.
  /*! \note Parallelization parameters are initialized here.
   */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::Init(ClassModel& Model)
  {
    // To be called even if there is no parallelization.
    BaseModuleParallel::Init(Model);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      BuildPartition_s();
    if (ParallelizedOnAerSpecies())
      BuildPartition_s_aer();
    if (!ParallelizedOnSpecies() || !ParallelizedOnAerSpecies())
      {
	// For advection along x.
	BuildPartition_y();
	// For advection along y and z.
	BuildPartition_x();
      }
#endif
  }


  //! Performs an integration over one time step.
  /*!
    \param Model (input/output) model to be updated. Its interface must
    contain:
    <ul>
    <li> GetDelta_t()
    <li> GetNs()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> HasDepositionVelocity(int)
    <li> A3("DepositionVelocity_i")
    <li> A3("DepositionVelocity_f")
    <li> HasSurfaceEmission(int)
    <li> A3("SurfaceEmission_i")
    <li> A3("SurfaceEmission_f")
    <li> D3("ZonalDiffusionCoefficient")
    <li> D3("MeridionalDiffusionCoefficient")
    <li> D3("VerticalDiffusionCoefficient_i")
    <li> D3("VerticalDiffusionCoefficient_f")
    <li> D3("AirDensity_i")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> CellCenterDistance_x
    <li> CellCenterDistance_y
    <li> CellCenterDistance_z
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::Forward(ClassModel& Model)
  {

    int Ns = Model.GetNs();
    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    /*** Mesh ***/

    Array<T, 1> ModifiedCellCenterDistance_x(Nx);
    Array<T, 1> ModifiedCellCenterDistance_y(Ny);
    Array<T, 1> ModifiedCellCenterDistance_z(Nz);

    for (int i = 1; i < Nx; i++)
      ModifiedCellCenterDistance_x(i) = Model.CellCenterDistance_x(i - 1);
    ModifiedCellCenterDistance_x(0) = Model.CellWidth_x(0);

    for (int i = 1; i < Ny; i++)
      ModifiedCellCenterDistance_y(i) = Model.CellCenterDistance_y(i - 1);
    ModifiedCellCenterDistance_y(0) = Model.CellWidth_y(0);

    for (int i = 1; i < Nz; i++)
      ModifiedCellCenterDistance_z(i) = Model.CellCenterDistance_z(i - 1);
    ModifiedCellCenterDistance_z(0) = Model.CellWidth_z(0);

    /*** Diffusion ***/

    zonal_diffusion_
      = (Model.D3("ZonalDiffusionCoefficient").GetMax() > T(0.));
    meridional_diffusion_
      = (Model.D3("MeridionalDiffusionCoefficient").GetMax() > T(0.));

    int first_index_along_s, last_index_along_s;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      {
	ScatterSlice_s_MPI(Model.GetConcentration());
	GetEdgePartition_s(first_index_along_s, last_index_along_s);
      }
    else
      {
	first_index_along_s = 0;
	last_index_along_s = Ns;
      }
#else
    first_index_along_s = 0;
    last_index_along_s = Ns;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)			\
  shared(Ns, Nz, Ny, Nx, ModifiedCellCenterDistance_x,			\
	 ModifiedCellCenterDistance_y, ModifiedCellCenterDistance_z)	\
  firstprivate(first_index_along_s, last_index_along_s)
#endif
    for(int s = first_index_along_s; s < last_index_along_s; s++)
      {

	/*** Deposition velocities ***/

	Array<T, 2> DepositionVelocity_i(shape(Ny, Nx));
	Array<T, 2> DepositionVelocity_f(shape(Ny, Nx));

	if (Model.HasDepositionVelocity(s))
	  {
	    int dep_s = Model.DepositionVelocityIndex(s);
	    DepositionVelocity_i
	      = Model.A3("DepositionVelocity_i")(dep_s, Range::all(),
						 Range::all()).copy();
	    DepositionVelocity_f
	      = Model.A3("DepositionVelocity_f")(dep_s, Range::all(),
						 Range::all()).copy();
	  }
	else
	  {
	    DepositionVelocity_i = 0;
	    DepositionVelocity_f = 0;
	  }

	/*** Surface emissions ***/

	Array<T, 2> SurfaceEmission_i(shape(Ny, Nx));
	Array<T, 2> SurfaceEmission_f(shape(Ny, Nx));

	if (Model.HasSurfaceEmission(s))
	  {
	    int dep_s = Model.SurfaceEmissionIndex(s);
	    SurfaceEmission_i
	      = Model.A3("SurfaceEmission_i")(dep_s, Range::all(),
					      Range::all()).copy();
	    SurfaceEmission_f
	      = Model.A3("SurfaceEmission_f")(dep_s, Range::all(),
					      Range::all()).copy();
	  }
	else
	  {
	    SurfaceEmission_i = 0;
	    SurfaceEmission_f = 0;
	  }

	/*** Concentrations ***/
	Array<T, 3> Concentration_s(&Model.GetConcentration()(s, 0, 0, 0),
				    shape(Nz, Ny, Nx));

	/*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	if (ParallelizedOnSpecies())
#endif
	  Forward(Model.GetCurrentTime(), Model.GetDelta_t(),
		  Model.D3("ZonalDiffusionCoefficient"),
		  Model.D3("MeridionalDiffusionCoefficient"),
		  Model.D3("VerticalDiffusionCoefficient_i"),
		  Model.D3("VerticalDiffusionCoefficient_f"),
		  DepositionVelocity_i, DepositionVelocity_f,
		  Model.CellWidth_x, Model.CellWidth_y, Model.CellWidth_z,
		  ModifiedCellCenterDistance_x,
		  ModifiedCellCenterDistance_y,
		  ModifiedCellCenterDistance_z,
		  SurfaceEmission_i, SurfaceEmission_f,
		  Nx, Ny, Nz, Model.D3("AirDensity_i"), Concentration_s);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	else
	  ForwardParallel(Model.GetCurrentTime(), Model.GetDelta_t(),
			  Model.D3("ZonalDiffusionCoefficient"),
			  Model.D3("MeridionalDiffusionCoefficient"),
			  Model.D3("VerticalDiffusionCoefficient_i"),
			  Model.D3("VerticalDiffusionCoefficient_f"),
			  DepositionVelocity_i, DepositionVelocity_f,
			  Model.CellWidth_x, Model.CellWidth_y,
			  Model.CellWidth_z,
			  ModifiedCellCenterDistance_x,
			  ModifiedCellCenterDistance_y,
			  ModifiedCellCenterDistance_z,
			  SurfaceEmission_i, SurfaceEmission_f,
			  Nx, Ny, Nz, Model.D3("AirDensity_i"),
			  Concentration_s);
#endif
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      GatherSlice_s_MPI(Model.GetConcentration());
#endif

    zonal_diffusion_ = true;
    meridional_diffusion_ = true;
    vertical_diffusion_ = true;
  }

  //! Performs an integration over one time step.
  /*!
    \param Model (input/output) model to be updated. Its interface must
    contain:
    <ul>
    <li> GetDelta_t()
    <li> GetNs()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> HasDepositionVelocity(int)
    <li> A3("DepositionVelocity_i")
    <li> A3("DepositionVelocity_f")
    <li> HasSurfaceEmission(int)
    <li> A3("SurfaceEmission_i")
    <li> A3("SurfaceEmission_f")
    <li> D3("ZonalDiffusionCoefficient")
    <li> D3("MeridionalDiffusionCoefficient")
    <li> D3("VerticalDiffusionCoefficient_i")
    <li> D3("VerticalDiffusionCoefficient_f")
    <li> D3("AirDensity_i")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> CellCenterDistance_x
    <li> CellCenterDistance_y
    <li> CellCenterDistance_z
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::ForwardXY(ClassModel& Model)
  {
    vertical_diffusion_ = false;

    int Ns = Model.GetNs();
    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    /*** Mesh ***/

    Array<T, 1> ModifiedCellCenterDistance_x(Nx);
    Array<T, 1> ModifiedCellCenterDistance_y(Ny);
    Array<T, 1> ModifiedCellCenterDistance_z(Nz);

    for (int i = 1; i < Nx; i++)
      ModifiedCellCenterDistance_x(i) = Model.CellCenterDistance_x(i - 1);
    ModifiedCellCenterDistance_x(0) = Model.CellWidth_x(0);

    for (int i = 1; i < Ny; i++)
      ModifiedCellCenterDistance_y(i) = Model.CellCenterDistance_y(i - 1);
    ModifiedCellCenterDistance_y(0) = Model.CellWidth_y(0);

    for (int i = 1; i < Nz; i++)
      ModifiedCellCenterDistance_z(i) = Model.CellCenterDistance_z(i - 1);
    ModifiedCellCenterDistance_z(0) = Model.CellWidth_z(0);

    /*** Diffusion ***/

    zonal_diffusion_
      = (Model.D3("ZonalDiffusionCoefficient").GetMax() > T(0.));
    meridional_diffusion_
      = (Model.D3("MeridionalDiffusionCoefficient").GetMax() > T(0.));

    int first_index_along_s, last_index_along_s;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      {
	ScatterSlice_s_MPI(Model.GetConcentration());
	GetEdgePartition_s(first_index_along_s, last_index_along_s);
      }
    else
      {
	first_index_along_s = 0;
	last_index_along_s = Ns;
      }
#else
    first_index_along_s = 0;
    last_index_along_s = Ns;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)			\
  shared(Ns, Nz, Ny, Nx, ModifiedCellCenterDistance_x,			\
	 ModifiedCellCenterDistance_y, ModifiedCellCenterDistance_z)	\
  firstprivate(first_index_along_s, last_index_along_s)
#endif
    for(int s = first_index_along_s; s < last_index_along_s; s++)
      {

	/*** Deposition velocities ***/

	Array<T, 2> DepositionVelocity_i(shape(Ny, Nx));
	Array<T, 2> DepositionVelocity_f(shape(Ny, Nx));

	if (Model.HasDepositionVelocity(s))
	  {
	    int dep_s = Model.DepositionVelocityIndex(s);
	    DepositionVelocity_i
	      = Model.A3("DepositionVelocity_i")(dep_s, Range::all(),
						 Range::all()).copy();
	    DepositionVelocity_f
	      = Model.A3("DepositionVelocity_f")(dep_s, Range::all(),
						 Range::all()).copy();
	  }
	else
	  {
	    DepositionVelocity_i = 0;
	    DepositionVelocity_f = 0;
	  }

	/*** Surface emissions ***/

	Array<T, 2> SurfaceEmission_i(shape(Ny, Nx));
	Array<T, 2> SurfaceEmission_f(shape(Ny, Nx));

	if (Model.HasSurfaceEmission(s))
	  {
	    int dep_s = Model.SurfaceEmissionIndex(s);
	    SurfaceEmission_i
	      = Model.A3("SurfaceEmission_i")(dep_s, Range::all(),
					      Range::all()).copy();
	    SurfaceEmission_f
	      = Model.A3("SurfaceEmission_f")(dep_s, Range::all(),
					      Range::all()).copy();
	  }
	else
	  {
	    SurfaceEmission_i = 0;
	    SurfaceEmission_f = 0;
	  }

	/*** Concentrations ***/
	Array<T, 3> Concentration_s(&Model.GetConcentration()(s, 0, 0, 0),
				    shape(Nz, Ny, Nx));

	/*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	if (ParallelizedOnSpecies())
#endif
	  Forward(Model.GetCurrentTime(), Model.GetDelta_t(),
		  Model.D3("ZonalDiffusionCoefficient"),
		  Model.D3("MeridionalDiffusionCoefficient"),
		  Model.D3("VerticalDiffusionCoefficient_i"),
		  Model.D3("VerticalDiffusionCoefficient_f"),
		  DepositionVelocity_i, DepositionVelocity_f,
		  Model.CellWidth_x, Model.CellWidth_y, Model.CellWidth_z,
		  ModifiedCellCenterDistance_x,
		  ModifiedCellCenterDistance_y,
		  ModifiedCellCenterDistance_z,
		  SurfaceEmission_i, SurfaceEmission_f,
		  Nx, Ny, Nz, Model.D3("AirDensity_i"), Concentration_s);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	else
	  ForwardParallel(Model.GetCurrentTime(), Model.GetDelta_t(),
			  Model.D3("ZonalDiffusionCoefficient"),
			  Model.D3("MeridionalDiffusionCoefficient"),
			  Model.D3("VerticalDiffusionCoefficient_i"),
			  Model.D3("VerticalDiffusionCoefficient_f"),
			  DepositionVelocity_i, DepositionVelocity_f,
			  Model.CellWidth_x, Model.CellWidth_y,
			  Model.CellWidth_z,
			  ModifiedCellCenterDistance_x,
			  ModifiedCellCenterDistance_y,
			  ModifiedCellCenterDistance_z,
			  SurfaceEmission_i, SurfaceEmission_f,
			  Nx, Ny, Nz, Model.D3("AirDensity_i"),
			  Concentration_s);
#endif
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      GatherSlice_s_MPI(Model.GetConcentration());
#endif

    zonal_diffusion_ = true;
    meridional_diffusion_ = true;
    vertical_diffusion_ = true;
  }

  //! Performs an integration over one time step.
  /*!
    \param Model (input/output) model to be updated. Its interface must
    contain:
    <ul>
    <li> GetDelta_t()
    <li> GetNs()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> HasDepositionVelocity(int)
    <li> A3("DepositionVelocity_i")
    <li> A3("DepositionVelocity_f")
    <li> HasSurfaceEmission(int)
    <li> A3("SurfaceEmission_i")
    <li> A3("SurfaceEmission_f")
    <li> D3("ZonalDiffusionCoefficient")
    <li> D3("MeridionalDiffusionCoefficient")
    <li> D3("VerticalDiffusionCoefficient_i")
    <li> D3("VerticalDiffusionCoefficient_f")
    <li> D3("AirDensity_i")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> CellCenterDistance_x
    <li> CellCenterDistance_y
    <li> CellCenterDistance_z
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::ForwardZ(ClassModel& Model)
  {

    zonal_diffusion_ = false;
    meridional_diffusion_ = false;
    vertical_diffusion_ = true;

    int Ns = Model.GetNs();
    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    /*** Mesh ***/

    Array<T, 1> ModifiedCellCenterDistance_x(Nx);
    Array<T, 1> ModifiedCellCenterDistance_y(Ny);
    Array<T, 1> ModifiedCellCenterDistance_z(Nz);

    for (int i = 1; i < Nx; i++)
      ModifiedCellCenterDistance_x(i) = Model.CellCenterDistance_x(i - 1);
    ModifiedCellCenterDistance_x(0) = Model.CellWidth_x(0);

    for (int i = 1; i < Ny; i++)
      ModifiedCellCenterDistance_y(i) = Model.CellCenterDistance_y(i - 1);
    ModifiedCellCenterDistance_y(0) = Model.CellWidth_y(0);

    for (int i = 1; i < Nz; i++)
      ModifiedCellCenterDistance_z(i) = Model.CellCenterDistance_z(i - 1);
    ModifiedCellCenterDistance_z(0) = Model.CellWidth_z(0);

    int first_index_along_s, last_index_along_s;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      {
	ScatterSlice_s_MPI(Model.GetConcentration());
	GetEdgePartition_s(first_index_along_s, last_index_along_s);
      }
    else
      {
	first_index_along_s = 0;
	last_index_along_s = Ns;
      }
#else
    first_index_along_s = 0;
    last_index_along_s = Ns;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)			\
  shared(Ns, Nz, Ny, Nx, ModifiedCellCenterDistance_x,			\
	 ModifiedCellCenterDistance_y, ModifiedCellCenterDistance_z)	\
  firstprivate(first_index_along_s, last_index_along_s)
#endif
    for(int s = first_index_along_s; s < last_index_along_s; s++)
      {

	/*** Deposition velocities ***/

	Array<T, 2> DepositionVelocity_i(shape(Ny, Nx));
	Array<T, 2> DepositionVelocity_f(shape(Ny, Nx));

	if (Model.HasDepositionVelocity(s))
	  {
	    int dep_s = Model.DepositionVelocityIndex(s);
	    DepositionVelocity_i
	      = Model.A3("DepositionVelocity_i")(dep_s, Range::all(),
						 Range::all()).copy();
	    DepositionVelocity_f
	      = Model.A3("DepositionVelocity_f")(dep_s, Range::all(),
						 Range::all()).copy();
	  }
	else
	  {
	    DepositionVelocity_i = 0;
	    DepositionVelocity_f = 0;
	  }

	/*** Surface emissions ***/

	Array<T, 2> SurfaceEmission_i(shape(Ny, Nx));
	Array<T, 2> SurfaceEmission_f(shape(Ny, Nx));

	if (Model.HasSurfaceEmission(s))
	  {
	    int dep_s = Model.SurfaceEmissionIndex(s);
	    SurfaceEmission_i
	      = Model.A3("SurfaceEmission_i")(dep_s, Range::all(),
					      Range::all()).copy();
	    SurfaceEmission_f
	      = Model.A3("SurfaceEmission_f")(dep_s, Range::all(),
					      Range::all()).copy();
	  }
	else
	  {
	    SurfaceEmission_i = 0;
	    SurfaceEmission_f = 0;
	  }

	/*** Concentrations ***/
	Array<T, 3> Concentration_s(&Model.GetConcentration()(s, 0, 0, 0),
				    shape(Nz, Ny, Nx));

	/*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	if (ParallelizedOnSpecies())
#endif
	  Forward(Model.GetCurrentTime(), Model.GetDelta_t(),
		  Model.D3("ZonalDiffusionCoefficient"),
		  Model.D3("MeridionalDiffusionCoefficient"),
		  Model.D3("VerticalDiffusionCoefficient_i"),
		  Model.D3("VerticalDiffusionCoefficient_f"),
		  DepositionVelocity_i, DepositionVelocity_f,
		  Model.CellWidth_x, Model.CellWidth_y, Model.CellWidth_z,
		  ModifiedCellCenterDistance_x,
		  ModifiedCellCenterDistance_y,
		  ModifiedCellCenterDistance_z,
		  SurfaceEmission_i, SurfaceEmission_f,
		  Nx, Ny, Nz, Model.D3("AirDensity_i"), Concentration_s);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	else
	  ForwardParallel(Model.GetCurrentTime(), Model.GetDelta_t(),
			  Model.D3("ZonalDiffusionCoefficient"),
			  Model.D3("MeridionalDiffusionCoefficient"),
			  Model.D3("VerticalDiffusionCoefficient_i"),
			  Model.D3("VerticalDiffusionCoefficient_f"),
			  DepositionVelocity_i, DepositionVelocity_f,
			  Model.CellWidth_x, Model.CellWidth_y,
			  Model.CellWidth_z,
			  ModifiedCellCenterDistance_x,
			  ModifiedCellCenterDistance_y,
			  ModifiedCellCenterDistance_z,
			  SurfaceEmission_i, SurfaceEmission_f,
			  Nx, Ny, Nz, Model.D3("AirDensity_i"),
			  Concentration_s);
#endif
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      GatherSlice_s_MPI(Model.GetConcentration());
#endif

    zonal_diffusion_ = true;
    meridional_diffusion_ = true;
    vertical_diffusion_ = true;
  }

  //! Performs one step of backward integration of adjoint model.
  /*!
    \param Model (input/output) model to be updated. Its interface must
    contain:
    <ul>
    <li> GetDelta_t()
    <li> GetNs()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> HasDepositionVelocity(int)
    <li> A3("DepositionVelocity_i")
    <li> A3("DepositionVelocity_f")
    <li> HasSurfaceEmission(int)
    <li> A3("SurfaceEmission_i")
    <li> A3("SurfaceEmission_f")
    <li> D3("ZonalDiffusionCoefficient")
    <li> D3("MeridionalDiffusionCoefficient")
    <li> D3("VerticalDiffusionCoefficient_i")
    <li> D3("VerticalDiffusionCoefficient_f")
    <li> D3("AirDensity_i")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> CellCenterDistance_x
    <li> CellCenterDistance_y
    <li> CellCenterDistance_z
    <li> GetConcentration()
    <li> GetConcentration_ccl()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::Backward(ClassModel& Model)
  {

    /*** Mesh ***/

    Array<T, 1> ModifiedCellCenterDistance_x(Model.GetNx());
    Array<T, 1> ModifiedCellCenterDistance_y(Model.GetNy());
    Array<T, 1> ModifiedCellCenterDistance_z(Model.GetNz());

    for (int i = 1; i < Model.GetNx(); i++)
      ModifiedCellCenterDistance_x(i) = Model.CellCenterDistance_x(i - 1);
    ModifiedCellCenterDistance_x(0) = Model.CellWidth_x(0);

    for (int i = 1; i < Model.GetNy(); i++)
      ModifiedCellCenterDistance_y(i) = Model.CellCenterDistance_y(i - 1);
    ModifiedCellCenterDistance_y(0) = Model.CellWidth_y(0);

    for (int i = 1; i < Model.GetNz(); i++)
      ModifiedCellCenterDistance_z(i) = Model.CellCenterDistance_z(i - 1);
    ModifiedCellCenterDistance_z(0) = Model.CellWidth_z(0);

    Data<T, 2> DepositionVelocity_i(shape(Model.GetNy(), Model.GetNx()));
    Data<T, 2> DepositionVelocity_f(shape(Model.GetNy(), Model.GetNx()));

    Data<T, 2> SurfaceEmission_i(shape(Model.GetNy(), Model.GetNx()));
    Data<T, 2> SurfaceEmission_f(shape(Model.GetNy(), Model.GetNx()));

    /*** Diffusion ***/

    zonal_diffusion_
      = (Model.D3("ZonalDiffusionCoefficient").GetMax() > T(0.));
    meridional_diffusion_
      = (Model.D3("MeridionalDiffusionCoefficient").GetMax() > T(0.));
    vertical_diffusion_
      = (Model.D3("VerticalDiffusionCoefficient_i").GetMax() > T(0.))
      && (Model.D3("VerticalDiffusionCoefficient_f").GetMax() > T(0.));

    for(int s = 0; s < Model.GetNs(); s++)
      {

	/*** Deposition velocities ***/

	if (Model.HasDepositionVelocity(s))
	  {
	    int dep_s = Model.DepositionVelocityIndex(s);
	    DepositionVelocity_i.GetArray()
	      = Model.A3("DepositionVelocity_i")(dep_s, Range::all(),
						 Range::all()).copy();
	    DepositionVelocity_f.GetArray()
	      = Model.A3("DepositionVelocity_f")(dep_s, Range::all(),
						 Range::all()).copy();
	  }
	else
	  {
	    DepositionVelocity_i.SetZero();
	    DepositionVelocity_f.SetZero();
	  }

	/*** Surface emissions ***/

	if (Model.HasSurfaceEmission(s))
	  {
	    int dep_s = Model.SurfaceEmissionIndex(s);
	    SurfaceEmission_i.GetArray()
	      = Model.A3("SurfaceEmission_i")(dep_s, Range::all(),
					      Range::all()).copy();
	    SurfaceEmission_f.GetArray()
	      = Model.A3("SurfaceEmission_f")(dep_s, Range::all(),
					      Range::all()).copy();
	  }
	else
	  {
	    SurfaceEmission_i.SetZero();
	    SurfaceEmission_f.SetZero();
	  }

	/*** Concentrations ***/

	Data<T, 3> Concentration(&Model.GetConcentration()(s, 0, 0, 0),
				 shape(Model.GetNz(), Model.GetNy(),
				       Model.GetNx()));
	Data<T, 3> Concentration_ccl(&Model.
				     GetConcentration_ccl()(s, 0, 0, 0),
				     shape(Model.GetNz(), Model.GetNy(),
					   Model.GetNx()));

	/*** Numerical integration ***/

	Backward(Model.GetCurrentTime(), Model.GetDelta_t(),
		 Model.D3("ZonalDiffusionCoefficient"),
		 Model.D3("MeridionalDiffusionCoefficient"),
		 Model.D3("VerticalDiffusionCoefficient_i"),
		 Model.D3("VerticalDiffusionCoefficient_f"),
		 DepositionVelocity_i, DepositionVelocity_f,
		 Model.CellWidth_x, Model.CellWidth_y, Model.CellWidth_z,
		 ModifiedCellCenterDistance_x,
		 ModifiedCellCenterDistance_y,
		 ModifiedCellCenterDistance_z,
		 SurfaceEmission_i, SurfaceEmission_f,
		 Model.GetNx(), Model.GetNy(), Model.GetNz(),
		 Model.D3("AirDensity_i"),
		 Concentration, Concentration_ccl);
      }

    zonal_diffusion_ = true;
    meridional_diffusion_ = true;
    vertical_diffusion_ = true;
  }



  //! Performs an integration over one time step for aerosols.
  /*!
    \param Model (input/output) model to be updated. Its interface must
    contain:
    <ul>
    <li> GetDelta_t()
    <li> GetNs_aer()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> HasDepositionVelocity_aer(int)
    <li> A4("DepositionVelocity_aer_i")
    <li> A4("DepositionVelocity_aer_f")
    <li> HasSurfaceEmission_aer(int)
    <li> A3("SurfaceEmission_aer_i")
    <li> A3("SurfaceEmission_aer_f")
    <li> D3("ZonalDiffusionCoefficient")
    <li> D3("MeridionalDiffusionCoefficient")
    <li> D3("VerticalDiffusionCoefficient_i")
    <li> D3("VerticalDiffusionCoefficient_f")
    <li> D3("AirDensity_i")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> CellCenterDistance_x
    <li> CellCenterDistance_y
    <li> CellCenterDistance_z
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void DiffusionROS2<T>::Forward_aer(ClassModel& Model)
  {

    /*** Mesh ***/

    Array<T, 1> ModifiedCellCenterDistance_x(Model.GetNx());
    Array<T, 1> ModifiedCellCenterDistance_y(Model.GetNy());
    Array<T, 1> ModifiedCellCenterDistance_z(Model.GetNz());

    for (int i = 1; i < Model.GetNx(); i++)
      ModifiedCellCenterDistance_x(i) = Model.CellCenterDistance_x(i - 1);
    ModifiedCellCenterDistance_x(0) = Model.CellWidth_x(0);

    for (int i = 1; i < Model.GetNy(); i++)
      ModifiedCellCenterDistance_y(i) = Model.CellCenterDistance_y(i - 1);
    ModifiedCellCenterDistance_y(0) = Model.CellWidth_y(0);

    for (int i = 1; i < Model.GetNz(); i++)
      ModifiedCellCenterDistance_z(i) = Model.CellCenterDistance_z(i - 1);
    ModifiedCellCenterDistance_z(0) = Model.CellWidth_z(0);

    /*** Diffusion ***/

    zonal_diffusion_
      = (Model.D3("ZonalDiffusionCoefficient").GetMax() > T(0.));
    meridional_diffusion_
      = (Model.D3("MeridionalDiffusionCoefficient").GetMax() > T(0.));
    vertical_diffusion_
      = (Model.D3("VerticalDiffusionCoefficient_i").GetMax() > T(0.))
      && (Model.D3("VerticalDiffusionCoefficient_f").GetMax() > T(0.));

    int first_index_along_s_aer, last_index_along_s_aer;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnAerSpecies())
      {
	ScatterSlice_s_aer_MPI(Model.GetConcentration_aer());
	GetEdgePartition_s_aer(first_index_along_s_aer,
			       last_index_along_s_aer);
      }
    else
      {
	first_index_along_s_aer = 0;
	last_index_along_s_aer = Model.GetNs_aer();
      }
#else
    first_index_along_s_aer = 0;
    last_index_along_s_aer = Model.GetNs_aer();
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)			\
  shared(ModifiedCellCenterDistance_x, ModifiedCellCenterDistance_y,	\
	 ModifiedCellCenterDistance_z)					\
  firstprivate(first_index_along_s_aer, last_index_along_s_aer)
#endif
    for(int s = first_index_along_s_aer; s < last_index_along_s_aer; s++)
      for(int b = 0; b < Model.GetNbin_aer(); b++)
	{

	  /*** Deposition velocities ***/

	  Array<T, 2> DepositionVelocity_i(shape(Model.GetNy(),
						 Model.GetNx()));
	  Array<T, 2> DepositionVelocity_f(shape(Model.GetNy(),
						 Model.GetNx()));
	  int dep_b;

	  if (Model.HasDepositionVelocity_aer(b))
	    {
	      dep_b = Model.DepositionVelocityIndex_aer(b);
	      DepositionVelocity_i
		= Model.A3("DepositionVelocity_aer_i")(dep_b, Range::all(),
						       Range::all()).copy();
	      DepositionVelocity_f
		= Model.A3("DepositionVelocity_aer_f")(dep_b, Range::all(),
						       Range::all()).copy();
	    }
	  else
	    {
	      DepositionVelocity_i = 0;
	      DepositionVelocity_f = 0;
	    }

	  /*** Surface emissions ***/

	  Array<T, 2> SurfaceEmission_i(shape(Model.GetNy(), Model.GetNx()));
	  Array<T, 2> SurfaceEmission_f(shape(Model.GetNy(), Model.GetNx()));
	  int emis_s, emis_b;
	  vector<int> index;

	  if (Model.HasSurfaceEmission_aer(s, b))
	    {
	      index = Model.SurfaceEmissionIndex_aer(s, b);
	      emis_s = index[0];
	      emis_b = index[1];
	      SurfaceEmission_i
		= Model.A4("SurfaceEmission_aer_i")(emis_s, emis_b,
						    Range::all(),
						    Range::all()).copy();
	      SurfaceEmission_f
		= Model.A4("SurfaceEmission_aer_f")(emis_s, emis_b,
						    Range::all(),
						    Range::all()).copy();
	    }
	  else
	    {
	      SurfaceEmission_i = 0;
	      SurfaceEmission_f = 0;
	    }

	  /*** Concentrations ***/

	  Array<T, 3> Concentration_aer_s(&Model.GetConcentration_aer()
					  (s, b, 0, 0, 0),
					  shape(Model.GetNz(), Model.GetNy(),
						Model.GetNx()));

	  /*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	  if (ParallelizedOnAerSpecies())
#endif
	    Forward(Model.GetCurrentTime(), Model.GetDelta_t(),
		    Model.D3("ZonalDiffusionCoefficient"),
		    Model.D3("MeridionalDiffusionCoefficient"),
		    Model.D3("VerticalDiffusionCoefficient_i"),
		    Model.D3("VerticalDiffusionCoefficient_f"),
		    DepositionVelocity_i, DepositionVelocity_f,
		    Model.CellWidth_x, Model.CellWidth_y, Model.CellWidth_z,
		    ModifiedCellCenterDistance_x,
		    ModifiedCellCenterDistance_y,
		    ModifiedCellCenterDistance_z,
		    SurfaceEmission_i, SurfaceEmission_f,
		    Model.GetNx(), Model.GetNy(), Model.GetNz(),
		    Model.D3("AirDensity_i"), Concentration_aer_s);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	  else
	    ForwardParallel(Model.GetCurrentTime(), Model.GetDelta_t(),
			    Model.D3("ZonalDiffusionCoefficient"),
			    Model.D3("MeridionalDiffusionCoefficient"),
			    Model.D3("VerticalDiffusionCoefficient_i"),
			    Model.D3("VerticalDiffusionCoefficient_f"),
			    DepositionVelocity_i, DepositionVelocity_f,
			    Model.CellWidth_x, Model.CellWidth_y,
			    Model.CellWidth_z,
			    ModifiedCellCenterDistance_x,
			    ModifiedCellCenterDistance_y,
			    ModifiedCellCenterDistance_z,
			    SurfaceEmission_i, SurfaceEmission_f,
			    Model.GetNx(), Model.GetNy(), Model.GetNz(),
			    Model.D3("AirDensity_i"), Concentration_aer_s);
#endif
	}

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnAerSpecies())
      GatherSlice_s_aer_MPI(Model.GetConcentration_aer());
#endif

    zonal_diffusion_ = true;
    meridional_diffusion_ = true;
    vertical_diffusion_ = true;
  }

  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Delta_t time step.
    \param ZonalDiffusionCoefficient zonal diffusion coefficients.
    \param MeridionalDiffusionCoefficient meridional diffusion coefficients.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the beginning of the time step.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the end of the time step.
    \param DepositionVelocity_i deposition velocities at the beginning of the
    time step.
    \param DepositionVelocity_f deposition velocities at the end of the time
    step.
    \param CellWidth_x cells widths along x in meters.
    \param CellWidth_y cells widths along y in meters.
    \param CellWidth_z cells widths along z in meters.
    \param ModifiedCellCenterDistance_x cells centers distances along x
    in meters.
    \param ModifiedCellCenterDistance_y cells centers distances along y
    in meters.
    \param ModifiedCellCenterDistance_z cells centers distances along z
    in meters.
    \param SurfaceEmission_i surface emissions at the beginning of the time
    step.
    \param SurfaceEmission_f surface emissions at the end of the time step.
    \param Nx number of cells along x.
    \param Ny number of cells along y.
    \param Nz number of cells along z.
    \param AirDensity air density.
    \param Concentration concentrations.
  */
  template<class T>
  void DiffusionROS2<T>
  ::Forward(T current_time, T Delta_t,
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
	    Data<T, 3>& AirDensity, Array<T, 3>& Concentration)
  {
    T final_time = current_time + Delta_t;

    // cout << "\t $$$$ before _diffX in DiffRos2: " << Concentration(0,20,5) << endl;

    // cout << "x, y, z: " << zonal_diffusion_ << " " << meridional_diffusion_ << " " << vertical_diffusion_ << endl;
    // cout << " current time: " << current_time << endl;

    if (zonal_diffusion_)
      _diffX(&Nx, &Ny, &Nz, &current_time, &final_time,
	     ZonalDiffusionCoefficient.GetData(),
	     ModifiedCellCenterDistance_x.data(), CellWidth_x.data(),
	     AirDensity.GetData(), Concentration.data());

    // cout << "\t $$$$ after _diffX in DiffRos2: " << Concentration(0,20,5) << endl;

    if (meridional_diffusion_)
      _diffY(&Nx, &Ny, &Nz, &current_time, &final_time,
	     MeridionalDiffusionCoefficient.GetData(),
	     ModifiedCellCenterDistance_y.data(), CellWidth_y.data(),
	     AirDensity.GetData(), Concentration.data());

    // cout << "\t $$$$ after _diffY in DiffRos2: " << Concentration(0,20,5) << endl;

    if (vertical_diffusion_)
      _diffZ(&Nx, &Ny, &Nz, &current_time, &final_time,
	     VerticalDiffusionCoefficient_i.GetData(),
	     DepositionVelocity_i.data(), SurfaceEmission_i.data(),
	     VerticalDiffusionCoefficient_f.GetData(),
	     DepositionVelocity_f.data(), SurfaceEmission_f.data(),
	     ModifiedCellCenterDistance_z.data(), CellWidth_z.data(),
	     AirDensity.GetData(), Concentration.data());

    // cout << "\t $$$$ after _diffZ in DiffRos2: " << Concentration(0,20,5) << endl;
  }

  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Delta_t time step.
    \param ZonalDiffusionCoefficient zonal diffusion coefficients.
    \param MeridionalDiffusionCoefficient meridional diffusion coefficients.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the beginning of the time step.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the end of the time step.
    \param DepositionVelocity_i deposition velocities at the beginning of the
    time step.
    \param DepositionVelocity_f deposition velocities at the end of the time
    step.
    \param CellWidth_x cells widths along x in meters.
    \param CellWidth_y cells widths along y in meters.
    \param CellWidth_z cells widths along z in meters.
    \param ModifiedCellCenterDistance_x cells centers distances along x
    in meters.
    \param ModifiedCellCenterDistance_y cells centers distances along y
    in meters.
    \param ModifiedCellCenterDistance_z cells centers distances along z
    in meters.
    \param SurfaceEmission_i surface emissions at the beginning of the time
    step.
    \param SurfaceEmission_f surface emissions at the end of the time step.
    \param Nx number of cells along x.
    \param Ny number of cells along y.
    \param Nz number of cells along z.
    \param AirDensity air density.
    \param Concentration concentrations.
  */
  template<class T>
  void DiffusionROS2<T>
  ::ForwardParallel(T current_time, T Delta_t,
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
		    Data<T, 3>& AirDensity, Array<T, 3>& Concentration)
  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    T final_time = current_time + Delta_t;
    const T factor = (1. + 1./sqrt(2.)) * Delta_t;

    int k, j, i;
    int first_index_along_x, first_index_along_y;
    int last_index_along_x, last_index_along_y;
    Array<T, 1> Rho(Nx);
    Array<T, 1> Concentration1D(Nx);
    Array<T, 1> DiffusionCoefficient(Nx + 1);

    // The jacobian matrix is never computed inside rosdiffX and rosdiffY.
    // It is computed once before the loop for a particular cell and
    // kept at this value thereafter.
    int literdiff = 0;

    if (zonal_diffusion_)
      {
	Array<T, 1> Matrix_x(Nx);
	Array<T, 1> Matrix_x_lower(Nx);
	Array<T, 1> Matrix_x_upper(Nx);

	// Computation of the jacobian matrix.
	Array<T, 1> Matrix_x_tmp(Nx);
	Array<T, 1> Matrix_xl_tmp(Nx);
	Array<T, 1> Matrix_xu_tmp(Nx);
	
	Rho = AirDensity()(0, 0, Range::all());
	DiffusionCoefficient =
	  ZonalDiffusionCoefficient()(0, 0, Range::all());

	_jacdiffX(&Nx, ModifiedCellCenterDistance_x.data(),
		  CellWidth_x.data(), DiffusionCoefficient.data(),
		  Rho.data(), Matrix_xl_tmp.data(), Matrix_x_tmp.data(),
		  Matrix_xu_tmp.data());

	for (i = 0; i < Nx; i++)
	  {
	    Matrix_x(i) = 1. - factor *  Matrix_x_tmp(i);
	    Matrix_x_lower(i) = - factor * Matrix_xl_tmp(i);
	    Matrix_x_upper(i)= - factor * Matrix_xu_tmp(i);
	  }

	if (Nx != 1)
	  {
	    // Buffers used to send and receive concentrations from process
	    // #0 to the other processes.
	    Array<T,3> OrdConc;
	    ScatterSlice_y_MPI(Concentration, OrdConc);
	    GetEdgePartition_y(first_index_along_y, last_index_along_y);
	    for (j = first_index_along_y; j < last_index_along_y; j++)
	      for (k = 0; k < Nz; k++)
		{
		  Rho = AirDensity()(k, j, Range::all());
		  DiffusionCoefficient =
		    ZonalDiffusionCoefficient()(k, j, Range::all());
		
		  // Rank of the indices of OrdConc has to match the one
		  // defined in ScatterSlice_y_MPI.
		  for (i = 0; i < Nx; i++)
		    Concentration1D(i) = OrdConc(j, k, i);

		  _rosdiffX(&Nx, ModifiedCellCenterDistance_x.data(),
			    CellWidth_x.data(), Matrix_x_lower.data(),
			    Matrix_x.data(), Matrix_x_upper.data(),
			    DiffusionCoefficient.data(), Rho.data(),
			    Concentration1D.data(), &current_time,
			    &final_time, &literdiff);

		  // Rank of the indices of OrdConc has to match the one
		  // defined in ScatterSlice_y_MPI.
		  OrdConc(j, k, Range::all()) = Concentration1D;
		}

	    GatherSlice_y_MPI(OrdConc, Concentration);
	  }
      }

    if (meridional_diffusion_)
      {
	Rho.resize(Ny);
	Concentration1D.resize(Ny);
	DiffusionCoefficient.resize(Ny + 1);

	Array<T, 1> Matrix_y(Ny);
	Array<T, 1> Matrix_y_lower(Ny);
	Array<T, 1> Matrix_y_upper(Ny);

	// Computation of the jacobian matrix.
	Array<T, 1> Matrix_y_tmp(Ny);
	Array<T, 1> Matrix_yl_tmp(Ny);
	Array<T, 1> Matrix_yu_tmp(Ny);

	Rho = AirDensity()(0,Range::all(), 0);
	DiffusionCoefficient =
	  MeridionalDiffusionCoefficient()(0, Range::all(), 0);

	_jacdiffY(&Ny, ModifiedCellCenterDistance_y.data(),
		  CellWidth_y.data(), DiffusionCoefficient.data(),
		  Rho.data(), Matrix_yl_tmp.data(), Matrix_y_tmp.data(),
		  Matrix_yu_tmp.data());

	for (j = 0; j < Ny; j++)
	  {
	    Matrix_y(j) = 1. - factor *  Matrix_y_tmp(j);
	    Matrix_y_lower(j) = - factor * Matrix_yl_tmp(j);
	    Matrix_y_upper(j)= - factor * Matrix_yu_tmp(j);
	  }

	if (Ny != 1)
	  {
	    // Buffers used to send and receive concentrations from process
	    // #0 to the other processes.
	    Array<T,3> OrdConc;
	    ScatterSlice_x_MPI(Concentration, OrdConc);
	    GetEdgePartition_x(first_index_along_x, last_index_along_x);
	    for (i = first_index_along_x; i < last_index_along_x; i++)
	      for (k = 0; k < Nz; k++)
		{
		  Rho = AirDensity()(k, Range::all(), i);
		  DiffusionCoefficient =
		    MeridionalDiffusionCoefficient()(k, Range::all(), i);

		  // Rank of the indices of OrdConc has to match the one
		  // defined in ScatterSlice_x_MPI.
		  for (j = 0; j < Ny; j++)
		    Concentration1D(j) = OrdConc(i, k, j);

		  _rosdiffY(&Ny, ModifiedCellCenterDistance_y.data(),
			    CellWidth_y.data(), Matrix_y_lower.data(),
			    Matrix_y.data(), Matrix_y_upper.data(),
			    DiffusionCoefficient.data(), Rho.data(),
			    Concentration1D.data(), &current_time,
			    &final_time, &literdiff);

		  // Rank of the indices of OrdConc has to match the one
		  // defined in ScatterSlice_x_MPI.
		  OrdConc(i, k, Range::all()) = Concentration1D;
		}
	    GatherSlice_x_MPI(OrdConc, Concentration);
	  }
      }

    if (vertical_diffusion_)
      {
	Rho.resize(Nz);
	Concentration1D.resize(Nz);
	DiffusionCoefficient.resize(Nz + 1);
	Array<T, 1> DiffusionCoefficient_f(Nz + 1);

	// Buffers used to send and receive concentrations from process
	// #0 to the other processes.
	Array<T,3> OrdConc;
	ScatterSlice_x_MPI(Concentration, OrdConc);
	GetEdgePartition_x(first_index_along_x, last_index_along_x);
	for (i = first_index_along_x; i < last_index_along_x; i++)
	  for (j = 0; j < Ny; j++)
	    {
	      Rho = AirDensity()(Range::all(), j, i);
	      DiffusionCoefficient =
		VerticalDiffusionCoefficient_i()(Range::all(), j, i);
	      DiffusionCoefficient_f =
		VerticalDiffusionCoefficient_f()(Range::all(), j, i);

	      // Rank of the indices of OrdConc has to match the one
	      // defined in ScatterSlice_x_MPI.
	      for (k = 0; k < Nz; k++)
		Concentration1D(k) = OrdConc(i, k, j);

	      _rosdiffZ(&Nz, Concentration1D.data(), &current_time,
			&final_time, DiffusionCoefficient.data(),
			DiffusionCoefficient_f.data(),
			CellWidth_z.data(),
			ModifiedCellCenterDistance_z.data(),
			&DepositionVelocity_i(j, i),
			&SurfaceEmission_i(j, i),
			&DepositionVelocity_f(j, i), &SurfaceEmission_f(j, i),
			Rho.data());


	      // Rank of the indices of OrdConc has to match the one
	      // defined in ScatterSlice_x_MPI.
	      OrdConc(i, Range::all(), j) = Concentration1D;
	    }
	GatherSlice_x_MPI(OrdConc, Concentration);
      }
#endif
  }


  //! Performs one step of backward integration of adjoint model.
  /*!
    \param current_time starting time in seconds.
    \param Delta_t time step.
    \param ZonalDiffusionCoefficient zonal diffusion coefficients.
    \param MeridionalDiffusionCoefficient meridional diffusion coefficients.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the beginning of the time step.
    \param VerticalDiffusionCoefficient_f vertical diffusion coefficients at
    the end of the time step.
    \param DepositionVelocity_i deposition velocities at the beginning of the
    time step.
    \param DepositionVelocity_f deposition velocities at the end of the time
    step.
    \param CellWidth_x cells widths along x in meters.
    \param CellWidth_y cells widths along y in meters.
    \param CellWidth_z cells widths along z in meters.
    \param ModifiedCellCenterDistance_x cells centers distances along x
    in meters.
    \param ModifiedCellCenterDistance_y cells centers distances along y
    in meters.
    \param ModifiedCellCenterDistance_z cells centers distances along z
    in meters.
    \param SurfaceEmission_i surface emissions at the beginning of the time
    step.
    \param SurfaceEmission_f surface emissions at the end of the time step.
    \param Nx number of cells along x.
    \param Ny number of cells along y.
    \param Nz number of cells along z.
    \param AirDensity air density.
    \param Concentration concentrations.
    \param Concentration_ccl adjoinot concentration data.
  */
  template<class T>
  void DiffusionROS2<T>
  ::Backward(T current_time, T Delta_t,
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
	     Data<T, 3>& Concentration_ccl)
  {
    T final_time = current_time + Delta_t;

    Array<T, 3> conc_sav_x(Concentration.GetArray().extent(0),
			   Concentration.GetArray().extent(1),
			   Concentration.GetArray().extent(2));
    Array<T, 3> conc_sav_y(Concentration.GetArray().extent(0),
			   Concentration.GetArray().extent(1),
			   Concentration.GetArray().extent(2));
    Array<T, 3> conc_sav_z(Concentration.GetArray().extent(0),
			   Concentration.GetArray().extent(1),
			   Concentration.GetArray().extent(2));
    Array<T, 3> res(Concentration.GetArray().extent(0),
		    Concentration.GetArray().extent(1),
		    Concentration.GetArray().extent(2));

    if (zonal_diffusion_)
      {
	conc_sav_x = Concentration.GetArray();
	_diffX(&Nx, &Ny, &Nz, &current_time, &final_time,
	       ZonalDiffusionCoefficient.GetData(),
	       ModifiedCellCenterDistance_x.data(), CellWidth_x.data(),
	       AirDensity.GetData(), Concentration.GetData());
      }

    if (meridional_diffusion_)
      {
	conc_sav_y = Concentration.GetArray();
	_diffY(&Nx, &Ny, &Nz, &current_time, &final_time,
	       MeridionalDiffusionCoefficient.GetData(),
	       ModifiedCellCenterDistance_y.data(), CellWidth_y.data(),
	       AirDensity.GetData(), Concentration.GetData());
      }

    if (vertical_diffusion_)
      {
	conc_sav_z = Concentration.GetArray();
	_diffZ(&Nx, &Ny, &Nz, &current_time, &final_time,
	       VerticalDiffusionCoefficient_i.GetData(),
	       DepositionVelocity_i.GetData(), SurfaceEmission_i.GetData(),
	       VerticalDiffusionCoefficient_f.GetData(),
	       DepositionVelocity_f.GetData(), SurfaceEmission_f.GetData(),
	       ModifiedCellCenterDistance_z.data(), CellWidth_z.data(),
	       AirDensity.GetData(), Concentration.GetData());
      }

    res = Concentration.GetArray();

    if (vertical_diffusion_)
      {
	Concentration.GetArray() = conc_sav_z;
	_diffZcl(&Nx, &Ny, &Nz, &current_time, &final_time,
		 VerticalDiffusionCoefficient_i.GetData(),
		 DepositionVelocity_i.GetData(), SurfaceEmission_i.GetData(),
		 VerticalDiffusionCoefficient_f.GetData(),
		 DepositionVelocity_f.GetData(), SurfaceEmission_f.GetData(),
		 ModifiedCellCenterDistance_z.data(), CellWidth_z.data(),
		 AirDensity.GetData(), Concentration.GetData(),
		 Concentration_ccl.GetData());
      }

    if (meridional_diffusion_)
      {
	Concentration.GetArray() = conc_sav_y;
	_diffYcl(&Nx, &Ny, &Nz, &current_time, &final_time,
		 MeridionalDiffusionCoefficient.GetData(),
		 ModifiedCellCenterDistance_y.data(), CellWidth_y.data(),
		 AirDensity.GetData(), Concentration.GetData(),
		 Concentration_ccl.GetData());
      }

    if (zonal_diffusion_)
      {
	Concentration.GetArray() = conc_sav_x;
	_diffXcl(&Nx, &Ny, &Nz, &current_time, &final_time,
		 ZonalDiffusionCoefficient.GetData(),
		 ModifiedCellCenterDistance_x.data(), CellWidth_x.data(),
		 AirDensity.GetData(), Concentration.GetData(),
		 Concentration_ccl.GetData());
      }

    Concentration.GetArray() = res;
  }

    
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_DIFFUSIONROS2_CXX
#endif
