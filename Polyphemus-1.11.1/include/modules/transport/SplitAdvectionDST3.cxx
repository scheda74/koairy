// Copyright (C) 2007-2008, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_SPLITADVECTIONDST3_CXX


#include "SplitAdvectionDST3.hxx"

#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"

namespace Polyphemus
{


  //! Initialization of the scheme.
  /*! \note Parallelization parameters are initialized here.
   */
  template<class T>
  template<class ClassModel>
  void SplitAdvectionDST3<T>::Init(ClassModel& Model)
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
    <li> HasBoundaryCondition(int)
    <li> A3("BoundaryCondition_z")
    <li> A4("BoundaryCondition_y")
    <li> A4("BoundaryCondition_x")
    <li> D3("ZonalWind")
    <li> D3("MeridionalWind")
    <li> D3("VerticalWind")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void SplitAdvectionDST3<T>::Forward(ClassModel& Model)
  {
    int Ns = Model.GetNs();
    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    // Computes the number of subcycles in each direction in order to satisfy
    // the CFL.
    int Nt_x, Nt_y, Nt_z;
    T limit_Delta_t;
    limit_Delta_t = min(Model.CellWidth_x)
      / (max(-Model.D3("ZonalWind").GetMin(),
             Model.D3("ZonalWind").GetMax()) + 1.e-10);
    Nt_x = int(ceil(Model.GetDelta_t() / limit_Delta_t));
    limit_Delta_t = min(Model.CellWidth_y)
      / (max(-Model.D3("MeridionalWind").GetMin(),
             Model.D3("MeridionalWind").GetMax()) + 1.e-10);
    Nt_y = int(ceil(Model.GetDelta_t() / limit_Delta_t));
    limit_Delta_t = min(Model.CellWidth_z)
      / (max(-Model.D3("VerticalWind").GetMin(),
             Model.D3("VerticalWind").GetMax()) + 1.e-10);
    Nt_z = int(ceil(Model.GetDelta_t() / limit_Delta_t));

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
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  shared(Ns, Nz, Ny, Nx, Nt_z, Nt_y, Nt_x)		\
  firstprivate(first_index_along_s, last_index_along_s)
#endif
    for (int s = first_index_along_s; s < last_index_along_s; s++)
      {
        int with_bc;

        /*** Boundary conditions ***/
        Array<T, 3> BoundaryCondition_x_i(shape(Nz, Ny, 2));
        Array<T, 3> BoundaryCondition_y_i(shape(Nz, 2, Nx));
        Array<T, 2> BoundaryCondition_z_i(shape(Ny, Nx));

        if (Model.HasBoundaryCondition(s))
          {
            int bc_s = Model.BoundaryConditionIndex(s);
            BoundaryCondition_x_i
              = Model.A4("BoundaryCondition_x")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_y_i
              = Model.A4("BoundaryCondition_y")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_z_i
              = Model.A3("BoundaryCondition_z")(bc_s, Range::all(),
                                                Range::all()).copy();
            with_bc = 1;
          }
        else
          {
            BoundaryCondition_x_i = 0;
            BoundaryCondition_y_i = 0;
            BoundaryCondition_z_i = 0;

            with_bc = 0;
          }

        /*** Concentrations ***/
        Array<T, 3> Concentration_s(&Model.GetConcentration()(s, 0, 0, 0),
                                    shape(Nz, Ny, Nx));

        /*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
        if (ParallelizedOnSpecies())
#endif
          Forward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                  Model.D3("VerticalWind"), with_bc, Model.CellWidth_x,
                  Model.CellWidth_y, Model.CellWidth_z,
                  Nx, Ny, Nz, BoundaryCondition_x_i,
                  BoundaryCondition_y_i, BoundaryCondition_z_i,
                  Model.GetDelta_t(), Nt_x, Nt_y, Nt_z, Concentration_s);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
        else
          ForwardParallel(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                          Model.D3("VerticalWind"), with_bc,
                          Model.CellWidth_x, Model.CellWidth_y,
                          Model.CellWidth_z,
                          Nx, Ny, Nz, BoundaryCondition_x_i,
                          BoundaryCondition_y_i, BoundaryCondition_z_i,
                          Model.GetDelta_t(), Nt_x, Nt_y, Nt_z,
                          Concentration_s);
#endif
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnSpecies())
      GatherSlice_s_MPI(Model.GetConcentration());
#endif
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
    <li> HasBoundaryCondition_aer(int, int)
    <li> A3("BoundaryCondition_aer_z")
    <li> A4("BoundaryCondition_aer_y")
    <li> A4("BoundaryCondition_aer_x")
    <li> D3("ZonalWind")
    <li> D3("MeridionalWind")
    <li> D3("VerticalWind")
    <li> CellWidth_x
    <li> CellWidth_y
    <li> CellWidth_z
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void SplitAdvectionDST3<T>::Forward_aer(ClassModel& Model)
  {
    // Computes the number of subcycles in each direction in order to satisfy
    // the CFL.
    int Nt_x, Nt_y, Nt_z;
    T limit_Delta_t;
    limit_Delta_t = min(Model.CellWidth_x)
      / (max(-Model.D3("ZonalWind").GetMin(),
             Model.D3("ZonalWind").GetMax()) + 1.e-10);
    Nt_x = int(ceil(Model.GetDelta_t() / limit_Delta_t));
    limit_Delta_t = min(Model.CellWidth_y)
      / (max(-Model.D3("MeridionalWind").GetMin(),
             Model.D3("MeridionalWind").GetMax()) + 1.e-10);
    Nt_y = int(ceil(Model.GetDelta_t() / limit_Delta_t));
    limit_Delta_t = min(Model.CellWidth_z)
      / (max(-Model.D3("VerticalWind").GetMin(),
             Model.D3("VerticalWind").GetMax()) + 1.e-10);
    Nt_z = int(ceil(Model.GetDelta_t() / limit_Delta_t));


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
#pragma omp parallel for num_threads(Nthreads_openmp)		\
  shared(Nt_z, Nt_y, Nt_x)					\
  firstprivate(first_index_along_s_aer, last_index_along_s_aer)
#endif
    for (int b = 0; b < Model.GetNbin_aer(); b++)
      {
        int bc_b, with_number_bc;

        for (int s = first_index_along_s_aer; s < last_index_along_s_aer; s++)
          {
            int with_bc;
            int bc_s;
            vector<int> index;

            /*** Boundary conditions ***/
            Array<T, 3> BoundaryCondition_x_i(shape(Model.GetNz(),
                                                    Model.GetNy(), 2));
            Array<T, 3> BoundaryCondition_y_i(shape(Model.GetNz(), 2,
                                                    Model.GetNx()));
            Array<T, 2> BoundaryCondition_z_i(shape(Model.GetNy(),
                                                    Model.GetNx()));

            if (Model.HasBoundaryCondition_aer(s, b))
              {
                index = Model.BoundaryConditionIndex_aer(s, b);
                bc_s = index[0];
                bc_b = index[1];
                BoundaryCondition_x_i
                  = Model.A5("BoundaryCondition_x_aer")(bc_s, bc_b,
                                                        Range::all(),
                                                        Range::all(),
                                                        Range::all()).copy();
                BoundaryCondition_y_i
                  = Model.A5("BoundaryCondition_y_aer")(bc_s, bc_b,
                                                        Range::all(),
                                                        Range::all(),
                                                        Range::all()).copy();
                BoundaryCondition_z_i
                  = Model.A4("BoundaryCondition_z_aer")(bc_s, bc_b,
                                                        Range::all(),
                                                        Range::all()).copy();
                with_bc = 1;
              }
            else
              {
                BoundaryCondition_x_i = 0;
                BoundaryCondition_y_i = 0;
                BoundaryCondition_z_i = 0;
                
                with_bc = 0;
              }

            /*** Concentrations ***/

            Array<T, 3>
              Concentration(&Model.GetConcentration_aer()(s, b, 0, 0, 0),
                            shape(Model.GetNz(), Model.GetNy(), Model.GetNx()));

            /*** Numerical integration ***/
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
            if (ParallelizedOnAerSpecies())
#endif
              Forward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                      Model.D3("VerticalWind"), with_bc, Model.CellWidth_x,
                      Model.CellWidth_y, Model.CellWidth_z, Model.GetNx(),
                      Model.GetNy(), Model.GetNz(), BoundaryCondition_x_i,
                      BoundaryCondition_y_i, BoundaryCondition_z_i,
                      Model.GetDelta_t(), Nt_x, Nt_y, Nt_z, Concentration);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
            else
              ForwardParallel(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                              Model.D3("VerticalWind"), with_bc,
                              Model.CellWidth_x, Model.CellWidth_y,
                              Model.CellWidth_z, Model.GetNx(),
                              Model.GetNy(), Model.GetNz(),
                              BoundaryCondition_x_i, BoundaryCondition_y_i,
                              BoundaryCondition_z_i, Model.GetDelta_t(),
                              Nt_x, Nt_y, Nt_z, Concentration);
#endif
          }

        if (Model.HasNumberConcentration_aer())
          {

	    /*** Number Boundary conditions ***/
            Array<T, 3> NumberBoundaryCondition_x_i(shape(Model.GetNz(),
                                                          Model.GetNy(), 2));
            Array<T, 3> NumberBoundaryCondition_y_i(shape(Model.GetNz(), 2,
							  Model.GetNx()));
            Array<T, 2> NumberBoundaryCondition_z_i(shape(Model.GetNy(),
                                                          Model.GetNx()));
		
            if (Model.HasNumberBoundaryCondition_aer(b))
              {
                bc_b = Model.NumberBoundaryConditionIndex_aer(b);
                NumberBoundaryCondition_x_i
                  = Model.A4("NumberBoundaryCondition_x_aer")(bc_b,
                                                              Range::all(),
                                                              Range::all(),
                                                              Range::all()).copy();
                NumberBoundaryCondition_y_i
                  = Model.A4("NumberBoundaryCondition_y_aer")(bc_b,
                                                              Range::all(),
                                                              Range::all(),
                                                              Range::all()).copy();
                NumberBoundaryCondition_z_i
                  = Model.A3("NumberBoundaryCondition_z_aer")(bc_b,
                                                              Range::all(),
                                                              Range::all()).copy();
                
                with_number_bc = 1;
              }
            else
              {
                NumberBoundaryCondition_x_i = 0;
                NumberBoundaryCondition_y_i = 0;
                NumberBoundaryCondition_z_i = 0;
		
                with_number_bc = 0;
              }

	    /*** Number Concentrations ***/

            Array<T, 3>
              NumberConcentration(&Model.GetNumberConcentration_aer()(b, 0, 0, 0),
                                  shape(Model.GetNz(), Model.GetNy(), Model.GetNx()));

	    /*** Numerical integration ***/

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	    if (ParallelizedOnAerSpecies())
#endif
              Forward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                      Model.D3("VerticalWind"), with_number_bc, Model.CellWidth_x,
                      Model.CellWidth_y, Model.CellWidth_z, Model.GetNx(),
                      Model.GetNy(), Model.GetNz(), NumberBoundaryCondition_x_i,
                      NumberBoundaryCondition_y_i, NumberBoundaryCondition_z_i,
                      Model.GetDelta_t(), Nt_x, Nt_y, Nt_z, NumberConcentration);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
	    else
	      ForwardParallel(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
			      Model.D3("VerticalWind"), with_number_bc,
			      Model.CellWidth_x, Model.CellWidth_y,
			      Model.CellWidth_z, Model.GetNx(),
			      Model.GetNy(), Model.GetNz(),
			      NumberBoundaryCondition_x_i, NumberBoundaryCondition_y_i,
			      NumberBoundaryCondition_z_i, Model.GetDelta_t(),
			      Nt_x, Nt_y, Nt_z, NumberConcentration);
            
#endif
          }
      }


#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (ParallelizedOnAerSpecies())
      GatherSlice_s_aer_MPI(Model.GetConcentration_aer());
#endif
  }

  //! Performs an integration over one time step.
  /*!
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param VerticalWind vertical wind.
    \param with_bc 1 if boundary conditions are included, 0 otherwise.
    \param CellWidth_x cells widths along x in meters.
    \param CellWidth_y cells widths along y in meters.
    \param CellWidth_z cells widths along z in meters.
    \param Nx number of cells along x.
    \param Ny number of cells along y.
    \param Nz number of cells along z.
    \param BoundaryCondition_x boundary conditions along x.
    \param BoundaryCondition_y boundary conditions along y.
    \param BoundaryCondition_z boundary conditions along z.
    \param Delta_t time step.
    \param Nt_x number of subcycles for the time integration along x.
    \param Nt_y number of subcycles for the time integration along y.
    \param Nt_z number of subcycles for the time integration along z.
    \param Concentration concentrations.
  */
  template<class T>
  void SplitAdvectionDST3<T>::Forward(Data<T, 3>& ZonalWind,
                                      Data<T, 3>& MeridionalWind,
                                      Data<T, 3>& VerticalWind,
                                      int with_bc, Array<T, 1>& CellWidth_x,
                                      Array<T, 1>& CellWidth_y,
                                      Array<T, 1>& CellWidth_z,
                                      int Nx, int Ny, int Nz,
                                      Array<T, 3>& BoundaryCondition_x,
                                      Array<T, 3>& BoundaryCondition_y,
                                      Array<T, 2>& BoundaryCondition_z,
                                      T Delta_t, int Nt_x, int Nt_y, int Nt_z,
                                      Array<T, 3>& Concentration)
  {
    int k, j, i, h;

    // Subcycle timestep.
    T Delta_t_sub;

    bool bc = with_bc == 1;

    Array<T, 1> Wind(Nx + 1);
    Array<T, 1> BoundaryCondition(2);
    Array<T, 1> Concentration1D(Nx);

    // Along x.
    Delta_t_sub = Delta_t / T(Nt_x);
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        {
          Wind = ZonalWind()(k, j, Range::all());
          if (bc)
            {
              BoundaryCondition(0) = BoundaryCondition_x(k, j, 0);
              BoundaryCondition(1) = BoundaryCondition_x(k, j, 1);
            }
          Concentration1D = Concentration(k, j, Range::all());
          for (h = 0; h < Nt_x; h++)
            _split_advection(&Nx, CellWidth_x.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());
          Concentration(k, j, Range::all()) = Concentration1D;
        }

    // Along y.
    Delta_t_sub = Delta_t / T(Nt_y);
    Wind.resize(Ny + 1);
    Concentration1D.resize(Ny);
    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        {
          Wind = MeridionalWind()(k, Range::all(), i);
          if (bc)
            {
              BoundaryCondition(0) = BoundaryCondition_y(k, 0, i);
              BoundaryCondition(1) = BoundaryCondition_y(k, 1, i);
            }
          Concentration1D = Concentration(k, Range::all(), i);
          for (h = 0; h < Nt_y; h++)
            _split_advection(&Ny, CellWidth_y.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());
          Concentration(k, Range::all(), i) = Concentration1D;
        }

    // Along z.
    Delta_t_sub = Delta_t / T(Nt_z);
    Wind.resize(Nz + 1);
    Concentration1D.resize(Nz);
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        {
          Wind = VerticalWind()(Range::all(), j, i);
          if (bc)
            {
              BoundaryCondition(0) = 0;
              BoundaryCondition(1) = BoundaryCondition_z(j, i);
            }
          Concentration1D = Concentration(Range::all(), j, i);
          for (h = 0; h < Nt_z; h++)
            _split_advection(&Nz, CellWidth_z.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());
          Concentration(Range::all(), j, i) = Concentration1D;
        }
  }

  //! Performs an integration over one time step.
  /*!
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param VerticalWind vertical wind.
    \param with_bc 1 if boundary conditions are included, 0 otherwise.
    \param CellWidth_x cells widths along x in meters.
    \param CellWidth_y cells widths along y in meters.
    \param CellWidth_z cells widths along z in meters.
    \param Nx number of cells along x.
    \param Ny number of cells along y.
    \param Nz number of cells along z.
    \param BoundaryCondition_x boundary conditions along x.
    \param BoundaryCondition_y boundary conditions along y.
    \param BoundaryCondition_z boundary conditions along z.
    \param Delta_t time step.
    \param Nt_x number of subcycles for the time integration along x.
    \param Nt_y number of subcycles for the time integration along y.
    \param Nt_z number of subcycles for the time integration along z.
    \param Concentration concentrations.
  */
  template<class T>
  void SplitAdvectionDST3<T>::ForwardParallel(Data<T, 3>& ZonalWind,
                                              Data<T, 3>& MeridionalWind,
                                              Data<T, 3>& VerticalWind,
                                              int with_bc,
                                              Array<T, 1>& CellWidth_x,
                                              Array<T, 1>& CellWidth_y,
                                              Array<T, 1>& CellWidth_z,
                                              int Nx, int Ny, int Nz,
                                              Array<T, 3>&
                                              BoundaryCondition_x,
                                              Array<T, 3>&
                                              BoundaryCondition_y,
                                              Array<T, 2>&
                                              BoundaryCondition_z,
                                              T Delta_t,
                                              int Nt_x, int Nt_y, int Nt_z,
                                              Array<T, 3>& Concentration)
  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    int k, j, i, h;
    int first_index_along_x, first_index_along_y;
    int last_index_along_x, last_index_along_y;

    // Subcycle timestep.
    T Delta_t_sub;

    bool bc = with_bc == 1;

    Array<T, 1> Wind(Nx + 1);
    Array<T, 1> BoundaryCondition(2);
    Array<T, 1> Concentration1D(Nx);

    // Advection along x.
    // Buffers used to send and receive concentrations from process #0 to the
    // other processes.
    Array<T, 3> OrdConc;
    ScatterSlice_y_MPI(Concentration, OrdConc);
    GetEdgePartition_y(first_index_along_y, last_index_along_y);
    Delta_t_sub = Delta_t / T(Nt_x);
    for (j = first_index_along_y; j < last_index_along_y; j++)
      for (k = 0; k < Nz; k++)
        {
          Wind = ZonalWind()(k, j, Range::all());
          if (bc)
            {
              BoundaryCondition(0) = BoundaryCondition_x(k, j, 0);
              BoundaryCondition(1) = BoundaryCondition_x(k, j, 1);
            }
          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_y_MPI.
          for (i = 0; i < Nx; i++)
            Concentration1D(i) = OrdConc(j, k, i);

          for (h = 0; h < Nt_x; h++)
            _split_advection(&Nx, CellWidth_x.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());

          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_y_MPI.
          OrdConc(j, k, Range::all()) = Concentration1D;
        }
    GatherSlice_y_MPI(OrdConc, Concentration);

    // Advection along y.
    ScatterSlice_x_MPI(Concentration, OrdConc);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
    Delta_t_sub = Delta_t / T(Nt_y);
    Wind.resize(Ny + 1);
    Concentration1D.resize(Ny);
    for (i = first_index_along_x; i < last_index_along_x; i++)
      for (k = 0; k < Nz; k++)
        {
          Wind = MeridionalWind()(k, Range::all(), i);
          if (bc)
            {
              BoundaryCondition(0) = BoundaryCondition_y(k, 0, i);
              BoundaryCondition(1) = BoundaryCondition_y(k, 1, i);
            }

          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_x_MPI.
          for (j = 0; j < Ny; j++)
            Concentration1D(j) = OrdConc(i, j, k);

          for (h = 0; h < Nt_y; h++)
            _split_advection(&Ny, CellWidth_y.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());

          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_x_MPI.
          OrdConc(i, Range::all(), k) = Concentration1D;
        }
    GatherSlice_x_MPI(OrdConc, Concentration);

    // Advection along z.
    // Nx_para, first, last, count and offset were already computed for
    // advection along y.
    ScatterSlice_x_MPI(Concentration, OrdConc);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
    Delta_t_sub = Delta_t / T(Nt_z);
    Wind.resize(Nz + 1);
    Concentration1D.resize(Nz);
    for (i = first_index_along_x; i < last_index_along_x; i++)
      for (j = 0; j < Ny; j++)
        {
          Wind = VerticalWind()(Range::all(), j, i);
          if (bc)
            {
              BoundaryCondition(0) = 0.;
              BoundaryCondition(1) = BoundaryCondition_z(j, i);
            }

          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_x_MPI.
          for (k = 0; k < Nz; k++)
            Concentration1D(k) = OrdConc(i, j, k);

          for (h = 0; h < Nt_z; h++)
            _split_advection(&Nz, CellWidth_z.data(), Wind.data(),
                             &with_bc, BoundaryCondition.data(),
                             &Delta_t_sub, Concentration1D.data());

          // Rank of the indices of OrdConc must match the one defined in
          // ScatterSlice_x_MPI.
          OrdConc(i, j, Range::all()) = Concentration1D;
        }
    GatherSlice_x_MPI(OrdConc, Concentration);
#endif
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_SPLITADVECTIONDST3_CXX
#endif
