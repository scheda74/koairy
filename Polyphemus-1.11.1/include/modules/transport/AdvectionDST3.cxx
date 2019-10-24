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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_ADVECTIONDST3_CXX


#include "AdvectionDST3.hxx"

namespace Polyphemus
{


  //! Initialization of the scheme.
  /*! \note Empty method.
   */
  template<class T>
  template<class ClassModel>
  void AdvectionDST3<T>::Init(ClassModel& Model)
  {
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
  void AdvectionDST3<T>::Forward(ClassModel& Model)
  {
    Data<T, 3> BoundaryCondition_x_i(shape(Model.GetNz(), Model.GetNy(), 2));
    Data<T, 3> BoundaryCondition_y_i(shape(Model.GetNz(), 2, Model.GetNx()));
    Data<T, 2> BoundaryCondition_z_i(shape(Model.GetNy(), Model.GetNx()));

    int with_bc;
    for (int s = 0; s < Model.GetNs(); s++)
      {

        /*** Boundary conditions ***/

        if (Model.HasBoundaryCondition(s))
          {
            int bc_s = Model.BoundaryConditionIndex(s);
            BoundaryCondition_x_i.GetArray()
              = Model.A4("BoundaryCondition_x")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_y_i.GetArray()
              = Model.A4("BoundaryCondition_y")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_z_i.GetArray()
              = Model.A3("BoundaryCondition_z")(bc_s, Range::all(),
                                                Range::all()).copy();
            with_bc = 1;
          }
        else
          {
            BoundaryCondition_x_i.SetZero();
            BoundaryCondition_y_i.SetZero();
            BoundaryCondition_z_i.SetZero();

            with_bc = 0;
          }

        /*** Concentrations ***/

        Data<T, 3> Concentration(&Model.GetConcentration()(s, 0, 0, 0),
                                 shape(Model.GetNz(), Model.GetNy(),
                                       Model.GetNx()));

        /*** Numerical integration ***/

        Forward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                Model.D3("VerticalWind"), with_bc, Model.CellWidth_x,
                Model.CellWidth_y, Model.CellWidth_z, Model.GetNx(),
                Model.GetNy(), Model.GetNz(), BoundaryCondition_x_i,
                BoundaryCondition_y_i, BoundaryCondition_z_i,
                Model.GetDelta_t(), Concentration);
      }
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
    <li> GetConcentratio_ccl()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void AdvectionDST3<T>::Backward(ClassModel& Model)
  {
    Data<T, 3> BoundaryCondition_x_i(shape(Model.GetNz(), Model.GetNy(), 2));
    Data<T, 3> BoundaryCondition_y_i(shape(Model.GetNz(), 2, Model.GetNx()));
    Data<T, 2> BoundaryCondition_z_i(shape(Model.GetNy(), Model.GetNx()));

    int with_bc;

    for (int s = 0; s < Model.GetNs(); s++)
      {

        /*** Boundary conditions ***/

        if (Model.HasBoundaryCondition(s))
          {
            int bc_s = Model.BoundaryConditionIndex(s);
            BoundaryCondition_x_i.GetArray()
              = Model.A4("BoundaryCondition_x")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_y_i.GetArray()
              = Model.A4("BoundaryCondition_y")(bc_s, Range::all(),
                                                Range::all(),
                                                Range::all()).copy();
            BoundaryCondition_z_i.GetArray()
              = Model.A3("BoundaryCondition_z")(bc_s, Range::all(),
                                                Range::all()).copy();

            with_bc = 1;
          }
        else
          {
            BoundaryCondition_x_i.SetZero();
            BoundaryCondition_y_i.SetZero();
            BoundaryCondition_z_i.SetZero();

            with_bc = 0;
          }


        /*** Temporary adjoint data for boundary conditions ***/

        Data<T, 3> BoundaryCondition_x_i_ccl(shape(Model.GetNz(),
                                                   Model.GetNy(), 2));
        Data<T, 3> BoundaryCondition_y_i_ccl(shape(Model.GetNz(),
                                                   2, Model.GetNx()));
        Data<T, 2> BoundaryCondition_z_i_ccl(shape(Model.GetNy(),
                                                   Model.GetNx()));

        BoundaryCondition_x_i_ccl.SetZero();
        BoundaryCondition_y_i_ccl.SetZero();
        BoundaryCondition_z_i_ccl.SetZero();

        /*** Concentrations ***/

        Data<T, 3> Concentration(&Model.GetConcentration()(s, 0, 0, 0),
                                 shape(Model.GetNz(), Model.GetNy(),
                                       Model.GetNx()));
        Data<T, 3> Concentration_ccl(&Model
                                     .GetConcentration_ccl()(s, 0, 0, 0),
                                     shape(Model.GetNz(), Model.GetNy(),
                                           Model.GetNx()));

        /*** Numerical integration ***/

        Backward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                 Model.D3("VerticalWind"), with_bc, Model.CellWidth_x,
                 Model.CellWidth_y, Model.CellWidth_z, Model.GetNx(),
                 Model.GetNy(), Model.GetNz(), BoundaryCondition_x_i,
                 BoundaryCondition_y_i, BoundaryCondition_z_i,
                 Model.GetDelta_t(), Concentration,
                 BoundaryCondition_x_i_ccl, BoundaryCondition_y_i_ccl,
                 BoundaryCondition_z_i_ccl, Concentration_ccl);
      }
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
  void AdvectionDST3<T>::Forward_aer(ClassModel& Model)
  {
    Data<T, 3> BoundaryCondition_x_i(shape(Model.GetNz(), Model.GetNy(), 2));
    Data<T, 3> BoundaryCondition_y_i(shape(Model.GetNz(), 2, Model.GetNx()));
    Data<T, 2> BoundaryCondition_z_i(shape(Model.GetNy(), Model.GetNx()));

    int with_bc;
    vector<int> index;
    int bc_s, bc_b;
    for (int s = 0; s < Model.GetNs_aer(); s++)
      for (int b = 0; b < Model.GetNbin_aer(); b++)
        {

          /*** Boundary conditions ***/

          if (Model.HasBoundaryCondition_aer(s, b))
            {
              index = Model.BoundaryConditionIndex_aer(s, b);
              bc_s = index[0];
              bc_b = index[1];
              BoundaryCondition_x_i.GetArray()
                = Model.A5("BoundaryCondition_x_aer")(bc_s, bc_b,
                                                      Range::all(),
                                                      Range::all(),
                                                      Range::all()).copy();
              BoundaryCondition_y_i.GetArray()
                = Model.A5("BoundaryCondition_y_aer")(bc_s, bc_b,
                                                      Range::all(),
                                                      Range::all(),
                                                      Range::all()).copy();
              BoundaryCondition_z_i.GetArray()
                = Model.A4("BoundaryCondition_z_aer")(bc_s, bc_b,
                                                      Range::all(),
                                                      Range::all()).copy();
              with_bc = 1;
            }
          else
            {
              BoundaryCondition_x_i.SetZero();
              BoundaryCondition_y_i.SetZero();
              BoundaryCondition_z_i.SetZero();

              with_bc = 0;
            }

          /*** Concentrations ***/

          Data<T, 3>
            Concentration(&Model.GetConcentration_aer()(s, b, 0, 0, 0),
                          shape(Model.GetNz(), Model.GetNy(), Model.GetNx()));

          /*** Numerical integration ***/

          Forward(Model.D3("ZonalWind"), Model.D3("MeridionalWind"),
                  Model.D3("VerticalWind"), with_bc, Model.CellWidth_x,
                  Model.CellWidth_y, Model.CellWidth_z, Model.GetNx(),
                  Model.GetNy(), Model.GetNz(), BoundaryCondition_x_i,
                  BoundaryCondition_y_i, BoundaryCondition_z_i,
                  Model.GetDelta_t(), Concentration);
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
    \param Concentration concentrations.
  */
  template<class T>
  void AdvectionDST3<T>::Forward(Data<T, 3>& ZonalWind,
                                 Data<T, 3>& MeridionalWind,
                                 Data<T, 3>& VerticalWind,
                                 int with_bc, Array<T, 1>& CellWidth_x,
                                 Array<T, 1>& CellWidth_y,
                                 Array<T, 1>& CellWidth_z,
                                 int Nx, int Ny, int Nz,
                                 Data<T, 3>& BoundaryCondition_x,
                                 Data<T, 3>& BoundaryCondition_y,
                                 Data<T, 2>& BoundaryCondition_z,
                                 T Delta_t, Data<T, 3>& Concentration)
  {
    _advection(&Nx, &Ny, &Nz, CellWidth_x.data(), CellWidth_y.data(),
               CellWidth_z.data(), ZonalWind.GetData(),
               MeridionalWind.GetData(), VerticalWind.GetData(),
               &with_bc, BoundaryCondition_x.GetData(),
               BoundaryCondition_y.GetData(),
               BoundaryCondition_z.GetData(),
               &Delta_t, Concentration.GetData());
  }


  //! Performs one step of backward integration of adjoint model.
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
    \param Concentration concentrations.
    \param Concentration_ccl adjoinot concentration data.
  */
  template<class T>
  void AdvectionDST3<T>::Backward(Data<T, 3>& ZonalWind,
                                  Data<T, 3>& MeridionalWind,
                                  Data<T, 3>& VerticalWind,
                                  int with_bc, Array<T, 1>& CellWidth_x,
                                  Array<T, 1>& CellWidth_y,
                                  Array<T, 1>& CellWidth_z,
                                  int Nx, int Ny, int Nz,
                                  Data<T, 3>& BoundaryCondition_x,
                                  Data<T, 3>& BoundaryCondition_y,
                                  Data<T, 2>& BoundaryCondition_z,
                                  T Delta_t,
                                  Data<T, 3>& Concentration,
                                  Data<T, 3>& BoundaryCondition_x_ccl,
                                  Data<T, 3>& BoundaryCondition_y_ccl,
                                  Data<T, 2>& BoundaryCondition_z_ccl,
                                  Data<T, 3>& Concentration_ccl)
  {
    _advectioncl(&Nx, &Ny, &Nz, CellWidth_x.data(), CellWidth_y.data(),
                 CellWidth_z.data(), ZonalWind.GetData(),
                 MeridionalWind.GetData(), VerticalWind.GetData(),
                 &with_bc, BoundaryCondition_x.GetData(),
                 BoundaryCondition_y.GetData(),
                 BoundaryCondition_z.GetData(),
                 &Delta_t, Concentration.GetData(),
                 BoundaryCondition_x_ccl.GetData(),
                 BoundaryCondition_y_ccl.GetData(),
                 BoundaryCondition_z_ccl.GetData(),
                 Concentration_ccl.GetData());
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_ADVECTIONDST3_CXX
#endif
