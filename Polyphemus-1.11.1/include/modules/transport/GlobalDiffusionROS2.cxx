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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_GLOBALDIFFUSIONROS2_CXX


#include "GlobalDiffusionROS2.hxx"


namespace Polyphemus
{


  //! Default constructor.
  /*!
    \note Zonal, meridional and vertical diffusion are taken into account in
    the default configuration.
  */
  template<class T>
  GlobalDiffusionROS2<T>::GlobalDiffusionROS2():
    zonal_diffusion_(true), meridional_diffusion_(true),
    vertical_diffusion_(true)
  {
  }


  //! Initialization of the scheme.
  /*! \note Empty method.
   */
  template<class T>
  template<class ClassModel>
  void GlobalDiffusionROS2<T>::Init(ClassModel& Model)
  {
    // Checks that the domain is global.
    if (abs(Model.GetDelta_x() * T(Model.GetNx()) - 360.)
        > Model.GetDelta_x() * 1.e-6)
      throw string("Error in 'GlobalDiffusionROS2': ")
        + "the domain is not global (along longitude).";
    if (abs(Model.GetDelta_y() * T(Model.GetNy()) - 180.)
        > Model.GetDelta_y() * 1.e-6)
      throw string("Error in 'GlobalDiffusionROS2': ")
        + "the domain is not global (along latitude).";
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
  void GlobalDiffusionROS2<T>::Forward(ClassModel& Model)
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

    for (int s = 0; s < Model.GetNs(); s++)
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

        /*** Numerical integration ***/

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
                Model.D3("AirDensity_i"), Concentration);
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
  void GlobalDiffusionROS2<T>::Forward_aer(ClassModel& Model)
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

    vector<int> index;
    int dep_b;
    int emis_s, emis_b;
    for (int s = 0; s < Model.GetNs_aer(); s++)
      for (int b = 0; b < Model.GetNbin_aer(); b++)
        {

          /*** Deposition velocities ***/

          if (Model.HasDepositionVelocity_aer(b))
            {
              dep_b = Model.DepositionVelocityIndex_aer(b);
              DepositionVelocity_i.GetArray()
                = Model.A3("DepositionVelocity_aer_i")(dep_b, Range::all(),
                                                       Range::all()).copy();
              DepositionVelocity_f.GetArray()
                = Model.A3("DepositionVelocity_aer_f")(dep_b, Range::all(),
                                                       Range::all()).copy();
            }
          else
            {
              DepositionVelocity_i.SetZero();
              DepositionVelocity_f.SetZero();
            }

          /*** Surface emissions ***/

          if (Model.HasSurfaceEmission_aer(s, b))
            {
              index = Model.SurfaceEmissionIndex_aer(s, b);
              emis_s = index[0];
              emis_b = index[1];
              SurfaceEmission_i.GetArray()
                = Model.A4("SurfaceEmission_aer_i")(emis_s, emis_b,
                                                    Range::all(),
                                                    Range::all()).copy();
              SurfaceEmission_f.GetArray()
                = Model.A4("SurfaceEmission_aer_f")(emis_s, emis_b,
                                                    Range::all(),
                                                    Range::all()).copy();
            }
          else
            {
              SurfaceEmission_i.SetZero();
              SurfaceEmission_f.SetZero();
            }

          /*** Concentrations ***/

          Data<T, 3> Concentration(&Model.GetConcentration_aer()(s, b,
                                                                 0, 0, 0),
                                   shape(Model.GetNz(), Model.GetNy(),
                                         Model.GetNx()));

          /*** Numerical integration ***/

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
                  Model.D3("AirDensity_i"), Concentration);
        }

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
  void GlobalDiffusionROS2<T>
  ::Forward(T current_time, T Delta_t,
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
            Data<T, 3>& AirDensity, Data<T, 3>& Concentration)
  {
    T final_time = current_time + Delta_t;

    if (zonal_diffusion_)
      _global_diffX(&Nx, &Ny, &Nz, &current_time, &final_time,
                    ZonalDiffusionCoefficient.GetData(),
                    ModifiedCellCenterDistance_x.data(), CellWidth_x.data(),
                    AirDensity.GetData(), Concentration.GetData());

    if (meridional_diffusion_)
      _global_diffY(&Nx, &Ny, &Nz, &current_time, &final_time,
                    MeridionalDiffusionCoefficient.GetData(),
                    ModifiedCellCenterDistance_y.data(), CellWidth_y.data(),
                    AirDensity.GetData(), Concentration.GetData());

    if (vertical_diffusion_)
      _global_diffZ(&Nx, &Ny, &Nz, &current_time, &final_time,
                    VerticalDiffusionCoefficient_i.GetData(),
                    DepositionVelocity_i.GetData(),
                    SurfaceEmission_i.GetData(),
                    VerticalDiffusionCoefficient_f.GetData(),
                    DepositionVelocity_f.GetData(),
                    SurfaceEmission_f.GetData(),
                    ModifiedCellCenterDistance_z.data(), CellWidth_z.data(),
                    AirDensity.GetData(), Concentration.GetData());
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_GLOBALDIFFUSIONROS2_CXX
#endif
