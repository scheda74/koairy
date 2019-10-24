// Copyright (C) 2009, ENPC - INRIA - EDF R&D
// Author(s): Pierre Tran
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


#ifndef POLYPHEMUS_FILE_MODELS_LAGRANGIANTRANSPORT_HXX


#include <vector>
#include <list>
#include <algorithm>
#include "newran.h"
#include "AtmoData.hxx"
#include "BaseModel.cxx"


namespace Polyphemus
{

  using namespace std;
  using namespace AtmoData;


  //////////////////////////
  // LAGRANGIAN TRANSPORT //
  //////////////////////////


  //! This class is a lagrangian solver for a transport equation.
  template<class T, class ClassParticle>
  class LagrangianTransport: public BaseModel<T>
  {

  protected:

    /*** Configuration ***/

    //! List of particles.
    list<ClassParticle> particles_list;

    /*** Domain ***/

    //! Domain end along x.
    T x_max;

    //! Domain end along y.
    T y_max;

    //! Cell widths along x in meters.
    Array<T, 1> CellWidth_x;
    //! Cell widths along y in meters.
    Array<T, 1> CellWidth_y;
    //! Cell widths along z in meters.
    Array<T, 1> CellWidth_z;

    //! Distances between cell centers along x in meters.
    Array<T, 1> CellCenterDistance_x;
    //! Distances between cell centers along y in meters.
    Array<T, 1> CellCenterDistance_y;
    //! Distances between cell centers along z in meters.
    Array<T, 1> CellCenterDistance_z;

    /*** Winds ***/

    //! Zonal wind at current date.
    Data<T, 3> ZonalWind_i;
    //! Meridional wind at current date.
    Data<T, 3> MeridionalWind_i;
    //! Vertical wind at current date.
    Data<T, 3> VerticalWind_i;

    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_i;
    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_f;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_i;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_f;
    //! Vertical wind buffer.
    Data<T, 3> FileVerticalWind_i;
    //! Vertical wind buffer.
    Data<T, 3> FileVerticalWind_f;

    /*** Other meteorological fields ***/

    //! Temperature at current date.
    Data<T, 3> Temperature_i;
    //! Temperature at next date.
    Data<T, 3> Temperature_f;
    //! Pressure at current date.
    Data<T, 3> Pressure_i;
    //! Pressure at next date.
    Data<T, 3> Pressure_f;
    //! Air density at current date.
    Data<T, 3> AirDensity_i;
    //! Air density at interfaces along z and at current date.
    Data<T, 3> AirDensity_interf_z_i;
    //! Air density at interfaces along y and at current date.
    Data<T, 3> AirDensity_interf_y_i;
    //! Air density at interfaces along x and at current date.
    Data<T, 3> AirDensity_interf_x_i;
    //! Air density at next date.
    Data<T, 3> AirDensity_f;
    //! Air density at interfaces along z and at next date.
    Data<T, 3> AirDensity_interf_z_f;
    //! Air density at interfaces along y and at next date.
    Data<T, 3> AirDensity_interf_y_f;
    //! Air density at interfaces along x and at next date.
    Data<T, 3> AirDensity_interf_x_f;

    //! Temperature buffer.
    Data<T, 3> FileTemperature_i;
    //! Temperature buffer.
    Data<T, 3> FileTemperature_f;
    //! Pressure buffer.
    Data<T, 3> FilePressure_i;
    //! Pressure buffer.
    Data<T, 3> FilePressure_f;

    /*** Diffusion ***/

    //! Horizontal diffusion coefficient.
    T horizontal_diffusion;
    //! Zonal diffusion coefficient at current date.
    Data<T, 3> ZonalDiffusionCoefficient_i;
    //! Meridional diffusion coefficient at current date.
    Data<T, 3> MeridionalDiffusionCoefficient_i;
    //! Vertical diffusion coefficient at current date.
    Data<T, 3> VerticalDiffusionCoefficient_i;
    //! Vertical diffusion coefficient at next date.
    Data<T, 3> VerticalDiffusionCoefficient_f;
    //! Vertical diffusion coefficient buffer.
    Data<T, 3> FileVerticalDiffusionCoefficient_i;
    //! Vertical diffusion coefficient buffer.
    Data<T, 3> FileVerticalDiffusionCoefficient_f;
    //! Horizontal diffusion coefficient used for the gaussian kernel.
    T gaussian_kernel_horizontal_diffusion;


    /*** Emission ***/

    //! Approximate time step between two particles release at a given source.
    T Delta_t_particle_emission;

    /*** Random terms ***/

    //! Seed of the random generator.
    NEWRAN::LGM_mixed *urng;
    //! Uniform random number generator.
    NEWRAN::Uniform RandomNumber;


  public:

    /*** Constructor and destructor ***/

    LagrangianTransport(string config_file);
    virtual ~LagrangianTransport();

    /*** Access methods ***/

    T GetHorizontalDiffusion()
    {
      return horizontal_diffusion;
    }
    T GetGaussianKernelHorizontalDiffusion()
    {
      return gaussian_kernel_horizontal_diffusion;
    }
    T GetX_max()
    {
      return x_max;
    }
    T GetY_max()
    {
      return y_max;
    }

    /*** Configuration ***/

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    /*** Initializations ***/

    void Allocate();
    void Init();
    void InitStep();

    virtual void SetDate(Date date);

    /*** Integration ***/

    void Transport();
    void PointEmission();
    void Forward();

    /*** Random number generation ***/

    T GetRandomNumber();

    /*** Other methods ***/

    void TransformMeridionalWind(Data<T, 3>& MeridionalWind);
    void TransformZonalWind(Data<T, 3>& ZonalWind);
    void TransformMeridionalDiffusion(Data<T, 3>& MeridionalDiffusion_);
    void TransformZonalDiffusion(Array<T, 1>& GridY_interf_,
                                 Data<T, 3>& ZonalDiffusion_);

    void ComputeVerticalWind(Array<T, 1>& CellWidth_x_,
                             Array<T, 1>& CellWidth_y_,
                             Array<T, 1>& CellWidth_z_,
                             Data<T, 3>& ZonalWind_,
                             Data<T, 3>& MeridionalWind_,
                             Data<T, 3>& VerticalWind_);
    void ComputeVerticalWind(Array<T, 1>& CellWidth_x_,
                             Array<T, 1>& CellWidth_y_,
                             Array<T, 1>& CellWidth_z_,
                             Data<T, 3>& AirDensity_interf_x_,
                             Data<T, 3>& ZonalWind_,
                             Data<T, 3>& AirDensity_interf_y_,
                             Data<T, 3>& MeridionalWind_,
                             Data<T, 3>& AirDensity_interf_z_,
                             Data<T, 3>& VerticalWind_);

    void ComputeAirDensity(Data<T, 3>& Temperature_, Data<T, 3>& Pressure_,
                           Data<T, 3>& AirDensity_);

    void InterpolateInterface_z(Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_z_);
    void InterpolateInterface_y(Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_y_);
    void InterpolateInterface_x(Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_x_);

    void InterpolateInterface_z(Array<T, 1>& CellCenterDistance_z_,
                                Array<T, 1>& CellWidth_z_, Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_z_);
    void InterpolateInterface_y(Array<T, 1>& CellCenterDistance_y_,
                                Array<T, 1>& CellWidth_y_, Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_y_);
    void InterpolateInterface_x(Array<T, 1>& CellCenterDistance_x_,
                                Array<T, 1>& CellWidth_x_, Data<T, 3>& Data_,
                                Data<T, 3>& Data_interf_x_);

    void ComputeCellWidth(T Delta_x_, Array<T, 1>& GridY_interf_,
                          Array<T, 1>& CellWidth_x_,
                          Array<T, 1>& CellWidth_y_);
    void ComputeCellCenterDistance(T Delta_x_, Array<T, 1>& GridY_interf_,
                                   Array<T, 1>& CellCenterDistance_x_,
                                   Array<T, 1>& CellCenterDistance_y_);

    virtual T GetConcentration(int species, T z, T y, T x);
    using BaseModel<T>::GetConcentration;
    virtual void ComputeConcentration();
    virtual void ComputeConcentration(const vector<int>& species_index,
                                      const vector<int>& levels);
    void GetCellIndices(T lon, T lat, T height,
                        int& index_z, int& index_y, int& index_x);

  protected:

    virtual void InitAllData();

  };


  /////////////////////////////
  // RELATED PREDICATE CLASS //
  /////////////////////////////


  //! Free-standing predicate class.
  template<class T>
  class IsOutsideOfDomain
  {

  private:

    //! Minimum abscissa in the curvilinear coordinate system.
    T x_min_;
    //! Maximum abscissa in the curvilinear coordinate system.
    T x_max_;
    //! Minimum ordinate in the curvilinear coordinate system.
    T y_min_;
    //! Maximum ordinate in the curvilinear coordinate system.
    T y_max_;

  public:

    /*** Constructor ***/

    template<typename ClassModel>
    IsOutsideOfDomain(ClassModel* Model);

    /*** Predicate ***/

    template<class ClassParticle>
    bool operator()(const ClassParticle& particle)
    {
      return particle.IsOutside(x_min_, x_max_, y_min_, y_max_);
    }
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_LAGRANGIANTRANSPORT_HXX
#endif
