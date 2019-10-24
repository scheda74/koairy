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

// This file is part of the Eulerian model Polair3D.


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DTRANSPORT_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModel.cxx"


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
#define _compute_scavenging_coefficient compute_scavenging_coefficient__
#define _compute_scavenging_coefficient_pudykiewicz	\
  compute_scavenging_coefficient_pudykiewicz__
#else
#define _compute_scavenging_coefficient compute_scavenging_coefficient_
#define _compute_scavenging_coefficient_pudykiewicz	\
  compute_scavenging_coefficient_pudykiewicz_
#endif

  extern "C"
  {
    void _compute_scavenging_coefficient(int*, int*, int*, int*,
					 int*, int*, int*,
					 double*, double*, double*,
					 double*, double*, double*,
					 double*, double*);
    void _compute_scavenging_coefficient_pudykiewicz(int*, int*, int*, int*,
						     double*, double*,
						     double*, double*,
						     double*);
  }


  ///////////////////////
  // POLAIR3DTRANSPORT //
  ///////////////////////
  

  //! This class is a solver for an advection-diffusion equation.
  template<class T, class ClassAdvection, class ClassDiffusion>
  class Polair3DTransport: public BaseModel<T>
  {

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >::iterator
    input_files_iterator;

  public:
    
    /*** Configuration ***/
    
    //! Are coordinates Cartesian or latitude/longitude?
    bool option_cartesian;

    //! Is diffusion isotropic?
    bool option_isotropic_diffusion;

    //! Scavenging model.
    string scavenging_model;

    //! List of species with initial conditions.
    vector<string> species_list_ic;
    //! List of species with boundary conditions.
    vector<string> species_list_bc;
    //! List of species with deposition velocities.
    vector<string> species_list_dep;
    //! List of species with scavenging.
    vector<string> species_list_scav;
    //! List of species with surface emissions.
    vector<string> species_list_surf_emis;
    //! List of species with additional surface emissions.
    vector<string> species_list_add_surf_emis;
    //! List of species with volume emissions.
    vector<string> species_list_vol_emis;

    /*** Domain ***/

    //! Coordinates of boundary conditions along x.
    RegularGrid<T> GridX4D_interf_bc;
    //! Coordinates of boundary conditions along y.
    RegularGrid<T> GridY4D_interf_bc;
    
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

    /*** Species ***/

    //! Henry constants (mol / L / atm).
    map<string, T> henry_constant;
    //! Gas phase diffusivities (cm^2 / s).
    map<string, T> gas_phase_diffusivity;
    //! Constant scavenging coefficient (s^{-1}).
    map<string, T> scavenging_constant;
    //! Belot coefficient a.
    map<string, T> belot_constant_a;
    //! Belot coefficient b.
    map<string, T> belot_constant_b;

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
    //! Specific humidity at current date.
    Data<T, 3> SpecificHumidity_i;
    //! Specific humidity at next date.
    Data<T, 3> SpecificHumidity_f;
 
    //! Temperature buffer.
    Data<T, 3> FileTemperature_i;
    //! Temperature buffer.
    Data<T, 3> FileTemperature_f;
    //! Pressure buffer.
    Data<T, 3> FilePressure_i;
    //! Pressure buffer.
    Data<T, 3> FilePressure_f;
    //! Specific humidity buffer.
    Data<T, 3> FileSpecificHumidity_i;
    //! Specific humidity buffer.
    Data<T, 3> FileSpecificHumidity_f;

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

    /*** Initial conditions ***/

    // Number of species with initial conditions.
    int Ns_ic;

    /*** Boundary conditions ***/

    //! Number of species with boundary conditions.
    int Ns_bc;
    //! Grid for species with boundary conditions.
    RegularGrid<T> GridS_bc;

    //! Boundary conditions along z at current date.
    Data<T, 3> BoundaryCondition_z_i;
    //! Boundary conditions buffer.
    Data<T, 3> FileBoundaryCondition_z_i;
    //! Boundary conditions buffer.
    Data<T, 3> FileBoundaryCondition_z_f;

    //! Boundary conditions along y at current date.
    Data<T, 4> BoundaryCondition_y_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_y_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_y_f;

    //! Boundary conditions along x at current date.
    Data<T, 4> BoundaryCondition_x_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_x_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_x_f;

    /*** Loss terms ***/

    //! Number of species with deposition velocities.
    int Ns_dep;
    //! Grid for species with deposition velocities.
    RegularGrid<T> GridS_dep;
    //! Deposition velocities at current date.
    Data<T, 3> DepositionVelocity_i;
    //! Deposition velocities at next date.
    Data<T, 3> DepositionVelocity_f;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_i;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_f;
    //! Dry deposition fluxes at current date.
    Data<T, 3> DryDepositionFlux;

    //! Rain at current date.
    Data<T, 2> Rain_i;
    //! Rain at next date.
    Data<T, 2> Rain_f;
    //! Rain buffer.
    Data<T, 2> FileRain_i;
    //! Rain buffer.
    Data<T, 2> FileRain_f;

    //! Cloud height at current date.
    Data<T, 2> CloudHeight_i;
    //! Cloud height at next date.
    Data<T, 2> CloudHeight_f;
    //! Cloud height buffer.
    Data<T, 2> FileCloudHeight_i;
    //! Cloud height buffer.
    Data<T, 2> FileCloudHeight_f;

    //! Number of species with scavenging.
    int Ns_scav;
    //! Grid for species with scavenging.
    RegularGrid<T> GridS_scav;
    //! Scavenging coefficients at current date.
    Data<T, 4> ScavengingCoefficient_i;
    //! Scavenging coefficients at next date.
    Data<T, 4> ScavengingCoefficient_f;
    //! Wet deposition fluxes at current date.
    Data<T, 3> WetDepositionFlux;
    //! In cloud wet deposition fluxes at current date.
    Data<T, 3> InCloudWetDepositionFlux;
    //!  Scavenging rain threshold (mm / h).
    T scavenging_rain_threshold;

    /*** Source terms ***/

    //! Number of species with surface emissions.
    int Ns_surf_emis;
    //! Grid for species with surface emissions.
    RegularGrid<T> GridS_surf_emis;
    //! Surface emissions at current date.
    Data<T, 3> SurfaceEmission_i;
    //! Surface emissions at next date.
    Data<T, 3> SurfaceEmission_f;
    //! Surface emissions buffer.
    Data<T, 3> FileSurfaceEmission_i;
    //! Surface emissions buffer.
    Data<T, 3> FileSurfaceEmission_f;

    //! Number of species with additional surface emissions.
    int Ns_add_surf_emis;
    //! Grid for species with additional surface emissions.
    RegularGrid<T> GridS_add_surf_emis;
    //! Additional surface emissions at current date.
    Data<T, 3> AdditionalSurfaceEmission_i;
    //! Additional surface emissions at next date.
    Data<T, 3> AdditionalSurfaceEmission_f;
    //! Additional surface emissions buffer.
    Data<T, 3> FileAdditionalSurfaceEmission_i;
    //! Additional surface emissions buffer.
    Data<T, 3> FileAdditionalSurfaceEmission_f;

    //! Number species with volume emissions.
    int Ns_vol_emis;
    //! Grid for species with volume emissions.
    RegularGrid<T> GridS_vol_emis;
    //! Number of species with volume emissions.
    int Nz_vol_emis;
    //! Grid for altitudes of volume emissions.
    RegularGrid<T> GridZ_vol_emis;
    //! Volume emissions at current date.
    Data<T, 4> VolumeEmission_i;
    //! Volume emissions at next date.
    Data<T, 4> VolumeEmission_f;
    //! Volume emissions buffer.
    Data<T, 4> FileVolumeEmission_i;
    //! Volume emissions buffer.
    Data<T, 4> FileVolumeEmission_f;

    /*** Numerical schemes ***/

    //! Advection numerical scheme.
    ClassAdvection Advection_;
    //! Diffusion numerical scheme.
    ClassDiffusion Diffusion_;

  public:
    
    /*** Constructor and destructor ***/
    
    Polair3DTransport(string config_file);
    virtual ~Polair3DTransport();
    
    /*** Configuration ***/

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    bool HasInitialCondition(int s) const;
    bool HasInitialCondition(string name) const;
    int InitialConditionIndex(int s) const;
    int InitialConditionIndex(string name) const;
    string InitialConditionName(int s) const;
    int InitialConditionGlobalIndex(int s) const;
    
    bool HasBoundaryCondition(int s) const;
    bool HasBoundaryCondition(string name) const;
    int BoundaryConditionIndex(int s) const;
    int BoundaryConditionIndex(string name) const;
    string BoundaryConditionName(int s) const;
    int BoundaryConditionGlobalIndex(int s) const;
    
    bool HasDepositionVelocity(int s) const;
    bool HasDepositionVelocity(string name) const;
    int DepositionVelocityIndex(int s) const;
    int DepositionVelocityIndex(string name) const;
    string DepositionVelocityName(int s) const;
    int DepositionVelocityGlobalIndex(int s) const;
    
    bool HasScavenging(int s) const;
    bool HasScavenging(string name) const;
    int ScavengingIndex(int s) const;
    int ScavengingIndex(string name) const;
    string ScavengingName(int s) const;
    int ScavengingGlobalIndex(int s) const;
    
    bool HasSurfaceEmission(int s) const;
    bool HasSurfaceEmission(string name) const;
    int SurfaceEmissionIndex(int s) const;
    int SurfaceEmissionIndex(string name) const;
    string SurfaceEmissionName(int s) const;
    int SurfaceEmissionGlobalIndex(int s) const;

    string AdditionalSurfaceEmissionName(int s) const;
    
    bool HasVolumeEmission(int s) const;
    bool HasVolumeEmission(string name) const;
    int VolumeEmissionIndex(int s) const;
    int VolumeEmissionIndex(string name) const;
    string VolumeEmissionName(int s) const;
    int VolumeEmissionGlobalIndex(int s) const;
    
    /*** Initializations ***/

    virtual void Allocate();
    void Init();
    void InitStep();
    void InitScavengingCoefficient(Data<T, 4>& ScavengingCoefficient_);
    void InitScavengingCoefficient(Data<T, 2>& Rain_,
				   Data<T, 4>& ScavengingCoefficient_);
    void InitScavengingCoefficient(Data<T, 3>& Temperature_,
				   Data<T, 3>& Pressure_,
				   Data<T, 3>& SpecificHumidity_,
				   Data<T, 4>& ScavengingCoefficient_);
    void InitScavengingCoefficient(Data<T, 3>& Temperature_,
				   Data<T, 3>& Pressure_,
				   Data<T, 2>& CloudHeight_,
				   Data<T, 2>& Rain_,
				   Data<T, 4>& ScavengingCoefficient_);

    virtual void SetDate(Date date);

    /*** Integration ***/

    void Advection();
    void Advection_b();
    void Diffusion();
    void DiffusionXY();
    void DiffusionZ();
    void Diffusion_b();
    void PointEmission();

    void Forward();
    void SetBackward(bool flag);
    void Backward();

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

    T ComputeCellVolume(T Delta_x_, T Delta_y_, T Delta_z_, T lat);
    void ComputeCellWidth(T Delta_x_, Array<T, 1>& GridY_interf_,
			  Array<T, 1>& CellWidth_x_,
			  Array<T, 1>& CellWidth_y_);
    void ComputeCellCenterDistance(T Delta_x_, Array<T, 1>& GridY_interf_,
				   Array<T, 1>& CellCenterDistance_x_,
				   Array<T, 1>& CellCenterDistance_y_);

    virtual T GetConcentration(int species, T z, T y, T x);
    using BaseModel<T>::GetConcentration;
    void GetCellIndices(T lon, T lat, T height,
			int& index_z, int& index_y, int& index_x);
  protected:

    virtual void InitAllData();
    
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DTRANSPORT_HXX
#endif
