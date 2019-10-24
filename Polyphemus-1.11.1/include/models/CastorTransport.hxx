// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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

// This file is part of the Eulerian model Castor.

// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).


#ifndef POLYPHEMUS_FILE_MODELS_CASTORTRANSPORT_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModel.cxx"


namespace Polyphemus
{

  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // CASTORTRANSPORT //
  /////////////////////


  //! This class is a solver for an advection-diffusion equation.
  /*!  CastorTransport is a clone of Chimere with respect to transport.
   */
  template<class T, class ClassTransport>
  class CastorTransport: public BaseModel<T>
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

    //! List of species with initial conditions.
    vector<string> species_list_ic;
    //! List of species with boundary conditions.
    vector<string> species_list_bc;
    //! List of species with deposition velocities.
    vector<string> species_list_dep;
    //! List of species with volume emissions.
    vector<string> species_list_vol_emis;

    string chem_dir;

    /*** Domain ***/

    //! Coordinates of boundary conditions along x.
    RegularGrid<T> GridX4D_interf_bc;
    //! Coordinates of boundary conditions along y.
    RegularGrid<T> GridY4D_interf_bc;

    //! Cell widths along x in meters.
    Array<T, 1> CellWidth_x;
    //! Cell widths along y in meters.
    Array<T, 1> CellWidth_y;

    /*** Coordinates ***/

    //! Altitudes of layer interfaces (meters).
    Data<T, 3> Altitude;
    //! Buffer for altitudes of layer interfaces.
    Data<T, 3> FileAltitude_i;
    //! Buffer for altitudes of layer interfaces.
    Data<T, 3> FileAltitude_f;

    //! Layer thicknesses (meters).
    Data<T, 3> Thickness;

    /*** Winds ***/

    //! Zonal wind at current date.
    Data<T, 3> ZonalWind;
    //! Meridional wind at current date.
    Data<T, 3> MeridionalWind;

    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_i;
    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_f;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_i;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_f;

    /*** Other meteorological fields ***/

    //! Temperature at current date.
    Data<T, 3> Temperature;
    //! Pressure at current date.
    Data<T, 3> Pressure;
    //! Air density at current date.
    Data<T, 3> AirDensity;
    //! Air density at interfaces along z and at current date.
    Data<T, 3> AirDensity_interf_z;
    //! Air density at interfaces along y and at current date.
    Data<T, 3> AirDensity_interf_y;
    //! Air density at interfaces along x and at current date.
    Data<T, 3> AirDensity_interf_x;

    //! Temperature buffer.
    Data<T, 3> FileTemperature_i;
    //! Temperature buffer.
    Data<T, 3> FileTemperature_f;
    //! Pressure buffer.
    Data<T, 3> FilePressure_i;
    //! Pressure buffer.
    Data<T, 3> FilePressure_f;
    //! Air density buffer.
    Data<T, 3> FileAirDensity_i;
    //! Air density buffer.
    Data<T, 3> FileAirDensity_f;

    /*** Diffusion ***/

    //! Vertical diffusion coefficient at current date.
    Data<T, 3> VerticalDiffusionCoefficient;
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
    Data<T, 3> BoundaryCondition_z;
    //! Boundary conditions buffer.
    Data<T, 3> FileBoundaryCondition_z_i;
    //! Boundary conditions buffer.
    Data<T, 3> FileBoundaryCondition_z_f;

    //! Boundary conditions along y at current date.
    Data<T, 4> BoundaryCondition_y;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_y_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_y_f;

    //! Boundary conditions along x at current date.
    Data<T, 4> BoundaryCondition_x;
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
    Data<T, 3> DepositionVelocity;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_i;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_f;

    /*** Source terms ***/

    //! Number species with volume emissions.
    int Ns_vol_emis;
    //! Grid for species with volume emissions.
    RegularGrid<T> GridS_vol_emis;
    //! Number of species with volume emissions.
    int Nz_vol_emis;
    //! Grid for altitudes of volume emissions.
    RegularGrid<T> GridZ_vol_emis;
    //! Volume emissions at current date.
    Data<T, 4> VolumeEmission;
    //! Volume emissions buffer.
    Data<T, 4> FileVolumeEmission_i;
    //! Volume emissions buffer.
    Data<T, 4> FileVolumeEmission_f;

    /*** Integration ***/

    //! Date between the current date and the next date.
    Date intermediate_date;

    //! Loss term (for two-step).
    Data<T, 4> Loss;
    //! Production term (for two-step).
    Data<T, 4> Production;

    /*** State ***/

    //! Concentrations at previous step.
    Data<T, 4> PreviousConcentration;

    /*** Numerical schemes ***/

    //! Transport numerical scheme.
    ClassTransport Transport_;

  public:

    /*** Constructor and destructor ***/

    CastorTransport(string config_file);
    virtual ~CastorTransport();

    /*** Configuration ***/

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    T GetMeteoDelta_t();

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

    virtual void SetDate(Date date);

    /*** Integration ***/

    void Forward();

    /*** Other methods ***/

    void ComputeAirDensity(Data<T, 3>& Temperature_, Data<T, 3>& Pressure_,
                           Data<T, 3>& AirDensity_);

  protected:

    virtual void AddTime(T time);
    virtual void SubtractTime(T time);

    virtual void InitAllData();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CASTORTRANSPORT_HXX
#endif
