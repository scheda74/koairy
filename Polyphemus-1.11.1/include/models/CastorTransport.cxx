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


#ifndef POLYPHEMUS_FILE_MODELS_CASTORTRANSPORT_CXX


#include "CastorTransport.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassTransport>
  CastorTransport<T, ClassTransport>::CastorTransport(string config_file):
    BaseModel<T>(config_file)
  {

    /*** Managed data ***/

    this->option_manage["initial_condition"] = true;
    this->option_manage["temperature"] = true;
    this->option_manage["pressure"] = true;
    this->option_manage["air_density"] = true;
    this->option_manage["horizontal_wind"] = true;
    this->option_manage["vertical_diffusion"] = true;
    this->option_manage["boundary_condition"] = true;
    this->option_manage["deposition_velocity"] = true;
    this->option_manage["volume_emission"] = true;
    this->option_manage["scavenging_coefficient"] = true;

    /*** Pointers to 3D data ***/

    this->D3_map["Altitude"] = &Altitude;
    this->D3_map["Altitude_i"] = &Altitude;

    this->D3_map["Thickness"] = &Thickness;
    this->D3_map["Thickness_i"] = &Thickness;

    this->D3_map["Temperature"] = &Temperature;
    this->D3_map["Temperature_i"] = &Temperature;

    this->D3_map["Pressure"] = &Pressure;
    this->D3_map["Pressure_i"] = &Pressure;

    this->D3_map["AirDensity"] = &AirDensity;
    this->D3_map["AirDensity_i"] = &AirDensity;

    this->D3_map["MeridionalWind"] = &MeridionalWind;
    this->D3_map["MeridionalWind_i"] = &MeridionalWind;

    this->D3_map["ZonalWind"] = &ZonalWind;
    this->D3_map["ZonalWind_i"] = &ZonalWind;

    this->D3_map["VerticalDiffusionCoefficient"]
      = &VerticalDiffusionCoefficient;
    this->D3_map["VerticalDiffusionCoefficient_i"]
      = &VerticalDiffusionCoefficient;

    this->D3_map["BoundaryCondition_z"] = &BoundaryCondition_z;
    this->D3_map["BoundaryCondition_z_i"] = &BoundaryCondition_z;

    this->D3_map["DepositionVelocity"] = &DepositionVelocity;
    this->D3_map["DepositionVelocity_i"] = &DepositionVelocity;

    /*** Pointers to 4D data ***/

    this->D4_map["BoundaryCondition_y"] = &BoundaryCondition_y;
    this->D4_map["BoundaryCondition_y_i"] = &BoundaryCondition_y;

    this->D4_map["BoundaryCondition_x"] = &BoundaryCondition_x;
    this->D4_map["BoundaryCondition_x_i"] = &BoundaryCondition_x;

    this->D4_map["VolumeEmission"] = &VolumeEmission;
    this->D4_map["VolumeEmission_i"] = &VolumeEmission;

    /*** Fields species lists ***/

    this->field_species["InitialCondition"] = &species_list_ic;
    this->field_species["BoundaryCondition_z"] = &species_list_bc;
    this->field_species["BoundaryCondition_y"] = &species_list_bc;
    this->field_species["BoundaryCondition_x"] = &species_list_bc;
    this->field_species["DepositionVelocity"] = &species_list_dep;
    this->field_species["VolumeEmission"] = &species_list_vol_emis;
  }


  //! Destructor.
  template<class T, class ClassTransport>
  CastorTransport<T, ClassTransport>::~CastorTransport()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[domain]");

    this->config.SetSection("[options]");
    this->config.PeekValue("With_transport",
                           this->option_process["with_transport"]);
    this->config.PeekValue("With_initial_condition",
                           this->option_process["with_initial_condition"]);
    this->config.PeekValue("Interpolated_initial_condition",
                           this->option_process["interpolated"]);
    this->config.PeekValue("With_boundary_condition",
                           this->option_process["with_boundary_condition"]);
    this->config.PeekValue("With_deposition",
                           this->option_process["with_deposition"]);
    this->config.PeekValue("With_volume_emission",
                           this->option_process["with_volume_emission"]);

    /*** Input files ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Meteorological files.
    this->input_files["meteo"].Read(data_description_file, "meteo");

    // Boundary conditions files.
    if (this->option_process["with_boundary_condition"])
      this->input_files["boundary_condition"].Read(data_description_file,
                                                   "boundary_condition");
    else
      this->input_files["boundary_condition"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["boundary_condition"].Begin();
         i != this->input_files["boundary_condition"].End(); i++)
      species_list_bc.push_back(i->first);
    Ns_bc = int(species_list_bc.size());

    // Initial conditions files.
    if (this->option_process["with_initial_condition"]
        && !this->option_process["interpolated"])
      this->input_files["initial_condition"].ReadFiles(data_description_file,
                                                       "initial_condition");
    else
      this->input_files["initial_condition"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["initial_condition"].Begin();
         i != this->input_files["initial_condition"].End(); i++)
      species_list_ic.push_back(i->first);
    if (this->option_process["with_initial_condition"]
        && this->option_process["interpolated"])
      Ns_ic = Ns_bc;
    else
      Ns_ic = int(species_list_ic.size());

    // Deposition velocities files.
    if (this->option_process["with_deposition"])
      this->input_files["deposition_velocity"].Read(data_description_file,
                                                    "deposition");
    else
      this->input_files["deposition_velocity"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["deposition_velocity"].Begin();
         i != this->input_files["deposition_velocity"].End(); i++)
      species_list_dep.push_back(i->first);
    Ns_dep = int(species_list_dep.size());

    // Volume emission files.
    if (this->option_process["with_volume_emission"])
      this->input_files["volume_emission"].Read(data_description_file,
                                                "volume_emission");
    else
      this->input_files["volume_emission"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["volume_emission"].Begin();
         i != this->input_files["volume_emission"].End(); i++)
      species_list_vol_emis.push_back(i->first);
    Ns_vol_emis = int(species_list_vol_emis.size());
    // Additional variable: number of levels.
    if (this->option_process["with_volume_emission"])
      {
        data_description_stream.SetSection("[volume_emission]");
        data_description_stream.PeekValue("Nz", "> 0", Nz_vol_emis);
      }
    else
      Nz_vol_emis = 0;
  }


  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::CheckConfiguration()
  {
    if (this->option_manage["temperature"]
        && this->input_files["meteo"]("Temperature").empty())
      throw "Temperature is needed but no input data file was provided.";
    if (this->option_manage["pressure"]
        && this->input_files["meteo"]("Pressure").empty())
      throw "Pressure is needed but no input data file was provided.";

    if (this->option_manage["horizontal_wind"]
        && this->input_files["meteo"]("ZonalWind").empty())
      throw "Zonal wind is needed but no input data file was provided.";
    if (this->option_manage["horizontal_wind"]
        && this->input_files["meteo"]("MeridionalWind").empty())
      throw "Meridional wind is needed but no input data file was provided.";

    if (this->option_manage["vertical_diffusion"]
        && this->input_files["meteo"]("VerticalDiffusion").empty())
      throw string("Vertical diffusion is needed but no input data file")
        + " was provided.";

    if (this->option_process["with_initial_condition"]
        && this->option_process["interpolated"]
        && !this->option_process["with_boundary_condition"])
      throw string("Initial conditions cannot be interpolated from")
        + " boundary conditions since boundary conditions are not activated.";
  }


  //! Returns the time step of meteorological data.
  /*!
    \return The time step (seconds) of meteorological data.
  */
  template<class T, class ClassTransport>
  T CastorTransport<T, ClassTransport>::GetMeteoDelta_t()
  {
    return this->input_files["meteo"].GetDelta_t();
  }


  //! Checks whether a species has initial conditions.
  /*!
    \param s species global index.
    \return True if the species has initial conditions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool CastorTransport<T, ClassTransport>::HasInitialCondition(int s) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(),
                this->GetSpeciesName(s)) != species_list_ic.end();
  }


  //! Checks whether a species has initial conditions.
  /*!
    \param name species name.
    \return True if the species has initial conditions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool
  CastorTransport<T, ClassTransport>::HasInitialCondition(string name) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(), name)
      != species_list_ic.end();
  }


  //! Returns the index in initial conditions of a given species.
  /*!
    \param s species global index.
    \return The species index in initial conditions.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>::InitialConditionIndex(int s) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(),
                this->GetSpeciesName(s)) - species_list_ic.begin();
  }


  //! Returns the index in initial conditions of a given species.
  /*!
    \param name species name.
    \return The species index in initial conditions.
  */
  template<class T, class ClassTransport>
  int
  CastorTransport<T, ClassTransport>::InitialConditionIndex(string name) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(), name)
      - species_list_ic.begin();
  }


  //! Returns the name of a species with initial conditions.
  /*!
    \param s species index in initial conditions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  string CastorTransport<T, ClassTransport>::InitialConditionName(int s) const
  {
    return species_list_ic.at(s);
  }


  //! Returns the name of a species with initial conditions.
  /*!
    \param s species index in initial conditions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  int
  CastorTransport<T, ClassTransport>::InitialConditionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(InitialConditionName(s));
  }


  //! Checks whether a species has boundary conditions.
  /*!
    \param s species global index.
    \return True if the species has boundary conditions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool CastorTransport<T, ClassTransport>::HasBoundaryCondition(int s) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(),
                this->GetSpeciesName(s)) != species_list_bc.end();
  }


  //! Checks whether a species has boundary conditions.
  /*!
    \param name species name.
    \return True if the species has boundary conditions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool
  CastorTransport<T, ClassTransport>::HasBoundaryCondition(string name) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(), name)
      != species_list_bc.end();
  }


  //! Returns the index in boundary conditions of a given species.
  /*!
    \param s species global index.
    \return The species index in boundary conditions.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>::BoundaryConditionIndex(int s) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(),
                this->GetSpeciesName(s)) - species_list_bc.begin();
  }


  //! Returns the index in boundary conditions of a given species.
  /*!
    \param name species name.
    \return The species index in boundary conditions.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>
  ::BoundaryConditionIndex(string name) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(), name)
      - species_list_bc.begin();
  }


  //! Returns the name of a species with boundary conditions.
  /*!
    \param s species index in boundary conditions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  string
  CastorTransport<T, ClassTransport>::BoundaryConditionName(int s) const
  {
    return species_list_bc.at(s);
  }


  //! Returns the name of a species with boundary conditions.
  /*!
    \param s species index in boundary conditions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>
  ::BoundaryConditionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(BoundaryConditionName(s));
  }


  //! Checks whether a species has deposition velocities.
  /*!
    \param s species global index.
    \return True if the species has deposition velocities, false otherwise.
  */
  template<class T, class ClassTransport>
  bool CastorTransport<T, ClassTransport>::HasDepositionVelocity(int s) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(),
                this->GetSpeciesName(s)) != species_list_dep.end();
  }


  //! Checks whether a species has deposition velocities.
  /*!
    \param name species name.
    \return True if the species has deposition velocities, false otherwise.
  */
  template<class T, class ClassTransport>
  bool CastorTransport<T, ClassTransport>
  ::HasDepositionVelocity(string name) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(), name)
      != species_list_dep.end();
  }


  //! Returns the index in deposition velocities of a given species.
  /*!
    \param s species global index.
    \return The species index in deposition velocities.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>::DepositionVelocityIndex(int s) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(),
                this->GetSpeciesName(s)) - species_list_dep.begin();
  }


  //! Returns the index in deposition velocities of a given species.
  /*!
    \param name species name.
    \return The species index in deposition velocities.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>
  ::DepositionVelocityIndex(string name) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(), name)
      - species_list_dep.begin();
  }


  //! Returns the name of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The species name.
  */
  template<class T, class ClassTransport>
  string
  CastorTransport<T, ClassTransport>::DepositionVelocityName(int s) const
  {
    return species_list_dep.at(s);
  }


  //! Returns the name of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The species name.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>
  ::DepositionVelocityGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(DepositionVelocityName(s));
  }


  //! Checks whether a species has volume emissions.
  /*!
    \param s species global index.
    \return True if the species has volume emissions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool CastorTransport<T, ClassTransport>::HasVolumeEmission(int s) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                this->GetSpeciesName(s)) != species_list_vol_emis.end();
  }


  //! Checks whether a species has volume emissions.
  /*!
    \param name species name.
    \return True if the species has volume emissions, false otherwise.
  */
  template<class T, class ClassTransport>
  bool
  CastorTransport<T, ClassTransport>::HasVolumeEmission(string name) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                name) != species_list_vol_emis.end();
  }


  //! Returns the index in volume emissions of a given species.
  /*!
    \param s species global index.
    \return The species index in volume emissions.
  */
  template<class T, class ClassTransport>
  int CastorTransport<T, ClassTransport>::VolumeEmissionIndex(int s) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                this->GetSpeciesName(s)) - species_list_vol_emis.begin();
  }


  //! Returns the index in volume emissions of a given species.
  /*!
    \param name species name.
    \return The species index in volume emissions.
  */
  template<class T, class ClassTransport>
  int
  CastorTransport<T, ClassTransport>::VolumeEmissionIndex(string name) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                name) - species_list_vol_emis.begin();
  }


  //! Returns the name of a species with volume emissions.
  /*!
    \param s species index in volume emissions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  string CastorTransport<T, ClassTransport>::VolumeEmissionName(int s) const
  {
    return species_list_vol_emis.at(s);
  }


  //! Returns the name of a species with volume emissions.
  /*!
    \param s species index in volume emissions.
    \return The species name.
  */
  template<class T, class ClassTransport>
  int
  CastorTransport<T, ClassTransport>::VolumeEmissionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(VolumeEmissionName(s));
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::Allocate()
  {
    BaseModel<T>::Allocate();

    /*** Mesh dimensions ***/

    CellWidth_x.resize(this->Nx);
    CellWidth_y.resize(this->Ny);

    /*** Boundary conditions ***/

    GridS_bc = RegularGrid<T>(Ns_bc);

    GridY4D_interf_bc = RegularGrid<T>(this->y_min - this->Delta_y,
                                       T(this->Ny + 1) * this->Delta_y, 2);
    GridX4D_interf_bc = RegularGrid<T>(this->x_min - this->Delta_x,
                                       T(this->Nx + 1) * this->Delta_x, 2);

    // Along Z.
    BoundaryCondition_z.Resize(this->GridS3D, this->GridY3D, this->GridX3D);
    FileBoundaryCondition_z_i.Resize(this->GridS3D, this->GridY3D,
                                     this->GridX3D);
    FileBoundaryCondition_z_f.Resize(this->GridS3D, this->GridY3D,
                                     this->GridX3D);

    // Along Y.
    BoundaryCondition_y.Resize(this->GridS4D, this->GridZ4D,
                               GridY4D_interf_bc, this->GridX4D);
    FileBoundaryCondition_y_i.Resize(this->GridS4D, this->GridZ4D,
                                     GridY4D_interf_bc, this->GridX4D);
    FileBoundaryCondition_y_f.Resize(this->GridS4D, this->GridZ4D,
                                     GridY4D_interf_bc, this->GridX4D);

    // Along X.
    BoundaryCondition_x.Resize(this->GridS4D, this->GridZ4D,
                               this->GridY4D, GridX4D_interf_bc);
    FileBoundaryCondition_x_i.Resize(this->GridS4D, this->GridZ4D,
                                     this->GridY4D, GridX4D_interf_bc);
    FileBoundaryCondition_x_f.Resize(this->GridS4D, this->GridZ4D,
                                     this->GridY4D, GridX4D_interf_bc);

    /*** Coordinates ***/

    Altitude.Resize(this->Nz + 1, this->Ny, this->Nx);
    FileAltitude_i.Resize(this->Nz + 1, this->Ny, this->Nx);
    FileAltitude_f.Resize(this->Nz + 1, this->Ny, this->Nx);

    Thickness.Resize(this->Nz, this->Ny, this->Nx);

    /*** Temperature ***/

    Temperature.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Pressure ***/

    Pressure.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FilePressure_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FilePressure_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Air density ***/

    AirDensity.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileAirDensity_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileAirDensity_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Winds ***/

    ZonalWind.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileZonalWind_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileZonalWind_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    MeridionalWind.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileMeridionalWind_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileMeridionalWind_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Diffusion coefficients ***/

    VerticalDiffusionCoefficient.Resize(this->GridZ3D_interf, this->GridY3D,
                                        this->GridX3D);
    FileVerticalDiffusionCoefficient_i.Resize(this->GridZ3D_interf,
                                              this->GridY3D, this->GridX3D);
    FileVerticalDiffusionCoefficient_f.Resize(this->GridZ3D_interf,
                                              this->GridY3D, this->GridX3D);

    /*** Deposition velocities ***/

    GridS_dep = RegularGrid<T>(Ns_dep);

    DepositionVelocity.Resize(GridS_dep, this->GridY3D, this->GridX3D);
    FileDepositionVelocity_i.Resize(GridS_dep, this->GridY3D, this->GridX3D);
    FileDepositionVelocity_f.Resize(GridS_dep, this->GridY3D, this->GridX3D);

    /*** Volume emissions ***/

    GridS_vol_emis = RegularGrid<T>(Ns_vol_emis);
    GridZ_vol_emis = RegularGrid<T>(Nz_vol_emis);

    VolumeEmission.Resize(GridS_vol_emis, GridZ_vol_emis,
                          this->GridY4D, this->GridX4D);
    FileVolumeEmission_i.Resize(GridS_vol_emis, GridZ_vol_emis,
                                this->GridY4D, this->GridX4D);
    FileVolumeEmission_f.Resize(GridS_vol_emis, GridZ_vol_emis,
                                this->GridY4D, this->GridX4D);

    /*** State ***/

    PreviousConcentration.Resize(this->GridS4D, this->GridZ4D,
                                 this->GridY4D, this->GridX4D);
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::Init()
  {
    BaseModel<T>::Init();

    /*** Coordinates ***/

    CellWidth_x.resize(this->Ny);
    const T pi = 3.14159265358979323846264;
    T ratio = pi / 180.;
    const T Earth_radius = 6371000.;
    // const T Earth_radius = 6371229.;
    for (int j = 0; j < this->Ny; j++)
      {
        CellWidth_x(j) = this->Delta_x * ratio * Earth_radius
          * cos(ratio * this->GridY3D(j));
        CellWidth_y(j) = this->Delta_y * ratio * Earth_radius;
      }

    /*** Species data ***/

    ConfigStream species_data_stream(this->GetSpeciesFile());
    string species;

    /*** Input data ***/

    intermediate_date = this->current_date;
    intermediate_date.AddSeconds(this->Delta_t / 2.);
    CastorTransport<T, ClassTransport>::InitAllData();

    /*** Initial conditions ***/

    if (this->option_manage["initial_condition"])
      {
        this->Concentration.SetZero();

        if (this->option_process["with_initial_condition"]
            && !this->option_process["interpolated"])
          for (int i = 0; i < Ns_ic; i++)
            {
              string filename
                = this->input_files["initial_condition"](species_list_ic[i]);
              int index = this->GetSpeciesIndex(species_list_ic[i]);
              Data<T, 3>
                Concentration_tmp(&(this->Concentration)(index, 0, 0, 0),
                                  shape(this->Nz, this->Ny, this->Nx));
              if (is_num(filename))
                Concentration_tmp.Fill(to_num<T>(filename));
              else
                FormatBinary<float>().Read(filename, Concentration_tmp);
            }

        if (this->option_process["with_initial_condition"]
            && this->option_process["interpolated"])
          for (int s = 0; s < this->Ns; s++)
            for (int k = 0; k < this->Nz; k++)
              for (int j = 0; j < this->Ny; j++)
                for (int i = 0; i < this->Nx; i++)
                  {
                    T conc_west = this->BoundaryCondition_x(s, k, j, 0)
                      / AirDensity(k, j, 0);
                    T conc_east = this->BoundaryCondition_x(s, k, j, 1)
                      / AirDensity(k, j, this->Nx - 1);
                    T conc_south = this->BoundaryCondition_y(s, k, 0, i)
                      / AirDensity(k, 0, i);
                    T conc_north = this->BoundaryCondition_y(s, k, 1, i)
                      / AirDensity(k, this->Ny - 1, i);
                    conc_west *= .5 * T(this->Nx - i) / T(this->Nx + 1.);
                    conc_east *= .5 * T(i + 1) / T(this->Nx + 1);
                    conc_south *= .5 * T(this->Ny - j) / T(this->Ny + 1.);
                    conc_north *= .5 * T(j + 1) / T(this->Ny + 1);
                    this->Concentration(s, k, j, i)
                      = conc_west + conc_east + conc_north + conc_south;
                    this->Concentration(s, k, j, i) *= AirDensity(k, j, i);
                  }
      }

    PreviousConcentration() = this->Concentration();

    /*** Numerical schemes ***/

    if (this->option_process["with_transport"])
      Transport_.Init(*this);
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::InitStep()
  {
    int k, j, i;

    /*** Altitude ***/

    this->InitData("meteo", "Altitude", FileAltitude_i,
                   FileAltitude_f, this->intermediate_date, Altitude);
    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
        for (i = 0; i < this->Nx; i++)
          Thickness(k, j, i) = Altitude(k + 1, j, i) - Altitude(k, j, i);

    /*** Temperature ***/

    if (this->option_manage["temperature"])
      this->InitData("meteo", "Temperature", FileTemperature_i,
                     FileTemperature_f, intermediate_date, Temperature);

    /*** Pressure ***/

    if (this->option_manage["pressure"])
      this->InitData("meteo", "Pressure", FilePressure_i,
                     FilePressure_f, intermediate_date, Pressure);

    /*** Air density ***/

    ComputeAirDensity(Temperature, Pressure, AirDensity);
    ComputeAirDensity(FileTemperature_i, FilePressure_i, FileAirDensity_i);
    ComputeAirDensity(FileTemperature_f, FilePressure_f, FileAirDensity_f);

    /*** Winds ***/

    if (this->option_manage["horizontal_wind"])
      {
        this->InitData("meteo", "ZonalWind", FileZonalWind_i,
                       FileZonalWind_f, intermediate_date, ZonalWind);
        this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
                       FileMeridionalWind_f, intermediate_date,
                       MeridionalWind);
      }

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      this->UpdateData("meteo", "VerticalDiffusion",
                       FileVerticalDiffusionCoefficient_i,
                       FileVerticalDiffusionCoefficient_f,
                       VerticalDiffusionCoefficient);

    /*** Boundary conditions ***/

    if (this->option_manage["boundary_condition"])
      {
        BoundaryCondition_z() = 1.e-10;
        BoundaryCondition_y() = 1.e-10;
        BoundaryCondition_x() = 1.e-10;
      }

    if (this->option_manage["boundary_condition"])
      for (int i = 0; i < Ns_bc; i++)
        {
          int j = BoundaryConditionGlobalIndex(i);
          string filename
            = this->input_files["boundary_condition"](species_list_bc[i]);
          Date date = this->input_files["boundary_condition"].GetDateMin();
          T Delta_t = this->input_files["boundary_condition"].GetDelta_t();

          this->InitData(find_replace(filename, "&c", "z"), date, Delta_t,
                         FileBoundaryCondition_z_i, FileBoundaryCondition_z_f,
                         this->intermediate_date, j, BoundaryCondition_z);
          this->InitData(find_replace(filename, "&c", "y"), date, Delta_t,
                         FileBoundaryCondition_y_i, FileBoundaryCondition_y_f,
                         this->intermediate_date, j, BoundaryCondition_y);
          this->InitData(find_replace(filename, "&c", "x"), date, Delta_t,
                         FileBoundaryCondition_x_i, FileBoundaryCondition_x_f,
                         this->intermediate_date, j, BoundaryCondition_x);
        }

    if (this->option_manage["boundary_condition"])
      {
        int s, k, j, i;
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              {
                BoundaryCondition_x(s, k, j, 0)
                  *= AirDensity(k, j, 0) * 1.e-9;
                BoundaryCondition_x(s, k, j, 1)
                  *= AirDensity(k, j, this->Nx - 1) * 1.e-9;
              }

        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (i = 0; i < this->Nx; i++)
              {
                BoundaryCondition_y(s, k, 0, i)
                  *= AirDensity(k, 0, i) * 1.e-9;
                BoundaryCondition_y(s, k, 1, i)
                  *= AirDensity(k, this->Ny - 1, i) * 1.e-9;
              }

        for (s = 0; s < this->Ns; s++)
          for (j = 0; j < this->Ny; j++)
            for (i = 0; i < this->Nx; i++)
              BoundaryCondition_z(s, j, i)
                *= AirDensity(this->Nz - 1, j, i) * 1.e-9;
      }

    /*** Deposition velocities ***/

    if (this->option_manage["deposition_velocity"])
      for (int i = 0; i < Ns_dep; i++)
        this->InitData("deposition_velocity", species_list_dep[i],
                       FileDepositionVelocity_i, FileDepositionVelocity_f,
                       this->intermediate_date, i, DepositionVelocity);

    /*** Volume emissions ***/

    if (this->option_manage["volume_emission"])
      for (int i = 0; i < Ns_vol_emis; i++)
        this->InitData("volume_emission", species_list_vol_emis[i],
                       FileVolumeEmission_i, FileVolumeEmission_f,
                       this->intermediate_date, i, VolumeEmission, false);

    /*** Transport initialization ***/

    if (this->option_process["with_transport"])
      Transport_.InitStep(*this);
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::SetDate(Date date)
  {
    BaseModel<T>::SetDate(date);
    intermediate_date = this->current_date;
    intermediate_date.AddSeconds(this->Delta_t / 2.);
    CastorTransport<T, ClassTransport>::InitAllData();
  }


  ///////////////////////////
  // NUMERICAL INTEGRATION //
  ///////////////////////////


  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    adds volume emissions. These three processes are split (operator
    splitting).
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::Forward()
  {
    int s, h, k, j, i;

    Data<T, 4> Concentration_tmp(this->GridS4D, this->GridZ4D,
                                 this->GridY4D, this->GridX4D);

    T conc;
    for (s = 0; s < this->Ns; s++)
      for (k = 0; k < this->Nz; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              Concentration_tmp(s, k, j, i)
                = (4. *  this->Concentration(s, k, j, i)
                   -  this->PreviousConcentration(s, k, j, i)) / 3.;
              conc = this->Concentration(s, k, j, i);
              this->Concentration(s, k, j, i) = 2. * conc
                - this->PreviousConcentration(s, k, j, i);
              this->PreviousConcentration(s, k, j, i) = conc;
            }
    this->Concentration.ThresholdMin(1.);
    Concentration_tmp.ThresholdMin(1.);

    int Niter = 1;
    // During the first hour (spin-up).
    if (T(this->step) * this->Delta_t < 3600.)
      Niter = 5;

    // To speed up tests.
    bool deposition_velocity(this->option_manage["deposition_velocity"]);
    bool with_transport(this->option_process["with_transport"]);

    T loss, production;
    T factor = 2. / 3. * this->Delta_t;
    for (h = 0; h < Niter; h++)
      for (k = 0; k < this->Nz; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            for (s = 0; s < this->Ns; s++)
              {
                loss = 0.;
                production = 0.;

                // Emissions.
                if (k < this->Nz_vol_emis && this->HasVolumeEmission(s))
                  {
                    int e = this->VolumeEmissionIndex(s);
                    production += this->VolumeEmission(e, k, j, i) / 100.
                      / (this->Altitude(k + 1, j, i)
                         - this->Altitude(k, j, i));
                  }

                // Deposition.
                if (deposition_velocity && k == 0
                    && this->HasDepositionVelocity(s))
                  {
                    int e = this->DepositionVelocityIndex(s);
                    loss += this->DepositionVelocity(e, j, i)
                      * this->Concentration(s, k, j, i);
                  }

                if (with_transport)
                  Transport_.LossProduction(*this, s, k, j, i,
                                            loss, production);

                this->Concentration(s, k, j, i) =
                  (Concentration_tmp(s, k, j, i) + factor * production)
                  / (1. + factor * loss
                     / this->Concentration(s, k, j, i));
              }

    this->AddTime(this->Delta_t);
    this->step++;
  }


  ///////////////////
  // OTHER METHODS //
  ///////////////////


  //! Computes air density on the basis of temperature and pressure.
  /*! Formula: AirDensity = Pressure / (287. * Temperature).
    \param Temperature_ temperature.
    \param Pressure_ pressure.
    \param AirDensity_ (output) air density.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>
  ::ComputeAirDensity(Data<T, 3>& Temperature_, Data<T, 3>& Pressure_,
                      Data<T, 3>& AirDensity_)
  {
    int k, j, i;

    int Nz = AirDensity_.GetLength(0);
    int Ny = AirDensity_.GetLength(1);
    int Nx = AirDensity_.GetLength(2);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          AirDensity_(k, j, i) = 7.2868e16 * Pressure_(k, j, i)
            / Temperature_(k, j, i);
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Adds time to current time.
  /*! Sets the previous date to the (old) current date, adds \a time seconds
    to the current date, sets the next date to the (new) current date plus
    \a time seconds, and sets the intermediate date accordingly.
    \param time seconds to be added.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::AddTime(T time)
  {
    BaseModel<T>::AddTime(time);
    intermediate_date = this->current_date;
    intermediate_date.AddSeconds(this->Delta_t / 2.);
  }


  //! Subtracts time to current time.
  /*! Sets the next date to the (old) current date, subtracts \a time seconds
    to the current date, sets the previous date to the (new) current date
    minus \a time seconds, and sets the intermediate date accordingly.
    \param time seconds to be subtracted.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::SubtractTime(T time)
  {
    BaseModel<T>::AddTime(time);
    intermediate_date = this->current_date;
    intermediate_date.AddSeconds(this->Delta_t / 2.);
  }


  //! Moves model input-data to the current date.
  /*! This method prepares the model for a time integration from the current
    date. It reads input data to be read before InitStep and Forward.
  */
  template<class T, class ClassTransport>
  void CastorTransport<T, ClassTransport>::InitAllData()
  {
    int k, j, i;

    /*** Altitude ***/

    this->InitData("meteo", "Altitude", FileAltitude_i,
                   FileAltitude_f, this->intermediate_date, Altitude);
    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
        for (i = 0; i < this->Nx; i++)
          Thickness(k, j, i) = Altitude(k + 1, j, i) - Altitude(k, j, i);

    /*** Pressure, temperature and air density ***/

    if (this->option_manage["temperature"])
      this->InitData("meteo", "Temperature", FileTemperature_i,
                     FileTemperature_f, this->intermediate_date, Temperature);
    if (this->option_manage["pressure"])
      this->InitData("meteo", "Pressure", FilePressure_i,
                     FilePressure_f, this->intermediate_date, Pressure);

    ComputeAirDensity(Temperature, Pressure, AirDensity);
    ComputeAirDensity(FileTemperature_i, FilePressure_i, FileAirDensity_i);
    ComputeAirDensity(FileTemperature_f, FilePressure_f, FileAirDensity_f);

    /*** Winds ***/

    if (this->option_manage["horizontal_wind"])
      {
        this->InitData("meteo", "ZonalWind", FileZonalWind_i,
                       FileZonalWind_f, intermediate_date, ZonalWind);
        this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
                       FileMeridionalWind_f, intermediate_date,
                       MeridionalWind);
      }

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      this->InitData("meteo", "VerticalDiffusion",
                     FileVerticalDiffusionCoefficient_i,
                     FileVerticalDiffusionCoefficient_f,
                     this->intermediate_date, VerticalDiffusionCoefficient);

    /*** Boundary conditions ***/


    if (this->option_manage["boundary_condition"])
      {
        BoundaryCondition_z() = 1.e-10;
        BoundaryCondition_y() = 1.e-10;
        BoundaryCondition_x() = 1.e-10;
      }

    if (this->option_manage["boundary_condition"])
      for (int i = 0; i < Ns_bc; i++)
        {
          int j = BoundaryConditionGlobalIndex(i);
          string filename
            = this->input_files["boundary_condition"](species_list_bc[i]);
          Date date = this->input_files["boundary_condition"].GetDateMin();
          T Delta_t = this->input_files["boundary_condition"].GetDelta_t();

          this->InitData(find_replace(filename, "&c", "z"), date, Delta_t,
                         FileBoundaryCondition_z_i, FileBoundaryCondition_z_f,
                         this->intermediate_date, j, BoundaryCondition_z);
          this->InitData(find_replace(filename, "&c", "y"), date, Delta_t,
                         FileBoundaryCondition_y_i, FileBoundaryCondition_y_f,
                         this->intermediate_date, j, BoundaryCondition_y);
          this->InitData(find_replace(filename, "&c", "x"), date, Delta_t,
                         FileBoundaryCondition_x_i, FileBoundaryCondition_x_f,
                         this->intermediate_date, j, BoundaryCondition_x);
        }

    if (this->option_manage["boundary_condition"])
      {
        int s, k, j, i;
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              {
                BoundaryCondition_x(s, k, j, 0)
                  *= AirDensity(k, j, 0) * 1.e-9;
                BoundaryCondition_x(s, k, j, 1)
                  *= AirDensity(k, j, this->Nx - 1) * 1.e-9;
              }

        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (i = 0; i < this->Nx; i++)
              {
                BoundaryCondition_y(s, k, 0, i)
                  *= AirDensity(k, 0, i) * 1.e-9;
                BoundaryCondition_y(s, k, 1, i)
                  *= AirDensity(k, this->Ny - 1, i) * 1.e-9;
              }

        for (s = 0; s < this->Ns; s++)
          for (j = 0; j < this->Ny; j++)
            for (i = 0; i < this->Nx; i++)
              BoundaryCondition_z(s, j, i)
                *= AirDensity(this->Nz - 1, j, i) * 1.e-9;
      }

    /*** Deposition velocities ***/

    if (this->option_manage["deposition_velocity"])
      for (int i = 0; i < Ns_dep; i++)
        this->InitData("deposition_velocity", species_list_dep[i],
                       FileDepositionVelocity_i, FileDepositionVelocity_f,
                       this->intermediate_date, i, DepositionVelocity);

    /*** Volume emissions ***/

    if (this->option_manage["volume_emission"])
      for (int i = 0; i < Ns_vol_emis; i++)
        this->InitData("volume_emission", species_list_vol_emis[i],
                       FileVolumeEmission_i, FileVolumeEmission_f,
                       this->intermediate_date, i, VolumeEmission, false);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CASTORTRANSPORT_CXX
#endif
