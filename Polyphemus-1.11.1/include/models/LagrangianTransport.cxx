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



#ifndef POLYPHEMUS_FILE_MODELS_LAGRANGIANTRANSPORT_CXX


#include "LagrangianTransport.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassParticle>
  LagrangianTransport<T, ClassParticle>
  ::LagrangianTransport(string config_file): BaseModel<T>(config_file)
  {

    /*** Random generator initialization ***/

    urng = new NEWRAN::LGM_mixed(0.46875);
    NEWRAN::Random::Set(*urng);


    /*** Managed data ***/

    this->option_manage["temperature"] = true;
    this->option_manage["pressure"] = true;
    this->option_manage["horizontal_wind"] = true;
    this->option_manage["vertical_wind"] = true;
    this->option_manage["horizontal_diffusion"] = true;
    this->option_manage["vertical_diffusion"] = true;

    /*** Pointers to 3D data ***/

    this->D3_map["Temperature"] = &Temperature_f;
    this->D3_map["Temperature_i"] = &Temperature_i;
    this->D3_map["Temperature_f"] = &Temperature_f;

    this->D3_map["Pressure"] = &Pressure_f;
    this->D3_map["Pressure_i"] = &Pressure_i;
    this->D3_map["Pressure_f"] = &Pressure_f;

    this->D3_map["AirDensity"] = &AirDensity_f;
    this->D3_map["AirDensity_i"] = &AirDensity_i;
    this->D3_map["AirDensity_f"] = &AirDensity_f;

    this->D3_map["VerticalWind"] = &VerticalWind_i;
    this->D3_map["VerticalWind_i"] = &VerticalWind_i;

    this->D3_map["MeridionalWind"] = &MeridionalWind_i;
    this->D3_map["MeridionalWind_i"] = &MeridionalWind_i;

    this->D3_map["ZonalWind"] = &ZonalWind_i;
    this->D3_map["ZonalWind_i"] = &ZonalWind_i;

    this->D3_map["VerticalDiffusionCoefficient"]
      = &VerticalDiffusionCoefficient_f;
    this->D3_map["VerticalDiffusionCoefficient_i"]
      = &VerticalDiffusionCoefficient_i;
    this->D3_map["VerticalDiffusionCoefficient_f"]
      = &VerticalDiffusionCoefficient_f;

    this->D3_map["MeridionalDiffusionCoefficient"]
      = &MeridionalDiffusionCoefficient_i;
    this->D3_map["MeridionalDiffusionCoefficient_i"]
      = &MeridionalDiffusionCoefficient_i;

    this->D3_map["ZonalDiffusionCoefficient"]
      = &ZonalDiffusionCoefficient_i;
    this->D3_map["ZonalDiffusionCoefficient_i"]
      = &ZonalDiffusionCoefficient_i;
  }


  //! Destructor.
  template<class T, class ClassParticle>
  LagrangianTransport<T, ClassParticle>::~LagrangianTransport()
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
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[domain]");

    this->config.SetSection("[options]");
    this->config.PeekValue("With_air_density",
                           this->option_process["with_air_density"]);
    this->config.PeekValue("With_point_emission",
                           this->option_process["with_point_emission"]);

    /*** Disables management of useless fields ***/

    // Pressure and temperature are only needed for air density
    // or microphysical model of scavenging.
    this->option_manage["temperature"] = this->option_manage["temperature"]
      && this->option_process["with_air_density"];

    this->option_manage["pressure"] = this->option_manage["pressure"]
      && this->option_process["with_air_density"];

    /*** Input files ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Meteorological files.
    this->input_files["meteo"].Read(data_description_file, "meteo");

    // Point emissions.
    if (this->option_process["with_point_emission"])
      {
        string point_emission_file;
        data_description_stream.SetSection("[point_emission]");
        data_description_stream.PeekValue("file", point_emission_file);
        this->PointEmissionManager = new BasePointEmission<T>();
        this->PointEmissionManager->Init(point_emission_file,
                                         this->species_list);
        data_description_stream.PeekValue("Delta_t_particle_emission",
                                          "positive",
                                          Delta_t_particle_emission);
      }

    /*** Rest of the configuration ***/

    // Horizontal diffusion coefficient.
    this->config.PeekValue("Horizontal_diffusion", "positive",
                           horizontal_diffusion);
    // Horizontal diffusion coefficient for the gaussian kernels.
    this->config.PeekValue("Gaussian_kernel_horizontal_diffusion", "positive",
                           gaussian_kernel_horizontal_diffusion);

    /*** Species data ***/

    ConfigStream species_data_stream(this->GetSpeciesFile());
    string species;
  }


  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::CheckConfiguration()
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
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::Allocate()
  {
    BaseModel<T>::Allocate();

    /*** Mesh dimensions ***/

    CellWidth_x.resize(this->Nx);
    CellWidth_y.resize(this->Ny);
    CellWidth_z.resize(this->Nz);

    CellCenterDistance_x.resize(this->Nx - 1);
    CellCenterDistance_y.resize(this->Ny - 1);
    CellCenterDistance_z.resize(this->Nz - 1);

    /*** Temperature ***/

    Temperature_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    Temperature_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Pressure ***/

    Pressure_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    Pressure_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FilePressure_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FilePressure_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Air density ***/

    AirDensity_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    AirDensity_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    AirDensity_interf_z_i.Resize(this->GridZ3D_interf, this->GridY3D,
                                 this->GridX3D);
    AirDensity_interf_z_f.Resize(this->GridZ3D_interf, this->GridY3D,
                                 this->GridX3D);
    AirDensity_interf_y_i.Resize(this->GridZ3D, this->GridY3D_interf,
                                 this->GridX3D);
    AirDensity_interf_y_f.Resize(this->GridZ3D, this->GridY3D_interf,
                                 this->GridX3D);
    AirDensity_interf_x_i.Resize(this->GridZ3D, this->GridY3D,
                                 this->GridX3D_interf);
    AirDensity_interf_x_f.Resize(this->GridZ3D, this->GridY3D,
                                 this->GridX3D_interf);

    /*** Winds ***/

    ZonalWind_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D_interf);
    FileZonalWind_i.Resize(this->GridZ3D, this->GridY3D,
                           this->GridX3D_interf);
    FileZonalWind_f.Resize(this->GridZ3D, this->GridY3D,
                           this->GridX3D_interf);

    MeridionalWind_i.Resize(this->GridZ3D, this->GridY3D_interf,
                            this->GridX3D);
    FileMeridionalWind_i.Resize(this->GridZ3D, this->GridY3D_interf,
                                this->GridX3D);
    FileMeridionalWind_f.Resize(this->GridZ3D, this->GridY3D_interf,
                                this->GridX3D);

    VerticalWind_i.Resize(this->GridZ3D_interf, this->GridY3D, this->GridX3D);

    /*** Diffusion coefficients ***/

    ZonalDiffusionCoefficient_i.Resize(this->GridZ3D, this->GridY3D,
                                       this->GridX3D_interf);
    MeridionalDiffusionCoefficient_i.Resize(this->GridZ3D,
                                            this->GridY3D_interf,
                                            this->GridX3D);

    VerticalDiffusionCoefficient_i.Resize(this->GridZ3D_interf, this->GridY3D,
                                          this->GridX3D);
    VerticalDiffusionCoefficient_f.Resize(this->GridZ3D_interf, this->GridY3D,
                                          this->GridX3D);
    FileVerticalDiffusionCoefficient_i.Resize(this->GridZ3D_interf,
                                              this->GridY3D, this->GridX3D);
    FileVerticalDiffusionCoefficient_f.Resize(this->GridZ3D_interf,
                                              this->GridY3D, this->GridX3D);
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::Init()
  {
    BaseModel<T>::Init();

    /*** Initialization ***/
    this->Concentration.SetZero();
    x_max = this->x_min + this->Delta_x * (this->Nx - 1);
    y_max = this->y_min + this->Delta_y * (this->Ny - 1);

    /*** Mesh dimensions ***/

    ComputeCellWidth(this->Delta_x, this->GridY3D_interf.GetArray(),
                     CellWidth_x, CellWidth_y);
    for (int k = 0; k < this->Nz; k++)
      CellWidth_z(k) = this->GridZ3D_interf(k + 1) - this->GridZ3D_interf(k);

    ComputeCellCenterDistance(this->Delta_x,
                              this->GridY3D_interf.GetArray(),
                              CellCenterDistance_x, CellCenterDistance_y);
    for (int k = 0; k < this->Nz - 1; k++)
      CellCenterDistance_z(k) = this->GridZ3D(k + 1) - this->GridZ3D(k);

    /*** Input data ***/

    LagrangianTransport<T, ClassParticle>
      ::SetDate(this->Date_min);
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are needed for the current step.
   */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::InitStep()
  {

    /*** Temperature ***/

    if (this->option_manage["temperature"])
      this->UpdateData("meteo", "Temperature", FileTemperature_i,
                       FileTemperature_f, Temperature_i, Temperature_f);

    /*** Pressure ***/

    if (this->option_manage["pressure"])
      this->UpdateData("meteo", "Pressure", FilePressure_i,
                       FilePressure_f, Pressure_i, Pressure_f);

    /*** Air density ***/

    if (this->option_process["with_air_density"])
      {
        AirDensity_i.GetArray() = AirDensity_f.GetArray();
        ComputeAirDensity(Temperature_f, Pressure_f, AirDensity_f);

        AirDensity_interf_z_i.GetArray() = AirDensity_interf_z_f.GetArray();
        AirDensity_interf_y_i.GetArray() = AirDensity_interf_y_f.GetArray();
        AirDensity_interf_x_i.GetArray() = AirDensity_interf_x_f.GetArray();
      }

    /*** Winds ***/

    if (this->option_manage["horizontal_wind"])
      {
        this->UpdateData("meteo", "ZonalWind", FileZonalWind_i,
                         FileZonalWind_f, ZonalWind_i);
        this->UpdateData("meteo", "MeridionalWind", FileMeridionalWind_i,
                         FileMeridionalWind_f, MeridionalWind_i);
      }

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      {
        this->UpdateData("meteo", "VerticalDiffusion",
                         FileVerticalDiffusionCoefficient_i,
                         FileVerticalDiffusionCoefficient_f,
                         VerticalDiffusionCoefficient_i,
                         VerticalDiffusionCoefficient_f);
      }


    /*** Particles state (age, height...) ***/
    for_each(particles_list.begin(), particles_list.end(),
             bind2nd(mem_fun_ref(&ClassParticle::Update), this->Delta_t));
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date Date.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::SetDate(Date date)
  {
    BaseModel<T>::SetDate(date);
    LagrangianTransport<T, ClassParticle>::InitAllData();
  }


  ///////////////////////////
  // NUMERICAL INTEGRATION //
  ///////////////////////////


  //! Performs one advection-diffusion step.
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::Transport()
  {
    bool ind_print = true;
    for (typename list<ClassParticle>::iterator
           particle = particles_list.begin();
         particle != particles_list.end(); particle++)
      particle->Transport(this);
  }


  //! Performs the emission of particles at point emissions.
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::PointEmission()
  {
    int Nemis = this->PointEmissionManager->GetNumberEmission();
    for (int emission = 0; emission < Nemis; emission++)
      {
        if (this->PointEmissionManager->
            IsEmitting(this->current_date, this->next_date, emission))
          {
            T lon_emission, lat_emission, z_emission;
            this->PointEmissionManager->
              GetEmissionCoordinates(lon_emission, lat_emission, z_emission,
                                     emission);

            T emission_timespan = this->PointEmissionManager->
              GetEmissionTimespan(this->current_date, this->next_date,
                                  emission);
            vector<int> emitted_species_index =
              this->PointEmissionManager->GetEmittedSpeciesIndex(emission);
            int Ns = emitted_species_index.size();
            vector<T> quantity(Ns);

            int N_emitted_particles =
              int(emission_timespan / Delta_t_particle_emission);

            T effective_Delta_t_particle_emission = emission_timespan /
              N_emitted_particles;

            // Computes the mass attributed to each particle.
            for (int species = 0; species < Ns; species++)
              {
                T rate = this->PointEmissionManager->GetRate(emission,
                                                             species);
                quantity[species] = rate
                  * effective_Delta_t_particle_emission;
              }

            // Creates the emitted particles.
            for (int count = 0; count < N_emitted_particles; count++)
              {
                // Particle age is initialized to a negative value in order
                // to get the right age at the end of the current timestep.
                ClassParticle
                  emitted_particle(z_emission, lat_emission, lon_emission,
                                   emitted_species_index, quantity,
                                   emission_timespan - count
                                   * effective_Delta_t_particle_emission);
                particles_list.push_back(emitted_particle);
              }
          }
      }
  }


  //! Performs one step forward.
  /*! It performs one advection-diffusion step and then
    adds emissions. These two processes are split (operator
    splitting).
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::Forward()
  {
    /*** Air density ***/

    if (this->option_process["with_air_density"])
      {
        InterpolateInterface_z(CellCenterDistance_z, CellWidth_z,
                               AirDensity_f, AirDensity_interf_z_f);
        InterpolateInterface_y(CellCenterDistance_y, CellWidth_y,
                               AirDensity_f, AirDensity_interf_y_f);
        InterpolateInterface_x(CellCenterDistance_x, CellWidth_x,
                               AirDensity_f, AirDensity_interf_x_f);
      }

    /*** Winds ***/

    if (this->option_manage["horizontal_wind"])
      {
        TransformZonalWind(ZonalWind_i);
        TransformMeridionalWind(MeridionalWind_i);
      }

    if (this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
        ComputeVerticalWind(CellWidth_x, CellWidth_y, CellWidth_z,
                            AirDensity_interf_x_i, ZonalWind_i,
                            AirDensity_interf_y_i, MeridionalWind_i,
                            AirDensity_interf_z_i, VerticalWind_i);
      else
        ComputeVerticalWind(CellWidth_x, CellWidth_y, CellWidth_z,
                            ZonalWind_i, MeridionalWind_i, VerticalWind_i);

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      {
        // Computes rho * Kz.
        if (this->option_process["with_air_density"])
          VerticalDiffusionCoefficient_f.GetArray() =
            AirDensity_interf_z_f.GetArray()
            * VerticalDiffusionCoefficient_f.GetArray();
      }

    /*** Emission ***/

    if (this->option_process["with_point_emission"])
      PointEmission();

    /*** Time integration ***/

    Transport();

    // Removes the particles that are considered outside of the domain.
    IsOutsideOfDomain<T> ParticleIsOutsideOfDomain(this);
    particles_list.remove_if(ParticleIsOutsideOfDomain);

    this->AddTime(this->Delta_t);
    this->step++;
  }


  //////////////////////////////
  // RANDOM NUMBER GENERATION //
  //////////////////////////////


  //! Returns a uniform random number from the range (-sqrt(3), sqrt(3)).
  template<class T, class ClassParticle>
  T LagrangianTransport<T, ClassParticle>::
  GetRandomNumber()
  {
    return (3.464102 * RandomNumber.Next() - 1.732051);
  }


  ///////////////////
  // OTHER METHODS //
  ///////////////////


  //! Computes cells widths (in meters).
  /*!
    \param Del ta_x_ step along x in degrees.
    \param GridY_interf_ interface coordinates along y in degrees.
    \param CellWidth_x_ (output) cells widths along x in meters.
    \param CellWidth_y_ (output) cells widths along y in meters.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeCellWidth(T Delta_x_, Array<T, 1>& GridY_interf_,
                     Array<T, 1>& CellWidth_x_, Array<T, 1>& CellWidth_y_)
  {
    const T earth_radius = 6371229.;
    const T pi(3.14159265358979323846264);

    for (int i = 0; i < this->Nx; i++)
      CellWidth_x_(i) = earth_radius * Delta_x_ * pi / 180.;

    Array<T, 1> GridY_interf_trans(this->Ny + 1);
    for (int j = 0; j < this->Ny + 1; j++)
      GridY_interf_trans(j)
        = earth_radius * sin(GridY_interf_(j) * pi / 180.);

    for (int j = 0; j < this->Ny; j++)
      CellWidth_y_(j) = GridY_interf_trans(j + 1) - GridY_interf_trans(j);
  }


  //! Computes the distances between cells centers (in meters).
  /*!
    \param Delta_x_ step along x in degrees.
    \param GridY_interf_ interface coordinates along y in degrees.
    \param CellCenterDistance_x_ (output) distances between cells centers
    along x in meters.
    \param CellCenterDistance_y_ (output) distances between cells centers
    along y in meters.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeCellCenterDistance(T Delta_x_, Array<T, 1>& GridY_interf_,
                              Array<T, 1>& CellCenterDistance_x_,
                              Array<T, 1>& CellCenterDistance_y_)
  {
    const T earth_radius = 6371229.;
    const T pi(3.14159265358979323846264);

    for (int i = 0; i < this->Nx - 1; i++)
      CellCenterDistance_x_(i) = earth_radius * Delta_x_ * pi / 180.;

    Array<T, 1> GridY_interf_trans(this->Ny + 1), GridY_trans(this->Ny);
    for (int j = 0; j < this->Ny + 1; j++)
      GridY_interf_trans(j)
        = earth_radius * sin(GridY_interf_(j) * pi / 180.);
    for (int j = 0; j < this->Ny; j++)
      GridY_trans(j)
        = (GridY_interf_trans(j) + GridY_interf_trans(j + 1)) / 2.;

    for (int j = 0; j < this->Ny - 1; j++)
      CellCenterDistance_y_(j) = GridY_trans(j + 1) - GridY_trans(j);
  }


  //! Computes concentration at the point of coordinates (x,y,z).
  /*!
    \param s index of species.
    \param z height.
    \param lat latitude.
    \param lon longitude.
  */
  template<class T, class ClassParticle>
  T LagrangianTransport<T, ClassParticle>
  ::GetConcentration(int s, T z, T lat, T lon)
  {
    T concentration = 0.;

    typename list<ClassParticle>::iterator current_particle;
    for (current_particle = particles_list.begin();
         current_particle != particles_list.end(); current_particle++)
      {
        concentration +=
          current_particle->GetConcentrationContributionOnPoint
          (s, z, lat, lon);
      }
    return concentration;
  }


  //! Computes concentrations on the whole grid.
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeConcentration()
  {
    int s, k, j, i;
    for (s = 0; s < this->GetNs(); s++)
      for (k = 0; k < this->GetNz(); k++)
        for (j = 0; j < this->GetNy(); j++)
          for (i = 0; i < this->GetNx(); i++)
            this->Concentration(s, k, j, i) =
              GetConcentration(s, GetConcentration()[1](k),
                               GetConcentration()[2](j),
                               GetConcentration()[3](i));
  }


  //! Computes concentrations for a list of species and levels.
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeConcentration(const vector<int>& species_index,
                         const vector<int>& levels)
  {
    unsigned int s, k;
    int j, i;
    for (s = 0; s < species_index.size(); s++)
      for (k = 0; k < levels.size(); k++)
        for (j = 0; j < this->GetNy(); j++)
          for (i = 0; i < this->GetNx(); i++)
            this->Concentration(species_index[s], levels[k], j, i) =
              GetConcentration(s, GetConcentration()[1](levels[k]),
                               GetConcentration()[2](j),
                               GetConcentration()[3](i));
  }


  /*! \brief Transforms the zonal wind to ease the numerical integration of
    advection. */
  /*! Formula: ZonalWind = ZonalWind / cos(latitude).
    \param ZonalWind (input/output) zonal wind.
    \note Coordinates associated with ZonalWind must be in degrees.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::TransformZonalWind(Data<T, 3>& ZonalWind)
  {
    int i, j, k;

    const T pi(3.14159265358979323846264);
    const T ratio = pi / 180.;

    int Nx = ZonalWind.GetLength(2);
    int Ny = ZonalWind.GetLength(1);
    int Nz = ZonalWind.GetLength(0);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          ZonalWind(k, j, i) /= cos(ZonalWind[1].Value(k, j, i)
                                    * ratio);
  }


  /*! \brief Transforms the meridional wind to ease the numerical integration
    of advection. */
  /*! Formula: MeridionalWind = MeridionalWind * cos(latitude).
    \param MeridionalWind (input/output) meridional wind.
    \note Coordinates associated with MeridionalWind must be in degrees.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::TransformMeridionalWind(Data<T, 3>& MeridionalWind)
  {

    int i, j, k;

    const T pi(3.14159265358979323846264);
    const T ratio = pi / 180.;

    int Nx = MeridionalWind.GetLength(2);
    int Ny = MeridionalWind.GetLength(1);
    int Nz = MeridionalWind.GetLength(0);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          MeridionalWind(k, j, i) *= cos(MeridionalWind[1].Value(k, j, i)
                                         * ratio);
  }


  /*! \brief Transforms the zonal diffusion coefficients to ease the numerical
    integration of diffusion. */
  /*!
    \param GridY_interf_ interfaces coordinates (in degrees) along y.
    \param ZonalDiffusion_ (input/output) zonal diffusion coefficients.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::TransformZonalDiffusion(Array<T, 1>& GridY_interf_,
                            Data<T, 3>& ZonalDiffusion_)
  {
    int i, j, k;

    const T earth_radius = 6371229.;
    const T pi(3.14159265358979323846264);

    int Nx = ZonalDiffusion_.GetLength(2) - 1;
    int Ny = ZonalDiffusion_.GetLength(1);
    int Nz = ZonalDiffusion_.GetLength(0);

    Array<T, 1> GridY_interf_trans(Ny + 1), GridY_trans(Ny);
    for (int j = 0; j < Ny + 1; j++)
      GridY_interf_trans(j)
        = earth_radius * sin(GridY_interf_(j) * pi / 180.);
    for (int j = 0; j < Ny; j++)
      GridY_trans(j)
        = (GridY_interf_trans(j) + GridY_interf_trans(j + 1)) / 2.;

    Array<T, 1> factor(Ny);
    for (int j = 0; j < Ny; j++)
      factor(j) = 1. - GridY_trans(j) * GridY_trans(j)
        / (earth_radius * earth_radius);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx + 1; i++)
          ZonalDiffusion_(k, j, i) /= factor(j);
  }


  /*! \brief Transforms the meridional diffusion coefficients to ease the
    numerical integration of diffusion. */
  /*!
    \param MeridionalDiffusion_ (input/output) meridional diffusion
    coefficients.
    \note Coordinates associated with MeridionalDiffusion_ must be in degrees.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::TransformMeridionalDiffusion(Data<T, 3>& MeridionalDiffusion_)
  {
    int i, j, k;

    const T pi(3.14159265358979323846264);
    const T ratio = pi / 180.;

    int Nx = MeridionalDiffusion_.GetLength(2);
    int Ny = MeridionalDiffusion_.GetLength(1) - 1;
    int Nz = MeridionalDiffusion_.GetLength(0);

    Array<T, 1> factor(Ny + 1);
    for (int j = 0; j < Ny + 1; j++)
      factor(j) = cos(ratio * MeridionalDiffusion_[1].Value(0, j, 0))
        * cos(ratio * MeridionalDiffusion_[1].Value(0, j, 0));

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny + 1; j++)
        for (i = 0; i < Nx; i++)
          MeridionalDiffusion_(k, j, i) *= factor(j);
  }


  //! Computes the vertical wind.
  /*! Computes the vertical wind so that div(wind) = 0.
    \param CellWidth_x_ cells widths in meters along x.
    \param CellWidth_y_ cells widths in meters along y.
    \param CellWidth_z_ cells widths in meters along z.
    \param ZonalWind_ zonal wind.
    \param MeridionalWind_ meridonal wind.
    \param VerticalWind_ (output) vertical wind.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeVerticalWind(Array<T, 1>& CellWidth_x_, Array<T, 1>& CellWidth_y_,
                        Array<T, 1>& CellWidth_z_,
                        Data<T, 3>& ZonalWind_, Data<T, 3>& MeridionalWind_,
                        Data<T, 3>& VerticalWind_)
  {
    int k, j, i;
    int Nz = ZonalWind_.GetLength(0);
    int Ny = ZonalWind_.GetLength(1);
    int Nx = MeridionalWind_.GetLength(2);

    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        VerticalWind_(0, j, i) = 0.;

    for (k = 1; k < Nz + 1; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          VerticalWind_(k, j, i) = VerticalWind_(k - 1, j, i)
            - (ZonalWind_(k - 1, j, i + 1) - ZonalWind_(k - 1, j, i))
            * CellWidth_z_(k - 1) / CellWidth_x_(i)
            - (MeridionalWind_(k - 1, j + 1, i) - MeridionalWind_(k - 1, j, i))
            * CellWidth_z_(k - 1) / CellWidth_y_(j);
  }


  //! Computes the vertical wind.
  /*! Computes the vertical wind so that div(air_density * wind) = 0.
    \param CellWidth_x_ cells widths in meters along x.
    \param CellWidth_y_ cells widths in meters along y.
    \param CellWidth_z_ cells widths in meters along z.
    \param AirDensity_interf_x_ air density at cells interfaces along x.
    \param ZonalWind_ zonal wind.
    \param AirDensity_interf_y_ air density at cells interfaces along y.
    \param MeridionalWind_ meridonal wind.
    \param AirDensity_interf_y_ air density at cells interfaces along z.
    \param VerticalWind_ (output) vertical wind.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeVerticalWind(Array<T, 1>& CellWidth_x_, Array<T, 1>& CellWidth_y_,
                        Array<T, 1>& CellWidth_z_,
                        Data<T, 3>& AirDensity_interf_x_,
                        Data<T, 3>& ZonalWind_,
                        Data<T, 3>& AirDensity_interf_y_,
                        Data<T, 3>& MeridionalWind_,
                        Data<T, 3>& AirDensity_interf_z_,
                        Data<T, 3>& VerticalWind_)
  {
    int k, j, i;
    int Nz = ZonalWind_.GetLength(0);
    int Ny = ZonalWind_.GetLength(1);
    int Nx = MeridionalWind_.GetLength(2);

    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        VerticalWind_(0, j, i) = 0.;

    for (k = 1; k < Nz + 1; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            VerticalWind_(k, j, i)
              = VerticalWind_(k - 1, j, i) * AirDensity_interf_z_(k - 1, j, i)
              - (ZonalWind_(k - 1, j, i + 1) * AirDensity_interf_x_(k - 1, j, i + 1)
                 - ZonalWind_(k - 1, j, i) * AirDensity_interf_x_(k - 1, j, i))
              * CellWidth_z_(k - 1) / CellWidth_x_(i)
              - (MeridionalWind_(k - 1, j + 1, i)
                 * AirDensity_interf_y_(k - 1, j + 1, i)
                 - MeridionalWind_(k - 1, j, i)
                 * AirDensity_interf_y_(k - 1, j, i))
              * CellWidth_z_(k - 1) / CellWidth_y_(j);
            VerticalWind_(k, j, i) /= AirDensity_interf_z_(k, j, i);
          }
  }


  //! Computes air density on the basis of temperature and pressure.
  /*! Formula: AirDensity = Pressure / (287. * Temperature).
    \param Temperature_ temperature.
    \param Pressure_ pressure.
    \param AirDensity_ (output) air density.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::ComputeAirDensity(Data<T, 3>& Temperature_, Data<T, 3>& Pressure_,
                      Data<T, 3>& AirDensity_)
  {
    const T r = 287.;

    int k, j, i;

    int Nz = AirDensity_.GetLength(0);
    int Ny = AirDensity_.GetLength(1);
    int Nx = AirDensity_.GetLength(2);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          AirDensity_(k, j, i) = Pressure_(k, j, i)
            / (r * Temperature_(k, j, i));
  }


  //! Interpolates data on cells interfaces along z.
  /*!
    \param Data_ input data.
    \param Data_interf_z_ (output) data interpolated on interfaces along z.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_z(Data<T, 3>& Data_, Data<T, 3>& Data_interf_z_)
  {
    int k, k_in, j, i;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    T weight_0, weight_1;

    for (k = 0; k < Nz + 1; k++)
      {
        if (k == 0)
          k_in = 1;
        else if (k == Nz)
          k_in = Nz - 1;
        else
          k_in = k;

        weight_1 = (Data_interf_z_[0].Value(k, 0, 0)
                    - Data_[0].Value(k_in - 1, 0, 0))
          / (Data_[0].Value(k_in, 0, 0) - Data_[0].Value(k_in - 1, 0, 0));
        weight_0 = 1. - weight_1;

        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Data_interf_z_(k, j, i) =  weight_0 * Data_(k_in - 1, j, i)
              + weight_1 * Data_(k_in, j, i);
      }
  }


  //! Interpolates data on cells interfaces along z.
  /*!
    \param CellCenterDistance_z_ distances between cells centers along z.
    \param CellWidth_z_ cells widths along z.
    \param Data_ input data on cell centers.
    \param Data_interf_z_ (output) data interpolated on interfaces along z.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_z(Array<T, 1>& CellCenterDistance_z_,
                           Array<T, 1>& CellWidth_z_,
                           Data<T, 3>& Data_, Data<T, 3>& Data_interf_z_)
  {
    int k, j, i;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    T weight;

    weight = CellWidth_z_(0) / (2. * CellCenterDistance_z_(0));
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        Data_interf_z_(0, j, i) = (1. + weight) * Data_(0, j, i)
          - weight * Data_(1, j, i);

    for (k = 1; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          Data_interf_z_(k, j, i) =
            0.5 * (Data_(k - 1, j, i) * CellWidth_z_(k)
                   + Data_(k, j, i) * CellWidth_z_(k - 1))
            / CellCenterDistance_z_(k - 1);

    weight = CellWidth_z_(Nz - 1) / (2. * CellCenterDistance_z_(Nz - 2));
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        Data_interf_z_(Nz, j, i) = (1. + weight) * Data_(Nz - 1, j, i)
          - weight * Data_(Nz - 2, j, i);
  }


  //! Interpolates data on cells interfaces along y.
  /*!
    \param Data_ input data.
    \param Data_interf_y_ (output) data interpolated on interfaces along y.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_y(Data<T, 3>& Data_, Data<T, 3>& Data_interf_y_)
  {
    int k, j, j_in, i;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    // Assumption...
    T weight_0, weight_1;

    for (j = 0; j < Ny + 1; j++)
      {
        if (j == 0)
          j_in = 1;
        else if (j == Ny)
          j_in = Ny - 1;
        else
          j_in = j;

        weight_1 = (Data_interf_y_[1].Value(0, j, 0)
                    - Data_[1].Value(0, j_in - 1, 0))
          / (Data_[1].Value(0, j_in, 0) - Data_[1].Value(0, j_in - 1, 0));
        weight_0 = 1. - weight_1;

        for (k = 0; k < Nz; k++)
          for (i = 0; i < Nx; i++)
            Data_interf_y_(k, j, i) =  weight_0 * Data_(k, j_in - 1, i)
              + weight_1 * Data_(k, j_in, i);
      }
  }


  //! Interpolates data on cells interfaces along y.
  /*!
    \param CellCenterDistance_y_ distances between cells centers along y.
    \param CellWidth_y_ cells widths along y.
    \param Data_ input data on cell centers.
    \param Data_interf_y_ (output) data interpolated on interfaces along y.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_y(Array<T, 1>& CellCenterDistance_y_,
                           Array<T, 1>& CellWidth_y_,
                           Data<T, 3>& Data_, Data<T, 3>& Data_interf_y_)
  {
    int k, j, i;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    T weight;

    weight = CellWidth_y_(0) / (2. * CellCenterDistance_y_(0));
    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        Data_interf_y_(k, 0, i) = (1. + weight) * Data_(k, 0, i)
          - weight * Data_(k, 1, i);

    for (j = 1; j < Ny; j++)
      for (k = 0; k < Nz; k++)
        for (i = 0; i < Nx; i++)
          Data_interf_y_(k, j, i) =
            0.5 * (Data_(k, j - 1, i) * CellWidth_y_(j)
                   + Data_(k, j, i) * CellWidth_y_(j - 1))
            / CellCenterDistance_y_(j - 1);

    weight = CellWidth_y_(Ny - 1) / (2. * CellCenterDistance_y_(Ny - 2));
    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        Data_interf_y_(k, Ny, i) = (1. + weight) * Data_(k, Ny - 1, i)
          - weight * Data_(k, Ny - 2, i);
  }


  //! Interpolates data on cells interfaces along x.
  /*!
    \param Data_ input data.
    \param Data_interf_x_ (output) data interpolated on interfaces along x.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_x(Data<T, 3>& Data_, Data<T, 3>& Data_interf_x_)
  {
    int k, j, i, i_in;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    // Assumption...
    T weight_0, weight_1;

    for (i = 0; i < Nx + 1; i++)
      {
        if (i == 0)
          i_in = 1;
        else if (i == Nx)
          i_in = Nx - 1;
        else
          i_in = i;

        weight_1 = (Data_interf_x_[2].Value(0, 0, i)
                    - Data_[2].Value(0, 0, i_in - 1))
          / (Data_[2].Value(0, 0, i_in) - Data_[2].Value(0, 0, i_in - 1));
        weight_0 = 1. - weight_1;

        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            Data_interf_x_(k, j, i) =  weight_0 * Data_(k, j, i_in - 1)
              + weight_1 * Data_(k, j, i_in);
      }
  }


  //! Interpolates data on cells interfaces along x.
  /*!
    \param CellCenterDistance_x_ distances between cells centers along x.
    \param CellWidth_x_ cells widths along x.
    \param Data_ input data on cell centers.
    \param Data_interf_x_ (output) data interpolated on interfaces along x.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::InterpolateInterface_x(Array<T, 1>& CellCenterDistance_x_,
                           Array<T, 1>& CellWidth_x_,
                           Data<T, 3>& Data_, Data<T, 3>& Data_interf_x_)
  {
    int k, j, i;

    int Nz = Data_.GetLength(0);
    int Ny = Data_.GetLength(1);
    int Nx = Data_.GetLength(2);

    T weight;

    weight = CellWidth_x_(0) / (2. * CellCenterDistance_x_(0));
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        Data_interf_x_(k, j, 0) = (1. + weight) * Data_(k, j, 0)
          - weight * Data_(k, j, 1);

    for (i = 1; i < Nx; i++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          Data_interf_x_(k, j, i) =
            0.5 * (Data_(k, j, i - 1) * CellWidth_x_(i)
                   + Data_(k, j, i) * CellWidth_x_(i - 1))
            / CellCenterDistance_x_(i - 1);

    weight = CellWidth_x_(Nx - 1) / (2. * CellCenterDistance_x_(Nx - 2));
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        Data_interf_x_(k, j, Nx) = (1. + weight) * Data_(k, j, Nx - 1)
          - weight * Data_(k, j, Nx - 2);
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Moves model input-data to the current date.
  /*! This method prepares the model for a time integration from the current
    date. It reads input data to be read before InitStep and Forward.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>::InitAllData()
  {

    /*** Pressure, temperature and air density ***/

    if (this->option_manage["temperature"])
      this->InitData("meteo", "Temperature", FileTemperature_i,
                     FileTemperature_f, this->current_date, Temperature_f);
    if (this->option_manage["pressure"])
      this->InitData("meteo", "Pressure", FilePressure_i,
                     FilePressure_f, this->current_date, Pressure_f);

    if (this->option_process["with_air_density"])
      {
        ComputeAirDensity(Temperature_f, Pressure_f, AirDensity_f);

        InterpolateInterface_z(CellCenterDistance_z, CellWidth_z,
                               AirDensity_f, AirDensity_interf_z_f);
        InterpolateInterface_y(CellCenterDistance_y, CellWidth_y,
                               AirDensity_f, AirDensity_interf_y_f);
        InterpolateInterface_x(CellCenterDistance_x, CellWidth_x,
                               AirDensity_f, AirDensity_interf_x_f);
      }
    else
      {
        AirDensity_i.Fill(1.);
        AirDensity_f.Fill(1.);

        AirDensity_interf_z_i.Fill(1.);
        AirDensity_interf_y_i.Fill(1.);
        AirDensity_interf_x_i.Fill(1.);

        AirDensity_interf_z_f.Fill(1.);
        AirDensity_interf_y_f.Fill(1.);
        AirDensity_interf_x_f.Fill(1.);
      }

    /*** Winds ***/

    if (this->option_manage["horizontal_wind"])
      {
        this->InitData("meteo", "ZonalWind", FileZonalWind_i,
                       FileZonalWind_f, this->current_date, ZonalWind_i);
        TransformZonalWind(ZonalWind_i);

        this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
                       FileMeridionalWind_f, this->current_date,
                       MeridionalWind_i);
        TransformMeridionalWind(MeridionalWind_i);
      }

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      this->InitData("meteo", "VerticalDiffusion",
                     FileVerticalDiffusionCoefficient_i,
                     FileVerticalDiffusionCoefficient_f,
                     this->current_date, VerticalDiffusionCoefficient_f);

    // Computes rho * Kz.
    if (this->option_process["with_air_density"]
        && this->option_manage["vertical_diffusion"])
      VerticalDiffusionCoefficient_f.GetArray() =
        AirDensity_interf_z_f.GetArray()
        * VerticalDiffusionCoefficient_f.GetArray();

    if (this->option_manage["horizontal_diffusion"])
      {
        ZonalDiffusionCoefficient_i.Fill(horizontal_diffusion);
        MeridionalDiffusionCoefficient_i.Fill(horizontal_diffusion);
      }

    if (this->option_manage["horizontal_diffusion"])
      {
        TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
                                ZonalDiffusionCoefficient_i);
        TransformMeridionalDiffusion(MeridionalDiffusionCoefficient_i);
      }
  }


  //! Computes indices of the cell containing a given point in Eulerian grid.
  /*!
    \param lon longitude of the point (degrees).
    \param lat latitude of the point (degrees).
    \param height of the point (meters).
    \param k (output) cell index along z.
    \param j (output) cell index along y.
    \param i (output) cell index along x.
  */
  template<class T, class ClassParticle>
  void LagrangianTransport<T, ClassParticle>
  ::GetCellIndices(T lon, T lat, T height,
                   int& index_z, int& index_y, int& index_x)
  {
    index_x = int((lon - this->x_min + this->Delta_x / 2.) / this->Delta_x);
    index_y = int((lat - this->y_min + this->Delta_y / 2.) / this->Delta_y);
    for (int k = 0; k < this->Nz; k++)
      {
        index_z = k;
        if (this->GridZ3D_interf(k + 1) > height)
          break;
      }
  }


  ////////////////////////////////////
  // RELATED PREDICATE CLASS METHOD //
  ////////////////////////////////////


  //! Constructor for boundary coordinates in the curvilinear coordinates
  //system.
  /*! This constructor should be extended if the related ClassParticle
    do not use the same curvilinear coordinates system.
    \param Model Model object.
  */
  template<class T>
  template<typename ClassModel>
  IsOutsideOfDomain<T>
  ::IsOutsideOfDomain(ClassModel* Model)
  {
    const T earth_radius(6371229.);
    const T pi(3.1416);

    x_min_ = Model->GetX_min() *  earth_radius * pi / 180.;
    x_max_ = Model->GetX_max() * earth_radius * pi / 180.;
    y_min_ = sin(Model->GetY_min() * pi / 180.) *  earth_radius;
    y_max_ = sin(Model->GetY_max() * pi / 180.) * earth_radius;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_LAGRANGIANTRANSPORT_CXX
#endif
