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


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DTRANSPORT_CXX


#include "Polair3DTransport.hxx"

#include "AtmoData.hxx"
#include "BaseModel.cxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Default constructor.
  template<class T, class ClassAdvection, class ClassDiffusion>
  Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::Polair3DTransport(): BaseModel<T>()
  {
  }


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::Polair3DTransport(string config_file)
  {
    Construct(config_file);

    /*** Input parameters ***/

    this->RegisterParameter("Rain_i");
    this->RegisterParameter("Rain_f");
    this->RegisterParameter("CloudBaseHeight_i");
    this->RegisterParameter("CloudBaseHeight_f");
    this->RegisterParameter("CloudTopHeight_i");
    this->RegisterParameter("CloudTopHeight_f");
    this->RegisterParameter("Temperature_i");
    this->RegisterParameter("Temperature_f");
    this->RegisterParameter("Pressure_i");
    this->RegisterParameter("Pressure_f");
    this->RegisterParameter("AirDensity_i");
    this->RegisterParameter("AirDensity_f");
    this->RegisterParameter("MeridionalWind_i");
    this->RegisterParameter("ZonalWind_i");
    this->RegisterParameter("SpecificHumidity_i");
    this->RegisterParameter("SpecificHumidity_f");
    this->RegisterParameter("VerticalDiffusionCoefficient_i");
    this->RegisterParameter("VerticalDiffusionCoefficient_f");
    this->RegisterParameter("BoundaryCondition_z_i");
    this->RegisterParameter("DepositionVelocity_i");
    this->RegisterParameter("DepositionVelocity_f");
    this->RegisterParameter("SurfaceEmission_i");
    this->RegisterParameter("SurfaceEmission_f");
    this->RegisterParameter("BoundaryCondition_y_i");
    this->RegisterParameter("BoundaryCondition_x_i");
    this->RegisterParameter("VolumeEmission_i");
    this->RegisterParameter("VolumeEmission_f");
    this->RegisterParameter("ScavengingBelowCloudCoefficient_i");
    this->RegisterParameter("ScavengingBelowCloudCoefficient_f");
    this->RegisterParameter("ScavengingInCloudCoefficient_i");
    this->RegisterParameter("ScavengingInCloudCoefficient_f");
  }


  //! Constructs the model.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::Construct(string config_file)
  {
    BaseModel<T>::Construct(config_file);

    /*** Managed data ***/

    this->option_manage["initial_condition"] = true;
    this->option_manage["temperature"] = true;
    this->option_manage["pressure"] = true;
    this->option_manage["horizontal_wind"] = true;
    this->option_manage["vertical_wind"] = true;
    this->option_manage["rain"] = true;
    this->option_manage["specific_humidity"] = true;
    this->option_manage["cloud_base_height"] = true;
    this->option_manage["cloud_top_height"] = true;
    this->option_manage["horizontal_diffusion"] = true;
    this->option_manage["vertical_diffusion"] = true;
    this->option_manage["boundary_condition"] = true;
    this->option_manage["deposition_velocity"] = true;
    this->option_manage["surface_emission"] = true;
    this->option_manage["additional_surface_emission"] = true;
    this->option_manage["volume_emission"] = true;
    this->option_manage["scavenging_below_cloud_coefficient"] = true;
    this->option_manage["scavenging_in_cloud_coefficient"] = true;

    /*** Pointers to 2D data ***/

    this->Register(Rain_i, Rain_f, "Rain");
    this->Register(CloudBaseHeight_i, CloudBaseHeight_f, "CloudBaseHeight");
    this->Register(CloudTopHeight_i, CloudTopHeight_f, "CloudTopHeight");

    /*** Pointers to 3D data ***/

    this->Register(Temperature_i, Temperature_f, "Temperature");
    this->Register(Pressure_i, Pressure_f, "Pressure");
    this->Register(AirDensity_i, AirDensity_f, "AirDensity");
    this->Register(VerticalWind_i, "VerticalWind");
    this->Register(MeridionalWind_i, "MeridionalWind");
    this->Register(ZonalWind_i, "ZonalWind");
    this->Register(SpecificHumidity_i, SpecificHumidity_f, "SpecificHumidity");
    this->Register(VerticalDiffusionCoefficient_i,
                   VerticalDiffusionCoefficient_f,
                   "VerticalDiffusionCoefficient");
    this->Register(MeridionalDiffusionCoefficient_i,
                   "MeridionalDiffusionCoefficient");
    this->Register(ZonalDiffusionCoefficient_i, "ZonalDiffusionCoefficient");
    this->Register(BoundaryCondition_z_i, "BoundaryCondition_z");
    this->Register(DepositionVelocity_i, DepositionVelocity_f,
                   "DepositionVelocity");
    this->Register(SurfaceEmission_i, SurfaceEmission_f, "SurfaceEmission");
    this->Register(DryDepositionFlux, "DryDepositionFlux");
    this->Register(WetDepositionFlux, "WetDepositionFlux");
    this->Register(InCloudWetDepositionFlux , "InCloudWetDepositionFlux");

    /*** Pointers to 4D data ***/

    this->Register(BoundaryCondition_y_i, "BoundaryCondition_y");
    this->Register(BoundaryCondition_x_i, "BoundaryCondition_x");
    this->Register(VolumeEmission_i, VolumeEmission_f, "VolumeEmission");
    this->Register(ScavengingBelowCloudCoefficient_i,
                   ScavengingBelowCloudCoefficient_f,
                   "ScavengingBelowCloudCoefficient");
    this->Register(ScavengingInCloudCoefficient_i,
                   ScavengingInCloudCoefficient_f,
                   "ScavengingInCloudCoefficient");

    /*** Fields species lists ***/

    this->field_species["InitialCondition"] = &species_list_ic;
    this->field_species["BoundaryCondition_z"] = &species_list_bc;
    this->field_species["BoundaryCondition_z_i"] = &species_list_bc;
    this->field_species["BoundaryCondition_y"] = &species_list_bc;
    this->field_species["BoundaryCondition_y_i"] = &species_list_bc;
    this->field_species["BoundaryCondition_x"] = &species_list_bc;
    this->field_species["BoundaryCondition_x_i"] = &species_list_bc;
    this->field_species["DepositionVelocity"] = &species_list_dep;
    this->field_species["DepositionVelocity_i"] = &species_list_dep;
    this->field_species["DepositionVelocity_f"] = &species_list_dep;
    this->field_species["ScavengingBelowCloudCoefficient"]
      = &species_list_scav;
    this->field_species["ScavengingBelowCloudCoefficient_i"]
      = &species_list_scav;
    this->field_species["ScavengingBelowCloudCoefficient_f"]
      = &species_list_scav;
    this->field_species["ScavengingInCloudCoefficient"]
      = &species_list_scav;
    this->field_species["ScavengingInCloudCoefficient_i"]
      = &species_list_scav;
    this->field_species["ScavengingInCloudCoefficient_f"]
      = &species_list_scav;
    this->field_species["SurfaceEmission"] = &species_list_surf_emis;
    this->field_species["SurfaceEmission_i"] = &species_list_surf_emis;
    this->field_species["SurfaceEmission_f"] = &species_list_surf_emis;
    this->field_species["VolumeEmission"] = &species_list_vol_emis;
    this->field_species["VolumeEmission_i"] = &species_list_vol_emis;
    this->field_species["VolumeEmission_f"] = &species_list_vol_emis;
  }


  //! Destructor.
  template<class T, class ClassAdvection, class ClassDiffusion>
  Polair3DTransport<T, ClassAdvection, ClassDiffusion>::~Polair3DTransport()
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ReadConfiguration()
  {

    BaseModel<T>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[domain]");
    this->config.PeekValue("Cartesian", option_cartesian);

    this->config.SetSection("[options]");
    this->config.PeekValue("With_advection",
                           this->option_process["with_advection"]);
    this->config.PeekValue("With_diffusion",
                           this->option_process["with_diffusion"]);
    this->config.PeekValue("With_air_density",
                           this->option_process["with_air_density"]);
    this->config.PeekValue("With_initial_condition",
                           this->option_process["with_initial_condition"]);
    this->config.PeekValue("With_boundary_condition",
                           this->option_process["with_boundary_condition"]);
    this->config.PeekValue("With_deposition",
                           this->option_process["with_deposition"]);
    this->config.PeekValue("With_point_emission",
                           this->option_process["with_point_emission"]);
    this->config.PeekValue("With_surface_emission",
                           this->option_process["with_surface_emission"]);
    this->config.PeekValue("With_additional_surface_emission",
                           this->option_process
                           ["with_additional_surface_emission"]);
    this->config.PeekValue("With_volume_emission",
                           this->option_process["with_volume_emission"]);
    this->config.PeekValue("Collect_dry_flux",
                           this->option_process["collect_dry_flux"]);
    this->config.PeekValue("Collect_wet_flux",
                           this->option_process["collect_wet_flux"]);

    this->config.PeekValue("Scavenging_below_cloud_model",
                           "none|constant|belot|microphysical|microphysical-ph",
                           this->scavenging_below_cloud_model);
    this->config.PeekValue("Scavenging_in_cloud_model",
                           "none|belot|pudykiewicz|aqueous",
                           this->scavenging_in_cloud_model);
    this->option_process["with_transport_in_cloud_scavenging"] =
      this->scavenging_in_cloud_model == "belot"
      || this->scavenging_in_cloud_model == "pudykiewicz";

    /*** Disables management of useless fields ***/

    // Pressure and temperature are only needed for air density
    // or two models of scavenging.
    this->option_manage["temperature"] = this->option_manage["temperature"]
      && (this->option_process["with_air_density"]
          || this->scavenging_below_cloud_model == "microphysical"
          || this->scavenging_in_cloud_model == "pudykiewicz");
    this->option_manage["pressure"] = this->option_manage["pressure"]
      && (this->option_process["with_air_density"]
          || this->scavenging_below_cloud_model == "microphysical"
          || this->scavenging_in_cloud_model == "pudykiewicz");

    // Winds are only needed for advection.
    this->option_manage["horizontal_wind"]
      = this->option_manage["horizontal_wind"]
      && this->option_process["with_advection"];
    this->option_manage["vertical_wind"]
      = this->option_manage["vertical_wind"]
      && this->option_manage["horizontal_wind"];

    // Diffusion coefficients are only needed for diffusion.
    this->option_manage["horizontal_diffusion"]
      = this->option_manage["horizontal_diffusion"]
      && this->option_process["with_diffusion"];
    this->option_manage["vertical_diffusion"]
      = this->option_manage["vertical_diffusion"]
      && this->option_process["with_diffusion"];

    // Rain and cloud height are only needed for scavenging.
    this->option_manage["rain"] = this->option_manage["rain"]
      && (this->scavenging_below_cloud_model != "none"
          || this->scavenging_in_cloud_model != "none");
    this->option_manage["cloud_base_height"]
      = this->option_manage["cloud_base_height"]
      && (this->scavenging_below_cloud_model != "none"
          || this->scavenging_in_cloud_model != "none");
    this->option_manage["cloud_top_height"] =
      this->option_manage["cloud_top_height"]
      && this->option_process["with_transport_in_cloud_scavenging"];
    // Specific humidity is only needed for scavenging in passive transport.
    this->option_manage["specific_humidity"]
      = this->option_manage["specific_humidity"]
      && this->scavenging_in_cloud_model == "pudykiewicz";

    /*** Input files ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Meteorological files.
    this->input_files["meteo"].Read(data_description_file, "meteo");

    // Initial conditions files.
    if (this->option_process["with_initial_condition"])
      this->input_files["initial_condition"].ReadFiles(data_description_file,
                                                       "initial_condition");
    else
      this->input_files["initial_condition"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["initial_condition"].Begin();
         i != this->input_files["initial_condition"].End(); i++)
      species_list_ic.push_back(i->first);
    Ns_ic = int(species_list_ic.size());

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

    // Scavenging.
    if (this->scavenging_below_cloud_model != "none"
        || this->option_process["with_transport_in_cloud_scavenging"])
      this->input_files["scavenging_coefficient"]
        .ReadFields(data_description_file, "scavenging");
    else
      this->input_files["scavenging_coefficient"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["scavenging_coefficient"].Begin();
         i != this->input_files["scavenging_coefficient"].End(); i++)
      species_list_scav.push_back(i->first);
    Ns_scav = int(species_list_scav.size());

    // Point emissions.
    if (this->option_process["with_point_emission"])
      {
        string point_emission_file;
        data_description_stream.SetSection("[point_emission]");
        data_description_stream.PeekValue("file", point_emission_file);
        this->PointEmissionManager = new BasePointEmission<T>();
        this->PointEmissionManager->Init(point_emission_file,
                                         this->species_list);
      }

    // Surface emission files.
    if (this->option_process["with_surface_emission"])
      this->input_files["surface_emission"].Read(data_description_file,
                                                 "surface_emission");
    else
      this->input_files["surface_emission"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["surface_emission"].Begin();
         i != this->input_files["surface_emission"].End(); i++)
      species_list_surf_emis.push_back(i->first);
    Ns_surf_emis = int(species_list_surf_emis.size());

    // Additional surface emission files.
    if (this->option_process["with_additional_surface_emission"])
      this->input_files["additional_surface_emission"]
        .Read(data_description_file, "additional_surface_emission");
    else
      this->input_files["additional_surface_emission"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["additional_surface_emission"].Begin();
         i != this->input_files["additional_surface_emission"].End(); i++)
      species_list_add_surf_emis.push_back(i->first);
    Ns_add_surf_emis = int(species_list_add_surf_emis.size());

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
        data_description_stream.PeekValue("Nz", "positive", Nz_vol_emis);
      }
    else
      Nz_vol_emis = 0;

    /*** Rest of the configuration ***/

    // Horizontal diffusion coefficient.
    this->config.PeekValue("Horizontal_diffusion", "positive",
                           horizontal_diffusion);
    // Isotropy.
    this->config.PeekValue("Isotropic_diffusion", option_isotropic_diffusion);

    /*** Species data ***/

    ConfigStream species_data_stream(this->GetSpeciesFile());
    string species;

    if (this->scavenging_below_cloud_model == "microphysical")
      {
        // Henry constants in mol / L / atm.
        species_data_stream.SetSection("[henry]");
        while (!species_data_stream.IsEmpty())
          {
            species = species_data_stream.GetElement();
            species_data_stream.GetNumber(henry_constant[species]);
          }

        // Gas phase diffusivity constants in cm^2 / s.
        species_data_stream.SetSection("[diffusivity]");
        while (!species_data_stream.IsEmpty())
          {
            species = species_data_stream.GetElement();
            species_data_stream.GetNumber(gas_phase_diffusivity[species]);
          }
      }
    else if (this->scavenging_below_cloud_model == "belot")
      {
        // Initializes a and b coefficients.
        species_data_stream.SetSection("[belot_below_cloud]");
        species = species_data_stream.PeekElement();
        if (species == "all")
          {
            T a, b;
            species_data_stream.GetNumber(a);
            species_data_stream.GetNumber(b);
            for (int i = 0; i < Ns_scav; i++)
              {
                belot_below_cloud_constant_a[ScavengingName(i)] = a;
                belot_below_cloud_constant_b[ScavengingName(i)] = b;
              }
          }
        else
          while (!species_data_stream.IsEmpty())
            {
              species = species_data_stream.GetElement();
              species_data_stream
                .GetNumber(belot_below_cloud_constant_a[species]);
              species_data_stream
                .GetNumber(belot_below_cloud_constant_b[species]);
            }
      }
    else if (this->scavenging_below_cloud_model == "constant")
      {
        species_data_stream.SetSection("[scavenging_coefficient]");
        // Initializes rain threshold.
        species_data_stream.GetValue("Scavenging_rain_threshold", "positive",
                                     scavenging_rain_threshold);
        // Constant scavenging coefficient in s^{-1}.
        species = species_data_stream.PeekElement();
        if (species == "all")
          {
            T scav_constant;
            species_data_stream.GetNumber(scav_constant);
            for (int i = 0; i < Ns_scav; i++)
              scavenging_constant[ScavengingName(i)] = scav_constant;
          }
        else
          while (!species_data_stream.IsEmpty())
            {
              species = species_data_stream.GetElement();
              if (!this->IsSpecies(species))
                throw string("ERROR! Species \"") + species
                  + string("\" found in section [scavenging_coefficient]")
                  + " is not a model species.";
              species_data_stream.GetNumber(scavenging_constant[species]);
            }
      }

    if (this->scavenging_in_cloud_model == "belot")
      {
        // Initializes a and b coefficients.
        species_data_stream.SetSection("[belot_in_cloud]");
        species = species_data_stream.PeekElement();
        if (species == "all")
          {
            T a, b;
            species_data_stream.GetNumber(a);
            species_data_stream.GetNumber(b);
            for (int i = 0; i < Ns_scav; i++)
              {
                belot_in_cloud_constant_a[ScavengingName(i)] = a;
                belot_in_cloud_constant_b[ScavengingName(i)] = b;
              }
          }
        else
          while (!species_data_stream.IsEmpty())
            {
              species = species_data_stream.GetElement();
              species_data_stream
                .GetNumber(belot_in_cloud_constant_a[species]);
              species_data_stream
                .GetNumber(belot_in_cloud_constant_b[species]);
            }
      }
  }


  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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

    if (this->option_manage["rain"]
        && this->input_files["meteo"]("Rain").empty())
      throw "Rain is needed but no input data file was provided.";
    if (this->option_manage["cloud_base_height"]
        && this->input_files["meteo"]("CloudBaseHeight").empty())
      throw "Cloud base height is needed but no input data file was "
        "provided.";
    if (this->option_manage["cloud_top_height"]
        && this->input_files["meteo"]("CloudTopHeight").empty())
      throw "Cloud top height is needed but no input data file was provided.";
    if (this->option_manage["specific_humidity"]
        && this->input_files["meteo"]("SpecificHumidity").empty())
      throw string("Specific humidity is needed but no input data file was")
        + " provided.";

    if (this->option_manage["vertical_diffusion"]
        && this->input_files["meteo"]("VerticalDiffusion").empty())
      throw string("Vertical diffusion is needed but no input data file")
        + " was provided.";

    if (!this->option_process["with_diffusion"])
      {
        if (this->option_process["with_deposition"])
          throw "Deposition velocities cannot be included without diffusion.";
        if (this->option_process["with_surface_emission"])
          throw "Surface emissions cannot be included without diffusion.";
      }

    if (!this->option_process["with_surface_emission"]
        && this->option_process["with_additional_surface_emission"])
      throw string("Additional surface emissions are added only ")
        + "if surface emissions are taken into account.";

    // Checks that all species with additional emissions are already in
    // species with surface emissions.
    if (this->option_process["with_additional_surface_emission"])
      for (int s = 0; s < Ns_add_surf_emis; s++)
        {
          string name = AdditionalSurfaceEmissionName(s);
          if (!HasSurfaceEmission(name))
            throw string("Species \"") + name
              + string("\" is in additional surface emissions (section")
              + string(" [additional_surface_emission]), but it is not in")
              + string(" surface emissions (section [surface_emission]).");
        }

    if (!this->option_process["with_deposition"]
        && this->option_process["collect_dry_flux"])
      throw string("Dry deposition fluxes cannot be collected") +
        " without deposition.";

    if (this->scavenging_below_cloud_model == "microphysical")
      {
        // Checks if Henry constant is well initialized for each scavenged
        // species.
        for (int i = 0; i < Ns_scav; i++)
          if (henry_constant.find(ScavengingName(i)) == henry_constant.end())
            throw string("ERROR! Henry constant for species \"")
              + ScavengingName(i)
              + string("\" not found in section [henry].");

        // Checks if gas phase diffusivity constant is well initialized for
        // all scavenged species.
        for (int i = 0; i < Ns_scav; i++)
          if (gas_phase_diffusivity.find(ScavengingName(i))
              == gas_phase_diffusivity.end())
            throw string("ERROR! Gas phase diffusivity constant for")
              + string(" species \"") + ScavengingName(i)
              + string("\" not found in section [diffusivity].");
      }
    else if (this->scavenging_below_cloud_model == "belot")
      {
        // Checks if the coefficients a and b are well initialized for each
        // scavenged species.
        for (int i = 0; i < Ns_scav; i++)
          if (belot_below_cloud_constant_a.find(ScavengingName(i))
              == belot_below_cloud_constant_a.end()
              || belot_below_cloud_constant_b.find(ScavengingName(i))
              == belot_below_cloud_constant_b.end())
            throw string("ERROR! Coefficient a or b for species \"")
              + ScavengingName(i)
              + string("\" not found in section [belot_below_cloud].");
      }
    else if (this->scavenging_below_cloud_model == "constant")
      {
        // Checks if the scavenging coefficient is well initialized for each
        // scavenged species.
        for (int i = 0; i < Ns_scav; i++)
          if (scavenging_constant.find(ScavengingName(i))
              == scavenging_constant.end())
            throw string("ERROR! Scavenging coefficient for species \"")
              + ScavengingName(i)
              + string("\" not found in section [scavenging_coefficient].");
      }
    else if ((this->scavenging_below_cloud_model == "none"
              || this->scavenging_in_cloud_model == "none")
             && this->option_process["collect_wet_flux"])
      throw "Wet deposition cannot be collected without scavenging.";

    if (this->scavenging_in_cloud_model == "belot")
      {
        // Checks if the coefficients a and b are well initialized for each
        // scavenged species.
        for (int i = 0; i < Ns_scav; i++)
          if (belot_in_cloud_constant_a.find(ScavengingName(i))
              == belot_in_cloud_constant_a.end()
              || belot_in_cloud_constant_b.find(ScavengingName(i))
              == belot_in_cloud_constant_b.end())
            throw string("ERROR! Coefficient a or b for species \"")
              + ScavengingName(i)
              + string("\" not found in section [belot_in_cloud].");
      }
  }


  //! Checks whether a species has initial conditions.
  /*!
    \param s species global index.
    \return True if the species has initial conditions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasInitialCondition(int s) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(),
                this->GetSpeciesName(s)) != species_list_ic.end();
  }


  //! Checks whether a species has initial conditions.
  /*!
    \param name species name.
    \return True if the species has initial conditions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasInitialCondition(string name) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(), name)
      != species_list_ic.end();
  }


  //! Returns the index in initial conditions of a given species.
  /*!
    \param s species global index.
    \return The species index in initial conditions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitialConditionIndex(int s) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(),
                this->GetSpeciesName(s)) - species_list_ic.begin();
  }


  //! Returns the index in initial conditions of a given species.
  /*!
    \param name species name.
    \return The species index in initial conditions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitialConditionIndex(string name) const
  {
    return find(species_list_ic.begin(), species_list_ic.end(), name)
      - species_list_ic.begin();
  }


  //! Returns the name of a species with initial conditions.
  /*!
    \param s species index in initial conditions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitialConditionName(int s) const
  {
    return species_list_ic.at(s);
  }


  //! Returns the name of a species with initial conditions.
  /*!
    \param s species index in initial conditions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitialConditionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(InitialConditionName(s));
  }


  //! Checks whether a species has boundary conditions.
  /*!
    \param s species global index.
    \return True if the species has boundary conditions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasBoundaryCondition(int s) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(),
                this->GetSpeciesName(s)) != species_list_bc.end();
  }


  //! Checks whether a species has boundary conditions.
  /*!
    \param name species name.
    \return True if the species has boundary conditions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasBoundaryCondition(string name) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(), name)
      != species_list_bc.end();
  }


  //! Returns the index in boundary conditions of a given species.
  /*!
    \param s species global index.
    \return The species index in boundary conditions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::BoundaryConditionIndex(int s) const
  {
    return find(species_list_bc.begin(), species_list_bc.end(),
                this->GetSpeciesName(s)) - species_list_bc.begin();
  }


  //! Returns the index in boundary conditions of a given species.
  /*!
    \param name species name.
    \return The species index in boundary conditions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::BoundaryConditionName(int s) const
  {
    return species_list_bc.at(s);
  }

  
  //! Returns the name of the model.
  /*!
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetModelName() const
  {
    return "Polair3D";
  }

  //! Returns the name of a species with boundary conditions.
  /*!
    \param s species index in boundary conditions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::BoundaryConditionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(BoundaryConditionName(s));
  }


  //! Checks whether a species has deposition velocities.
  /*!
    \param s species global index.
    \return True if the species has deposition velocities, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasDepositionVelocity(int s) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(),
                this->GetSpeciesName(s)) != species_list_dep.end();
  }


  //! Checks whether a species has deposition velocities.
  /*!
    \param name species name.
    \return True if the species has deposition velocities, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::DepositionVelocityIndex(int s) const
  {
    return find(species_list_dep.begin(), species_list_dep.end(),
                this->GetSpeciesName(s)) - species_list_dep.begin();
  }


  //! Returns the index in deposition velocities of a given species.
  /*!
    \param name species name.
    \return The species index in deposition velocities.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::DepositionVelocityName(int s) const
  {
    return species_list_dep.at(s);
  }


  //! Returns the index of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The index of a species in the species list.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::DepositionVelocityGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(DepositionVelocityName(s));
  }


  //! Checks whether a species is scavenged.
  /*!
    \param s species global index.
    \return True if the species is scavenged, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasScavenging(int s) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(),
                this->GetSpeciesName(s)) != species_list_scav.end();
  }


  //! Checks whether a species is scavenged.
  /*!
    \param name species name.
    \return True if the species is scavenged, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasScavenging(string name) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(), name)
      != species_list_scav.end();
  }


  //! Returns the index in scavenged species of a given species.
  /*!
    \param s species global index.
    \return The species index in scavenged species.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ScavengingIndex(int s) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(),
                this->GetSpeciesName(s)) - species_list_scav.begin();
  }


  //! Returns the index in scavenged species of a given species.
  /*!
    \param name species name.
    \return The species index in scavenged species.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ScavengingIndex(string name) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(), name)
      - species_list_scav.begin();
  }


  //! Returns the name of a species that is scavenged.
  /*!
    \param s species index in scavenged species.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ScavengingName(int s) const
  {
    return species_list_scav.at(s);
  }


  //! Returns the name of a species that is scavenged.
  /*!
    \param s species index in scavenged species.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ScavengingGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(ScavengingName(s));
  }


  //! Checks whether a species has surface emissions.
  /*!
    \param s species global index.
    \return True if the species has surface emissions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasSurfaceEmission(int s) const
  {
    return find(species_list_surf_emis.begin(), species_list_surf_emis.end(),
                this->GetSpeciesName(s)) != species_list_surf_emis.end();
  }


  //! Checks whether a species has surface emissions.
  /*!
    \param name species name.
    \return True if the species has surface emissions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasSurfaceEmission(string name) const
  {
    return find(species_list_surf_emis.begin(), species_list_surf_emis.end(),
                name) != species_list_surf_emis.end();
  }


  //! Returns the index in surface emissions of a given species.
  /*!
    \param s species global index.
    \return The species index in surface emissions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SurfaceEmissionIndex(int s) const
  {
    return find(species_list_surf_emis.begin(), species_list_surf_emis.end(),
                this->GetSpeciesName(s)) - species_list_surf_emis.begin();
  }


  //! Returns the index in surface emissions of a given species.
  /*!
    \param name species name.
    \return The species index in surface emissions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SurfaceEmissionIndex(string name) const
  {
    return find(species_list_surf_emis.begin(), species_list_surf_emis.end(),
                name) - species_list_surf_emis.begin();
  }


  //! Returns the name of a species with surface emissions.
  /*!
    \param s species index in surface emissions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SurfaceEmissionName(int s) const
  {
    return species_list_surf_emis.at(s);
  }


  //! Returns the name of a species with surface emissions.
  /*!
    \param s species index in surface emissions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SurfaceEmissionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(SurfaceEmissionName(s));
  }


  //! Returns the name of a species with additional surface emissions.
  /*!
    \param s species index in additional surface emissions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::AdditionalSurfaceEmissionName(int s) const
  {
    return species_list_add_surf_emis.at(s);
  }


  //! Checks whether a species has volume emissions.
  /*!
    \param s species global index.
    \return True if the species has volume emissions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasVolumeEmission(int s) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                this->GetSpeciesName(s)) != species_list_vol_emis.end();
  }


  //! Checks whether a species has volume emissions.
  /*!
    \param name species name.
    \return True if the species has volume emissions, false otherwise.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  bool Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::HasVolumeEmission(string name) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                name) != species_list_vol_emis.end();
  }


  //! Returns the index in volume emissions of a given species.
  /*!
    \param s species global index.
    \return The species index in volume emissions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::VolumeEmissionIndex(int s) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                this->GetSpeciesName(s)) - species_list_vol_emis.begin();
  }


  //! Returns the index in volume emissions of a given species.
  /*!
    \param name species name.
    \return The species index in volume emissions.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::VolumeEmissionIndex(string name) const
  {
    return find(species_list_vol_emis.begin(), species_list_vol_emis.end(),
                name) - species_list_vol_emis.begin();
  }


  //! Returns the name of a species with volume emissions.
  /*!
    \param s species index in volume emissions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  string Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::VolumeEmissionName(int s) const
  {
    return species_list_vol_emis.at(s);
  }


  //! Returns the name of a species with volume emissions.
  /*!
    \param s species index in volume emissions.
    \return The species name.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  int Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::VolumeEmissionGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(VolumeEmissionName(s));
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Allocate()
  {
    BaseModel<T>::Allocate();

    /*** Mesh dimensions ***/

    CellWidth_x.resize(this->Nx);
    CellWidth_y.resize(this->Ny);
    CellWidth_z.resize(this->Nz);

    CellCenterDistance_x.resize(this->Nx - 1);
    CellCenterDistance_y.resize(this->Ny - 1);
    CellCenterDistance_z.resize(this->Nz - 1);

    /*** Boundary conditions ***/

    GridS_bc = RegularGrid<T>(Ns_bc);
    GridY4D_interf_bc = RegularGrid<T>(this->y_min - this->Delta_y,
                                       T(this->Ny + 1) * this->Delta_y, 2);
    GridX4D_interf_bc = RegularGrid<T>(this->x_min - this->Delta_x,
                                       T(this->Nx + 1) * this->Delta_x, 2);

    // Along Z.
    BoundaryCondition_z_i.Resize(GridS_bc, this->GridY3D, this->GridX3D);
    FileBoundaryCondition_z_i.Resize(GridS_bc, this->GridY3D, this->GridX3D);
    FileBoundaryCondition_z_f.Resize(GridS_bc, this->GridY3D, this->GridX3D);

    // Along Y.
    BoundaryCondition_y_i.Resize(GridS_bc, this->GridZ4D,
                                 GridY4D_interf_bc, this->GridX4D);
    FileBoundaryCondition_y_i.Resize(GridS_bc, this->GridZ4D,
                                     GridY4D_interf_bc, this->GridX4D);
    FileBoundaryCondition_y_f.Resize(GridS_bc, this->GridZ4D,
                                     GridY4D_interf_bc, this->GridX4D);

    // Along X.
    BoundaryCondition_x_i.Resize(GridS_bc, this->GridZ4D,
                                 this->GridY4D, GridX4D_interf_bc);
    FileBoundaryCondition_x_i.Resize(GridS_bc, this->GridZ4D,
                                     this->GridY4D, GridX4D_interf_bc);
    FileBoundaryCondition_x_f.Resize(GridS_bc, this->GridZ4D,
                                     this->GridY4D, GridX4D_interf_bc);

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

    /*** Rain ***/

    if (this->option_manage["rain"])
      {
        Rain_i.Resize(this->GridY2D, this->GridX2D);
        Rain_f.Resize(this->GridY2D, this->GridX2D);
        FileRain_i.Resize(this->GridY2D, this->GridX2D);
        FileRain_f.Resize(this->GridY2D, this->GridX2D);
      }

    /*** Cloud height ***/

    if (this->option_manage["cloud_base_height"])
      {
        CloudBaseHeight_i.Resize(this->GridY2D, this->GridX2D);
        CloudBaseHeight_f.Resize(this->GridY2D, this->GridX2D);
        FileCloudBaseHeight_i.Resize(this->GridY2D, this->GridX2D);
        FileCloudBaseHeight_f.Resize(this->GridY2D, this->GridX2D);
      }

    if (this->option_manage["cloud_top_height"])
      {
        CloudTopHeight_i.Resize(this->GridY2D, this->GridX2D);
        CloudTopHeight_f.Resize(this->GridY2D, this->GridX2D);
        FileCloudTopHeight_i.Resize(this->GridY2D, this->GridX2D);
        FileCloudTopHeight_f.Resize(this->GridY2D, this->GridX2D);
      }

    /*** Specific humidity ***/

    if (this->option_manage["specific_humidity"])
      {
        SpecificHumidity_i.Resize(this->GridZ3D, this->GridY3D,
                                  this->GridX3D);
        SpecificHumidity_f.Resize(this->GridZ3D, this->GridY3D,
                                  this->GridX3D);
        FileSpecificHumidity_i.Resize(this->GridZ3D, this->GridY3D,
                                      this->GridX3D);
        FileSpecificHumidity_f.Resize(this->GridZ3D, this->GridY3D,
                                      this->GridX3D);
      }

    /*** Scavenging coefficients ***/

    GridS_scav = RegularGrid<T>(Ns_scav);

    ScavengingBelowCloudCoefficient_i.Resize(GridS_scav, this->GridZ4D,
                                             this->GridY4D, this->GridX4D);
    ScavengingBelowCloudCoefficient_f.Resize(GridS_scav, this->GridZ4D,
                                             this->GridY4D, this->GridX4D);

    ScavengingInCloudCoefficient_i.Resize(GridS_scav, this->GridZ4D,
                                          this->GridY4D, this->GridX4D);
    ScavengingInCloudCoefficient_f.Resize(GridS_scav, this->GridZ4D,
                                          this->GridY4D, this->GridX4D);

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

    /*** Deposition velocities ***/

    GridS_dep = RegularGrid<T>(Ns_dep);

    DepositionVelocity_i.Resize(GridS_dep, this->GridY3D, this->GridX3D);
    DepositionVelocity_f.Resize(GridS_dep, this->GridY3D, this->GridX3D);
    FileDepositionVelocity_i.Resize(GridS_dep, this->GridY3D, this->GridX3D);
    FileDepositionVelocity_f.Resize(GridS_dep, this->GridY3D, this->GridX3D);

    /*** Surface emissions ***/

    GridS_surf_emis = RegularGrid<T>(Ns_surf_emis);

    SurfaceEmission_i.Resize(GridS_surf_emis, this->GridY3D, this->GridX3D);
    SurfaceEmission_f.Resize(GridS_surf_emis, this->GridY3D, this->GridX3D);
    FileSurfaceEmission_i.Resize(GridS_surf_emis,
                                 this->GridY3D, this->GridX3D);
    FileSurfaceEmission_f.Resize(GridS_surf_emis,
                                 this->GridY3D, this->GridX3D);

    /*** Additional surface emissions ***/

    GridS_add_surf_emis = RegularGrid<T>(Ns_add_surf_emis);

    AdditionalSurfaceEmission_i.Resize(GridS_add_surf_emis,
                                       this->GridY3D, this->GridX3D);
    AdditionalSurfaceEmission_f.Resize(GridS_add_surf_emis,
                                       this->GridY3D, this->GridX3D);
    FileAdditionalSurfaceEmission_i.Resize(GridS_add_surf_emis,
                                           this->GridY3D, this->GridX3D);
    FileAdditionalSurfaceEmission_f.Resize(GridS_add_surf_emis,
                                           this->GridY3D, this->GridX3D);

    /*** Volume emissions ***/

    GridS_vol_emis = RegularGrid<T>(Ns_vol_emis);
    GridZ_vol_emis = RegularGrid<T>(Nz_vol_emis);

    VolumeEmission_i.Resize(GridS_vol_emis, GridZ_vol_emis,
                            this->GridY4D, this->GridX4D);
    VolumeEmission_f.Resize(GridS_vol_emis, GridZ_vol_emis,
                            this->GridY4D, this->GridX4D);
    FileVolumeEmission_i.Resize(GridS_vol_emis, GridZ_vol_emis,
                                this->GridY4D, this->GridX4D);
    FileVolumeEmission_f.Resize(GridS_vol_emis, GridZ_vol_emis,
                                this->GridY4D, this->GridX4D);

    /*** Deposition fluxes ***/

    if (this->option_process["collect_dry_flux"])
      DryDepositionFlux.Resize(GridS_dep, this->GridY3D, this->GridX3D);

    if (this->option_process["collect_wet_flux"])
      {
        if (this->scavenging_below_cloud_model != "none")
          WetDepositionFlux.Resize(GridS_scav, this->GridY3D, this->GridX3D);
        if (this->scavenging_in_cloud_model != "none")
          InCloudWetDepositionFlux.Resize(this->GridS_scav, this->GridY3D,
                                          this->GridX3D);
      }
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Init()
  {
    BaseModel<T>::Init();

    /*** Mesh dimensions ***/

    if (option_cartesian)
      {
        CellWidth_x = this->Delta_x;
        CellWidth_y = this->Delta_y;
      }
    else
      ComputeCellWidth(this->Delta_x, this->GridY3D_interf.GetArray(),
                       CellWidth_x, CellWidth_y);
    for (int k = 0; k < this->Nz; k++)
      CellWidth_z(k) = this->GridZ3D_interf(k + 1) - this->GridZ3D_interf(k);

    if (option_cartesian)
      {
        CellCenterDistance_x = this->Delta_x;
        CellCenterDistance_y = this->Delta_y;
      }
    else
      ComputeCellCenterDistance(this->Delta_x,
                                this->GridY3D_interf.GetArray(),
                                CellCenterDistance_x, CellCenterDistance_y);
    for (int k = 0; k < this->Nz - 1; k++)
      CellCenterDistance_z(k) = this->GridZ3D(k + 1) - this->GridZ3D(k);

    /*** Input data ***/

    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitAllData();

    /*** Initial conditions ***/

    if (this->option_manage["initial_condition"])
      {
        this->Concentration.SetZero();

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
      }

    /*** Deposition ***/

    if (this->option_process["collect_dry_flux"])
      for (int s = 0; s < Ns_dep; s++)
        for (int j = 0; j < this->Ny; j++)
          for (int i = 0; i < this->Nx; i++)
            DryDepositionFlux(s, j, i) = DepositionVelocity_f(s, j, i)
              * this->Concentration(DepositionVelocityGlobalIndex(s),
                                    0, j, i);

    if (this->option_process["collect_wet_flux"]
        && this->scavenging_below_cloud_model != "none")
      {
        this->WetDepositionFlux.SetZero();
        // Depth of the layer that is below the cloud.
        T depth_below;
        for (int s = 0; s < Ns_scav; s++)
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_f(j, i) > 0.)
                for (int k = 0; k < this->Nz; k++)
                  if (this->CloudBaseHeight_f(j, i) > this->GridZ4D_interf(k))
                    {
                      depth_below = min(this->CloudBaseHeight_f(j, i),
                                        this->GridZ4D_interf(k + 1))
                        - this->GridZ4D_interf(k);
                      WetDepositionFlux(s, j, i) +=
                        this->Concentration(ScavengingGlobalIndex(s), k, j, i)
                        / this->Delta_t
                        * (1. - exp(- this->Delta_t *
                                    ScavengingBelowCloudCoefficient_f
                                    (s, k, j, i))) * depth_below;
                    }
      }

    if (this->option_process["collect_wet_flux"]
        && this->scavenging_in_cloud_model != "none")
      {
        this->InCloudWetDepositionFlux.SetZero();
        if (this->option_process["with_transport_in_cloud_scavenging"])
          {
            // Depth of the layer that is inside the cloud.
            T depth_inside;
            for (int s = 0; s < Ns_scav; s++)
              for (int j = 0; j < this->Ny; j++)
                for (int i = 0; i < this->Nx; i++)
                  for (int k = 0; k < this->Nz; k++)
                    if (this->CloudBaseHeight_f(j, i)
                        < this->GridZ4D_interf(k + 1)
                        && this->CloudTopHeight_f(j, i)
                        > this->GridZ4D_interf(k))
                      {
                        depth_inside = min(this->CloudTopHeight_f(j, i),
                                           this->GridZ4D_interf(k + 1))
                          - max(this->CloudBaseHeight_f(j, i),
                                this->GridZ4D_interf(k));
                        InCloudWetDepositionFlux(s, j, i) +=
                          this->Concentration(ScavengingGlobalIndex(s), k, j, i)
                          / this->Delta_t
                          * (1. - exp(- this->Delta_t *
                                      ScavengingBelowCloudCoefficient_f
                                      (s, k, j, i))) * depth_inside;
                      }
          }
      }

    /*** Numerical schemes ***/

    Advection_.Init(*this);
    Diffusion_.Init(*this);
  }

  //! Model initialization for each step.
  /*! It reads on file the data that are needed for the current step.
   */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitStep()
  {
    if (this->data_gone_through_initstep
        && this->data_date == this->current_date)
      return;

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

    /*** Specific humidity, Rain, cloud height and scavenging
         coefficients ***/

    if (this->option_manage["specific_humidity"])
      this->UpdateData("meteo", "SpecificHumidity", FileSpecificHumidity_i,
                       FileSpecificHumidity_f, SpecificHumidity_i,
                       SpecificHumidity_f);

    if (this->option_manage["rain"])
      this->UpdateData("meteo", "Rain", FileRain_i, FileRain_f,
                       Rain_i, Rain_f, false);

    if (this->option_manage["cloud_base_height"])
      this->UpdateData("meteo", "CloudBaseHeight", FileCloudBaseHeight_i,
                       FileCloudBaseHeight_f, CloudBaseHeight_i,
                       CloudBaseHeight_f);

    if (this->option_manage["cloud_top_height"])
      this->UpdateData("meteo", "CloudTopHeight", FileCloudTopHeight_i,
                       FileCloudTopHeight_f, CloudTopHeight_i,
                       CloudTopHeight_f);

    if (this->option_manage["scavenging_below_cloud_coefficient"]
        && this->scavenging_below_cloud_model != "none")
      {
        ScavengingBelowCloudCoefficient_i.GetArray()
          = ScavengingBelowCloudCoefficient_f.GetArray();
        if (this->scavenging_below_cloud_model == "microphysical")
          InitScavengingCoefficient(Temperature_f, Pressure_f,
                                    CloudBaseHeight_f, Rain_f,
                                    ScavengingBelowCloudCoefficient_f);
        else if (this->scavenging_below_cloud_model == "belot")
          InitScavengingCoefficient(Rain_f, belot_below_cloud_constant_a,
                                    belot_below_cloud_constant_b,
                                    ScavengingBelowCloudCoefficient_f);
        else if (this->scavenging_below_cloud_model == "constant")
          InitScavengingCoefficient(ScavengingBelowCloudCoefficient_f);
      }

    if (this->option_manage["scavenging_in_cloud_coefficient"]
        && this->option_process["with_transport_in_cloud_scavenging"])
      {
        ScavengingInCloudCoefficient_i.GetArray()
          = ScavengingInCloudCoefficient_f.GetArray();
        if (this->scavenging_in_cloud_model == "pudykiewicz")
          InitScavengingCoefficient(Temperature_f, Pressure_f,
                                    SpecificHumidity_f,
                                    ScavengingInCloudCoefficient_f);
        else if (this->scavenging_in_cloud_model == "belot")
          InitScavengingCoefficient(Rain_f, belot_in_cloud_constant_a,
                                    belot_in_cloud_constant_b,
                                    ScavengingInCloudCoefficient_f);
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

    /*** Boundary conditions ***/

    if (this->option_manage["boundary_condition"])
      for (int i = 0; i < Ns_bc; i++)
        {
          string filename
            = this->input_files["boundary_condition"](species_list_bc[i]);
          Date date = this->input_files["boundary_condition"].GetDateMin();
          T Delta_t = this->input_files["boundary_condition"].GetDelta_t();
          this->UpdateData(find_replace(filename, "&c", "z"), date, Delta_t,
                           FileBoundaryCondition_z_i,
                           FileBoundaryCondition_z_f,
                           i, BoundaryCondition_z_i);
          this->UpdateData(find_replace(filename, "&c", "y"), date, Delta_t,
                           FileBoundaryCondition_y_i,
                           FileBoundaryCondition_y_f,
                           i, BoundaryCondition_y_i);
          this->UpdateData(find_replace(filename, "&c", "x"), date, Delta_t,
                           FileBoundaryCondition_x_i,
                           FileBoundaryCondition_x_f,
                           i, BoundaryCondition_x_i);
        }

    /*** Deposition velocities ***/

    if (this->option_manage["deposition_velocity"])
      for (int i = 0; i < Ns_dep; i++)
        this->UpdateData("deposition_velocity", species_list_dep[i],
                         FileDepositionVelocity_i, FileDepositionVelocity_f,
                         i, DepositionVelocity_i, DepositionVelocity_f);

    /*** Surface emissions ***/

    if (this->option_manage["surface_emission"])
      for (int i = 0; i < Ns_surf_emis; i++)
        this->UpdateData("surface_emission", species_list_surf_emis[i],
                         FileSurfaceEmission_i, FileSurfaceEmission_f, i,
                         SurfaceEmission_i, SurfaceEmission_f);

    /*** Additional surface emissions ***/

    if (this->option_manage["additional_surface_emission"])
      for (int i = 0; i < Ns_add_surf_emis; i++)
        this->UpdateData("additional_surface_emission",
                         species_list_add_surf_emis[i],
                         FileAdditionalSurfaceEmission_i,
                         FileAdditionalSurfaceEmission_f, i,
                         AdditionalSurfaceEmission_i,
                         AdditionalSurfaceEmission_f);

    if (this->option_manage["additional_surface_emission"])
      for (int s = 0; s < Ns_add_surf_emis; s++)
        {
          int gs = SurfaceEmissionIndex(AdditionalSurfaceEmissionName(s));
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              SurfaceEmission_f(gs, j, i)
                += AdditionalSurfaceEmission_f(s, j, i);
        }

    /*** Volume emissions ***/

    if (this->option_manage["volume_emission"])
      for (int i = 0; i < Ns_vol_emis; i++)
        this->UpdateData("volume_emission", species_list_vol_emis[i],
                         FileVolumeEmission_i, FileVolumeEmission_f, i,
                         VolumeEmission_i, VolumeEmission_f);

    this->data_date = this->current_date;
    this->data_gone_through_initstep = true;
    this->data_gone_through_forward = false;
  }


  //! Initializes scavenging coefficients.
  /*! Assigns a constant value to below-cloud scavenging coefficients
    for each species.
    \param ScavengingCoefficient_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitScavengingCoefficient(Data<T, 4>& ScavengingCoefficient_)
  {
    int Ns_scav = ScavengingCoefficient_.GetLength(0);
    int Nz = ScavengingCoefficient_.GetLength(1);
    int Ny = ScavengingCoefficient_.GetLength(2);
    int Nx = ScavengingCoefficient_.GetLength(3);

    for (int s = 0; s < Ns_scav; s++)
      {
        T scavenging_value = this->scavenging_constant[ScavengingName(s)];
        for (int k = 0; k < Nz; k++)
          for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
              ScavengingCoefficient_(s, k, j, i) = scavenging_value;
      }
  }

  //! Initializes scavenging coefficients.
  /*! Computes scavenging coefficients following Belot (1988).
    \param Rain_ rain intensity (mm/h).
    \param[in] belot_constant_a Belot coefficient a.
    \param[in] belot_constant_b Belot coefficient b.
    \param ScavengingCoefficient_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitScavengingCoefficient(Data<T, 2>& Rain_,
                              map<string, T>& belot_constant_a,
                              map<string, T>& belot_constant_b,
                              Data<T, 4>& ScavengingCoefficient_)
  {
    int Ns_scav = ScavengingCoefficient_.GetLength(0);
    int Nz = ScavengingCoefficient_.GetLength(1);
    int Ny = ScavengingCoefficient_.GetLength(2);
    int Nx = ScavengingCoefficient_.GetLength(3);

    ScavengingCoefficient_.Fill(0.);
    for (int s = 0; s < Ns_scav; s++)
      {
        T cste_a = belot_constant_a[ScavengingName(s)];
        T cste_b = belot_constant_b[ScavengingName(s)];
        for (int j = 0; j < Ny; j++)
          for (int i = 0; i < Nx; i++)
            if (Rain_(j, i) > 0.)
              {
                T scavenging_value = cste_a * pow(Rain_(j, i), cste_b);
                for (int k = 0; k < Nz; k++)
                  ScavengingCoefficient_(s, k, j, i) = scavenging_value;
              }
      }
  }


  //! Initializes scavenging coefficients.
  /*! Computes scavenging coefficients following Pudykiewicz (1989).
    \param Temperature_ temperature (K).
    \param Pressure_ pressure (Pa).
    \param ScavengingCoefficient_ (output) the scavenging
    coefficients (1/s).
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitScavengingCoefficient(Data<T, 3>& Temperature_,
                              Data<T, 3>& Pressure_,
                              Data<T, 3>& SpecificHumidity_,
                              Data<T, 4>& ScavengingCoefficient_)
  {
    int Nz = ScavengingCoefficient_.GetLength(1);
    int Ny = ScavengingCoefficient_.GetLength(2);
    int Nx = ScavengingCoefficient_.GetLength(3);

    _compute_scavenging_coefficient_pudykiewicz(&Nx, &Ny, &Nz, &Ns_scav,
                                                this->GridZ4D
                                                .GetArray().data(),
                                                Temperature_.GetData(),
                                                Pressure_.GetData(),
                                                SpecificHumidity_.GetData(),
                                                ScavengingCoefficient_
                                                .GetData());
  }


  //! Initializes scavenging coefficients.
  /*! Computes below-cloud scavenging coefficients of a reversibly soluble
    gas.
    \param Temperature_ temperature (K).
    \param Pressure_ pressure (Pa).
    \param CloudBaseHeight_ cloud basis height (m).
    \param Rain_ rain intensity (mm/h).
    \param ScavengingCoefficient_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::InitScavengingCoefficient(Data<T, 3>& Temperature_,
                              Data<T, 3>& Pressure_,
                              Data<T, 2>& CloudBaseHeight_,
                              Data<T, 2>& Rain_,
                              Data<T, 4>& ScavengingCoefficient_)
  {
    int Nz = ScavengingCoefficient_.GetLength(1);
    int Ny = ScavengingCoefficient_.GetLength(2);
    int Nx = ScavengingCoefficient_.GetLength(3);

    // Henry constants of scavenged species.
    Array<T, 1> henry_constant_(Ns_scav);
    // Gas phase diffusivities of scavenged species.
    Array<T, 1> gas_phase_diffusivity_(Ns_scav);

    // Extracts species data for scavenged species only.
    for (int s = 0; s < this->Ns_scav; s++)
      {
        henry_constant_(s) = henry_constant[ScavengingName(s)];
        gas_phase_diffusivity_(s) = gas_phase_diffusivity[ScavengingName(s)];
      }

    int one = 1;
    _compute_scavenging_coefficient(&Nx, &Ny, &Nz, &Ns_scav,
                                    &one, &one, &one,
                                    this->GridZ4D.GetArray().data(),
                                    gas_phase_diffusivity_.data(),
                                    henry_constant_.data(),
                                    Temperature_.GetData(),
                                    Pressure_.GetData(),
                                    Rain_.GetData(),
                                    CloudBaseHeight_.GetData(),
                                    ScavengingCoefficient_.GetData());
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SetDate(Date date)
  {
    BaseModel<T>::SetDate(date);

    if (date != this->data_date)
      Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitAllData();
  }


  ///////////////////////////
  // NUMERICAL INTEGRATION //
  ///////////////////////////


  //! Performs one advection step.
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Advection()
  {
    Advection_.Forward(*this);
  }


  //! Performs one advection step backward.
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Advection_b()
  {
    Advection_.Backward(*this);
  }


  //! Performs one diffusion step.
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Diffusion()
  {
    Diffusion_.Forward(*this);
  }

  //! Performs one diffusion step backward.
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Diffusion_b()
  {
    Diffusion_.Backward(*this);
  }


  //! Performs the time integration of point emissions.
  template<class T, class ClassAdvection, class ClassDiffusion>
  void
  Polair3DTransport<T, ClassAdvection, ClassDiffusion>::PointEmission()
  {
    int k, j, i, s, Ns, Npoint;
    int Nemis = this->PointEmissionManager->GetNumberEmission();

    for (int emission = 0; emission < Nemis; emission++)
      {
        vector<int> emitted_species_index =
          this->PointEmissionManager->GetEmittedSpeciesIndex(emission);
        Ns = emitted_species_index.size();

        for (int species = 0; species < Ns; species++)
          if (this->PointEmissionManager->
              IsEmitting(this->current_date, this->next_date, emission))
            {
              s = emitted_species_index[species];
              Array<T, 2> point_emission;
              this->PointEmissionManager->
                GetEmission(this->current_date, this->next_date, species,
                            emission, point_emission);
              Npoint = point_emission.extent(0);

              for (int index = 0; index < Npoint; index++)
                {
                  GetCellIndices(point_emission(index, 0),
                                 point_emission(index, 1),
                                 point_emission(index, 2),
                                 k, j, i);
                  T cell_volume;
                  if (option_cartesian)
                    cell_volume = CellWidth_z(k) * CellWidth_y(j)
                      * CellWidth_x(i);
                  else
                    cell_volume
                      = ComputeCellVolume(this->Delta_x, this->Delta_y,
                                          CellWidth_z(k),
                                          this->y_min + j * this->Delta_y);
                  this->Concentration(s, k, j, i) += point_emission(index, 3)
                    / cell_volume;
                }
            }
      }
  }


  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    adds volume emissions. These three processes are split (operator
    splitting).
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Forward()
  {

    /*** Air density ***/

    if (!this->data_gone_through_forward
        && this->option_process["with_air_density"])
      {
        InterpolateInterface_z(CellCenterDistance_z, CellWidth_z,
                               AirDensity_f, AirDensity_interf_z_f);
        InterpolateInterface_y(CellCenterDistance_y, CellWidth_y,
                               AirDensity_f, AirDensity_interf_y_f);
        InterpolateInterface_x(CellCenterDistance_x, CellWidth_x,
                               AirDensity_f, AirDensity_interf_x_f);
      }

    /*** Winds ***/

    if (!this->data_gone_through_forward
        && this->option_manage["horizontal_wind"])
      if (!option_cartesian)
        {
          TransformZonalWind(ZonalWind_i);
          TransformMeridionalWind(MeridionalWind_i);
        }

    if (!this->data_gone_through_forward
        && this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
        ComputeVerticalWind(CellWidth_x, CellWidth_y, CellWidth_z,
                            AirDensity_interf_x_i, ZonalWind_i,
                            AirDensity_interf_y_i, MeridionalWind_i,
                            AirDensity_interf_z_i, VerticalWind_i);
      else
        ComputeVerticalWind(CellWidth_x, CellWidth_y, CellWidth_z,
                            ZonalWind_i, MeridionalWind_i, VerticalWind_i);

    /*** Diffusion coefficients ***/

    if (!this->data_gone_through_forward
        && this->option_manage["vertical_diffusion"])
      {
        // Computes rho * Kz.
        if (this->option_process["with_air_density"])
          VerticalDiffusionCoefficient_f.GetArray() =
            AirDensity_interf_z_f.GetArray()
            * VerticalDiffusionCoefficient_f.GetArray();
      }

    if (!this->data_gone_through_forward
        && this->option_manage["horizontal_diffusion"]
        && option_isotropic_diffusion)
      {
        LinearInterpolationRegular(VerticalDiffusionCoefficient_i,
                                   ZonalDiffusionCoefficient_i);
        ZonalDiffusionCoefficient_i.ThresholdMin(0.);
        LinearInterpolationRegular(VerticalDiffusionCoefficient_i,
                                   MeridionalDiffusionCoefficient_i);
        MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

        if (!option_cartesian)
          {
            TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
                                    ZonalDiffusionCoefficient_i);
            TransformMeridionalDiffusion(MeridionalDiffusionCoefficient_i);
          }
      }

    this->data_gone_through_forward = true;

    /*** Time integration ***/

    if (this->option_process["with_advection"])
      Advection();
    //! Dry fluxes are collected before dry deposition is performed
    //! in diffusion.
    if (this->option_process["collect_dry_flux"])
      for (int s = 0; s < Ns_dep; s++)
        for (int j = 0; j < this->Ny; j++)
          for (int i = 0; i < this->Nx; i++)
            DryDepositionFlux(s, j, i) = 0.5
              * (DepositionVelocity_i(s, j, i)
                 + DepositionVelocity_f(s, j, i))
              * this->Concentration(DepositionVelocityGlobalIndex(s),
                                    0, j, i);
    if (this->option_process["with_diffusion"])
      Diffusion();
    if (this->option_process["with_point_emission"])
      PointEmission();
    if (this->option_process["with_volume_emission"])
      for (int s = 0; s < this->Ns; s++)
        {
          int k, j, i, emis_s;
          if (HasVolumeEmission(s))
            {
              emis_s = VolumeEmissionIndex(s);
              for (k = 0; k < Nz_vol_emis; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration(s, k, j, i) +=
                      this->Delta_t * VolumeEmission_i(emis_s, k, j, i);
            }
        }

    if (this->option_process["collect_wet_flux"])
      {
        if (this->scavenging_below_cloud_model != "none")
          this->WetDepositionFlux.SetZero();
        if (this->scavenging_in_cloud_model != "none")
          this->InCloudWetDepositionFlux.SetZero();
      }
    if (this->scavenging_below_cloud_model != "none")
      for (int s = 0; s < Ns_scav; s++)
        {
          T scavenging_ratio;
          T cloud_height_mean;
          // Depth of the layer that is below the cloud.
          T depth_below;
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  cloud_height_mean = 0.5
                    * (this->CloudBaseHeight_i(j, i)
                       + this->CloudBaseHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (cloud_height_mean > this->GridZ4D_interf(k))
                      {
                        depth_below = min(cloud_height_mean,
                                          this->GridZ4D_interf(k + 1))
                          - this->GridZ4D_interf(k);
                        scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (ScavengingBelowCloudCoefficient_i(s, k, j, i)
                                 + ScavengingBelowCloudCoefficient_f(s, k, j, i)));
                        if (this->option_process["collect_wet_flux"])
                          WetDepositionFlux(s, j, i) +=
                            this->Concentration(ScavengingGlobalIndex(s),
                                                k, j, i)
                            / this->Delta_t * (1. - scavenging_ratio)
                            * depth_below;
                        this->Concentration(ScavengingGlobalIndex(s), k, j, i)
                          *= 1. - (1. - scavenging_ratio)
                          * depth_below / (this->GridZ4D_interf(k + 1)
                                           - this->GridZ4D_interf(k));
                      }
                }
        }

    if (this->option_process["with_transport_in_cloud_scavenging"])
      for (int s = 0; s < Ns_scav; s++)
        {
          T scavenging_ratio;
          T cloud_base_height_mean, cloud_top_height_mean;
          // Depth of the layer that is inside the cloud.
          T depth_inside;
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  cloud_base_height_mean = 0.5
                    * (this->CloudBaseHeight_i(j, i)
                       + this->CloudBaseHeight_f(j, i));
                  cloud_top_height_mean = 0.5
                    * (this->CloudTopHeight_i(j, i)
                       + this->CloudTopHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (cloud_base_height_mean < this->GridZ4D_interf(k + 1)
                        && cloud_top_height_mean > this->GridZ4D_interf(k))
                      {
                        depth_inside = min(cloud_top_height_mean,
                                           this->GridZ4D_interf(k + 1))
                          - max(cloud_base_height_mean,
                                this->GridZ4D_interf(k));
                        scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (ScavengingInCloudCoefficient_i(s, k, j, i)
                                 + ScavengingInCloudCoefficient_f
                                 (s, k, j, i)));
                        if (this->option_process["collect_wet_flux"])
                          InCloudWetDepositionFlux(s, j, i) +=
                            this->Concentration(ScavengingGlobalIndex(s),
                                                k, j, i)
                            / this->Delta_t * (1. - scavenging_ratio)
                            * depth_inside;
                        this->Concentration(ScavengingGlobalIndex(s), k, j, i)
                          *= 1. - (1. - scavenging_ratio)
                          * depth_inside / (this->GridZ4D_interf(k + 1)
                                            - this->GridZ4D_interf(k));
                      }
                }
        }

    this->AddTime(this->Delta_t);
    this->step++;
  }


  //! Prepares for backward integration.
  /*! It sets flag for backward integration, then allocates memories for
    adjoint concentration data and adjoint field data if necessary. The
    adjoint concentration data and adjoint field data are intialized to zero.
    \param flag true for backward integration; false for forward integration.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::SetBackward(bool flag)
  {
    BaseModel<T>::SetBackward(flag);
  }


  //! Performs one step of backward integration of adjoint model.
  /*! It performs one advection step, then one diffusion step and finally
    adds volume emissions. These three processes are split (operator
    splitting). Then adjoint model of the above processes is integrated
    one step backward.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Backward()
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
      if (!option_cartesian)
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

    if (this->option_manage["horizontal_diffusion"]
        && option_isotropic_diffusion)
      {
        LinearInterpolationRegular(VerticalDiffusionCoefficient_i,
                                   ZonalDiffusionCoefficient_i);
        ZonalDiffusionCoefficient_i.ThresholdMin(0.);
        LinearInterpolationRegular(VerticalDiffusionCoefficient_i,
                                   MeridionalDiffusionCoefficient_i);
        MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

        if (!option_cartesian)
          {
            TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
                                    ZonalDiffusionCoefficient_i);
            TransformMeridionalDiffusion(MeridionalDiffusionCoefficient_i);
          }
      }

    /*** Forward integrations and trajectory generations ***/

    int n_tap = 2;
    Array<T, 5> conc_tap(n_tap, this->Ns, this->Nz, this->Ny, this->Nx);

    if (this->option_process["with_advection"])
      {
        Array<T, 4> conc(&conc_tap(0, 0, 0, 0, 0),
                         shape(this->Ns, this->Nz, this->Ny, this->Nx));
        conc = this->Concentration.GetArray();

        Advection();
      }

    if (this->option_process["with_diffusion"])
      {
        Array<T, 4> conc(&conc_tap(1, 0, 0, 0, 0),
                         shape(this->Ns, this->Nz, this->Ny, this->Nx));
        conc = this->Concentration.GetArray();

        Diffusion();
      }

    if (this->option_process["with_point_emission"])
      PointEmission();

    if (this->option_process["with_volume_emission"])
      for (int s = 0; s < this->Ns; s++)
        {
          int k, j, i, emis_s;
          if (HasVolumeEmission(s))
            {
              emis_s = VolumeEmissionIndex(s);
              for (k = 0; k < Nz_vol_emis; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration(s, k, j, i) +=
                      this->Delta_t * VolumeEmission_i(emis_s, k, j, i);
            }
        }

    for (int s = 0; s < this->Ns; s++)
      if (HasScavenging(s) && this->scavenging_below_cloud_model != "none")
        {
          // Depth of the layer that is below the cloud.
          T depth_below;
          int scav_s = this->ScavengingIndex(s);
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  T CloudBaseHeight_mean = 0.5 *
                    (this->CloudBaseHeight_i(j, i)
                     + this->CloudBaseHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (CloudBaseHeight_mean > this->GridZ4D_interf(k))
                      {
                        depth_below = min(CloudBaseHeight_mean,
                                          this->GridZ4D_interf(k + 1))
                          - this->GridZ4D_interf(k);
                        T scavenging_ratio
                          = exp(-0.5 * this->Delta_t
                                * (this->ScavengingBelowCloudCoefficient_i
                                   (scav_s, k, j, i)
                                   + this->ScavengingBelowCloudCoefficient_f
                                   (scav_s, k, j, i)));
                        this->Concentration(s, k, j, i) *= 1.
                          - (1. - scavenging_ratio)
                          * depth_below / (this->GridZ4D_interf(k + 1)
                                           - this->GridZ4D_interf(k));
                      }
                }
        }

    for (int s = 0; s < this->Ns; s++)
      if (HasScavenging(s) && this->option_process["with_transport_in_cloud_scavenging"])
        {
          // Depth of the layer that is inside the cloud.
          T depth_inside;
          int scav_s = this->ScavengingIndex(s);
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  T CloudBaseHeight_mean = 0.5 *
                    (this->CloudBaseHeight_i(j, i)
                     + this->CloudBaseHeight_f(j, i));
                  T CloudTopHeight_mean = 0.5 *
                    (this->CloudTopHeight_i(j, i)
                     + this->CloudTopHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (CloudBaseHeight_mean < this->GridZ4D_interf(k + 1)
                        && CloudTopHeight_mean > this->GridZ4D_interf(k))
                      {
                        depth_inside = min(CloudTopHeight_mean,
                                           this->GridZ4D_interf(k + 1))
                          - max(CloudBaseHeight_mean,
                                this->GridZ4D_interf(k));
                        T scavenging_ratio
                          = exp(-0.5 * this->Delta_t
                                * (this->ScavengingInCloudCoefficient_i
                                   (scav_s, k, j, i)
                                   + this->ScavengingInCloudCoefficient_f
                                   (scav_s, k, j, i)));
                        this->Concentration(s, k, j, i) *= 1.
                          - (1. - scavenging_ratio)
                          * depth_inside / (this->GridZ4D_interf(k + 1)
                                            - this->GridZ4D_interf(k));
                      }
                }
        }

    Array<T, 4> res(this->Ns, this->Nz, this->Ny, this->Nx);
    res = this->Concentration.GetArray();


    /*** Backward integrations ***/

    for (int s = 0; s < this->Ns; s++)
      if (HasScavenging(s)
          && this->option_process["with_transport_in_cloud_scavenging"])
        {
          // Depth of the layer that is inside the cloud.
          T depth_inside;
          int scav_s = this->ScavengingIndex(s);
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  T CloudBaseHeight_mean = 0.5 *
                    (this->CloudBaseHeight_i(j, i)
                     + this->CloudBaseHeight_f(j, i));
                  T CloudTopHeight_mean = 0.5 *
                    (this->CloudTopHeight_i(j, i)
                     + this->CloudTopHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (CloudBaseHeight_mean < this->GridZ4D_interf(k + 1)
                        && CloudTopHeight_mean > this->GridZ4D_interf(k))
                      {
                        depth_inside = min(CloudTopHeight_mean,
                                           this->GridZ4D_interf(k + 1))
                          - max(CloudBaseHeight_mean,
                                this->GridZ4D_interf(k));
                        T scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (this->ScavengingInCloudCoefficient_i
                                 (scav_s, k, j, i)
                                 + this->ScavengingInCloudCoefficient_f
                                 (scav_s, k, j, i)));
                        this->Concentration_ccl(s, k, j, i) *= 1.
                          - (1. - scavenging_ratio)
                          * depth_inside / (this->GridZ4D_interf(k + 1)
                                            - this->GridZ4D_interf(k));
                      }
                }
        }

    for (int s = 0; s < this->Ns; s++)
      if (HasScavenging(s) && this->scavenging_below_cloud_model != "none")
        {
          // Depth of the layer that is below the cloud.
          T depth_below;
          int scav_s = this->ScavengingIndex(s);
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  T CloudBaseHeight_mean = 0.5 *
                    (this->CloudBaseHeight_i(j, i)
                     + this->CloudBaseHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (CloudBaseHeight_mean > this->GridZ4D_interf(k))
                      {
                        depth_below = min(CloudBaseHeight_mean,
                                          this->GridZ4D_interf(k + 1))
                          - this->GridZ4D_interf(k);
                        T scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (this->ScavengingBelowCloudCoefficient_i
                                 (scav_s, k, j, i)
                                 + this->ScavengingBelowCloudCoefficient_f
                                 (scav_s, k, j, i)));
                        this->Concentration_ccl(s, k, j, i) *= 1.
                          - (1. - scavenging_ratio)
                          * depth_below / (this->GridZ4D_interf(k + 1)
                                           - this->GridZ4D_interf(k));
                      }
                }
        }

    if (this->option_process["with_diffusion"])
      {
        Array<T, 4> conc(&conc_tap(1, 0, 0, 0, 0),
                         shape(this->Ns, this->Nz, this->Ny, this->Nx));
        this->Concentration.GetArray() = conc;

        Diffusion_b();
      }

    if (this->option_process["with_advection"])
      {
        Array<T, 4> conc(&conc_tap(0, 0, 0, 0, 0),
                         shape(this->Ns, this->Nz, this->Ny, this->Nx));
        this->Concentration.GetArray() = conc;

        Advection_b();
      }

    this->AddTime(this->Delta_t);
    this->step++;


    /*** Resets concentration results ***/

    this->Concentration.GetArray() = res;
  }


  ///////////////////
  // OTHER METHODS //
  ///////////////////


  //! Computes the volume (in cubic meters) of a cell at a given latitude.
  /*!
    \param Delta_x_ step along x in degrees.
    \param Delta_y_ step along y in degrees.
    \param Delta_z_ cell width along z in meters.
    \param lat latitude of the cell in degrees.
    \return The cell volume in cubic meters.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  T Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::ComputeCellVolume(T Delta_x_, T Delta_y_, T Delta_z_, T lat)
  {
    const T earth_radius_2 = 6371229. * 6371229.;
    const T pi(3.14159265358979323846264);

    return earth_radius_2 * cos(lat * pi / 180.)
      * Delta_x_ * (pi / 180.) * Delta_y_ * (pi / 180.) * Delta_z_;
  }


  //! Computes cells widths (in meters).
  /*!
    \param Delta_x_ step along x in degrees.
    \param GridY_interf_ interface coordinates along y in degrees.
    \param CellWidth_x_ (output) cells widths along x in meters.
    \param CellWidth_y_ (output) cells widths along y in meters.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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


  /*! \brief Transforms the zonal wind to ease the numerical integration of
    advection. */
  /*! Formula: ZonalWind = ZonalWind / cos(latitude).
    \param ZonalWind (input/output) zonal wind.
    \note Coordinates associated with ZonalWind must be in degrees.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
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
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitAllData()
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
        if (!option_cartesian)
          TransformZonalWind(ZonalWind_i);

        this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
                       FileMeridionalWind_f, this->current_date,
                       MeridionalWind_i);
        if (!option_cartesian)
          TransformMeridionalWind(MeridionalWind_i);
      }

    /*** Specific humidity, Rain, cloud height and scavenging
         coefficients ***/

    if (this->option_manage["specific_humidity"])
      this->InitData("meteo", "SpecificHumidity", FileSpecificHumidity_i,
                     FileSpecificHumidity_f, this->current_date,
                     SpecificHumidity_f);

    if (this->option_manage["rain"])
      this->InitData("meteo", "Rain", FileRain_i, FileRain_f,
                     this->current_date, Rain_f, false);

    if (this->option_manage["cloud_base_height"])
      this->InitData("meteo", "CloudBaseHeight", FileCloudBaseHeight_i,
                     FileCloudBaseHeight_f, this->current_date,
                     CloudBaseHeight_f);

    if (this->option_manage["cloud_top_height"])
      this->InitData("meteo", "CloudTopHeight", FileCloudTopHeight_i,
                     FileCloudTopHeight_f, this->current_date,
                     CloudTopHeight_f);

    if (this->option_manage["scavenging_below_cloud_coefficient"]
        && this->scavenging_below_cloud_model != "none")
      {
        this->scavenging_rain_threshold = 0.;
        if (this->scavenging_below_cloud_model == "microphysical")
          InitScavengingCoefficient(Temperature_f, Pressure_f,
                                    CloudBaseHeight_f, Rain_f,
                                    ScavengingBelowCloudCoefficient_f);
        else if (this->scavenging_below_cloud_model == "belot")
          InitScavengingCoefficient(Rain_f, belot_below_cloud_constant_a,
                                    belot_below_cloud_constant_b,
                                    ScavengingBelowCloudCoefficient_f);
        else if (this->scavenging_below_cloud_model == "constant")
          InitScavengingCoefficient(ScavengingBelowCloudCoefficient_f);
      }

    if (this->option_manage["scavenging_in_cloud_coefficient"]
        && this->option_process["with_transport_in_cloud_scavenging"])
      {
        if (this->scavenging_in_cloud_model == "pudykiewicz")
          InitScavengingCoefficient(Temperature_f, Pressure_f,
                                    SpecificHumidity_f,
                                    ScavengingInCloudCoefficient_f);
        if (this->scavenging_in_cloud_model == "belot")
          InitScavengingCoefficient(Rain_f, belot_in_cloud_constant_a,
                                    belot_in_cloud_constant_b,
                                    ScavengingInCloudCoefficient_f);
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
      if (!option_isotropic_diffusion)
        {
          ZonalDiffusionCoefficient_i.Fill(horizontal_diffusion);
          MeridionalDiffusionCoefficient_i.Fill(horizontal_diffusion);
        }
      else
        {
          LinearInterpolationRegular(VerticalDiffusionCoefficient_f,
                                     ZonalDiffusionCoefficient_i);
          ZonalDiffusionCoefficient_i.ThresholdMin(0.);
          LinearInterpolationRegular(VerticalDiffusionCoefficient_f,
                                     MeridionalDiffusionCoefficient_i);
          MeridionalDiffusionCoefficient_i.ThresholdMin(0.);
        }

    if (this->option_manage["horizontal_diffusion"] && !option_cartesian)
      {
        TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
                                ZonalDiffusionCoefficient_i);
        TransformMeridionalDiffusion(MeridionalDiffusionCoefficient_i);
      }

    /*** Boundary conditions ***/

    if (this->option_manage["boundary_condition"])
      for (int i = 0; i < Ns_bc; i++)
        {
          string filename
            = this->input_files["boundary_condition"](species_list_bc[i]);
          Date date = this->input_files["boundary_condition"].GetDateMin();
          T Delta_t = this->input_files["boundary_condition"].GetDelta_t();

          this->InitData(find_replace(filename, "&c", "z"), date, Delta_t,
                         FileBoundaryCondition_z_i, FileBoundaryCondition_z_f,
                         this->current_date, i, BoundaryCondition_z_i);
          this->InitData(find_replace(filename, "&c", "y"), date, Delta_t,
                         FileBoundaryCondition_y_i, FileBoundaryCondition_y_f,
                         this->current_date, i, BoundaryCondition_y_i);
          this->InitData(find_replace(filename, "&c", "x"), date, Delta_t,
                         FileBoundaryCondition_x_i, FileBoundaryCondition_x_f,
                         this->current_date, i, BoundaryCondition_x_i);
        }

    /*** Deposition velocities ***/

    if (this->option_manage["deposition_velocity"])
      for (int i = 0; i < Ns_dep; i++)
        this->InitData("deposition_velocity", species_list_dep[i],
                       FileDepositionVelocity_i, FileDepositionVelocity_f,
                       this->current_date, i, DepositionVelocity_f);

    /*** Surface emissions ***/

    if (this->option_manage["surface_emission"])
      for (int i = 0; i < Ns_surf_emis; i++)
        this->InitData("surface_emission", species_list_surf_emis[i],
                       FileSurfaceEmission_i, FileSurfaceEmission_f,
                       this->current_date, i, SurfaceEmission_f);

    /*** Additional surface emissions ***/

    if (this->option_manage["additional_surface_emission"])
      for (int i = 0; i < Ns_add_surf_emis; i++)
        this->InitData("additional_surface_emission",
                       species_list_add_surf_emis[i],
                       FileAdditionalSurfaceEmission_i,
                       FileAdditionalSurfaceEmission_f,
                       this->current_date, i, AdditionalSurfaceEmission_f);

    if (this->option_manage["additional_surface_emission"])
      for (int s = 0; s < Ns_add_surf_emis; s++)
        {
          int gs = SurfaceEmissionIndex(AdditionalSurfaceEmissionName(s));
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              SurfaceEmission_f(gs, j, i)
                += AdditionalSurfaceEmission_f(s, j, i);
        }

    /*** Volume emissions ***/

    if (this->option_manage["volume_emission"])
      for (int i = 0; i < Ns_vol_emis; i++)
        this->InitData("volume_emission", species_list_vol_emis[i],
                       FileVolumeEmission_i, FileVolumeEmission_f,
                       this->current_date, i, VolumeEmission_f);

    this->data_date = this->current_date;
    this->data_gone_through_initstep = false;
    this->data_gone_through_forward = false;
  }

  //! Returns the DryDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  Data<T, 3>& Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetDryDepositionFlux()
  {
    return this->DryDepositionFlux;
  }

  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  Data<T, 3>& Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetWetDepositionFlux()
  {
    return this->WetDepositionFlux;
  }

  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  Data<T, 3>& Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetInCloudWetDepositionFlux()
  {
    return this->InCloudWetDepositionFlux;
  }
  //! Computes concentration at a given point.
  /*!
    \return The concentration at the point.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  T Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetConcentration(int species, T z, T y, T x)
  {
    T concentration;
    Data<T, 3>
      Concentration_tmp(&this->Concentration
                        (species, 0, 0, 0),
                        shape(this->Nz, this->Ny, this->Nx));
    Concentration_tmp.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    Array<T, 1> Coord(3);
    Coord(0) = z;
    Coord(1) = y;
    Coord(2) = x;
    LinearInterpolationPoint(Concentration_tmp, Coord, concentration);
    return concentration;

  }


  /*! Computes indices of the cell containing a given point in Eulerian grid
    regular along x and y.
  */
  /*!
    \param lon longitude of the point (degrees).
    \param lat latitude of the point (degrees).
    \param height of the point (meters).
    \param k (output) cell index along z.
    \param j (output) cell index along y.
    \param i (output) cell index along x.
  */
  template<class T, class ClassAdvection, class ClassDiffusion>
  void Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  ::GetCellIndices(T lon, T lat, T height,
                   int& index_z, int& index_y, int& index_x)
  {
    index_x = int((lon - this->x_min + this->Delta_x / 2.) / this->Delta_x);
    index_y = int((lat - this->y_min + this->Delta_y / 2.) / this->Delta_y);
    index_z = 0;
    for (int k = 0; k < this->Nz; k++)
      {
        index_z = k;
        if (this->GridZ3D_interf(k + 1) > height)
          break;
      }
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DTRANSPORT_CXX
#endif
