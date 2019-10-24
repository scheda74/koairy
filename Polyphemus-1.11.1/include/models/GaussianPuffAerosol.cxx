// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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

// This file is part of a Gaussian puff model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFAEROSOL_CXX


#include "GaussianPuffAerosol.hxx"
#define PI 3.141592653589793115997963468544185161590576171875

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration filename.
  */
  template<class T, class ClassChemistry>
  GaussianPuffAerosol<T, ClassChemistry>
  ::GaussianPuffAerosol(string config_file):
    GaussianPuffChemistry<T, ClassChemistry>(config_file)
  {
  }


  //! Destructor.
  template<class T, class ClassChemistry>
  GaussianPuffAerosol<T, ClassChemistry>::~GaussianPuffAerosol()
  {
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::ReadConfiguration()
  {
    GaussianPuffChemistry<T, ClassChemistry>::ReadConfiguration();

    /*** Species and bins ***/

    // Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    // Section "[aerosol_species]" contains all aerosol species names.
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());

   // Reads mass density aerosol
    map<string, double> parameter;
    map<string, double>::iterator iter;
    string species;

    ConfigStream config_species(this->GetSpeciesFile());
    config_species.SetSection("[mass_density_aer]");
    Mass_Density_aer.resize(this->Ns_aer);
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    
    for (int i = 0; i < this->Ns_aer; i++)
      {
	iter = parameter.find(this->species_list_aer[i]);
	if(iter != parameter.end())
	  Mass_Density_aer[i] = iter->second;
	else
	  throw string("Module: no mass density for ")
	    + string("aerosol species \"") + this->species_list_aer[i]
	    + string("\". Please provide one.");
      }
    parameter.clear();

    // Reads bin bounds.
    this->config.SetSection("[domain]");
    this->config.Find("Bin_bounds");
    bin_list = split(this->config.GetLine());
    this->Nbin_aer = int(bin_list.size()) - 1;
    this->Nsize_section_aer = int(bin_list.size()) - 1;

    // vector<string> species_bin;
    // string species;
    // int bin_index
    //SZ
    if (this->option_process["with_external_composition"])
    {//read aerosol groups
      // cout<<"with_external_composition"<<endl;
      species_stream.SetSection("[aerosol_groups]");
      while (!species_stream.IsEmpty())
	this->groups_list_aer.push_back(species_stream.GetElement());

      this->Ngroup_aer = int(this->groups_list_aer.size());
      //read relations between aerosol species and groups
      species_stream.SetSection("[aerosol_species_group_relations]");
      map<string, string> parameter2;
      map<string, string>::iterator iter2;
      while (!species_stream.IsEmpty())
      {
	species = species_stream.GetElement();
	parameter2[species] = species_stream.GetElement();
      }
      this->aerosol_species_group_relation.resize(this->Ns_aer);
      for (int i = 0; i < this->Ns_aer; i++)
      {
	iter2 = parameter2.find(this->species_list_aer[i]);
	if(iter2 != parameter2.end())
	  {
	    int group_index = this->GetGroupIndex_aer(iter2->second);
	    this->aerosol_species_group_relation(i) = group_index;
	  }
	else
	  this->aerosol_species_group_relation(i) = -1;
      }
      parameter2.clear();
    }
    else//in case of internal mixing
    {
      this->groups_list_aer.push_back("None");
      this->Ngroup_aer=int(this->groups_list_aer.size());
      this->aerosol_species_group_relation.resize(this->Ns_aer);
      for (int i = 0; i < this->Ns_aer; i++)
      {
	this->aerosol_species_group_relation(i)=0;
      }
    }

    this->Ncomposition_aer=1;
    this->Nbin_aer=this->Nsize_section_aer;
    this->composition_bounds.Resize(1,1,2);
    this->composition_bounds(0,0,0)=0.0;
    this->composition_bounds(0,0,1)=1.0;

    //For number computation.
    this->config.SetSection("[options]");
    this->config.PeekValue("With_number_concentration",
			   this->option_process["with_number_concentration"]);

    /*** Output options ***/

    this->config.SetSection("[output]");

    if (this->option_process["with_chemistry"])
      if (this->option_process["with_output_plume_mass"])
        {
          this->config.PeekValue("Directory_output_puff_mass", directory_puff_mass);
          this->config.PeekValue("Delta_t_output", delta_t_output);
        }

    this->config.SetSection("[output]");
    if (this->option_process["with_chemistry"]){
      if (this->option_process["with_output_plume_concentration"]){
	this->config.PeekValue("Directory_output_puff_concentration", directory_puff_concentration);
	this->config.PeekValue("Directory_output_puff_background", directory_puff_background_concentration);
	this->config.PeekValue("Delta_t_output", delta_t_output);

      }
    }

    this->config.SetSection("[output]");
    if (this->option_process["with_chemistry"]){
      if (this->option_process["with_output_plume_coordinate"]){
	this->config.PeekValue("Directory_output_puff_coordinate", directory_puff_coordinate);
	this->config.PeekValue("Delta_t_output", delta_t_output);
      }
    }

    this->config.SetSection("[output]");
    if (this->option_process["with_chemistry"]){
      if (this->option_process["with_number_concentration"])
	if (this->option_process["with_output_plume_number"]){
	  this->config.PeekValue("Directory_output_puff_number", directory_puff_number);
	  this->config.PeekValue("Directory_output_puff_number_concentration", directory_puff_number_concentration);
	  this->config.PeekValue("Delta_t_output", delta_t_output);
	}
    }
    this->config.SetSection("[output]");
    this->config.PeekValue("Save_plume", this->option_process["save_plume"]);
  }


  //! Allocates memory.
  /*! Allocates the grids and the concentration Data for aerosol species.
   */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::Allocate()
  {
    GaussianPuffChemistry<T, ClassChemistry>::Allocate();

    BinBound_aer.resize(this->Nbin_aer + 1);
    // Reads bin bounds in micrometers and converts it to meters.
    for (int i = 0; i < this->Nbin_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]);
    wet_diameter_aer.resize(this->Nbin_aer);

    /*** Additional grids ***/

    GridS5D_aer = RegularGrid<T>(this->Ns_aer);
    GridB5D_aer = RegularGrid<T>(this->Nbin_aer);
    GridZ5D = RegularGrid<T>(this->Nz);
    GridY5D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridX5D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);

    GridZ4D = RegularGrid<T>(this->Nz);
    GridY4D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridX4D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);

    GridZ4D.SetVariable(1);
    GridZ4D.SetDuplicate(false);

    GridY4D.SetVariable(2);
    GridY4D.SetDuplicate(false);

    GridX4D.SetVariable(3);
    GridX4D.SetDuplicate(false);
    GridS4D = RegularGrid<T>(this->Ns);

    GridS4D_number = RegularGrid<T>(this->Nbin_aer);

    /*** State ***/

    this->Concentration_aer.Resize(GridS5D_aer, GridB5D_aer,
                                   GridZ5D, GridY5D, GridX5D);

    if (this->option_process["with_number_concentration"])
      {
	this->NumberConcentration_aer.Resize(GridS4D_number,
					     GridZ4D, GridY4D, GridX4D);
      }

  }


  // //! Model initialization.
  // /*! It reads the configuration.
  //  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::Init()
  {
    GaussianPuffChemistry<T, ClassChemistry>::Init();
  }

  //! Source initialization.
  /*! It sets all sources from a text file. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species: species name,
    <li> Date beg: the release date,
    <li> Quantity: the released quantity (mass unit),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees),
    </ul>
    \param file_puff The file containing all the sources parameters.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::InitSource(string file_puff, T delta_t_puff)
  {
    // List of emissions.
    this->PointEmissionManager = new BasePointEmission<T>();
    this->PointEmissionManager->Init(file_puff, this->species_list,
                                     this->species_list_aer, this->Nbin_aer);
  }


  //! Puffs initialization.
  /*!
    For each source, it creates the corresponding puff (in case of a puff
    source) or series of puffs (for a continuous source). The continuous
    source is discretized with the model time step.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::InitPuff()
  {
    int Nemis = this->PointEmissionManager->GetNumberEmission();
    for (int emission = 0; emission < Nemis; emission++)
      {
        vector<int> emitted_species_index =
          this->PointEmissionManager->GetEmittedSpeciesIndex(emission);
        int Ns_emitted = emitted_species_index.size();
        int s, Npoint;
        T time_puff;
        Array<T, 2> point_emission;
        T velocity = 0.;
        T temperature = 0.;
        T diameter = 0.;
        T width = 0.;
        T length = 0.;
    	T source_water = 0.;
        T volume_prev = 0.;
        string source_id;
        bool is_volume_source = false;

        // Getting source plume rise parameters.
        if (this->PointEmissionManager->HasPlumeRise(emission))
          this->PointEmissionManager->
            GetPlumeRiseParam(velocity, temperature, diameter, emission);

        if (this->PointEmissionManager->IsVolumeSource(emission))
          this->PointEmissionManager->
            GetVolumeSource(width, length, emission);

        is_volume_source = this->PointEmissionManager->
          IsVolumeSource(emission);

        source_id = this->PointEmissionManager->
          GetEmissionSourceId(emission);
        source_water = this->PointEmissionManager->
          GetEmissionSourceWater(emission);
        time_puff = this->current_date.GetSecondsFrom(this->Date_min);
        int time_step = int(time_puff / this->Delta_t);
        int Nt_puff = max(int(this->Delta_t_puff / this->Delta_t), 1);
        Date puff_next_date = this->current_date;
        puff_next_date.AddSeconds(this->Delta_t * Nt_puff);
        // If the source is emitting at time step t.
        if (this->PointEmissionManager->
            IsEmitting(this->current_date, this->next_date, emission)
            && time_step % Nt_puff == 0)
          {
            vector<T> quantity;

            // Getting source emission parameters for all species.
            for (int species = 0; species < this->Ns; species++)
              {
                vector<int>::iterator iter;
                iter = find(emitted_species_index.begin(),
                            emitted_species_index.end(), species);
                if (iter == emitted_species_index.end())
                  quantity.push_back(0.);
                else
                  {
                    for (int s = 0; s < Ns_emitted; s++)
                      if (emitted_species_index[s] == species)
                        this->PointEmissionManager->
                          GetEmission(this->current_date,
                                      puff_next_date, s, emission,
                                      point_emission);
                    quantity.push_back(point_emission(0, 3));
                  }
              }

            // Getting source emission parameters for aerosol species.
            vector<int> emitted_species_index_aer =
              this->PointEmissionManager->GetEmittedSpeciesIndex_aer(emission);
            int Ns_emitted_aer = emitted_species_index_aer.size();
            vector<map <string, string> > emitted_species_list_aer_bin =
              this->PointEmissionManager->GetEmittedAerosolSpecies(emission);

            Array<T, 2> quantity_aer;
            int species_aer, bin_aer;
     	    string species_name_tmp;
            quantity_aer.resize(this->Ns_aer, this->Nbin_aer);
            quantity_aer = 0.0;

            for (s = 0; s < Ns_emitted_aer; s++)
              {
                species_name_tmp = emitted_species_list_aer_bin[s]["species"];
                species_aer = this->GetSpeciesIndex_aer(species_name_tmp);

                int bin_aer = convert<int>(emitted_species_list_aer_bin[s]["bin"]);
                T quantity_aer_bin;
                this->PointEmissionManager->
                  GetEmission_aer(this->current_date,
                                  puff_next_date, s, emission,
                                  quantity_aer_bin);
                quantity_aer(species_aer, bin_aer) = quantity_aer_bin;
              }

	    //Getting source emissions for aerosol in number.
	    Array<T, 1> quantity_number_aer;
	    if (this->option_process["with_number_concentration"])
	      {
		ComputePuffNumberEmission_aer(emission, quantity_number_aer);
	      }
	    else
              {
                quantity_number_aer.resize(this->Nbin_aer);
                quantity_number_aer = 0.0;
              }
            
            // Creating one puff for each location with all species.
            Npoint = point_emission.extent(0);

            for (int index = 0; index < Npoint; index++)
              {
                Puff<T>* puff
                  = new PuffAerosol<T>(time_puff, velocity,
                                       temperature, diameter, width, length,
                                       quantity,
                                       point_emission(index, 0),
                                       point_emission(index, 1),
                                       point_emission(index, 2),
                                       source_water,
                                       volume_prev,
                                       this->species_list,
                                       this->photolysis_reaction_list,
                                       this->species_list_aer,
                                       this->bin_list,
                                       quantity_aer,      
                                       quantity_number_aer,
                                       is_volume_source, source_id);
                puff->InitPuff();
                this->PuffList.push_back(puff);
                this->Npuff++;
              }
          }
      } // for
    this->puff_index = 0;
    this->current_puff = this->PuffList.begin();
  }


  ////////////////////////////////////
  // BACKGROUND METEOROLOGICAL DATA //
  ///////////////////////////////////


  //! Initializes meteorological conditions.
  /*! It sets the meteorological data from a configuration file.  The
    situation is described in a dedicated section "[situation]" in which one
    finds the following entries:
    <ul>
    <li> Attenuation:
    <li> Pressure:
    <li> Specific_humidity:
    <li> Liquid_water_content:
    >li> pH:
    </ul>
    \param meteo ConfigStream instance through which all entries may be read to set the meteorological situation.
    \param show_meteo indicates whether the meteorological data is to be
    displayed on screen.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {
    GaussianPuffTransport<T>::InitMeteo(meteo, show_meteo);

    /*** For chemistry ***/

    if (this->option_process["with_chemistry"])
      {
        meteo.PeekValue("Attenuation", this->background_attenuation_);
        meteo.PeekValue("Pressure", this->background_pressure_);
        meteo.PeekValue("Specific_humidity", this->background_specific_humidity_);
        meteo.PeekValue("Liquid_water_content", background_liquid_water_content_);
        meteo.PeekValue("pH", background_ph_);
        InitChemistry(meteo);
      }
  }


  /*! Initializes scavenging coefficients.
    \param meteo ConfigStream instance through which scavenging coefficients
    may be read.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::InitScavenging(ConfigStream& meteo)
  {
    GaussianPuffTransport<T>::InitScavenging(meteo);
    scavenging_coefficient_aer.resize(this->Ns_aer, this->Nbin_aer);

    meteo.Find("Scavenging_coefficient_aerosol");
    vector<string> scav_coefficient =  split(meteo.GetLine());
    vector<string>::iterator iter;
    for (int i = 0; i < this->Ns_aer; i++)
      for (int b = 0; b < this->Nbin_aer; b++)
        {
          iter = find(scav_coefficient.begin(), scav_coefficient.end(),
                      this->species_list_aer[i]);
          if (iter++ == scav_coefficient.end() ||
              iter == scav_coefficient.end())
            throw string("Unable to find scavenging coefficient for")
              + string(" species \"") + this->species_list_aer[i] + "\".";
          scavenging_coefficient_aer(i, b) = to_num<T>(*iter);
        }
  }


  /*! Initializes deposition velocities.
    \param meteo ConfigStream instance through which deposition velocities
    may be read.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::InitDeposition(ConfigStream& meteo)
  {
    GaussianPuffTransport<T>::InitDeposition(meteo);
    deposition_velocity_aer.resize(this->Ns_aer, this->Nbin_aer);

    meteo.Find("Deposition_velocity_aerosol");
    vector<string> dep_velocity = split(meteo.GetLine());
    vector<string>::iterator iter;
    for (int  i = 0; i < this->Ns_aer; i++)
      for (int  b = 0; i < this->Nbin_aer; b++)
        {
          iter = find(dep_velocity.begin(), dep_velocity.end(),
                      this->species_list_aer[i]);
          if (iter++ == dep_velocity.end() ||
              iter == dep_velocity.end())
            throw string("Unable to find deposition velocity for")
              + string(" species \"") + this->species_list_aer[i] + "\".";
          deposition_velocity_aer(i, b) = to_num<T>(*iter);
        }
  }

  /*! Initializes photolysis parameters.
    \param Nr The number of photolysis reactions.
    \param photolysis_list  List of species with photolysis reactions.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::InitPhotolysis(int Nr, vector<string> photolysis_list)
  {
    GaussianPuffChemistry<T, ClassChemistry>
      ::InitPhotolysis(Nr, photolysis_list);

    this->Nr_photolysis = Nr;
    if (this->option_process["with_chemistry"]
        && this->option_process["with_photolysis"])
      {
        background_concentration_aer_.resize(this->Ns_aer, this->Nbin_aer);
        background_concentration_aer_ = 0.;

	background_concentration_number_.resize(this->Nbin_aer);
        background_concentration_number_ = 0.;


        this->Chemistry_.Init(*this);
      }
  }

  /*! Initializes background concentrations and photolysis parameters.
    \param meteo ConfigStream instance through which background concentrations
    and photolysis rates may be read.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::InitChemistry(ConfigStream& meteo)
  {
    // Longitude and latitude (for chemistry only).
    this->config.SetSection("[domain]");
    this->config.PeekValue("Longitude", this->background_longitude);
    this->config.PeekValue("Latitude", this->background_latitude);

    vector<string>::iterator iter;

    // Reads the list of species with photolysis.
    if (this->option_process["with_photolysis"])
      {
        ConfigStream species_stream(this->file_species);
        species_stream.SetSection("[photolysis]");
        while (!species_stream.IsEmpty())
          this->photolysis_reaction_list.push_back(species_stream.GetElement());
        this->Nr_photolysis = int(this->photolysis_reaction_list.size());
        this->photolysis_rate_.resize(this->Nr_photolysis);

        meteo.Find("Photolysis_rate");
        vector<string> photolysis =  split(meteo.GetLine());
        for (int i = 0; i < this->Nr_photolysis; i++)
          {
            iter = find(photolysis.begin(), photolysis.end(),
                        this->photolysis_reaction_list[i]);
            if (iter++ == photolysis.end() || iter == photolysis.end())
              throw string("Unable to find photolysis rate for")
                + string(" species \"") + this->photolysis_reaction_list[i] + "\".";
            this->photolysis_rate_(i) = to_num<T>(*iter);
          }
      }

    this->background_concentration_.resize(this->Ns);
    this->background_concentration_ = 0.;

    meteo.Find("Background_concentration");
    vector<string> concentration =  split(meteo.GetLine());
    for (iter = concentration.begin(); iter != concentration.end();
         iter++)
      for (int i = 0; i < this->Ns; i++)
        if (*iter == this->species_list[i])
          {
            iter++;
            this->background_concentration_(i) = to_num<T>(*iter);
          }

    background_concentration_aer_.resize(this->Ns_aer, this->Nbin_aer);
    background_concentration_aer_ = 0.;

    meteo.Find("Background_concentration_aer");
    vector<string> concentration_aer =  split(meteo.GetLine());
    for (iter = concentration_aer.begin(); iter != concentration_aer.end();
         iter++)
      for (int i = 0; i < this->Ns_aer; i++)
        if (*iter == this->species_list_aer[i])
          {
            iter++;
            for (int b = 0; b < this->Nbin_aer; b++)
              this->background_concentration_aer_(i, b) = to_num<T>(*iter);
          }

    this->Chemistry_.Init(*this);

  }


  //! Sets the current meteorological data to puff data, if available.
  /*! It sets the current meteorological data to puff data, if available,
    and to background meteorological data otherwise.
    \param puff the puff.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffCurrentMeteo(T interaction_coefficient, Puff<T>* puff)
  {
    GaussianPuffChemistry<T, ClassChemistry>::SetPuffCurrentMeteo(interaction_coefficient, puff);

    if (this->option_process["with_chemistry"])
      {
   
        if (puff->HasMeteo()){
            // puff->GetAdditionalMeteo(liquid_water_content_);
	  liquid_water_content_ = ComputePuffAdditionalLWC(this->temperature_,
							   this->pressure_,
							   this->specific_humidity_);

	}
        else
            liquid_water_content_ = background_liquid_water_content_;
      }
  }


  //! Sets the current meteorological data to puff data, if available.
  /*! It sets the current meteorological data to puff data, if available,
    and to background meteorological data otherwise.
    \param puff the puff.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetCurrentMeteo(Puff<T>* puff)
  {
    GaussianPuffChemistry<T, ClassChemistry>::SetCurrentMeteo(puff);

    if (this->option_process["with_chemistry"])
      {

        if (puff->HasMeteo())
          puff->GetAdditionalMeteo(liquid_water_content_);
        else
          liquid_water_content_ = background_liquid_water_content_;
      }
  }

  ////////////////////////////////////////
  // ACCESS METHODS FOR PUFF ATTRIBUTES //
  ///////////////////////////////////////


  //! Sets values of additional meteorological data for a given puff.
  /*! It sets the additionnal meteorological data (needed for chemistry)
    for a given puff.
    \param index: puff index.
    \param liquid_water_content: the liquid water content.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffAdditionalMeteo_aer(int index,
                               T liquid_water_content)
  {
    this->SetCurrentPuff(index);
    (*this->current_puff)->SetAdditionalMeteo(liquid_water_content);
  }


  //! Sets values of species background concentrations for a given puff.
  /*! It sets the background concentrations (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffBackgroundConcentration_aer(int index, Array<T, 2> concentration_aer)
  {
    this->SetCurrentPuff(index);

    int size_ns_aer = concentration_aer.size() / this->Nbin_aer;
    if (size_ns_aer != this->Ns_aer)
      throw string("The number of species_aer is set to ")
        + to_str<int>(this->Ns_aer) + " but the array is of size "
        + to_str<int>(size_ns_aer);

    for (int s = 0; s < this->Ns_aer; s++)
      for (int b = 0; b < this->Nbin_aer; b++)
        (*this->current_puff)->SetBackgroundConcentration(concentration_aer(s, b), s, b);
  }

  //! Gets values of species background concentrations for a given puff.
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::GetPuffBackgroundConcentration_aer(int index, int s, int b)
  {
    this->SetCurrentPuff(index);
    
    return (*this->current_puff)->GetBackgroundConcentration(s, b);
  }


  //! Sets values of species background concentrations for a given puff.
  /*! It sets the background concentrations (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffBackgroundConcentration_number(int index, Array<T, 1> concentration_number)
  {
    this->SetCurrentPuff(index);
    
    int size_ns_number = concentration_number.size();
    if (size_ns_number != this->Nbin_aer)
      throw string("The number of species_aer is set to ")
	+ to_str<int>(this->Nbin_aer) + " but the array is of size "
	+ to_str<int>(size_ns_number);

    for (int b = 0; b < this->Nbin_aer; b++) 
      (*this->current_puff)->SetBackgroundNumberConcentration(concentration_number(b), b);
  }

  //! Gets values of species background concentrations for a given puff.
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::GetPuffBackgroundConcentration_number(int index, int b)
  {
    this->SetCurrentPuff(index);
    
    return (*this->current_puff)->GetBackgroundNumberConcentration(b);
  }


  //! Sets species quantity for a given puff.
  /*! It sets the species quantity (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffQuantity_aer(int index, Array<T, 2> quantity_aer)
  {
    this->SetCurrentPuff(index);
    for(int s = 0; s < this->Ns_aer; s++)
      for (int b = 0; b < this->Nbin_aer; b++) 
	{
	  (*this->current_puff)->SetQuantity(quantity_aer(s, b), s, b);
	}
  }

  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::SetPuffQuantity_number(int index, Array<T, 1> quantity_number)
  {
    this->SetCurrentPuff(index);
    for (int b = 0; b < this->Nbin_aer; b++) 
      {
	(*this->current_puff)->SetNumberQuantity(quantity_number(b), b);
      }
  }



  //! Compute density
  /*
    return density (ug.um-3)
  */
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputeDensity(Data<T, 1> Conc_aer_tmp,
		   vector<float> Rho_species, float TotalMass, int Ns)
  {
    int s;
    float rho, subrho;

    rho = 0.0;
    subrho = 0.0;

    for (s = 0; s < Ns; s++)
      {
	subrho += Conc_aer_tmp(s) / Rho_species[s];
      }
    if (TotalMass == 0. or subrho == 0.)
      rho = 1.;
    else
      rho = 1.e-6 * TotalMass/subrho;
    
   
    return rho;
	
  }


  //! Compute number emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffNumberEmission_aer(int emission, Array<T, 1>& puff_number_emissions)
  {
    puff_number_emissions.resize(this->Nbin_aer);
    puff_number_emissions = 0.0;

    vector<int> emitted_species_index_number_aer =
              this->PointEmissionManager->GetEmittedSpeciesIndex_aer(emission);
    int Ns_emitted_number_aer = emitted_species_index_number_aer.size();

    vector<map <string, string> > emitted_species_list_number_aer_bin = 
              this->PointEmissionManager->GetEmittedAerosolSpecies(emission);
    int bin_aer;
    int species_aer;
    string species_name_tmp;

    T quantity_aer_number_bin_tmp;
    Array<T, 2> quantity_aer_number_tmp;
    quantity_aer_number_tmp.resize(this->Ns_aer, this->Nbin_aer);
    quantity_aer_number_tmp = 0.0;
 
    int i, j, k;
    float TotalMass; // ug.m-3
    float Rho_aer; 
    float MeanDiameter; //= 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um

    vector<float> Rho_species; // g.cm-3 -> 10-6 ug.um-3
	
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_aer);
    Conc_aer_tmp.SetZero();
	
    int Nt_puff = max(int(this->Delta_t_puff / this->Delta_t), 1);
    Date puff_next_date = this->current_date;
    puff_next_date.AddSeconds(this->Delta_t * Nt_puff);
    for (int s = 0; s < Ns_emitted_number_aer; s++)
      {

	species_name_tmp = emitted_species_list_number_aer_bin[s]["species"];
	species_aer = this->GetSpeciesIndex_aer(species_name_tmp);
	
	bin_aer = convert<int>(emitted_species_list_number_aer_bin[s]["bin"]);
	
	this->PointEmissionManager->
                      GetEmission_aer(this->current_date,
                                      puff_next_date, s, emission,
                                      quantity_aer_number_bin_tmp); 

	quantity_aer_number_tmp(species_aer, bin_aer) = quantity_aer_number_bin_tmp;
      }
    
    for (int b = 0; b < this->Nbin_aer; b++){
      TotalMass = 0.0;
      Rho_aer = 0.0;
      for (int s = 0; s < this->Ns_aer - 1; s++){
      	TotalMass = TotalMass + quantity_aer_number_tmp(s, b);
      	Conc_aer_tmp(s) = quantity_aer_number_tmp(s, b);
      }
      MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));
      // puff_number_emissions(b) = 1. / PI * 6. 
      // 	/ (MeanDiameter*MeanDiameter*MeanDiameter) / (1.e-6);
      Rho_aer = ComputeDensity(Conc_aer_tmp, Mass_Density_aer, TotalMass, this->Ns_aer - 1);

      puff_number_emissions(b) = TotalMass / Rho_aer / PI*6.
      	/(MeanDiameter*MeanDiameter*MeanDiameter);
    }

  }

  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::GetMassDensity_aer(Array<T,1>& MassDensity_aer)
  {
    MassDensity_aer = 0.;
    for (int s = 0; s < this->Nbin_aer; s++)
      MassDensity_aer(s) = Mass_Density_aer[s];
  }

 //! Compute number emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputeNumberPuff(int b, Array<T,2> concentration_list_aer)
  {
    float MeanDiameter; //= 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
    float Rho_aer; 

    Data<T, 1> Conc_aer_bin_tmp;
    int s;
    Conc_aer_bin_tmp.Resize(this->Ns_aer-1);
    Conc_aer_bin_tmp.SetZero();
    T tot_mass = 0.0;
    T number_concentration;
    MeanDiameter= 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); 

    for (s=0; s<this->Ns_aer-1; s++)
      {
	tot_mass += concentration_list_aer(s,b);
	Conc_aer_bin_tmp(s) = concentration_list_aer(s,b);
      }
    Rho_aer = ComputeDensity(Conc_aer_bin_tmp, Mass_Density_aer,
			     tot_mass, this->Ns_aer - 1);
    number_concentration = tot_mass / Rho_aer / PI*6.
      /(MeanDiameter*MeanDiameter*MeanDiameter);

    return number_concentration;

  }

  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputeOutputNumber(int b,
  			Array<T, 2> concentration_list_aer,
  			Array<T, 1> concentration_list_number,
  			Array<T, 2> background_concentration_aer,
  			Array<T, 1> background_concentration_number,
  			vector<float> Rho_species)
  {
    //Compute Rho_tot, Rho_back
    T rho_tot, rho_background;
    T conc_tot, conc_background;
    T Rho_0;
    T m_0;
    Rho_0 = 1.e-6;
    m_0 = 1.;
    conc_tot = 0.;
    conc_background = 0.;

    Data<T, 1> conc_tot_tmp(this->Ns_aer);
    Data<T, 1> conc_background_tmp(this->Ns_aer);
    conc_tot_tmp.SetZero();
    conc_background_tmp.SetZero();
    
    for (int s = 0; s < this->Ns_aer-1; s++)
      {
	conc_tot += concentration_list_aer(s, b);
	conc_tot_tmp(s) = concentration_list_aer(s, b);

	conc_background += background_concentration_aer(s, b);
	conc_background_tmp(s) = background_concentration_aer(s, b);
      }

    rho_tot = ComputeDensity(conc_tot_tmp, Rho_species, conc_tot, this->Ns_aer-1);
    rho_background = ComputeDensity(conc_background_tmp, Rho_species, conc_background, this->Ns_aer-1);
   
    if (background_concentration_number(b) == 0. && concentration_list_number(b) == 0.)
      return 0.;

    //Compute Dn
    T delta_number;
    T delta_ratio, background_ratio;   

    if (conc_tot > 0. && conc_background > 0.)
      {
	delta_ratio = (conc_tot - conc_background) / conc_tot;
	background_ratio = conc_background / conc_tot;

	delta_number = (concentration_list_number(b) * (rho_tot / conc_tot) 
			- background_ratio * background_concentration_number(b) 
			* (rho_background / conc_background));
	delta_number *= ((m_0 / Rho_0) / delta_ratio);
      }
    else if (conc_tot > 0. && conc_background == 0.)
      delta_number = concentration_list_number(b) * (rho_tot / conc_tot)
	* (m_0 / Rho_0);
    else if (conc_tot == 0. && conc_background > 0.)
      delta_number = -1. * (background_concentration_number(b) * (rho_background / conc_background)
			    * (m_0 / Rho_0));
    else
      delta_number = 0.;

    if (isnan(rho_background) || isnan(rho_tot) || isnan(delta_number) || rho_tot > 0.)
      {
    	cout << "In compute Output number: " << endl;
    	cout << "rho_back: " << rho_background << " rho_tot: " << rho_tot << endl;
    	cout << "conc_back: " << conc_background << " conc_tot: " << conc_tot << endl;
    	cout << "Ntot: " << concentration_list_number(b) << " Nback: " << background_concentration_number(b);
    	cout << " Dn: " << delta_number << endl;
      }


    return delta_number;       
  } 

  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputeInputNumber(int b,
		       Array<T, 2> concentration_list_aer,
		       Array<T, 1> delta_concentration_list_number,
		       Array<T, 2> background_concentration_aer,
		       Array<T, 1> background_concentration_number,
		       vector<float> Rho_species,
		       int write_in)
  {
    T number_tot;

    Data<T, 1> conc_background_tmp(this->Ns_aer);
    Data<T, 1> conc_tot_tmp(this->Ns_aer);
    T conc_background, conc_tot;
    T Rho_0;
    T m_0;
    Rho_0 = 1.e-6;
    m_0 = 1.;

    conc_background = 0.;
    conc_background_tmp.SetZero();
    conc_tot = 0.;
    conc_tot_tmp.SetZero();

    if (delta_concentration_list_number(b) <= 0. && background_concentration_number(b) == 0.)
      {
	number_tot = ComputeNumberPuff(b, concentration_list_aer);
	return number_tot;
      }

    for (int s = 0; s < this->Ns_aer-1; s++)
      {	
	conc_tot += concentration_list_aer(s, b);
	conc_background_tmp(s) = background_concentration_aer(s, b);
	conc_tot_tmp(s) = concentration_list_aer(s, b);
	conc_background += background_concentration_aer(s, b);	
      }
    
    if (conc_tot == 0.)
      return 0.;
    else if (conc_tot == conc_background)
      return background_concentration_number(b);
   
    //Compute rho
    T rho_background, rho_tot;

    rho_background = ComputeDensity(conc_background_tmp, Rho_species, conc_background, this->Ns_aer-1);
    rho_tot = ComputeDensity(conc_tot_tmp, Rho_species, conc_tot, this->Ns_aer-1);

    //    delta_concentration_list_number(b) *= 1.e-6; //Delta N normalized with rho0 = 1.e-6; 
    T delta_ratio, background_ratio;

    if (conc_background == 0 && conc_tot > 0.)
      number_tot = delta_concentration_list_number(b) 
	* (Rho_0 / m_0) * (conc_tot / rho_tot);
    else
      {

	delta_ratio = (conc_tot - conc_background) / conc_tot;
	background_ratio = conc_background / conc_tot;
	number_tot = (delta_ratio * delta_concentration_list_number(b) * (Rho_0 / m_0)
		      + background_ratio * background_concentration_number(b) 
		      * (rho_background / conc_background));
	number_tot *= (conc_tot / rho_tot);
      }
    write_in = 1;
    if (number_tot < 0.)
      {
    	//Compute number corresponding to puff mass and diameter
    	number_tot = ComputeNumberPuff(b, concentration_list_aer);
      }

    if (isnan(rho_background) || isnan(rho_tot) || write_in == 1)
      {
    	cout << "Dmean bin: " << 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    	cout << " rho_back: " << rho_background << " rho_tot: " << rho_tot ;

    	cout << " conc_back: " << conc_background << " conc_tot: " << conc_tot;

    	cout << " dc: " << conc_tot - conc_background << endl;

    	cout << " Nback: " << background_concentration_number(b);

    	cout << " Dn: " << delta_concentration_list_number(b) <<  " number_out: " << number_tot << endl;

      }

   

    return number_tot;        

  }





//! Compute number emission from mass concentration, diameter and

  //! density

  /*!

    \ param b bin number

    \ return number concentration (m-3)

  */

  template<class T, class ClassChemistry>

  T GaussianPuffAerosol<T, ClassChemistry>

  ::ComputeNumberPuffDiameter(int b, T diameter, Array<T,2> concentration_list_aer)

  {

    float Rho_aer; 



    Data<T, 1> Conc_aer_bin_tmp;

    int s;

    Conc_aer_bin_tmp.Resize(this->Ns_aer-1);

    Conc_aer_bin_tmp.SetZero();

    T tot_mass = 0.0;

    T number_concentration;

    diameter *= 1.e6; //um



    for (s=0; s<this->Ns_aer-1; s++)

      {

	tot_mass += concentration_list_aer(s,b);

	Conc_aer_bin_tmp(s) = concentration_list_aer(s,b);

      }

    Rho_aer = ComputeDensity(Conc_aer_bin_tmp, Mass_Density_aer,

			     tot_mass, this->Ns_aer-1); //ug/um3

    number_concentration = tot_mass / Rho_aer / PI*6.

      /(diameter*diameter*diameter);



    // cout << " in compute number puff: " << endl;

    // cout << " mass: " << tot_mass << " Rho: " << Rho_aer << endl;



    return number_concentration;



  }




  ///////////////////////////
  // COMPUTATIONAL METHODS //
  //////////////////////////

  //!  Performs the chemistry.
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::Chemistry()
  {
    // Definition of arrays.
    Array<T, 2> quantity_list(this->Npuff, this->Ns);
    Array<T, 2> quantity_list_number(this->Npuff, this->Nbin_aer);
    Array<T, 3> quantity_list_aer(this->Npuff, this->Ns_aer, this->Nbin_aer);
    Array<T, 2> interaction_coefficient(this->Npuff, this->Npuff);
    Array<T, 1> concentration_list(this->Ns);
    Array<T, 1> concentration_list_number(this->Nbin_aer);
    Array<T, 2> concentration_list_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> overlap_volume(this->Ns);
    Array<T, 1> overlap_volume_number(this->Nbin_aer);
    Array<T, 2> overlap_volume_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> overlap_volume2(this->Ns);
    Array<T, 1> overlap_volume_number2(this->Nbin_aer);
    Array<T, 2> overlap_volume_aer2(this->Ns_aer, this->Nbin_aer);

    Array<T, 1> background_concentration(this->Ns);
    Array<T, 1> background_concentration_number(this->Nbin_aer);
    Array<T, 2> background_concentration_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> photolysis_rate(this->Nr_photolysis);
    Array<T, 1> source(this->Ns);
    Data<T, 1> species_quantity(this->Ns);
    Data<T, 1> species_quantity_number(this->Nbin_aer);
    Data<T, 2> species_quantity_aer(this->Ns_aer, this->Nbin_aer);

    // Initialization of arrays.
    source = 0.;
    quantity_list = 0.;
    quantity_list_number = 0.;
    quantity_list_aer = 0.;
    interaction_coefficient = 0.;
    species_quantity.SetZero();
    species_quantity_number.SetZero();
    species_quantity_aer.SetZero();
    concentration_list = 0.;
    concentration_list_number = 0.;
    concentration_list_aer = 0.;

    // Other variables.
    int s, r, alpha, beta, i;
    typename list<Puff<T>* >::iterator iter, iter2;
    Puff<T>* puff;
    int b;

    // gas/particle compensation ..
    Array<T, 1> aer_compensation;
    aer_compensation.resize(this->Ns_aer);

    T concentration_list_aer_tot;
    int species_comp, species_comp_aer;
    
    // For interactions.
    T tmp, tmp2, c_init;
    Array<vector<int>, 1 > PuffInteractionList(this->Npuff);
    vector<int> PuffList_tmp;
    int Ninteraction;
    vector<T> write_puff;
    T time_puff_tmp;

    //To compensate gas/particle mass modification ...
    ConfigStream config_species(this->GetSpeciesFile());
    map<string, string> gas_aer_interaction;
    map<string, string>::iterator iter_interact_tmp;
    string species_interact;
    Array<int, 1> species_index_aerosol_interact_tmp;

    species_index_aerosol_interact_tmp.resize(this->Ns_aer);
    config_species.SetSection("[gas_species_aerosol_interact]");
    while (!config_species.IsEmpty())
      {
    	species_interact = config_species.GetElement();
    	gas_aer_interaction[species_interact] = config_species.GetElement();
      }

    for (int i = 0; i < this->Ns_aer; i++)
      {
    	iter_interact_tmp = gas_aer_interaction.find(this->species_list_aer[i]);
    	if (iter_interact_tmp != gas_aer_interaction.end())
    	  {
    	    int gas_index = this->GetSpeciesIndex(iter_interact_tmp->second);
    	    species_index_aerosol_interact_tmp(i) = gas_index;
    	  }
    	else
    	  species_index_aerosol_interact_tmp(i) = -1;
      }
    gas_aer_interaction.clear();

    //Check which species is in AEC
    Array<int, 1> is_organic(this->Ns_aer);
    is_organic = 0;
    config_species.SetSection("[aec_species]");
    string aec_specie;
    int aec_specie_index;
    while(!config_species.IsEmpty())
      {
	aec_specie = config_species.GetElement();
	aec_specie_index = this->GetSpeciesIndex_aer(aec_specie);
	is_organic(aec_specie_index) = 1;
      }
    config_species.SetSection("[poa_species]");
    while(!config_species.IsEmpty())
      {
	aec_specie = config_species.GetElement();
	aec_specie_index = this->GetSpeciesIndex_aer(aec_specie);
	is_organic(aec_specie_index) = 1;
      }
    config_species.SetSection("[pankow_species]");
    while(!config_species.IsEmpty())
      {
	aec_specie = config_species.GetElement();
	aec_specie_index = this->GetSpeciesIndex_aer(aec_specie);
	is_organic(aec_specie_index) = 1;
      }

    // Getting quantities, computing interaction coefficients for all puffs.
    alpha = 0;
    for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
      {
        // List of quantities for all puffs.
        for (s = 0; s < this->Ns; s++)
          quantity_list(alpha, s) = (*iter)->GetQuantity(s);

        for (s = 0; s < this->Ns_aer; s++)
          for (b = 0; b < this->Nbin_aer; b++)
            quantity_list_aer(alpha, s, b) = (*iter)->GetQuantity(s, b);

	T tot_number_tmp;
	tot_number_tmp = 0.;
    	for (b = 0; b < this->Nbin_aer; b++){
    	  quantity_list_number(alpha, b) = (*iter)->GetNumberQuantity(b);	    	   
    	}

        // Interaction coefficient for all puff pairs.
        if (this->option_process["with_puff_interaction"])
          {
            beta = 0;
            PuffList_tmp.clear();
            for (iter2 = this->PuffList.begin();
                 iter2 != this->PuffList.end(); iter2++)
              {
                tmp = this->ComputeInteractionCoefficient(*iter, *iter2);
                if (tmp != 0.)
                  {
                    interaction_coefficient(alpha, beta) = tmp;
                    PuffList_tmp.push_back(beta);
                  }
                beta++;
              }
            PuffInteractionList(alpha) = PuffList_tmp;
          }
        else
          interaction_coefficient(alpha, alpha) =
            this->ComputeInteractionCoefficient(*iter, *iter);
        alpha++;
      }

    // Chemistry for all puffs.
    alpha = 0;
    for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
      {
        puff = *iter;
    	time_puff_tmp = puff->GetReleaseTime();
    	SetPuffCurrentMeteo(interaction_coefficient(alpha, alpha), puff);
    	puff->SetAdditionalMeteo(this->attenuation_,this->pressure_,
    					this->specific_humidity_);

    	puff->SetAdditionalMeteo(liquid_water_content_);

//        SetCurrentMeteo(puff);
        overlap_volume = 0.;
        PuffList_tmp = PuffInteractionList(alpha);
        Ninteraction = PuffList_tmp.size();

        for (s = 0; s < this->Ns; s++)
          {
            // Background concentrations.
            if (puff->HasMeteo())
              background_concentration(s) =
                puff->GetBackgroundConcentration(s);
            else
              background_concentration(s) = this->background_concentration_(s);

            // List of concentrations to be modified by chemistry.
            tmp = 0.;
            tmp2 = 0.;
            if (this->option_process["with_puff_interaction"])
              for (i = 0; i < Ninteraction; i++)
                {
                  beta = PuffList_tmp[i];
                  tmp += interaction_coefficient(alpha, beta)
                    * quantity_list(beta, s);
                  tmp2 += interaction_coefficient(alpha, beta)
                    / interaction_coefficient(beta, beta);
                }
            else
              tmp = quantity_list(alpha, s)
                * interaction_coefficient(alpha, alpha);

            concentration_list(s) = tmp;

            if (concentration_list(s) + background_concentration(s) < 0.)
              {
                concentration_list(s) = 0.;
                background_concentration(s) = 0.;
              }

            // Overlap volume.
            if (this->option_process["with_puff_interaction"])
	      {
		if (concentration_list(s) != 0.)
		  overlap_volume(s) = quantity_list(alpha, s)
		    / concentration_list(s);
		else
		  overlap_volume(s) = 1. /
		    (interaction_coefficient(alpha, alpha) * tmp2);
		overlap_volume2(s) = 1. / 
		  (interaction_coefficient(alpha, alpha) * tmp2);
	      }
            else
              overlap_volume(s) = 1. / interaction_coefficient(alpha, alpha);

            // Total concentrations (puff and background).
            concentration_list(s) += background_concentration(s);

	    if (overlap_volume(s) > (1. / interaction_coefficient(alpha, alpha)))
	      overlap_volume(s) = 1./ interaction_coefficient(alpha, alpha);
	    if (overlap_volume2(s) > (1. / interaction_coefficient(alpha, alpha)))
	      overlap_volume2(s) = 1./ interaction_coefficient(alpha, alpha);
          }

        for (r = 0; r < this->Nr_photolysis; r++)
          if (puff->HasMeteo())
            photolysis_rate(r) = puff->GetPhotolysisRate(r);
          else
            photolysis_rate(r) = this->photolysis_rate_(r);


        /*** Aerosol chemistry ***/

        // Relative humidity
        _compute_relative_humidity(&(this->specific_humidity_),
                                   &(this->temperature_),
                                   &(this->pressure_),
                                   &relative_humidity_);
        relative_humidity_ =
          min(max(relative_humidity_, 0.05), 0.95);

        // Wet diameter
        this->InitWetDiameter_aer(relative_humidity_, this->temperature_,
                                  wet_diameter_aer);

        // In-cloud wet deposition
        double rain = 0.0;
        Array<T, 1> CurrentVerticalInterface(2);
        Array<T, 1> InCloudWetDepositionFlux1D(this->Ns);
        Array<T, 1> InCloudWetDepositionFlux_number1D(this->Nbin_aer);
        Array<T, 1> concentration_number_test(this->Nbin_aer);
        Array<T, 2> InCloudWetDepositionFlux_aer2D(this->Ns_aer, this->Nbin_aer);
        double lwc_avg = 0.0;
        double heightfog = 0.0;
        int ifog = 0;
        int b;
    	double tot_in_tmp, tot_0_tmp;
    	double coeff_out;

	Array<T, 1> deltaMass_aer(this->Ns_aer);
	deltaMass_aer = 0.;

        for (s = 0; s < this->Ns; s++)
          InCloudWetDepositionFlux1D(s) = 0.0;
        for (b = 0; b < this->Nbin_aer; b++)
            InCloudWetDepositionFlux_number1D(b) = 0.0;


	T Ctot_tmp;
	T Cback_tmp;
	T Qpuff_tmp;
	deltaMass_aer = 0.;
	aer_compensation = 0.;
        for (s = 0; s < this->Ns_aer ; s++)
          {
	    if (this->option_process["with_puff_interaction"] && Ninteraction > 1)
	      {
		tmp = 0.;
		tmp2 = 0.;
		Ctot_tmp = 0.;
		Cback_tmp = 0.;
		Qpuff_tmp = 0.;
		for (i = 0; i < Ninteraction; i++)
		  {
		    for (b = 0; b < this->Nbin_aer ; b++)
		      {		    		      
			beta = PuffList_tmp[i];
			tmp += interaction_coefficient(alpha, beta)
			  * quantity_list_aer(beta, s, b);
		      }
		    tmp2 += interaction_coefficient(alpha, beta)
		      / interaction_coefficient(beta, beta);
		  }

		for (b = 0; b < this->Nbin_aer ; b++)
		  {
		    if (puff->HasMeteo())
		      background_concentration_aer(s, b) =
			puff->GetBackgroundConcentration(s, b);
		    else
		      background_concentration_aer(s, b) = background_concentration_aer_(s, b);
		    
		    Cback_tmp += background_concentration_aer(s, b);
		    Qpuff_tmp += quantity_list_aer(alpha, s, b);
		  }

		Ctot_tmp = tmp;
		if (Ctot_tmp + Cback_tmp < 0.)
		  Ctot_tmp = 0.;

		
		if (Ctot_tmp > 0.)
		  for (b = 0; b < this->Nbin_aer ; b++)
		    {
		      overlap_volume_aer(s, b) = Qpuff_tmp / Ctot_tmp;
		      overlap_volume_aer2(s, b) = 1. / 
			(interaction_coefficient(alpha, alpha) * tmp2);
		  
		      if (overlap_volume_aer(s, b) > 1. / interaction_coefficient(alpha, alpha))
			overlap_volume_aer(s, b) = 1. / interaction_coefficient(alpha, alpha);
		      if (overlap_volume_aer2(s, b) > 1. / interaction_coefficient(alpha, alpha))
			overlap_volume_aer2(s, b) = 1. / interaction_coefficient(alpha, alpha);
		    }
		else
		  for (b = 0; b < this->Nbin_aer ; b++)
		    {
		      overlap_volume_aer(s, b) = 1. /
			(interaction_coefficient(alpha, alpha) * tmp2);
		      overlap_volume_aer2(s, b) = 1. / 
			(interaction_coefficient(alpha, alpha) * tmp2);
		      
		      if (overlap_volume_aer(s, b) > 1. / interaction_coefficient(alpha, alpha))
			overlap_volume_aer(s, b) = 1. / interaction_coefficient(alpha, alpha);
		      if (overlap_volume_aer2(s, b) > 1. / interaction_coefficient(alpha, alpha))
			overlap_volume_aer2(s, b) = 1. / interaction_coefficient(alpha, alpha);
		    }
	      }
	    else
	      for (b = 0; b < this->Nbin_aer ; b++)
		overlap_volume_aer(s, b) = 1./ interaction_coefficient(alpha, alpha);


	    for (b = 0; b < this->Nbin_aer ; b++)
	      {
		
		InCloudWetDepositionFlux_aer2D(s, b) = 0.0;
		//Background concentrations.
		if (puff->HasMeteo())
		  background_concentration_aer(s, b) =
		    puff->GetBackgroundConcentration(s, b);
		else
		  background_concentration_aer(s, b) = background_concentration_aer_(s, b);


                // List of concentrations to be modified by chemistry.
                tmp = 0.;
                tmp2 = 0.;
                if (this->option_process["with_puff_interaction"])
                  for (i = 0; i < Ninteraction; i++)
                    {
                      beta = PuffList_tmp[i];
                      tmp += interaction_coefficient(alpha, beta)
                        * quantity_list_aer(beta, s, b);
                      tmp2 += interaction_coefficient(alpha, beta)
                        / interaction_coefficient(beta, beta);
                    }
                else
                  tmp = quantity_list_aer(alpha, s, b)
                    * interaction_coefficient(alpha, alpha);
                concentration_list_aer(s, b) = tmp;

                if (concentration_list_aer(s, b) + background_concentration_aer(s, b) < 0.)
                  {
		    T delta_mass = concentration_list_aer(s, b) + background_concentration_aer(s, b);
		    deltaMass_aer(s) += delta_mass;
		    aer_compensation(s) += delta_mass;
		    concentration_list_aer(s, b) = 0.;
                  }
		else
		  concentration_list_aer(s, b) += background_concentration_aer(s, b);		  
	      }
	  } 

	T Ntot_tmp;
	T Cnback_tmp;
	T Cntot_tmp;
	
	if (this->option_process["with_puff_interaction"] && Ninteraction > 1)
	  {
	    tmp = 0.;
	    tmp2 = 0.;
	    Ntot_tmp = 0.;
	    Cnback_tmp = 0.;
	    Cntot_tmp = 0.;
	    for (i = 0; i < Ninteraction; i++)
	      {
		for (b = 0; b < this->Nbin_aer; b++)
		  {
		    beta = PuffList_tmp[i];
		    tmp += interaction_coefficient(alpha, beta)
		      * quantity_list_number(beta, b);
		  }
		tmp2 += interaction_coefficient(alpha, beta)
		  / interaction_coefficient(beta, beta);
	      }

	    for (b = 0; b < this->Nbin_aer; b++)
	      {
		if (puff->HasMeteo())
		  background_concentration_number(b) =
		    puff->GetBackgroundNumberConcentration(b);
		else
		  background_concentration_number(b) = this->background_concentration_number_(b);
		Cnback_tmp += background_concentration_number(b);
		Ntot_tmp += quantity_list_number(alpha, b);
	      }
	    Cntot_tmp = tmp;
	    if (Cntot_tmp + Cnback_tmp < 0.)
	      Cntot_tmp = 0.;
	    
	    if (Cntot_tmp > 0.)
	      for (b = 0; b < this->Nbin_aer; b++)
		{
		  overlap_volume_number(b) = Ntot_tmp / Cntot_tmp;
		  overlap_volume_number2(b) = 1. / 
		    (interaction_coefficient(alpha, alpha) * tmp2);

		  if (overlap_volume_number(b) > 1. / interaction_coefficient(alpha, alpha))
		    overlap_volume_number(b) = 1. / interaction_coefficient(alpha, alpha);
		  if (overlap_volume_number2(b) > 1. / interaction_coefficient(alpha, alpha))
		    overlap_volume_number2(b) = 1. / interaction_coefficient(alpha, alpha);
		}
	    else
	      for (b = 0; b < this->Nbin_aer; b++)
		{
		  overlap_volume_number(b) = 1. /
		    (interaction_coefficient(alpha, alpha) * tmp2);
		  overlap_volume_number2(b) = 1. / 
		    (interaction_coefficient(alpha, alpha) * tmp2);
		  
		  if (overlap_volume_number(b) > 1. / interaction_coefficient(alpha, alpha))
		    overlap_volume_number(b) = 1. / interaction_coefficient(alpha, alpha);
		  if (overlap_volume_number2(b) > 1. / interaction_coefficient(alpha, alpha))
		    overlap_volume_number2(b) = 1. / interaction_coefficient(alpha, alpha);
		}
	  }
	else
	  for (b = 0; b < this->Nbin_aer; b++)
	    overlap_volume_number(b) = 1./ interaction_coefficient(alpha, alpha);

	for (b = 0; b < this->Nbin_aer; b++)
	  {
	    if (puff->HasMeteo())
	      background_concentration_number(b) =
	    	puff->GetBackgroundNumberConcentration(b);
	    else
	      background_concentration_number(b) = this->background_concentration_number_(b);

	    tmp = 0.;
	    tmp2 = 0.;
	    if (this->option_process["with_puff_interaction"])
	      for (i = 0; i < Ninteraction; i++)
		{
		  beta = PuffList_tmp[i];
		  tmp += interaction_coefficient(alpha, beta)
		    * quantity_list_number(beta, b);
		  tmp2 += interaction_coefficient(alpha, beta)
		    / interaction_coefficient(beta, beta);
		}
	    else
	      tmp = quantity_list_number(alpha, b)
		* interaction_coefficient(alpha, alpha);

	    concentration_list_number(b) = tmp;
	    if (concentration_list_number(b) + background_concentration_number(b) > 0.)
	      concentration_list_number(b) += background_concentration_number(b);
	    else
	      concentration_list_number(b) = 0.;
	  }


	//redistribute DM
	Array<T, 2> concentration_list_aer0(this->Ns_aer, this->Nbin_aer);
	concentration_list_aer0 = 0.;
	T conc_tot_aer;
	conc_tot_aer = 0.;

	for (s = 0; s < this->Ns_aer; s++)
	  {
	    conc_tot_aer = 0.;
	    for (b = 0; b < this->Nbin_aer; b++)
	      {
		conc_tot_aer += concentration_list_aer(s, b);
		concentration_list_aer0(s, b) = concentration_list_aer(s, b);
	      }

	    if (conc_tot_aer <= 0.)
	      for (b = 0; b < this->Nbin_aer; b++)
		concentration_list_aer(s, b) = 0.;
	    else	     
	      for (b = 0; b < this->Nbin_aer; b++)
		{
		  T ratio_mass;			    
		  ratio_mass = concentration_list_aer(s, b) / conc_tot_aer;
		  if (deltaMass_aer(s) < 0.)
		    concentration_list_aer(s, b) = max(concentration_list_aer(s, b) 
						       + deltaMass_aer(s, b) * ratio_mass, 0.);

		}
	    //Redistributes aer_compensation on gas phase to solve equilibirum unstability
	    //for organic species only
	    if (is_organic(s) == 1 && aer_compensation(s) < 0.)
	      {
		int gas_index;
		gas_index = species_index_aerosol_interact_tmp(s);
		if (concentration_list(gas_index) + aer_compensation(s) > 0.)
		  concentration_list(gas_index) += aer_compensation(s);
		else
		  concentration_list(gas_index) = 0.;
	      }
	  }

	//Recompute number to keep diameter after DM redist
	for (b = 0; b < this->Nbin_aer; b++)
	  {
	    T number_new;
	    T rho_p1, rho_p0;
	    Data<T, 1> conc0_tmp(this->Ns_aer);
	    Data<T, 1> conc1_tmp(this->Ns_aer);
	    T conc0_tot, conc1_tot;
	    conc0_tmp.SetZero();
	    conc1_tmp.SetZero();
	    conc0_tot = 0.;
	    conc1_tot = 0.;
	    for (s = 0; s < this->Ns_aer-1; s++)
	      {
		conc0_tmp(s) = concentration_list_aer0(s, b);
		conc1_tmp(s) = concentration_list_aer(s, b);
		conc0_tot += concentration_list_aer0(s, b);
		conc1_tot += concentration_list_aer(s, b);
	      }
	    rho_p0 = ComputeDensity(conc0_tmp, Mass_Density_aer, conc0_tot, this->Ns_aer-1);
	    rho_p1 = ComputeDensity(conc1_tmp, Mass_Density_aer, conc1_tot, this->Ns_aer-1);
	    
	    if (conc1_tot == 0.)
	      number_new = 0.;
	    else
	      {
		if (conc0_tot == 0. || concentration_list_number(b) == 0.)
		  number_new = ComputeNumberPuff(b, concentration_list_aer);
		else
		  {
		    T ratio_number;
		    ratio_number = (conc1_tot * rho_p0) / (conc0_tot * rho_p1);
		    number_new = concentration_list_number(b) * ratio_number;
		  }
	      }
	    concentration_list_number(b) = number_new;
	    if (concentration_list_number(b) < 0.)
	      concentration_list_number(b) = 0.;	    
	  }


	for (b = 0; b < this->Nbin_aer; b++)
	  {
	    T tot_puff, tot_back;
	    tot_puff = 0.;
	    tot_back = 0.;
	    for (s = 0; s < this->Ns_aer-1; s++)
	      {
		tot_puff += concentration_list_aer(s, b);
		tot_back += background_concentration_aer(s, b);
	      }

	    if (concentration_list_number(b) == 0. && tot_puff > 0.)
	      concentration_list_number(b) = ComputeNumberPuff(b, concentration_list_aer);
	      
	    if (background_concentration_number(b) == 0. && tot_back > 0.)
	      background_concentration_number(b) = ComputeNumberPuff(b, background_concentration_aer);
	  }

        CurrentVerticalInterface(0) = 0.0;
        CurrentVerticalInterface(1) = 0.0;

        // Chemistry for background concentrations only (gaseous species).
        this->Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                                 this->attenuation_, this->specific_humidity_,
                                 this->temperature_, this->pressure_, source,
                                 photolysis_rate,
                                 T(this->next_date.
                                   GetSecondsFrom(this->current_date)),
                                 this->attenuation_, this->specific_humidity_,
                                 this->temperature_, this->pressure_, source,
                                 photolysis_rate, this->longitude_,
                                 this->latitude_, background_concentration,
                                 liquid_water_content_, wet_diameter_aer,
                                 background_concentration_aer, ph_,
				 background_concentration_number);

        // Chemistry for background concentrations only (particulate matters).
        this->Chemistry_.Forward_aer(T(this->current_date.GetNumberOfSeconds()),
                                     this->specific_humidity_,
                                     this->temperature_, this->pressure_,
                                     T(this->next_date.
                                       GetSecondsFrom(this->current_date)),
                                     background_concentration,
                                     liquid_water_content_, rain,
                                     CurrentVerticalInterface,
                                     background_concentration_aer,
                                     InCloudWetDepositionFlux1D,
                                     InCloudWetDepositionFlux_aer2D,
                                     ph_, lwc_avg, heightfog, ifog,
				     background_concentration_number,
    				     InCloudWetDepositionFlux_number1D);

        // Chemistry for total concentrations (gas species).
        this->Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                                 this->attenuation_, this->specific_humidity_,
                                 this->temperature_, this->pressure_, source,
                                 photolysis_rate,
                                 T(this->next_date.
                                   GetSecondsFrom(this->current_date)),
                                 this->attenuation_, this->specific_humidity_,
                                 this->temperature_, this->pressure_, source,
                                 photolysis_rate, this->longitude_,
                                 this->latitude_, concentration_list,
                                 liquid_water_content_, wet_diameter_aer,
                                 concentration_list_aer, ph_,
				 concentration_list_number);

        // Chemistry for total concentrations (particulate matters).
        this->Chemistry_.Forward_aer(T(this->current_date.GetNumberOfSeconds()),
                                     this->specific_humidity_,
                                     this->temperature_, this->pressure_,
                                     T(this->next_date.
                                       GetSecondsFrom(this->current_date)),
                                     concentration_list,
                                     liquid_water_content_, rain,
                                     CurrentVerticalInterface,
                                     concentration_list_aer,
                                     InCloudWetDepositionFlux1D,
                                     InCloudWetDepositionFlux_aer2D,
                                     ph_, lwc_avg, heightfog, ifog,
				     concentration_list_number,
    				     InCloudWetDepositionFlux_number1D);

	puff->SetPuffVolume(1. / interaction_coefficient(alpha, alpha));

        // New puff quantities.
        for (s = 0; s < this->Ns; s++)
          {
            T quantity = 0.;
            // Substracting background concentrations.
            concentration_list(s) -= background_concentration(s);
	    if (concentration_list(s) < 0.)
	      quantity = concentration_list(s) * overlap_volume2(s);
	    else
	      quantity = concentration_list(s) * overlap_volume(s);

            // Setting puff quantity.
            puff->SetQuantity(quantity, s);
            species_quantity(s) += quantity;

            // Setting puff background concentration.
            puff->SetBackgroundConcentration(background_concentration(s), s);
            this->background_concentration_(s) = background_concentration(s);
          }

	for (b = 0; b < this->Nbin_aer; b++)
	  {
	    T number_conc_delta;
	    number_conc_delta = concentration_list_number(b) - background_concentration_number(b);
	    if (number_conc_delta < 0.)
	      number_conc_delta *= overlap_volume_number2(b);
	    else
	      number_conc_delta *= overlap_volume_number(b);

	    puff->SetBackgroundNumberConcentration(background_concentration_number(b), b);
	    puff->SetNumberQuantity(number_conc_delta, b);
	  }


        for (s = 0; s < this->Ns_aer ; s++)
          for (b = 0; b < this->Nbin_aer ; b++)
            {
              T quantity = 0.;
	      quantity = concentration_list_aer(s, b) - background_concentration_aer(s, b);
	      if (quantity < 0.)
		quantity *= overlap_volume_aer2(s, b);
	      else
		quantity *= overlap_volume_aer(s, b);

              // Setting puff quantity.
              puff->SetQuantity(quantity, s, b);
              species_quantity_aer(s, b) += quantity;

              // Setting puff background concentration.
              puff->SetBackgroundConcentration(background_concentration_aer(s, b), s, b);
              background_concentration_aer_(s, b) = background_concentration_aer(s, b);

            }

        alpha++;
      } // iteration for puff

    if (this->option_process["with_output_plume_mass"])
      FormatBinary<float>().Append(species_quantity, this->file_mass);
  }


  // //!  Performs the chemistry when the feedback of the species concentration is taken into account from the Gaussian puff to the background.
  // template<class T, class ClassChemistry>
  // void GaussianPuffAerosol<T, ClassChemistry>
  // ::Chemistry(vector<list<int> > PuffCellList,
  //             vector<T> PuffCellVolume,
  //             Array<vector<T>, 1 >& PuffCellConcentration,
  //             Array<vector<T>, 2 >& PuffCellConcentration_aer)
  // {
  //   int Ncell = PuffCellVolume.size();
  //   int Ntot = this->Npuff + Ncell;

  //   // Definition of arrays.
  //   Array<T, 2> quantity_list(Ntot, this->Ns);
  //   Array<T, 3> quantity_list_aer(Ntot, this->Ns_aer, this->Nbin_aer);
  //   Array<T, 2> interaction_coefficient(Ntot, Ntot);
  //   Array<T, 1> concentration_list(this->Ns);
  //   Array<T, 2> concentration_list_aer(this->Ns_aer, this->Nbin_aer);
  //   Array<T, 1> overlap_volume(this->Ns);
  //   Array<T, 2> overlap_volume_aer(this->Ns_aer, this->Nbin_aer);
  //   Array<T, 1> background_concentration(this->Ns);
  //   Array<T, 2> background_concentration_aer(this->Ns_aer, this->Nbin_aer);
  //   Array<T, 1> photolysis_rate(this->Nr_photolysis);
  //   Array<T, 1> source(this->Ns);
  //   Data<T, 1> species_quantity(this->Ns);
  //   Data<T, 2> species_quantity_aer(this->Ns_aer, this->Nbin_aer);

  //   // Initialization of arrays.
  //   source = 0.;
  //   quantity_list = 0.;
  //   quantity_list_aer = 0.;
  //   interaction_coefficient = 0.;
  //   species_quantity.SetZero();
  //   species_quantity_aer.SetZero();

  //   // Other variables.
  //   int s, r, alpha, beta, i, b;
  //   typename list<Puff<T>* >::iterator iter, iter2;
  //   Puff<T>* puff;

  //   // For interactions.
  //   T tmp, tmp2;
  //   Array<vector<int>, 1 > PuffInteractionList(Ntot);
  //   vector<int> PuffList_tmp;
  //   int Ninteraction;

  //   // Getting quantities, computing interaction coefficients for all puffs.
  //   alpha = 0;
  //   for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
  //     {
  //       // List of quantities for all puffs.
  //       for (s = 0; s < this->Ns; s++)
  //         quantity_list(alpha, s) = (*iter)->GetQuantity(s);

  //       for (s = 0; s < this->Ns_aer; s++)
  //         for (b = 0; b < this->Nbin_aer; b++)
  //           quantity_list_aer(alpha, s, b) = (*iter)->GetQuantity(s, b);

  //       // Interaction coefficient for all puff pairs.
  //       if (this->option_process["with_puff_interaction"])
  //         {
  //           beta = 0;
  //           PuffList_tmp.clear();
  //           for (iter2 = this->PuffList.begin();
  //                iter2 != this->PuffList.end(); iter2++)
  //             {
  //               tmp = this->ComputeInteractionCoefficient(*iter, *iter2);
  //               if (tmp != 0.)
  //                 {
  //                   interaction_coefficient(alpha, beta) = tmp;
  //                   PuffList_tmp.push_back(beta);
  //                 }
  //               beta++;
  //             }
  //           PuffInteractionList(alpha) = PuffList_tmp;
  //         }
  //       else
  //         interaction_coefficient(alpha, alpha) =
  //           this->ComputeInteractionCoefficient(*iter, *iter);
  //       alpha++;
  //     }

  //   // Getting quantities, computing interaction coefficients for cells.
  //   typename list<int>::iterator iter3;
  //   list<int> pufflist_tmp;
  //   int index;
  //   for (i = 0; i < Ncell; i++)
  //     {
  //       index = this->Npuff + i;
  //       pufflist_tmp.clear();
  //       pufflist_tmp = PuffCellList[i];

  //       // Background quantities.
  //       for (s = 0; s < this->Ns; s++)
  //         quantity_list(index, s) = PuffCellVolume[i]
  //           * PuffCellConcentration(s)[i];

  //       for (s = 0; s < this->Ns_aer; s++)
  //         for (b = 0; b < this->Nbin_aer; b++)
  //           quantity_list_aer(index, s, b) = PuffCellVolume[i]
  //             * PuffCellConcentration_aer(s, b)[i];

  //       // Interaction coefficient between puff and cell.
  //       for (alpha = 0; alpha < this->Npuff; alpha++)
  //         {
  //           iter3 = find(pufflist_tmp.begin(), pufflist_tmp.end(), alpha);
  //           if (iter3 != pufflist_tmp.end())
  //             {
  //               interaction_coefficient(alpha, index) = 1. / PuffCellVolume[i];
  //               interaction_coefficient(index, alpha) = 1. / PuffCellVolume[i];
  //               PuffInteractionList(alpha).push_back(index);
  //               PuffInteractionList(index).push_back(alpha);
  //             }
  //         }
  //       // Interaction coefficient between two indentical cells.
  //       interaction_coefficient(index, index) = 1. / PuffCellVolume[i];
  //       PuffInteractionList(index).push_back(index);

  //     }

  //   // Chemistry for all puffs and cells.
  //   int puff_index;
  //   T quantity;
  //   for (alpha = 0; alpha < Ntot; alpha++)
  //     {
  //       // Looking for the appropriate puff to take the meteorological data.
  //       puff_index = 0;
  //       if (alpha < this->Npuff)
  //         puff_index = alpha;
  //       else
  //         puff_index = PuffInteractionList(alpha)[0];

  //       this->SetCurrentPuff(puff_index);
  //       puff = *this->current_puff;
  //       SetCurrentMeteo(puff);
  //       overlap_volume = 0.;
  //       PuffList_tmp = PuffInteractionList(alpha);
  //       Ninteraction = PuffList_tmp.size();

  //       // Computing overlap concentrations and overlap volumes.
  //       for (s = 0; s < this->Ns; s++)
  //         {
  //           // List of concentrations to be modified by chemistry.
  //           tmp = 0.;
  //           tmp2 = 0.;
  //           if (this->option_process["with_puff_interaction"])
  //             for (i = 0; i < Ninteraction; i++)
  //               {
  //                 beta = PuffList_tmp[i];
  //                 tmp += interaction_coefficient(alpha, beta)
  //                   * quantity_list(beta, s);
  //                 tmp2 += interaction_coefficient(alpha, beta)
  //                   / interaction_coefficient(beta, beta);
  //               }
  //           else
  //             tmp = quantity_list(alpha, s)
  //               * interaction_coefficient(alpha, alpha);
  //           concentration_list(s) = tmp;

  //           if (concentration_list(s) < 0.)
  //             concentration_list(s) = 0.;

  //           // Overlap volume.
  //           if (this->option_process["with_puff_interaction"])
  //             if (concentration_list(s) != 0.)
  //               overlap_volume(s) = quantity_list(alpha, s)
  //                 / concentration_list(s);
  //             else
  //               overlap_volume(s) = 1. /
  //                 (interaction_coefficient(alpha, alpha) * tmp2);
  //           else
  //             overlap_volume(s) = 1. / interaction_coefficient(alpha, alpha);
  //         }
  //       for (r = 0; r < this->Nr_photolysis; r++)
  //         if (puff->HasMeteo())
  //           photolysis_rate(r) = puff->GetPhotolysisRate(r);
  //         else
  //           photolysis_rate(r) = this->photolysis_rate_(r);

  //       /*** Relative humidity ***/

  //       _compute_relative_humidity(&(this->specific_humidity_),
  //                                  &(this->temperature_),
  //                                  &(this->pressure_),
  //                                  &relative_humidity_);
  //       relative_humidity_ =
  //         min(max(relative_humidity_, 0.), 1.);

  //       // Wet diameter
  //       this->InitWetDiameter_aer(relative_humidity_, this->temperature_,
  //                                 wet_diameter_aer);

  //       // Chemistry for total concentrations.
  //       this->Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
  //                                this->attenuation_, this->specific_humidity_,
  //                                this->temperature_, this->pressure_, source,
  //                                photolysis_rate,
  //                                T(this->next_date.GetNumberOfSeconds()),
  //                                this->attenuation_, this->specific_humidity_,
  //                                this->temperature_, this->pressure_, source,
  //                                photolysis_rate, this->longitude_,
  //                                this->latitude_, concentration_list,
  //                                liquid_water_content_, wet_diameter_aer,
  //                                concentration_list_aer, ph_);


  //       // Aerosol chemistry for background concentrations only.
  //       double rain = 0.0;
  //       Array<T, 1> CurrentVerticalInterface(2);
  //       Array<T, 1> InCloudWetDepositionFlux1D(this->Ns);
  //       Array<T, 2> InCloudWetDepositionFlux_aer2D(this->Ns_aer, this->Nbin_aer);
  //       double lwc_avg = 0.0;
  //       double heightfog = 0.0;
  //       int ifog = 0;
  //       int b;
  //       for (s = 0; s < this->Ns; s++)
  //         {
  //           InCloudWetDepositionFlux1D(s) = 0.0;
  //         }


  //       for (s = 0; s < this->Ns_aer ; s++)
  //         for (b = 0; b < this->Nbin_aer ; b++)
  //           {
  //             InCloudWetDepositionFlux_aer2D(s, b) = 0.0;

  //             // Background concentrations.
  //             if (puff->HasMeteo())
  //               background_concentration_aer(s, b) =
  //                 puff->GetBackgroundConcentration(s, b);
  //             else
  //               background_concentration_aer(s, b) = this->background_concentration_aer_(s, b);

  //             // List of concentrations to be modified by chemistry.
  //             tmp = 0.;
  //             tmp2 = 0.;
  //             if (this->option_process["with_puff_interaction"])
  //               for (i = 0; i < Ninteraction; i++)
  //                 {
  //                   beta = PuffList_tmp[i];
  //                   tmp += interaction_coefficient(alpha, beta)
  //                     * quantity_list_aer(beta, s, b);
  //                   tmp2 += interaction_coefficient(alpha, beta)
  //                     / interaction_coefficient(beta, beta);
  //                 }
  //             else
  //               tmp = quantity_list_aer(alpha, s, b)
  //                 * interaction_coefficient(alpha, alpha);
  //             concentration_list_aer(s, b) = tmp;

  //             if (concentration_list_aer(s, b) + background_concentration_aer(s, b) < 0.)
  //               {
  //                 concentration_list_aer(s, b) = 0.;
  //                 background_concentration_aer(s, b) = 0.;
  //               }

  //             // Overlap volume.
  //             if (this->option_process["with_puff_interaction"])
  //               if (concentration_list_aer(s, b) != 0.)
  //                 overlap_volume_aer(s, b) = quantity_list_aer(alpha, s, b)
  //                   / concentration_list_aer(s, b);
  //               else
  //                 overlap_volume_aer(s, b) = 1. /
  //                   (interaction_coefficient(alpha, alpha) * tmp2);
  //             else
  //               overlap_volume_aer(s, b) = 1. / interaction_coefficient(alpha, alpha);

  //             // Total concentrations (puff and background).
  //             concentration_list_aer(s, b) += background_concentration_aer(s, b);
  //           }

  //       CurrentVerticalInterface(0) = 0.0;
  //       CurrentVerticalInterface(1) = 0.0;

  //       this->Chemistry_.Forward_aer(T(this->current_date.GetNumberOfSeconds()),
  //                                    this->specific_humidity_,
  //                                    this->temperature_, this->pressure_,
  //                                    T(this->next_date.
  //                                      GetSecondsFrom(this->current_date)),
  //                                    concentration_list,
  //                                    liquid_water_content_, rain,
  //                                    CurrentVerticalInterface,
  //                                    concentration_list_aer,
  //                                    InCloudWetDepositionFlux1D,
  //                                    InCloudWetDepositionFlux_aer2D,
  //                                    ph_, lwc_avg, heightfog, ifog);


  //       // New puff quantities.
  //       if (alpha < this->Npuff)
  //         {
  //           for (s = 0; s < this->Ns; s++)
  //             {
  //               quantity = concentration_list(s) * overlap_volume(s);
  //               species_quantity(s) += quantity;
  //               puff->SetQuantity(quantity, s);
  //             }

  //           for (s = 0; s < this->Ns_aer ; s++)
  //             for (b = 0; b < this->Nbin_aer ; b++)
  //               {
  //                 quantity = concentration_list_aer(s, b) * overlap_volume_aer(s, b);
  //                 species_quantity_aer(s, b) += quantity;
  //                 // Setting puff quantity.
  //                 puff->SetQuantity(quantity, s, b);
  //               }
  //         }
  //       // New cell quantities.
  //       else
  //         {
  //           i = alpha - this->Npuff;

  //           // Getting background concentrations (without puffs).
  //           for (s = 0; s < this->Ns; s++)
  //             background_concentration(s) = PuffCellConcentration(s)[i];

  //           // Chemistry for background concentrations.
  //           this->Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
  //                                    this->attenuation_, this->specific_humidity_,
  //                                    this->temperature_, this->pressure_, source,
  //                                    photolysis_rate,
  //                                    T(this->next_date.GetNumberOfSeconds()),
  //                                    this->attenuation_, this->specific_humidity_,
  //                                    this->temperature_, this->pressure_, source,
  //                                    photolysis_rate, this->longitude_,
  //                                    this->latitude_, background_concentration,
  //                                    liquid_water_content_, wet_diameter_aer,
  //                                    background_concentration_aer, ph_);

  //           this->Chemistry_.Forward_aer(T(this->current_date.GetNumberOfSeconds()),
  //                                        this->specific_humidity_,
  //                                        this->temperature_, this->pressure_,
  //                                        T(this->next_date.
  //                                          GetSecondsFrom(this->current_date)),
  //                                        background_concentration,
  //                                        liquid_water_content_, rain,
  //                                        CurrentVerticalInterface,
  //                                        background_concentration_aer,
  //                                        InCloudWetDepositionFlux1D,
  //                                        InCloudWetDepositionFlux_aer2D,
  //                                        ph_, lwc_avg, heightfog, ifog);

  //           // Getting the concentration perturbation due to puffs.
  //           for (s = 0; s < this->Ns; s++)
  //             {
  //               quantity = concentration_list(s) * overlap_volume(s);
  //               PuffCellConcentration(s)[i] =
  //                 quantity * interaction_coefficient(alpha, alpha)
  //                 - background_concentration(s);
  //             }
  //           for (s = 0; s < this->Ns_aer ; s++)
  //             for (b = 0; b < this->Nbin_aer ; b++)
  //               {
  //                 quantity = concentration_list_aer(s, b) * overlap_volume_aer(s, b);
  //                 PuffCellConcentration(s, b)[i] =
  //                   quantity * interaction_coefficient(alpha, alpha)
  //                   - background_concentration_aer(s, b);
  //               }
  //         }
  //     }
  //   if (this->option_process["with_output_plume_mass"])
  //     FormatBinary<float>().Append(species_quantity, this->file_mass);
  // }


  //! Performs one step forward.
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::Forward()
  {
    this->Advection();
    this->Diffusion();

    int s, i, j, k;
    int b;
    this->Concentration.SetZero();
    this->Concentration_aer.SetZero();

    if (this->option_process["with_chemistry"])
      {
        Chemistry();
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              for (i = 0; i < this->Nx; i++)
                this->Concentration(s, k, j, i)
                  = this->background_concentration_(s);

        for (s = 0; s < this->Ns_aer; s++)
          for (b = 0; b < this->Nbin_aer ; b++)
            for (k = 0; k < this->Nz; k++)
              for (j = 0; j < this->Ny; j++)
                for (i = 0; i < this->Nx; i++)
                  this->Concentration_aer(s, b, k, j, i)
                    = this->background_concentration_aer_(s, b);
      }

    for (typename list<Puff<T>* >::iterator iter = this->PuffList.begin();
         iter != this->PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        SetCurrentMeteo(puff);
        for (s = 0; s < this->Ns; s++)
          if (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == s)
            {
              // Computing loss factors.
              this->ComputeLossFactor(puff, s);
              // Loop on all points to compute concentration.
              for (k = 0; k < this->Nz; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration(s, k, j, i) +=
                      this->ComputePuffConcentration(puff, s,
                                                     this->GridX4D(i),
                                                     this->GridY4D(j),
                                                     this->GridZ4D(k));
            }

        for (s = 0; s < this->Ns_aer; s++)
          for (b = 0; b < this->Nbin_aer; b++)
            {
              // Computing loss factors.
              ComputeLossFactorAerosol(puff, s, b);
              // Loop on all points to compute concentration.
              for (k = 0; k < this->Nz; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration_aer(s, b, k, j, i) +=
                      this->ComputePuffConcentrationAerosol(puff, s, b,
                                                            this->GridX5D(i),
                                                            this->GridY5D(j),
                                                            this->GridZ5D(k));
            }
      }
    this->AddTime(this->Delta_t);
    this->step++;
  }

  //! Computes the loss factor for a given puff, and all species.
  /*! Computes the loss factor for a given puff, and all species.
    \param index puff index.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::InitWetDiameter_aer(T& RelativeHumidity_,
                        T& Temperature_,
                        Array<T, 1>& WetDiameter_aer_)
  {
    int b;

    Array<T, 1> MeanDiameter(this->Nbin_aer);
    for (b = 0; b < this->Nbin_aer; b++)
      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    for (b = 0; b < this->Nbin_aer; b++)
      {
        _gerber_wet_diameter(&RelativeHumidity_,
                             &Temperature_, &MeanDiameter(b),
                             &WetDiameter_aer_(b));
        // Back to meters.
        WetDiameter_aer_(b) = WetDiameter_aer_(b) * 1.e-6;
      }
  }

  //! Computes the loss factor.
  /*! Computes the loss factor.
    \param puff the puff.
    \param s species index.
  */  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::ComputeLossFactorAerosol(Puff<T>* puff, int s, int b)
  {
    T loss_factor = 1.;
    T overcamp_factor = 1.;
    T distance = puff->GetDistance();
    T transfer_time = puff->GetPuffTime();
    T quantity = puff->GetQuantity(s, b);
    T z = puff->GetHeight();

    T rad, bio, scav, dep;
    rad = 0.;
    bio = 0.;
    if (puff->HasMeteo())
      {
        if (this->option_process["with_scavenging"])
          scav = puff->GetScavengingCoefficient(s, b);
        else
          scav = 0.;
        if (this->option_process["with_dry_deposition"])
          dep = puff->GetDepositionVelocity(s, b);
        else
          dep = 0.;
      }
    else
      {
        if (this->option_process["with_scavenging"])
          scav = scavenging_coefficient_aer(s, b);
        else
          scav = 0.;
        if (this->option_process["with_dry_deposition"])
          dep = deposition_velocity_aer(s, b);
        else
          dep = 0.;
      }
    this->ComputeLossFactor(distance, transfer_time, z, rad, bio, scav, dep,
                            loss_factor, overcamp_factor);
    puff->SetQuantity(loss_factor * quantity, s, b);
    puff->SetReflectionFactor(overcamp_factor, s, b);
  }

  //! Computes the species concentration due to the puff at a given point.
  /*!
    \param puff the puff.
    \param x abscissa of the point where concentration is computed (m).
    \param y ordinate of the point where concentration is computed (m).
    \param z height of the point where concentration is computed (m).
    \return The concentration at (\a x, \a y, \a z).
  */
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffConcentrationAerosol(Puff<T>* puff, int s, int b,
                                    T x, T y, T z)
  {

    // Puff quantity.
    T quantity = puff->GetQuantity(s, b);

    // Puff position.
    T z_c = puff->GetZ();
    T z_above = puff->GetHeightAboveBL();
    T penetration_factor = puff->GetPenetrationFactor();
    T overcamp_factor = puff->GetReflectionFactor(s, b);
    T y_c = puff->GetY();
    T x_c = puff->GetX();
    T time = puff->GetPuffTime();

    // Standard deviations.
    T sigma_x = puff->GetSigma_x();
    T sigma_y = puff->GetSigma_y();
    T sigma_z = puff->GetSigma_z();
    T initial_sigma_y_2 = puff->GetInitialSigma_y();
    T initial_sigma_z_2 = puff->GetInitialSigma_z();
    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
    sigma_z = sqrt(sigma_z * sigma_z + initial_sigma_z_2);

    const T pi = 3.14159265358979323846264;
    // (2. * pi) ^ (3/2)
    const T pi__1_5 = 15.749609945722419;

    // Downwind distance from puff center.
    T dx = (x - x_c) * this->cos_angle_ + (y - y_c) * this->sin_angle_;
    // Crosswind distance from puff center.
    T dy = (x_c - x) * this->sin_angle_ + (y - y_c) * this->cos_angle_;

    // Downwind contribution.
    T Fx = exp(-dx * dx  / (2. * sigma_x * sigma_x));

    // Crosswind contribution.
    T Fy = exp(-dy * dy  / (2. * sigma_y * sigma_y));

    // Ground reflection: plume is below the inversion layer.
    T FzG = 0.;
    T sigma_z_2 = sigma_z * sigma_z;
    if (this->inversion_height_ == 0. || z_c != 0.)
      if (z_c - sigma_z <= 0.)
        {
          T dzg = z + z_c;
          FzG = overcamp_factor * exp(-dzg * dzg / (2. * sigma_z_2));
        }
    T Fz = 0.;
    T concentration = 0.;

    // If there is an inversion layer.
    if (this->inversion_height_ > 0.)
      {
        T FzI = 0.;
        // Part of the puff that is below inversion layer.
        if (z_c != 0.)
          {
            // Vertical contribution.
            T dzr = z - z_c;
            Fz = exp(-dzr * dzr / (2. * sigma_z_2));

            // Reflection occurs only when puff touches the inversion layer.
            if (z_c + sigma_z >= this->inversion_height_)
              {
                T dzi = z + z_c - 2. * this->inversion_height_;
                T dzig = z - z_c + 2. * this->inversion_height_;
                T dzgi = z - z_c - 2. * this->inversion_height_;

                FzI += exp(-dzi * dzi / (2. * sigma_z_2))
                  + exp(-dzgi * dzgi / (2. * sigma_z_2))
                  + exp(-dzig * dzig / (2. * sigma_z_2));
              }

            // Plume only impacts below the inversion layer.
            if (z <= this->inversion_height_)
              {
                // Far field.
                if (sigma_z > 1.5 * this->inversion_height_)
                  concentration +=  quantity * Fx * Fy
                    / (2 * pi * sigma_x * sigma_y * this->inversion_height_);
                else
                  concentration += quantity
                    * (1. - penetration_factor)
                    * Fx * Fy * (Fz + FzG + FzI)
                    / (pi__1_5 * sigma_x * sigma_y * sigma_z);
              }
          }

        // Part of the puff that is above inversion layer.
        if (z_above != 0.)
          {
            // Vertical contribution.
            T dzr = z - z_above;
            Fz = exp(-dzr * dzr / (2. * sigma_z_2));

            if (this->option_gillani)
              {
                initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
                sigma_z = sqrt(initial_sigma_z_2 * (1. + 2.3 * sqrt(time)));
                sigma_z_2 = sigma_z * sigma_z;
              }
            // Reflection on the inversion layer.
            T dzi = z + z_above - 2. * this->inversion_height_;
            FzI += exp(-dzi * dzi / (2. * sigma_z_2));
            // Plume only impacts above the inversion layer.
            if (z >= this->inversion_height_)
              concentration +=  quantity
                * penetration_factor * Fx * Fy * (Fz + FzI)
                / (pi__1_5 * sigma_x * sigma_y * sigma_z);
          }
      }
    else
      {
        // Vertical contribution.
        T dzr = z - z_c;
        Fz = exp(-dzr * dzr / (2. * sigma_z_2));

        if (this->option_gillani && z_c >= this->boundary_height_)
          {
            initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
            sigma_z = sqrt(initial_sigma_z_2 * (1. + 2.3 * sqrt(time)));
            sigma_z_2 = sigma_z * sigma_z;
          }
        concentration =  quantity
          * Fx * Fy * (Fz + FzG)
          / (pi__1_5 * sigma_x * sigma_y * sigma_z);
      }
    return concentration;
  }

  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>::GetPuffQuantity_aer(int index, int s, int b)
  {
    this->SetCurrentPuff(index);
    if ((*this->current_puff)->GetNs_aer() == this->Ns_aer)
      return (*this->current_puff)->GetQuantity(s, b);
    else
      return 0.;
  }

  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>::GetPuffQuantity_number(int index, int b)
  {
    this->SetCurrentPuff(index);
    return (*this->current_puff)->GetNumberQuantity(b);
  }


  //! Computes puff additional liquid water content.
  // Compute additional liquid water content due to water emission
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffAdditionalLWC(T temperature,
                             T pressure,
                             T specific_humidity)
  {
    T Pvsat;
    T Qsat;
    T lwc_puff;
    
    Pvsat = 610.78 * exp(17.2694 * (temperature - 273.15) / (temperature - 35.86));
    // cout << "Pvsat : " << Pvsat << endl;

    Qsat = (0.62197 * Pvsat) / (pressure - (1 - 0.62197) * Pvsat);
    // cout << "Qsat : " << Qsat << endl;

    if (specific_humidity - Qsat > 0)
      lwc_puff = specific_humidity - Qsat;
    else
      lwc_puff = 0;
    // cout << "lwc1 : " << lwc_puff << endl;
    return lwc_puff;
  }

  //! Computes the integral of a puff over a given volume.
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffIntegral_aer(int index, int s, int b,
                            T x, T y, T z,
                            T lx, T ly, T lz)
  {
    this->SetCurrentPuff(index);
    this->SetCurrentMeteo(*this->current_puff);
    if ((*this->current_puff)->GetNs_aer() == this->Ns_aer)
      return ComputePuffIntegral_aer(*this->current_puff, s, b, x, y, z, lx, ly, lz);
    else
      return 0.;
  }

  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffIntegral_number(int index, int b,
                               T x, T y, T z,
                               T lx, T ly, T lz)
  {
    this->SetCurrentPuff(index);
    this->SetCurrentMeteo(*this->current_puff);
    return ComputePuffIntegral_number(*this->current_puff, b, 
                                      x, y, z, lx, ly, lz);
  }


  //!  Combine two overlapping puff.
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>
  ::CombineOverlappingPuff(int alpha, int beta, 
                           T volume_alpha, T volume_beta)
  {
    Puff<T>* puff_alpha;
    this->SetCurrentPuff(alpha);
    puff_alpha = (*this->current_puff);
    Puff<T>* puff_beta;
    this->SetCurrentPuff(beta);
    puff_beta = (*this->current_puff);
    T interact_alpha_beta;
    T overlap_tmp;

    interact_alpha_beta = this->ComputeInteractionCoefficient(puff_alpha, puff_beta);
    overlap_tmp = interact_alpha_beta * (volume_alpha * volume_beta);

    T z_tmp, z_above;
    T sigmax_alpha, sigmay_alpha, sigmaz_alpha, x_alpha, y_alpha, z_alpha;
    x_alpha = puff_alpha->GetX();
    y_alpha = puff_alpha->GetY();
    z_tmp = puff_alpha->GetZ();
    z_above = puff_alpha->GetHeightAboveBL();
    if (z_tmp != 0.)
      z_alpha = z_tmp;
    else
      z_alpha = z_above;

    sigmax_alpha = puff_alpha->GetSigma_x();
    sigmay_alpha = puff_alpha->GetSigma_y();
    sigmaz_alpha = puff_alpha->GetSigma_z();

    T sigmax_beta, sigmay_beta, sigmaz_beta, x_beta, y_beta, z_beta;
    x_beta = puff_beta->GetX();
    y_beta = puff_beta->GetY();
    z_tmp = puff_beta->GetZ();
    z_above = puff_beta->GetHeightAboveBL();
    if (z_tmp != 0.)
      z_beta = z_tmp;
    else
      z_beta = z_above;

    sigmax_beta = puff_beta->GetSigma_x();
    sigmay_beta = puff_beta->GetSigma_y();
    sigmaz_beta = puff_beta->GetSigma_z();

    T x_max, y_max, z_max, x_min, y_min, z_min;
    if (x_alpha > x_beta)
      {
	x_max = x_alpha;
	x_min = x_beta;
      }
    else
      {
	x_max = x_beta;
	x_min = x_alpha;
      }
    if (y_alpha > y_beta)
      {
	y_max = y_alpha;
	y_min = y_beta;
      }
    else
      {
	y_max = y_beta;
	y_min = y_alpha;
      }
    if (z_alpha > z_beta)
      {
	z_max = z_alpha;
	z_min = z_beta;
      }
    else
      {
	z_max = z_beta;
	z_min = z_alpha;
      }
    T delta_x, delta_y, delta_z;
    delta_x = x_max - x_min;
    delta_y = y_max - y_min;
    delta_z = z_max - z_min;

    T dx_new, dy_new, dz_new, x_out, y_out, z_out;
    T distance_alpha, distance_beta, distance_out, time_tmp;
    if (x_min == x_alpha)
      dx_new = volume_beta * delta_x / (volume_alpha + volume_beta);
    else
      dx_new = volume_alpha * delta_x / (volume_alpha + volume_beta);
    x_out = x_min + dx_new;
 
   if (y_min == y_alpha)
      dy_new = volume_beta * delta_y / (volume_alpha + volume_beta);
    else
      dy_new = volume_alpha * delta_y / (volume_alpha + volume_beta);
    y_out = y_min + dy_new;

    if (z_min == z_alpha)
      dz_new = volume_beta * delta_z / (volume_alpha + volume_beta);
    else
      dz_new = volume_alpha * delta_z / (volume_alpha + volume_beta);
    z_out = z_min + dz_new;

    distance_alpha = puff_alpha->GetDistance();
    distance_beta = puff_beta->GetDistance();
    distance_out = (distance_alpha + distance_beta) / 2.;
    time_tmp = puff_alpha->GetPuffTime();

    puff_alpha->SetPuffPosition(x_out, y_out, z_out, distance_out, time_tmp);

    // T coeff_sigma3, coeff_sigma;
    T sigmax_out, sigmay_out, sigmaz_out;

    T coeff_sy, coeff_sz, coeff_sx;
    T overlap_ratio;

    if (overlap_tmp < volume_alpha)
      {
	overlap_ratio = overlap_tmp / volume_alpha;
	coeff_sy = sigmay_beta / sigmay_alpha * (1-overlap_ratio);
	coeff_sx = sigmax_beta / sigmax_alpha * (1-overlap_ratio);
	coeff_sz = sigmaz_beta / sigmaz_alpha * (1-overlap_ratio);
      }
    else
      {
	coeff_sy = 0.;
	coeff_sx = 0.;
	coeff_sz = 0.;
      }

    sigmax_out = coeff_sx * sigmax_alpha + sigmax_alpha;
    sigmay_out = coeff_sy * sigmay_alpha + sigmay_alpha;
    sigmaz_out = coeff_sz * sigmaz_alpha + sigmaz_alpha;

    puff_alpha->SetSigma_x(sigmax_out);
    puff_alpha->SetSigma_y(sigmay_out);
    puff_alpha->SetSigma_z(sigmaz_out);

    T quantity_alpha, quantity_alpha_aer, quantity_alpha_number;
    T quantity_beta, quantity_beta_aer, quantity_beta_number;

    Array<T, 1> quantity_out(this->Ns);
    Array<T, 1> quantity_out_number(this->Nbin_aer);
    Array<T, 2> quantity_out_aer(this->Ns_aer, this->Nbin_aer);

    int s, b;
    for (s = 0; s < this->Ns; s++)
      {
	quantity_alpha = puff_alpha->GetQuantity(s);
	quantity_beta = puff_beta->GetQuantity(s);
	quantity_out(s) = quantity_alpha + quantity_beta;
	puff_alpha->SetQuantity(quantity_out(s), s);
      }

    for (s = 0; s < this->Ns_aer; s++)
      for (b = 0; b < this->Nbin_aer; b++)
	{
	  quantity_alpha_aer = puff_alpha->GetQuantity(s,b);
	  quantity_beta_aer = puff_beta->GetQuantity(s,b);
	  quantity_out_aer(s,b) = quantity_alpha_aer + quantity_beta_aer;
	  puff_alpha->SetQuantity(quantity_out_aer(s,b),s, b);
	}

    for (b = 0; b < this->Nbin_aer; b++)
      {
	quantity_alpha_number = puff_alpha->GetNumberQuantity(b);
	quantity_beta_number = puff_beta->GetNumberQuantity(b);
	quantity_out_number(b) = quantity_alpha_number + quantity_beta_number;
	puff_alpha->SetNumberQuantity(quantity_out_number(b), b);
      }
  }

  //! Computes the integral of a puff over a given volume.
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
  ::ComputePuffIntegral_aer(Puff<T>* puff, int s, int b,
                            T x, T y, T z,
                            T lx, T ly, T lz)
  {
    // sqrt(2.)
    const T sqrt_2 = 1.41421356237;

    // Puff quantity.
    T quantity = puff->GetQuantity(s, b);
    T overcamp_factor = puff->GetReflectionFactor(s, b);
    T concentration = 0.;

    // Puff center coordinates.
    T z_c;
    T y_;
    T x_;

    T z_bl = puff->GetZ();
    T z_above = puff->GetHeightAboveBL();
    if (z_bl != 0.)
      z_c = z_bl;
    else
      z_c = z_above;

    x_ = puff->GetX();
    y_ = puff->GetY();
    

    // Conversion.
    T x_c = x_ * this->cos_angle_ + y_ * this->sin_angle_;
    T y_c = y_ * this->cos_angle_ - x_ * this->sin_angle_;

    // Standard deviations.
    T sigma_x = puff->GetSigma_x();
    T sigma_y = puff->GetSigma_y();
    T sigma_z = puff->GetSigma_z();
    T initial_sigma_x_2;
    T initial_sigma_y_2;
    T initial_sigma_z_2;
    
    if (this->option_process["with_adms_dispersion_formula"])
      {
	initial_sigma_x_2 = 0.;
	initial_sigma_y_2 = 0.;
	initial_sigma_z_2 = 0.;
      }
    else
      {
	initial_sigma_x_2 = puff->GetInitialSigma_x();
	initial_sigma_y_2 = puff->GetInitialSigma_y();
	initial_sigma_z_2 = puff->GetInitialSigma_z();
      }
    
    sigma_x = sqrt(sigma_x * sigma_x + initial_sigma_x_2);
    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
    sigma_z = sqrt(sigma_z * sigma_z + initial_sigma_z_2);



    // Integration bounds.
    T x1, x2, y1, y2, z1, z2;
    x1 = x * this->cos_angle_ + y * this->sin_angle_ - lx / 2.;
    x2 = x * this->cos_angle_ + y * this->sin_angle_ + lx / 2.;

    y1 = y * this->cos_angle_ - x * this->sin_angle_ - ly / 2.;
    y2 = y * this->cos_angle_ - x * this->sin_angle_ + ly / 2.;

    z1 = max(z - lz / 2., 0.);
    z2 = z + lz / 2.;

    // Downwind contribution.
    T Fx1 = 0.5 * erf((x1 - x_c) / (sqrt_2 * sigma_x));
    T Fx2 = 0.5 * erf((x2 - x_c) / (sqrt_2 * sigma_x));

    // Crosswind contribution.
    T Fy1 = 0.5 * erf((y1 - y_c) / (sqrt_2 * sigma_y));
    T Fy2 = 0.5 * erf((y2 - y_c) / (sqrt_2 * sigma_y));

    // Vertical contribution.
    T Fz1 = 0.5 * erf((z1 - z_c)  / (sqrt_2 * sigma_z));
    T Fz2 = 0.5 * erf((z2 - z_c)  / (sqrt_2 * sigma_z));

    // Ground reflection: plume is below the inversion layer.
    T FzG1 = 0.;
    T FzG2 = 0.;

    T inversion_height_ = this->inversion_height_;
    if (inversion_height_ == 0. || z_c <= inversion_height_)
      if (z_c - sigma_z <= 0.)
	{
	  FzG1 = 0.5 * overcamp_factor * erf((z1 + z_c)
					     / (sqrt_2 * sigma_z));
	  FzG2 = 0.5 * overcamp_factor * erf((z2 + z_c)
					     / (sqrt_2 * sigma_z));
	}
    // If there is an inversion layer.
    if (inversion_height_ > 0.)
      {
	T FzI2 = 0.;
	T FzI1 = 0.;
	// When puff is below inversion layer.
	if (z_c <= inversion_height_)
	  {
	    // Reflection occurs only when puff touches the inversion layer.
	    if (z_c + sigma_z >= inversion_height_)
	      {
		T dzi = z2 + z_c - 2. * inversion_height_;
		T dzig = z2 - z_c + 2. * inversion_height_;
		T dzgi = z2 - z_c - 2. * inversion_height_;
		
		FzI2 += 0.5 * erf(dzi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzgi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzig / (sqrt_2 * sigma_z));

		dzi = z1 + z_c - 2. * inversion_height_;
		dzig = z1 - z_c + 2. * inversion_height_;
		dzgi = z1 - z_c - 2. * inversion_height_;
	
		FzI1 += 0.5 * erf(dzi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzgi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzig / (sqrt_2 * sigma_z));
	      }

	    // Plume only impacts below the inversion layer.
	    if (z <= inversion_height_)
	      {
		// Far field.
		if (sigma_z > 1.5 * inversion_height_)
		  concentration =  quantity * (Fx2 - Fx1)
		    * (Fy2 - Fy1) / (lx * ly * inversion_height_);
		else
		  concentration =  quantity
		    * (Fx2 - Fx1) * (Fy2 - Fy1)
		    * ((Fz2 - Fz1) + (FzG2 - FzG1) + (FzI2 - FzI1))
		    / (lx * ly * lz);
	      }
	    else
	      concentration = 0.;
	  }

	// When puff is above inversion layer.
	if (z_c > inversion_height_)
	  {
	    // Reflection on the inversion layer.
	    T dzi = z1 + z_c - 2. * inversion_height_;
	    FzI1 += 0.5 * erf(dzi / (sqrt_2 * sigma_z));
	    dzi = z2 + z_c - 2. * inversion_height_;
	    FzI2 += 0.5 * erf(dzi / (sqrt_2 * sigma_z));

	    // Plume only impacts above the inversion layer.
	    if (z >= inversion_height_)
	      concentration =  quantity
		* (Fx2 - Fx1) * (Fy2 - Fy1)
		* ((Fz2 - Fz1) + (FzI2 - FzI1)) / (lx * ly * lz);
	    else
	      concentration = 0.;
	  }
      }
    else
      {
	concentration =  quantity
	  * (Fx2 - Fx1) * (Fy2 - Fy1)
	  * ((Fz2 - Fz1) + (FzG2 - FzG1)) / (lx * ly * lz);
      }
    return concentration;
  }

  //! Computes the integral of a puff over a given volume.
  template<class T, class ClassChemistry>
  T GaussianPuffAerosol<T, ClassChemistry>
    ::ComputePuffIntegral_number(Puff<T>* puff, int b, 
                              T x, T y, T z,
                              T lx, T ly, T lz)
  {
    // sqrt(2.)
    const T sqrt_2 = 1.41421356237;

    // Puff quantity.
    T quantity = puff->GetNumberQuantity(b);
    T overcamp_factor = 1.;
    T concentration = 0.;

    // Puff center coordinates.
    T z_c;
    T z_bl = puff->GetZ();
    T z_above = puff->GetHeightAboveBL();
    if (z_bl != 0.)
      z_c = z_bl;
    else
      z_c = z_above;

    T y_ = puff->GetY();
    T x_ = puff->GetX();

    // Conversion.
    T x_c = x_ * this->cos_angle_ + y_ * this->sin_angle_;
    T y_c = y_ * this->cos_angle_ - x_ * this->sin_angle_;

    // Standard deviations.
    T sigma_x = puff->GetSigma_x();
    T sigma_y = puff->GetSigma_y();
    T sigma_z = puff->GetSigma_z();
    T initial_sigma_x_2;
    T initial_sigma_y_2;
    T initial_sigma_z_2;

    if (this->option_process["with_adms_dispersion_formula"])
      {
	initial_sigma_x_2 = 0.;
	initial_sigma_y_2 = 0.;
	initial_sigma_z_2 = 0.;
      }
    else
      {
	initial_sigma_x_2 = puff->GetInitialSigma_x();
	initial_sigma_y_2 = puff->GetInitialSigma_y();
	initial_sigma_z_2 = puff->GetInitialSigma_z();
      }

    sigma_x = sqrt(sigma_x * sigma_x + initial_sigma_x_2);
    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
    sigma_z = sqrt(sigma_z * sigma_z + initial_sigma_z_2);

    // Integration bounds.
    T x1, x2, y1, y2, z1, z2;
    x1 = x * this->cos_angle_ + y * this->sin_angle_ - lx / 2.;
    x2 = x * this->cos_angle_ + y * this->sin_angle_ + lx / 2.;

    y1 = y * this->cos_angle_ - x * this->sin_angle_ - ly / 2.;
    y2 = y * this->cos_angle_ - x * this->sin_angle_ + ly / 2.;

    z1 = max(z - lz / 2., 0.);
    z2 = z + lz / 2.;

    // Downwind contribution.
    T Fx1 = 0.5 * erf((x1 - x_c) / (sqrt_2 * sigma_x));
    T Fx2 = 0.5 * erf((x2 - x_c) / (sqrt_2 * sigma_x));

    // Crosswind contribution.
    T Fy1 = 0.5 * erf((y1 - y_c) / (sqrt_2 * sigma_y));
    T Fy2 = 0.5 * erf((y2 - y_c) / (sqrt_2 * sigma_y));

    // Vertical contribution.
    T Fz1 = 0.5 * erf((z1 - z_c)  / (sqrt_2 * sigma_z));
    T Fz2 = 0.5 * erf((z2 - z_c)  / (sqrt_2 * sigma_z));

    // Ground reflection: plume is below the inversion layer.
    T FzG1 = 0.;
    T FzG2 = 0.;

    T inversion_height_ = this->inversion_height_;
    if (inversion_height_ == 0. || z_c <= inversion_height_)
      if (z_c - sigma_z <= 0.)
	{
	  FzG1 = 0.5 * overcamp_factor * erf((z1 + z_c)
					     / (sqrt_2 * sigma_z));
	  FzG2 = 0.5 * overcamp_factor * erf((z2 + z_c)
					     / (sqrt_2 * sigma_z));
	}
    // If there is an inversion layer.
    if (inversion_height_ > 0.)
      {
	T FzI2 = 0.;
	T FzI1 = 0.;
	// When puff is below inversion layer.
	if (z_c <= inversion_height_)
	  {
	    // Reflection occurs only when puff touches the inversion layer.
	    if (z_c + sigma_z >= inversion_height_)
	      {
		T dzi = z2 + z_c - 2. * inversion_height_;
		T dzig = z2 - z_c + 2. * inversion_height_;
		T dzgi = z2 - z_c - 2. * inversion_height_;
		
		FzI2 += 0.5 * erf(dzi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzgi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzig / (sqrt_2 * sigma_z));

		dzi = z1 + z_c - 2. * inversion_height_;
		dzig = z1 - z_c + 2. * inversion_height_;
		dzgi = z1 - z_c - 2. * inversion_height_;
		
		FzI1 += 0.5 * erf(dzi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzgi / (sqrt_2 * sigma_z))
		  + 0.5 * erf(dzig / (sqrt_2 * sigma_z));
	      }

	    // Plume only impacts below the inversion layer.
	    if (z <= inversion_height_)
	      {
		// Far field.
		if (sigma_z > 1.5 * inversion_height_)
		  concentration =  quantity * (Fx2 - Fx1)
		    * (Fy2 - Fy1) / (lx * ly * inversion_height_);
		else
		  concentration =  quantity
		    * (Fx2 - Fx1) * (Fy2 - Fy1)
		    * ((Fz2 - Fz1) + (FzG2 - FzG1) + (FzI2 - FzI1))
		    / (lx * ly * lz);
	      }
	    else
	      concentration = 0.;
	  }

	// When puff is above inversion layer.
	if (z_c > inversion_height_)

	  {
	    // Reflection on the inversion layer.
	    T dzi = z1 + z_c - 2. * inversion_height_;
	    FzI1 += 0.5 * erf(dzi / (sqrt_2 * sigma_z));
	    dzi = z2 + z_c - 2. * inversion_height_;
	    FzI2 += 0.5 * erf(dzi / (sqrt_2 * sigma_z));

	    // Plume only impacts above the inversion layer.
	    if (z >= inversion_height_)
	      concentration =  quantity
		* (Fx2 - Fx1) * (Fy2 - Fy1)
		* ((Fz2 - Fz1) + (FzI2 - FzI1)) / (lx * ly * lz);
	    else
	      concentration = 0.;
	  }
      }
    else
      {
	concentration =  quantity
	  * (Fx2 - Fx1) * (Fy2 - Fy1)
	  * ((Fz2 - Fz1) + (FzG2 - FzG1)) / (lx * ly * lz);
      }
    return concentration;
  }

  //! Save plume quantity
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::PuffOutputSaver()
  {
    if (this->option_process["save_plume"])
      {
	if (this->option_process["with_output_plume_mass"])
	  {
	    SavePlumeQuantity();
	    if (this->option_process["with_number_concentration"] &&
		this->option_process["with_output_plume_number"])
	      SavePlumeNumberQuantity();
	  }
      }
    else
      {
	SavePuffQuantity();
	if (this->option_process["with_number_concentration"])
	  SavePuffNumberQuantity();
      }
  }

  //! Save plume quantity, one output file per source
   template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::SavePlumeQuantity()
  {
    typename list<Puff<T>* >::iterator iter;
    typename list<Puff<T>* >::reverse_iterator iter2;
    Puff<T>* puff;
    int Npuff, ipuff;
    Npuff = 100;
    
    vector<string> id_list;

    T x_p, y_p, z_p_c, z_p_above, z_p;
    T lat_p, lon_p;
    string puff_id;
    const T pi = 3.14159265358979323846264;
    const T earth_radius = 6371229.;
    
    if (this->PuffList.size() > 0)
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    int Nid_tmp = id_list.size();
	    int l = 0;
	    
	    puff_id = puff->GetSourceId();
	    if (Nid_tmp == 0)
	      id_list.push_back(puff_id);
	    else
	      while (l < Nid_tmp && puff_id != id_list[l])
		l++;
	    if (l == Nid_tmp && Nid_tmp != 0)
	      id_list.push_back(puff_id);		    
	  }
	    
	for (int l = 0; l < int(id_list.size()); l++)
	  {
	    Data<T, 2> puff_coordinate_out(Npuff, 6);
	    Data<T, 2> puff_mass_out(Npuff, this->Ns);
	    Data<T, 3> puff_mass_out_aer(Npuff, this->Ns_aer, this->Nbin_aer);
	    
	    puff_coordinate_out.Fill(-1.0);
	    puff_mass_out.Fill(-1.0);
	    puff_mass_out_aer.Fill(-1.0);

	    puff_id = id_list[l];
	    	
	    ipuff = 0;
	    for (iter2 = this->PuffList.rbegin(); iter2 != this->PuffList.rend(); iter2++)
	      {
		puff = *iter2;
		if (puff->GetSourceId() == puff_id)
		  {
		    //! Get puff coordinate
		    x_p=puff->GetX();
		    y_p=puff->GetY();
		    z_p_c=puff->GetZ();
		    z_p_above=puff->GetHeightAboveBL();
		    z_p;
		    
		    if (z_p_c !=0.)
		      z_p=z_p_c;
		    else
		      z_p = z_p_above;
		    
		    lat_p= (y_p / earth_radius) * (180. / pi);
		    lon_p = x_p / (earth_radius * cos(lat_p * pi / 180.)) * (180. / pi);
		    
		    puff_coordinate_out(ipuff, 0) = z_p;
		    puff_coordinate_out(ipuff, 1) = lat_p;
		    puff_coordinate_out(ipuff, 2) = lon_p; 
		    puff_coordinate_out(ipuff, 3) = puff->GetSigma_z();
		    puff_coordinate_out(ipuff, 4) = puff->GetSigma_y();
		    puff_coordinate_out(ipuff, 5) = puff->GetSigma_x();
		    
		    for (int s = 0; s < this->Ns; s++)	     
		      puff_mass_out(ipuff, s) = puff->GetQuantity(s);
		    
		    for (int s = 0; s < this->Ns_aer ; s++)
		      for (int b = 0; b < this->Nbin_aer ; b++)
			puff_mass_out_aer(ipuff, s, b) = puff->GetQuantity(s, b);
		    ipuff++;
		  }
	      }
	    FormatBinary<float>().Append(puff_coordinate_out, directory_puff_coordinate
					 + puff_id + "-coordinate.bin");
	    FormatBinary<float>().Append(puff_mass_out, directory_puff_mass
					 + puff_id + "-gas.bin");
	    FormatBinary<float>().Append(puff_mass_out_aer, directory_puff_mass
					 + puff_id + "-aer.bin");

	  }
      }
  }

  //! Save plume quantity, one output file per source
   template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::SavePlumeNumberQuantity()
  {
    typename list<Puff<T>* >::iterator iter;
    typename list<Puff<T>* >::reverse_iterator iter2;
    Puff<T>* puff;
    int Npuff, ipuff;
    Npuff = 100;
    
    vector<string> id_list;

    T x_p, y_p, z_p_c, z_p_above, z_p;
    T lat_p, lon_p;
    string puff_id;
    const T pi = 3.14159265358979323846264;
    const T earth_radius = 6371229.;
    
    if (this->PuffList.size() > 0)
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    int Nid_tmp = id_list.size();
	    int l = 0;
	    
	    puff_id = puff->GetSourceId();
	    if (Nid_tmp == 0)
	      id_list.push_back(puff_id);
	    else
	      while (l < Nid_tmp && puff_id != id_list[l])
		l++;
	    if (l == Nid_tmp && Nid_tmp != 0)
	      id_list.push_back(puff_id);		    
	  }
	    
	for (int l = 0; l < int(id_list.size()); l++)
	  {
	    Data<T, 2> puff_number_out(Npuff, this->Nbin_aer);
	    
	    puff_number_out.Fill(-1.0);
	    puff_id = id_list[l];
	    	
	    ipuff = 0;
	    for (iter2 = this->PuffList.rbegin(); iter2 != this->PuffList.rend(); iter2++)
	      {
		puff = *iter2;
		if (puff->GetSourceId() == puff_id)
		  {		    
		    for (int b = 0; b < this->Nbin_aer ; b++)
			puff_number_out(ipuff, b) = puff->GetNumberQuantity(b);
		    ipuff++;
		  }
	      }
	    FormatBinary<float>().Append(puff_number_out, directory_puff_number
					 + puff_id + "-number.bin");

	  }
      }
  }

  //!  Saves the quantities of species in Puffs.
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::SavePuffQuantity()
  {
    typename list<Puff<T>* >::iterator iter;
    Puff<T>* puff;
    bool savingpuff = true;
    Data<T, 1> species_quantity_out(this->Ns);
    Data<T, 2> species_quantity_aer_out(this->Ns_aer, this->Nbin_aer);
    Data<T, 1> species_concentration_out_tmp(this->Ns);
    Data<T, 2> species_concentration_aer_out_tmp(this->Ns_aer, this->Nbin_aer);
    Data<T, 1> species_concentration_out(this->Ns);
    Data<T, 2> species_concentration_aer_out(this->Ns_aer, this->Nbin_aer);
    Data<T, 1> species_concentration_background_out(this->Ns);
    Data<T, 2> species_concentration_aer_background_out(this->Ns_aer, this->Nbin_aer);
    Data<T, 1> puff_coordinate_out(8);
    species_quantity_out.SetZero();
    species_quantity_aer_out.SetZero();
    puff_coordinate_out.SetZero();

    if (this->option_process["with_output_plume_concentration"])
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    T release_time = puff->GetReleaseTime();
	    string puff_id = puff->GetSourceId();
	    string file_puff_concentration = directory_puff_concentration + to_str(release_time) +
	      "-" + puff_id;
	    string file_puff_concentration_aer = directory_puff_concentration + to_str(release_time) + 
	      "-aer-" + puff_id;

	    string file_background = directory_puff_background_concentration + 
	      to_str(release_time) + "-" + puff_id;
	    string file_background_aer = directory_puff_background_concentration + 
	      to_str(release_time) + "-aer-" + puff_id;

	    T volume_puff_tmp = 1. / this->ComputeInteractionCoefficient(puff, puff);


	    if (fmod(release_time,delta_t_output) == 0){
	      for (int s = 0; s < this->Ns; s++){
		species_concentration_out_tmp(s) = puff->GetQuantity(s);
		species_concentration_background_out(s) = puff->GetBackgroundConcentration(s);
		species_concentration_out(s) = species_concentration_out_tmp(s) /
		  volume_puff_tmp + species_concentration_background_out(s);
	      }
	      FormatBinary<float>().Append(species_concentration_out, file_puff_concentration);
	      FormatBinary<float>().Append(species_concentration_background_out, file_background);
	      
	      for (int s = 0; s < this->Ns_aer ; s++)
		for (int b = 0; b < this->Nbin_aer ; b++){
		  species_concentration_aer_out_tmp(s, b) = puff->GetQuantity(s, b);
		  species_concentration_aer_background_out(s, b) = puff->GetBackgroundConcentration(s, b);
		  species_concentration_aer_out(s,b) = species_concentration_aer_out_tmp(s, b) /
		    volume_puff_tmp + species_concentration_aer_background_out(s, b);
		}
	      FormatBinary<float>().Append(species_concentration_aer_out, file_puff_concentration_aer);
	      FormatBinary<float>().Append(species_concentration_aer_background_out, file_background_aer);

	    }
	  }
      }

    if (this->option_process["with_output_plume_mass"])
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    T release_time = puff->GetReleaseTime();
	    string puff_id = puff->GetSourceId();
	    string file_puff = directory_puff_mass + to_str(release_time) + "-" + puff_id;
	    string file_puff_aer = directory_puff_mass + to_str(release_time) + "-aer-" + puff_id;
	    string file_puff_coordinate = directory_puff_coordinate + 
	      to_str(release_time) + "-" + puff_id;

	    if (fmod(release_time,delta_t_output) == 0)
	      {
		for (int s = 0; s < this->Ns; s++)
		  species_quantity_out(s) = puff->GetQuantity(s);
		FormatBinary<float>().Append(species_quantity_out, file_puff);
		
		for (int s = 0; s < this->Ns_aer ; s++)
		  for (int b = 0; b < this->Nbin_aer ; b++)
		    species_quantity_aer_out (s, b) = puff->GetQuantity(s, b);
		FormatBinary<float>().Append(species_quantity_aer_out, file_puff_aer);	
	      }
	  }
      }

    if (this->option_process["with_output_plume_coordinate"])
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    T release_time = puff->GetReleaseTime();
	    string puff_id = puff->GetSourceId();
	    string file_puff_coordinate = directory_puff_coordinate + 
	      to_str(release_time) + "-" + puff_id;

	    //! Get puff coordinate
	    T x_p=puff->GetX();
	    T y_p=puff->GetY();
	    T z_p_c=puff->GetZ();
	    T z_p_above=puff->GetHeightAboveBL();
	    T z_p;
	    
	    if (z_p_c !=0.)
	      z_p=z_p_c;
	    else
	      z_p = z_p_above;


	    if (fmod(release_time,delta_t_output) == 0){

	      //! Convert cartesian coordinate to lat/lon
	      const T pi = 3.14159265358979323846264;
	      const T earth_radius = 6371229.;
	      T lat_p= (y_p / earth_radius) * (180. / pi);
	      T lon_p = x_p / (earth_radius * cos(lat_p * pi / 180.)) * (180. / pi); 

	      //! Get Puff LWC
	      T LWC; 
	      puff->GetAdditionalMeteo(LWC);

	      puff_coordinate_out(0) = z_p;
	      puff_coordinate_out(1) = lat_p;
	      puff_coordinate_out(2) = lon_p; 
	      puff_coordinate_out(3) = puff->GetSigma_z();
	      puff_coordinate_out(4) = puff->GetSigma_y();
	      puff_coordinate_out(5) = puff->GetSigma_x();
	      puff_coordinate_out(6) = puff->GetSpecificHumidity();
	      puff_coordinate_out(7) = LWC;
	      FormatBinary<float>().Append(puff_coordinate_out, file_puff_coordinate);
	    }
	  }
      } 
  }

  //!  Saves the quantities of species in Puffs.
  template<class T, class ClassChemistry>
  void GaussianPuffAerosol<T, ClassChemistry>::SavePuffNumberQuantity()
  {
    typename list<Puff<T>* >::iterator iter;
    Puff<T>* puff;
    bool savingpuff = true; 
    Data<T, 1> species_quantity_out(this->Nbin_aer);
    Data<T, 1> species_concentration_out_tmp(this->Nbin_aer);
    Data<T, 1> species_concentration_out(this->Nbin_aer);
    Data<T, 1> species_concentration_background_out(this->Nbin_aer);
    species_quantity_out.SetZero();

 
    if (this->option_process["with_output_plume_number"])
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    T release_time = puff->GetReleaseTime();
	    string puff_id = puff->GetSourceId();
	    string file_puff = directory_puff_number + to_str(release_time) + "-" + puff_id;
	      
	    for (int b = 0; b < this->Nbin_aer ; b++)
	      species_quantity_out(b) = puff->GetNumberQuantity(b);
	    FormatBinary<float>().Append(species_quantity_out, file_puff);

	  }
      }

    if (this->option_process["with_output_plume_concentration"])
      {
	for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
	  {
	    puff = *iter;
	    T release_time = puff->GetReleaseTime();
	    string puff_id = puff->GetSourceId();
	    string file_puff_concentration = directory_puff_number_concentration + to_str(release_time) +
	      "-" + puff_id;

	    T volume_puff_tmp = 1. / this->ComputeInteractionCoefficient(puff, puff);

	    if (fmod(release_time,delta_t_output) == 0){
	      for (int b = 0; b < this->Nbin_aer; b++){
		species_concentration_out_tmp(b) = puff->GetNumberQuantity(b);
		species_concentration_background_out(b) = puff->GetBackgroundNumberConcentration(b);
		species_concentration_out(b) = species_concentration_out_tmp(b) /
		  volume_puff_tmp + species_concentration_background_out(b);
	      }
	      FormatBinary<float>().Append(species_concentration_out, file_puff_concentration);
	    }
	  }
      }

  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFAEROSOL_CXX
#endif
