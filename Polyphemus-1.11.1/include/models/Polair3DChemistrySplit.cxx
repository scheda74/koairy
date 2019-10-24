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


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRY_CXX


#include "Polair3DChemistrySplit.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Polair3DChemistry(string config_file):
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>(config_file)
  {

    /*** Managed data ***/

    this->option_manage["photolysis_rate"] = true;
    this->option_manage["attenuation"] = true;
    this->option_manage["forced_concentration"] = true;

    /*** Pointers to 3D data ***/

    this->D3_map["Attenuation"] =   &Attenuation_f;
    this->D3_map["Attenuation_i"] = &Attenuation_i;
    this->D3_map["Attenuation_f"] = &Attenuation_f;

    /*** Pointers to 4D data ***/

    this->D4_map["PhotolysisRate"] = &PhotolysisRate_f;
    this->D4_map["PhotolysisRate_i"] = &PhotolysisRate_i;
    this->D4_map["PhotolysisRate_f"] = &PhotolysisRate_f;

    this->D4_map["ForcedConcentration"] = &ForcedConcentration_f;
    this->D4_map["ForcedConcentration_i"] = &ForcedConcentration_i;
    this->D4_map["ForcedConcentration_f"] = &ForcedConcentration_f;

    this->D4_map["Source"] = &Source_f;
    this->D4_map["Source_i"] = &Source_i;
    this->D4_map["Source_f"] = &Source_f;

    /*** Pointers to 4D adjoint data ***/

    this->D4_map["Source_i_ccl"] = &Source_i_ccl;
    this->D4_map["Source_f_ccl"] = &Source_f_ccl;

    /*** Fields species lists ***/

    this->field_species["PhotolysisRate"] = &photolysis_reaction_list;
    this->field_species["PhotolysisRate_i"] = &photolysis_reaction_list;
    this->field_species["PhotolysisRate_f"] = &photolysis_reaction_list;
    this->field_species["ForcedConcentration"] = &species_list_forced;
    this->field_species["ForcedConcentration_i"] = &species_list_forced;
    this->field_species["ForcedConcentration_f"] = &species_list_forced;
  }

  
  //! Destructor.
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::~Polair3DChemistry()
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
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ReadConfiguration()
  {

    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[options]");
    this->config.PeekValue("With_chemistry",
			   this->option_process["with_chemistry"]);
    this->config.PeekValue("With_forced_concentration",
			   this->option_process["with_forced_concentration"]);
    this->config.PeekValue("With_photolysis",
			   this->option_process["with_photolysis"]);
    this->config.PeekValue("Source_splitting", source_splitting);
    
    /*** Disables management of useless fields ***/

    // Pressure and temperature are not necessary without chemistry and air
    // density.
    this->option_manage["temperature"]
      = this->option_manage["temperature"]
      || this->option_process["with_air_density"]
      || this->option_process["with_chemistry"];
    this->option_manage["pressure"]
      = this->option_manage["pressure"]
      || this->option_process["with_air_density"]
      || this->option_process["with_chemistry"];

    // Specific humidity could be needed with scavenging.
    this->option_manage["specific_humidity"]
      = this->option_manage["specific_humidity"]
      || this->option_process["with_chemistry"];

    // Photolysis rates and attenuation coefficients are
    // not needed without chemistry.
    this->option_manage["photolysis_rate"]
      = this->option_manage["photolysis_rate"]
      && this->option_process["with_photolysis"]
      && this->option_process["with_chemistry"];
    this->option_manage["attenuation"]
      = this->option_manage["attenuation"]
      && this->option_process["with_chemistry"];

    /*** Input files ***/
    
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Forced concentrations.
    if (this->option_process["with_forced_concentration"])
      this->input_files["forced_concentration"].Read(data_description_file,
						     "forced_concentration");
    else
      this->input_files["forced_concentration"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["forced_concentration"].Begin();
	 i != this->input_files["forced_concentration"].End(); i++)
      species_list_forced.push_back(i->first);
    Ns_forced = int(species_list_forced.size());
    
    // Photolysis rates.
    if (this->option_process["with_photolysis"])
      this->input_files["photolysis_rates"].ReadFiles(data_description_file,
						      "photolysis_rates");
    else
      this->input_files["photolysis_rates"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["photolysis_rates"].Begin();
	 i != this->input_files["photolysis_rates"].End(); i++)
      photolysis_reaction_list.push_back(i->first);
    Nr_photolysis = int(photolysis_reaction_list.size());

    // Reads data description.
    if (this->option_process["with_photolysis"])
      {
	data_description_stream.SetSection("[photolysis_rates]");
	photolysis_date_min = data_description_stream.PeekValue("Date_min");
	data_description_stream.PeekValue("Delta_t", "> 0",
					  photolysis_delta_t);
	data_description_stream.PeekValue("Ndays", "> 0",
					  Nphotolysis_days);
	data_description_stream.PeekValue("Time_angle_min",
					  photolysis_time_angle_min);
	data_description_stream.PeekValue("Delta_time_angle", "> 0",
					  photolysis_delta_time_angle);
	data_description_stream.PeekValue("Ntime_angle", "> 0",
					  Nphotolysis_time_angle);
	data_description_stream.PeekValue("Latitude_min",
					  photolysis_latitude_min);
	data_description_stream.PeekValue("Delta_latitude", "> 0",
					  photolysis_delta_latitude);
	data_description_stream.PeekValue("Nlatitude", "> 0",
					  Nphotolysis_latitude);
	data_description_stream.Find("Altitudes");
	split(data_description_stream.GetLine(), altitudes_photolysis);
      }
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::CheckConfiguration()
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>
      ::CheckConfiguration();

    if (this->option_manage["specific_humidity"]
	&& this->input_files["meteo"]("SpecificHumidity").empty())
      throw string("Specific humidity is needed but no input data file was")
	+ " provided.";
    if (this->option_manage["attenuation"]
	&& this->input_files["meteo"]("Attenuation").empty())
      throw "Attenuation is needed but no input data file was provided.";

    if (this->option_cartesian && this->option_process["with_chemistry"])
      throw "Cartesian coordinates are not supported with chemistry.";
  }


  //! Checks whether a species has forced concentrations.
  /*!
    \param s species global index.
    \return True if the species has forced concentrations, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasForcedConcentration(int s) const
  {
    return find(species_list_forced.begin(), species_list_forced.end(),
		this->GetSpeciesName(s)) != species_list_forced.end();
  }


  //! Checks whether a species has forced concentrations.
  /*!
    \param name species name.
    \return True if the species has forced concentrations, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasForcedConcentration(string name) const
  {
    return find(species_list_forced.begin(), species_list_forced.end(), name)
      != species_list_forced.end();
  }


  //! Returns the index in forced concentrations of a given species.
  /*!
    \param s species global index.
    \return The species index in forced concentrations.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ForcedConcentrationIndex(int s) const
  {
    return find(species_list_forced.begin(), species_list_forced.end(),
		this->GetSpeciesName(s)) - species_list_forced.begin();
  }


  //! Returns the index in forced concentrations of a given species.
  /*!
    \param name species name.
    \return The species index in forced concentrations.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ForcedConcentrationIndex(string name) const
  {
    return find(species_list_forced.begin(), species_list_forced.end(), name)
      - species_list_forced.begin();
  }


  //! Returns the name of a species with forced concentrations.
  /*!
    \param s species index in forced concentrations.
    \return The species name.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  string Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ForcedConcentrationName(int s) const
  {
    return species_list_forced.at(s);
  }


  //! Returns the name of a species with forced concentrations.
  /*!
    \param s species index in forced concentrations.
    \return The species name.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ForcedConcentrationGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(ForcedConcentrationName(s));
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Allocate()
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Allocate();

    /*** Photolysis rates ***/
    
    if (this->option_process["with_photolysis"])
      {
	GridR_photolysis = RegularGrid<T>(Nr_photolysis);
    
	PhotolysisRate_i.Resize(GridR_photolysis, this->GridZ4D,
				this->GridY4D, this->GridX4D);
	PhotolysisRate_f.Resize(GridR_photolysis, this->GridZ4D,
				this->GridY4D, this->GridX4D);

	Grid_time_angle_photolysis =
	  RegularGrid<T>(photolysis_time_angle_min,
			 photolysis_delta_time_angle,
			 Nphotolysis_time_angle);
	Grid_latitude_photolysis = RegularGrid<T>(photolysis_latitude_min,
						  photolysis_delta_latitude,
						  Nphotolysis_latitude);

	Nphotolysis_z = int(altitudes_photolysis.size());
	GridZ_photolysis = RegularGrid<T>(Nphotolysis_z);
	for (unsigned int i = 0; i < altitudes_photolysis.size(); i++)
	  GridZ_photolysis(i) = to_num<T>(altitudes_photolysis[i]);
      }

    /*** Attenuation ***/

    Attenuation_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    Attenuation_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileAttenuation_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileAttenuation_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Forced concentrations ***/

    GridS_forced = RegularGrid<T>(Ns_forced);

    ForcedConcentration_i.Resize(GridS_forced, this->GridZ4D,
				 this->GridY4D, this->GridX4D);
    ForcedConcentration_f.Resize(GridS_forced, this->GridZ4D,
				 this->GridY4D, this->GridX4D);
    FileForcedConcentration_i.Resize(GridS_forced, this->GridZ4D,
				     this->GridY4D, this->GridX4D);
    FileForcedConcentration_f.Resize(GridS_forced, this->GridZ4D,
				     this->GridY4D, this->GridX4D);

    /*** Sources (source splitting) ***/

    if (source_splitting && this->option_process["with_chemistry"])
      {
	Source_i.Resize(this->GridS4D, this->GridZ4D,
			this->GridY4D, this->GridX4D);
	Source_f.Resize(this->GridS4D, this->GridZ4D,
			this->GridY4D, this->GridX4D);
      }
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Init()
  {

    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Init();

    /*** Input fields ***/

    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::InitAllData();

    /*** Chemical mechanism ***/

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
    if (this->option_process["with_chemistry"])
      Chemistry_.Init(*this);
#endif
  }
  

  //! Initializes photolysis rates.
  /*! The rates are computed on the basis of raw photolysis-rates read in
    files. The raw photolysis-rates depend upon the day, the time angle, the
    latitude and the altitude.
    \param date the date at which photolysis rates are needed.
    \param Rates (output) the photolysis rates.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitPhotolysis(Date date, Data<T, 4>& Rates)
  {
    int i, j, k;
    T time_angle;
    int angle_in, j_in, k_in;
    T alpha_angle, alpha_y, alpha_z;
    T one_alpha_angle, one_alpha_y, one_alpha_z;
    int nb_days;

    // Relevant step in input files.
    int day = int(T(date.GetSecondsFrom(photolysis_date_min))
		  / 86400. / photolysis_delta_t + .5);
    if (day >= Nphotolysis_days)
      throw string("There are not enough available data for photolysis ")
	+ string("rates. Missing days.");

    Data<T, 3> FileRates(Grid_time_angle_photolysis,
			 Grid_latitude_photolysis, GridZ_photolysis);

    // Loop over photolysis reactions.
    for (int r = 0; r < Nr_photolysis; r++)
      {
	string filename =
	  this->input_files["photolysis_rates"](photolysis_reaction_list[r]);

	FormatBinary<float> format;
	format.ReadRecord(filename, day, FileRates);

	// Interpolation.
	for (k = 0; k < this->Nz; k++)
	  for (j = 0; j < this->Ny; j++)
	    for (i = 0; i < this->Nx; i++)
	      {
		// Along z.
		k_in = 0;
		while (k_in < Nphotolysis_z - 1 && GridZ_photolysis(k_in)
		       < this->GridZ3D.Value(k, j, i))
		  k_in++;
		if (k_in == Nphotolysis_z - 1
		    && GridZ_photolysis(k_in) < this->GridZ3D.Value(k, j, i))
		  throw string("There are not enough available data for ")
		    + string("photolysis rates. Missing levels.");
		if (k_in > 0)
		  k_in--;
		alpha_z = (this->GridZ3D.Value(k, j, i) -
			   GridZ_photolysis(k_in))
		  / (GridZ_photolysis(k_in + 1) - GridZ_photolysis(k_in));

		// Along y (latitude).
		j_in = int((this->GridY3D.Value(k, j, i)
			    - photolysis_latitude_min)
			   / photolysis_delta_latitude);
		alpha_y = (this->GridY3D.Value(k, j, i)
			   - photolysis_latitude_min
			   - T(j_in) * photolysis_delta_latitude)
		  / photolysis_delta_latitude;

		// Time angle.
		time_angle = T(date.GetNumberOfSeconds()) / 3600.
		  - 12. + this->GridX3D.Value(k, j, i) / 15.;
		nb_days = int(time_angle / 24.);
		time_angle = abs(time_angle - 24. * T(nb_days));
		if (time_angle > 12.)
		  time_angle = 24. - time_angle;
		
		angle_in = int((time_angle - photolysis_time_angle_min)
			       / photolysis_delta_time_angle);
		alpha_angle = (time_angle - photolysis_time_angle_min
			       - T(angle_in) * photolysis_delta_time_angle)
		  / photolysis_delta_time_angle;

		one_alpha_angle = 1. - alpha_angle;
		one_alpha_y = 1. - alpha_y;
		one_alpha_z = 1. - alpha_z;

		if (angle_in >= Nphotolysis_time_angle - 1)
		  Rates(r, k, j, i) = 0.;
		else
		  Rates(r, k, j, i) =
		    one_alpha_z * one_alpha_y * one_alpha_angle
		    * FileRates(angle_in, j_in, k_in)
		    + alpha_z * one_alpha_y * one_alpha_angle
		    * FileRates(angle_in, j_in, k_in + 1)
		    + one_alpha_z * alpha_y * one_alpha_angle
		    * FileRates(angle_in, j_in + 1, k_in)
		    + alpha_z * alpha_y * one_alpha_angle
		    * FileRates(angle_in, j_in + 1, k_in + 1)
		    + one_alpha_z * one_alpha_y * alpha_angle
		    * FileRates(angle_in + 1, j_in, k_in)
		    + alpha_z * one_alpha_y * alpha_angle
		    * FileRates(angle_in + 1, j_in, k_in + 1)
		    + one_alpha_z * alpha_y * alpha_angle
		    * FileRates(angle_in + 1, j_in + 1, k_in)
		    + alpha_z * alpha_y * alpha_angle
		    * FileRates(angle_in + 1, j_in + 1, k_in + 1);
	      }
      }
    
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitStep()
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitStep();

    /*** Photolysis rates ***/
    
    PhotolysisRate_i.GetArray() = PhotolysisRate_f.GetArray();
    if (this->option_manage["photolysis_rate"])
      InitPhotolysis(this->next_date, PhotolysisRate_f);

    /*** Attenuation ***/
    
    if (this->option_manage["attenuation"])
      this->UpdateData("meteo", "Attenuation", FileAttenuation_i,
		       FileAttenuation_f, Attenuation_i, Attenuation_f);
    else
      Attenuation_i.GetArray() = Attenuation_f.GetArray();

    /*** Forced concentrations ***/

    if (this->option_manage["forced_concentration"])
      for (int i = 0; i < Ns_forced; i++)
	this->UpdateData("forced_concentration", species_list_forced[i],
			 FileForcedConcentration_i, FileForcedConcentration_f,
			 i, ForcedConcentration_i, ForcedConcentration_f);
    else
      ForcedConcentration_i.GetArray() = ForcedConcentration_f.GetArray();
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SetDate(Date date)
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::SetDate(date);
    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::InitAllData();
  }


  /////////////////
  // INTEGRATION //
  /////////////////


  //! Performs one step of chemistry integration.
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Chemistry()
  {
#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
    Chemistry_.Forward(*this);
#endif
  }


  //! Performs one step of chemistry integration backward.
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Chemistry_b()
  {
#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
    Chemistry_.Backward(*this);
#endif
  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SetForward()
  {

    /*** Air density ***/

    if (this->option_process["with_air_density"])
      {
	this->InterpolateInterface_z(this->CellCenterDistance_z,
                                     this->CellWidth_z,
				     this->AirDensity_f,
                                     this->AirDensity_interf_z_f);
	this->InterpolateInterface_y(this->CellCenterDistance_y,
                                     this->CellWidth_y,
				     this->AirDensity_f,
                                     this->AirDensity_interf_y_f);
	this->InterpolateInterface_x(this->CellCenterDistance_x,
                                     this->CellWidth_x,
				     this->AirDensity_f,
                                     this->AirDensity_interf_x_f);
      }
    
    /*** Wind ***/
    
    if (this->option_manage["horizontal_wind"])
      if (!this->option_cartesian)
	{
	  this->TransformZonalWind(this->ZonalWind_i);
	  this->TransformMeridionalWind(this->MeridionalWind_i);
	}

    if (this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z,
                                  this->AirDensity_interf_x_i,
				  this->ZonalWind_i,
                                  this->AirDensity_interf_y_i,
				  this->MeridionalWind_i,
                                  this->AirDensity_interf_z_i,
				  this->VerticalWind_i);
      else
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z, this->ZonalWind_i,
				  this->MeridionalWind_i,
                                  this->VerticalWind_i);

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      {
	// Computes rho * Kz.
	if (this->option_process["with_air_density"])
	  this->VerticalDiffusionCoefficient_f.GetArray() =
	    this->AirDensity_interf_z_f.GetArray()
	    * this->VerticalDiffusionCoefficient_f.GetArray();
      }
    
    if (this->option_manage["horizontal_diffusion"]
	&& this->option_isotropic_diffusion)
      {
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->ZonalDiffusionCoefficient_i);
	this->ZonalDiffusionCoefficient_i.ThresholdMin(0.);
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->MeridionalDiffusionCoefficient_i);
	this->MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

	if (!this->option_cartesian)
	  {
	    this->TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
					  this->ZonalDiffusionCoefficient_i);
	    this->TransformMeridionalDiffusion
              (this->MeridionalDiffusionCoefficient_i);
	  }
      }

    /*** Time integration ***/

    if (source_splitting && this->option_process["with_chemistry"])
      Source_f.GetArray() = this->Concentration.GetArray();
  }


  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::AdvectionSplit()
  {
    if (this->option_process["with_advection"])
      this->Advection();

    if (this->option_process["collect_dry_flux"])
      for (int s = 0; s < this->Ns_dep; s++)
	for (int j = 0; j < this->Ny; j++)
	  for (int i = 0; i < this->Nx; i++)
	    this->DryDepositionFlux(s, j, i) = 0.5
	      * (this->DepositionVelocity_i(s, j, i)
		 + this->DepositionVelocity_f(s, j, i))
	      * this->Concentration(this->DepositionVelocityGlobalIndex(s),
				    0, j, i);
  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::DiffusionSplit()
  {
    if (this->option_process["with_diffusion"])
      this->Diffusion();

    if (this->option_process["with_point_emission"])
      this->PointEmission();
  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::DiffusionSplitXY()
  {
    if (this->option_process["with_diffusion"])
      this->DiffusionXY();
  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::DiffusionSplitZ()
  {
    if (this->option_process["with_diffusion"])
      this->DiffusionZ();

    if (this->option_process["with_point_emission"])
      this->PointEmission();

  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ChemistrySplit()
  {
    int i, j, k;

    if (source_splitting && this->option_process["with_chemistry"])
      {
	Source_i.GetArray() =
	  (this->Concentration.GetArray() - Source_f.GetArray())
	  / this->Delta_t;
	this->Concentration.GetArray() = Source_f.GetArray();
	Source_f.GetArray() = Source_i.GetArray();
	
	if (this->option_process["with_volume_emission"])
	  for (int s = 0; s < this->Ns_vol_emis; s++)
	    {
	      int gs = this->VolumeEmissionGlobalIndex(s);
	      for (k = 0; k < this->Nz_vol_emis; k++)
		for (j = 0; j < this->Ny; j++)
		  for (i = 0; i < this->Nx; i++)
		    {
		      Source_i(gs, k, j, i)
			+= this->VolumeEmission_i(s, k, j, i);
		      Source_f(gs, k, j, i)
			+= this->VolumeEmission_f(s, k, j, i);
		    }
	    }
      }

    if (this->option_process["with_chemistry"])
      Chemistry();
    else if (this->option_process["with_volume_emission"])
      for (int s = 0; s < this->Ns_vol_emis; s++)
	{
	  int gs = this->VolumeEmissionGlobalIndex(s);
	  for (k = 0; k < this->Nz_vol_emis; k++)
	    for (j = 0; j < this->Ny; j++)
	      for (i = 0; i < this->Nx; i++)
		this->Concentration(gs, k, j, i)
		  += this->Delta_t * this->VolumeEmission_i(s, k, j, i);
	}

    if(this->option_process["collect_wet_flux"])
      this->WetDepositionFlux.SetZero();
    if (this->scavenging_model != "none")
      for (int s = 0; s < this->Ns_scav; s++)
	{
	  T scavenging_ratio;
	  T cloud_height_mean;
	  for (int j = 0; j < this->Ny; j++)
	    for (int i = 0; i < this->Nx; i++)
	      if (this->Rain_i(j, i) > 0.)
		{
		  cloud_height_mean = 0.5
		    * (this->CloudHeight_i(j, i) + this->CloudHeight_f(j, i));
		  for (int k = 0; k < this->Nz; k++)
		    if (cloud_height_mean > this->GridZ4D(k))
		      {
			scavenging_ratio =
			  exp(- 0.5 *
			      this->Delta_t *
			      (this->ScavengingCoefficient_i(s, k, j, i)
			       + this->ScavengingCoefficient_f(s, k, j, i)));
			if(this->option_process["collect_wet_flux"])
			  this->WetDepositionFlux(s, j, i) +=
			    this->Concentration
			    (this->ScavengingGlobalIndex(s), k, j, i)
			    / this->Delta_t * (1.- scavenging_ratio)
			    * (this->GridZ4D_interf(k+1)
			       - this->GridZ4D_interf(k));
			this->Concentration
			  (this->ScavengingGlobalIndex(s), k, j, i)
			  *= scavenging_ratio;
		      }
		}
	}

    this->AddTime(this->Delta_t);
    this->step++;

  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ForwardSplit()
  {

    SetForward();

    AdvectionSplit();

    DiffusionSplit();

    ChemistrySplit();

  }

  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Forward()
  {
    int i, j, k;

    /*** Air density ***/

    if (this->option_process["with_air_density"])
      {
	this->InterpolateInterface_z(this->CellCenterDistance_z,
                                     this->CellWidth_z,
				     this->AirDensity_f,
                                     this->AirDensity_interf_z_f);
	this->InterpolateInterface_y(this->CellCenterDistance_y,
                                     this->CellWidth_y,
				     this->AirDensity_f,
                                     this->AirDensity_interf_y_f);
	this->InterpolateInterface_x(this->CellCenterDistance_x,
                                     this->CellWidth_x,
				     this->AirDensity_f,
                                     this->AirDensity_interf_x_f);
      }
    
    /*** Wind ***/
    
    if (this->option_manage["horizontal_wind"])
      if (!this->option_cartesian)
	{
	  this->TransformZonalWind(this->ZonalWind_i);
	  this->TransformMeridionalWind(this->MeridionalWind_i);
	}

    if (this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z,
                                  this->AirDensity_interf_x_i,
				  this->ZonalWind_i,
                                  this->AirDensity_interf_y_i,
				  this->MeridionalWind_i,
                                  this->AirDensity_interf_z_i,
				  this->VerticalWind_i);
      else
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z, this->ZonalWind_i,
				  this->MeridionalWind_i,
                                  this->VerticalWind_i);

    /*** Diffusion coefficients ***/

    if (this->option_manage["vertical_diffusion"])
      {
	// Computes rho * Kz.
	if (this->option_process["with_air_density"])
	  this->VerticalDiffusionCoefficient_f.GetArray() =
	    this->AirDensity_interf_z_f.GetArray()
	    * this->VerticalDiffusionCoefficient_f.GetArray();
      }
    
    if (this->option_manage["horizontal_diffusion"]
	&& this->option_isotropic_diffusion)
      {
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->ZonalDiffusionCoefficient_i);
	this->ZonalDiffusionCoefficient_i.ThresholdMin(0.);
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->MeridionalDiffusionCoefficient_i);
	this->MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

	if (!this->option_cartesian)
	  {
	    this->TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
					  this->ZonalDiffusionCoefficient_i);
	    this->TransformMeridionalDiffusion
              (this->MeridionalDiffusionCoefficient_i);
	  }
      }

    /*** Time integration ***/

    if (source_splitting && this->option_process["with_chemistry"])
      Source_f.GetArray() = this->Concentration.GetArray();

    if (this->option_process["with_advection"])
      this->Advection();

    if (this->option_process["collect_dry_flux"])
      for (int s = 0; s < this->Ns_dep; s++)
	for (int j = 0; j < this->Ny; j++)
	  for (int i = 0; i < this->Nx; i++)
	    this->DryDepositionFlux(s, j, i) = 0.5
	      * (this->DepositionVelocity_i(s, j, i)
		 + this->DepositionVelocity_f(s, j, i))
	      * this->Concentration(this->DepositionVelocityGlobalIndex(s),
				    0, j, i);

    if (this->option_process["with_diffusion"])
      this->Diffusion();

    if (this->option_process["with_point_emission"])
      this->PointEmission();

    if (source_splitting && this->option_process["with_chemistry"])
      {
	Source_i.GetArray() =
	  (this->Concentration.GetArray() - Source_f.GetArray())
	  / this->Delta_t;
	this->Concentration.GetArray() = Source_f.GetArray();
	Source_f.GetArray() = Source_i.GetArray();
	
	if (this->option_process["with_volume_emission"])
	  for (int s = 0; s < this->Ns_vol_emis; s++)
	    {
	      int gs = this->VolumeEmissionGlobalIndex(s);
	      for (k = 0; k < this->Nz_vol_emis; k++)
		for (j = 0; j < this->Ny; j++)
		  for (i = 0; i < this->Nx; i++)
		    {
		      Source_i(gs, k, j, i)
			+= this->VolumeEmission_i(s, k, j, i);
		      Source_f(gs, k, j, i)
			+= this->VolumeEmission_f(s, k, j, i);
		    }
	    }
      }

    if (this->option_process["with_chemistry"])
      Chemistry();
    else if (this->option_process["with_volume_emission"])
      for (int s = 0; s < this->Ns_vol_emis; s++)
	{
	  int gs = this->VolumeEmissionGlobalIndex(s);
	  for (k = 0; k < this->Nz_vol_emis; k++)
	    for (j = 0; j < this->Ny; j++)
	      for (i = 0; i < this->Nx; i++)
		this->Concentration(gs, k, j, i)
		  += this->Delta_t * this->VolumeEmission_i(s, k, j, i);
	}

    if(this->option_process["collect_wet_flux"])
      this->WetDepositionFlux.SetZero();
    if (this->scavenging_model != "none")
      for (int s = 0; s < this->Ns_scav; s++)
	{
	  T scavenging_ratio;
	  T cloud_height_mean;
	  for (int j = 0; j < this->Ny; j++)
	    for (int i = 0; i < this->Nx; i++)
	      if (this->Rain_i(j, i) > 0.)
		{
		  cloud_height_mean = 0.5
		    * (this->CloudHeight_i(j, i) + this->CloudHeight_f(j, i));
		  for (int k = 0; k < this->Nz; k++)
		    if (cloud_height_mean > this->GridZ4D(k))
		      {
			scavenging_ratio =
			  exp(- 0.5 *
			      this->Delta_t *
			      (this->ScavengingCoefficient_i(s, k, j, i)
			       + this->ScavengingCoefficient_f(s, k, j, i)));
			if(this->option_process["collect_wet_flux"])
			  this->WetDepositionFlux(s, j, i) +=
			    this->Concentration
			    (this->ScavengingGlobalIndex(s), k, j, i)
			    / this->Delta_t * (1.- scavenging_ratio)
			    * (this->GridZ4D_interf(k+1)
			       - this->GridZ4D_interf(k));
			this->Concentration
			  (this->ScavengingGlobalIndex(s), k, j, i)
			  *= scavenging_ratio;
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
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SetBackward(bool flag)
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::SetBackward(flag);

    if (this->backward && this->option_process["with_chemistry"])
      {
	if (source_splitting)
	  {
	    Source_i_ccl.Resize(this->GridS4D, this->GridZ4D,
				this->GridY4D, this->GridX4D);
	    Source_f_ccl.Resize(this->GridS4D, this->GridZ4D,
				this->GridY4D, this->GridX4D);
	    Source_i_ccl.GetArray() = T(0.);
	    Source_f_ccl.GetArray() = T(0.);
	  }
	else
	  {
	    Source_i_ccl.Resize(this->GridS_vol_emis, this->GridZ4D,
				this->GridY4D, this->GridX4D);
	    Source_f_ccl.Resize(this->GridS_vol_emis, this->GridZ4D,
				this->GridY4D, this->GridX4D);
	    Source_i_ccl.GetArray() = T(0.);
	    Source_f_ccl.GetArray() = T(0.);
	  }
      }
  }


  //! Performs one step of backward integration of adjoint model.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting. Then adjoint model of the above processes is
    integrated one step backward.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Backward()
  {
    int i, j, k;

    /*** Air density ***/

    if (this->option_process["with_air_density"])
      {
	this->InterpolateInterface_z(this->CellCenterDistance_z,
                                     this->CellWidth_z,
				     this->AirDensity_f,
                                     this->AirDensity_interf_z_f);
	this->InterpolateInterface_y(this->CellCenterDistance_y,
                                     this->CellWidth_y,
				     this->AirDensity_f,
                                     this->AirDensity_interf_y_f);
	this->InterpolateInterface_x(this->CellCenterDistance_x,
                                     this->CellWidth_x,
				     this->AirDensity_f,
                                     this->AirDensity_interf_x_f);
      }
    
    /*** Wind ***/

    if (this->option_manage["horizontal_wind"])
      if (!this->option_cartesian)
	{
	  this->TransformZonalWind(this->ZonalWind_i);
	  this->TransformMeridionalWind(this->MeridionalWind_i);
	}

    if (this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z,
                                  this->AirDensity_interf_x_i,
				  this->ZonalWind_i,
                                  this->AirDensity_interf_y_i,
				  this->MeridionalWind_i,
                                  this->AirDensity_interf_z_i,
				  this->VerticalWind_i);
      else
	this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
				  this->CellWidth_z, this->ZonalWind_i,
				  this->MeridionalWind_i,
                                  this->VerticalWind_i);

    /*** Diffusion coefficients ***/
    
    if (this->option_manage["vertical_diffusion"])
      {
	// Computes rho * Kz.
	if (this->option_process["with_air_density"])
	  this->VerticalDiffusionCoefficient_f.GetArray() =
	    this->AirDensity_interf_z_f.GetArray()
	    * this->VerticalDiffusionCoefficient_f.GetArray();
      }
    
    if (this->option_manage["horizontal_diffusion"]
	&& this->option_isotropic_diffusion)
      {
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->ZonalDiffusionCoefficient_i);
	this->ZonalDiffusionCoefficient_i.ThresholdMin(0.);
	LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->MeridionalDiffusionCoefficient_i);
	this->MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

	if (!this->option_cartesian)
	  {
	    this->TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
					  this->ZonalDiffusionCoefficient_i);
	    this->TransformMeridionalDiffusion
              (this->MeridionalDiffusionCoefficient_i);
	  }
      }

    /*** Forward integrations and trajectory generations ***/

    if (source_splitting && this->option_process["with_chemistry"])
      Source_f.GetArray() = this->Concentration.GetArray();

    // Temporal storages for adjoint coding.
    int n_tap = 3;
    Array<T, 5> conc_tap(n_tap, this->Ns, this->Nz, this->Ny, this->Nx);
    Array<T, 4> S_i(this->Ns, this->Nz, this->Ny, this->Nx);
    Array<T, 4> S_f(this->Ns, this->Nz, this->Ny, this->Nx);

    // Initialization of adjoint data.
    Source_i_ccl.GetArray() = T(0.);
    Source_f_ccl.GetArray() = T(0.);

    if (this->option_process["with_advection"])
      {
	Array<T, 4> conc(&conc_tap(0, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	conc = this->Concentration.GetArray();

	this->Advection();
      }
   
    if (this->option_process["with_diffusion"])
      {
	Array<T, 4> conc(&conc_tap(1, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	conc = this->Concentration.GetArray();

	this->Diffusion();
      }

    if (this->option_process["with_point_emission"])
      this->PointEmission();

    if (source_splitting && this->option_process["with_chemistry"])
      {
	Source_i.GetArray() =
	  (this->Concentration.GetArray() - Source_f.GetArray())
	  / this->Delta_t;
	this->Concentration.GetArray() = Source_f.GetArray();
	Source_f.GetArray() = Source_i.GetArray();
	
	if (this->option_process["with_volume_emission"])
	  for (int s = 0; s < this->Ns_vol_emis; s++)
	    {
	      int gs = this->VolumeEmissionGlobalIndex(s);
	      for (k = 0; k < this->Nz_vol_emis; k++)
		for (j = 0; j < this->Ny; j++)
		  for (i = 0; i < this->Nx; i++)
		    {
		      Source_i(gs, k, j, i)
			+= this->VolumeEmission_i(s, k, j, i);
		      Source_f(gs, k, j, i)
			+= this->VolumeEmission_f(s, k, j, i);
		    }
	    }
      }
    
    if (this->option_process["with_chemistry"])
      {
	Array<T, 4> conc(&conc_tap(2, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	conc = this->Concentration.GetArray();

	if (source_splitting)
	  {
	    S_f = Source_f.GetArray();
	    S_i = Source_i.GetArray();
	  }

	Chemistry();
      }
    else if (this->option_process["with_volume_emission"])
      for (int s = 0; s < this->Ns_vol_emis; s++)
	{
	  int gs = this->VolumeEmissionGlobalIndex(s);
	  for (k = 0; k < this->Nz_vol_emis; k++)
	    for (j = 0; j < this->Ny; j++)
	      for (i = 0; i < this->Nx; i++)
		this->Concentration(gs, k, j, i)
		  += this->Delta_t * this->VolumeEmission_i(s, k, j, i);
	}

    for (int s = 0; s < this->Ns; s++)
      if (this->HasScavenging(s) && this->scavenging_model != "none")
	{
	  int scav_s = this->ScavengingIndex(s);
	  for (int j = 0; j < this->Ny; j++)
	    for (int i = 0; i < this->Nx; i++)
	      if (this->Rain_i(j, i) > 0.)
		{
		  T CloudHeight_mean = 0.5 *
		    (this->CloudHeight_i(j, i) + this->CloudHeight_f(j, i));
		  for (int k = 0; k < this->Nz; k++)
		    if (CloudHeight_mean < this->GridZ4D(k))
		      this->Concentration(s, k, j, i) *=
			exp(- 0.5 *
			    (this->ScavengingCoefficient_i(scav_s, k, j, i)
			     + this->ScavengingCoefficient_f(scav_s, k, j, i))
			    * this->Delta_t);
		}
	}

    Array<T, 4> res(this->Ns, this->Nz, this->Ny, this->Nx);
    res = this->Concentration.GetArray();

    /*** Backward integerations ***/

    for (int s = 0; s < this->Ns; s++)
      if (this->HasScavenging(s) && this->scavenging_model != "none")
	{
	  int scav_s = this->ScavengingIndex(s);
	  for (int j = 0; j < this->Ny; j++)
	    for (int i = 0; i < this->Nx; i++)
	      if (this->Rain_i(j, i) > 0.)
		{
		  T CloudHeight_mean = 0.5 *
		    (this->CloudHeight_i(j, i) + this->CloudHeight_f(j, i));
		  for (int k = 0; k < this->Nz; k++)
		    if (CloudHeight_mean < this->GridZ4D(k))
		      this->Concentration_ccl(s, k, j, i) *=
			exp(- 0.5 *
			    (this->ScavengingCoefficient_i(scav_s, k, j, i)
			     + this->ScavengingCoefficient_f(scav_s, k, j, i))
			    * this->Delta_t);
		}
	}

    if (this->option_process["with_chemistry"])
      {
	Array<T, 4> conc(&conc_tap(2, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	this->Concentration.GetArray() = conc;

	if (source_splitting)
	  {
	    Source_f.GetArray() = S_f;
	    Source_i.GetArray() = S_i;
	  }

	Chemistry_b();
      }

    if (source_splitting && this->option_process["with_chemistry"])
      this->Concentration_ccl.GetArray() +=
	(Source_f_ccl.GetArray() + Source_i_ccl.GetArray()) / this->Delta_t;

    if (this->option_process["with_diffusion"])
      {
	Array<T, 4> conc(&conc_tap(1, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	this->Concentration.GetArray() = conc;

	this->Diffusion_b();
      }
    
    if (this->option_process["with_advection"])
      {
	Array<T, 4> conc(&conc_tap(0, 0, 0, 0, 0),
			 shape(this->Ns, this->Nz, this->Ny, this->Nx));
	this->Concentration.GetArray() = conc;

	this->Advection_b();
      }

    if (source_splitting && this->option_process["with_chemistry"])
      this->Concentration_ccl.GetArray() = this->Concentration_ccl.GetArray()
	- (Source_f_ccl.GetArray() + Source_i_ccl.GetArray()) / this->Delta_t;

    this->AddTime(this->Delta_t);
    this->step++;

    /*** Resets concentration results ***/

    this->Concentration.GetArray() = res;
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the list of photolysis reactions.
  /*!
    \return The list of photolysis reactions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  vector<string>
  Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetPhotolysisReactionList() const
  {
    return photolysis_reaction_list;
  }


  //! Returns the number of photolysis reactions.
  /*!
    \return The number of photolysis reactions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNr_photolysis() const
  {
    return Nr_photolysis;
  }


  //! Returns the number of species with volume sources.
  /*! If source splitting is used, the number of species with volume sources
    is the total number of species; otherwise, it is the number of volume
    emissions.
    \return The number of species with volume sources.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNs_source()
  {
    if (source_splitting && this->option_process["with_chemistry"])
      return this->Ns;
    else
      return this->Ns_vol_emis;
  }


  //! Returns the number of levels of volume sources.
  /*! If source splitting is used, the number of levels of volume sources
    is the total number of vertical layers; otherwise, it is the number of
    levels in volume emissions.
    \return The number of levels of volume sources.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNz_source()
  {
    if (source_splitting && this->option_process["with_chemistry"])
      return this->Nz;
    else
      return this->Nz_vol_emis;
  }


  //! Returns the global index of a species with volume sources.
  /*!
    \param s species index in volume sources.
    \return The species global index.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SourceGlobalIndex(int s) const
  {
    if (source_splitting)
      return s;
    else
      {
	string species_name = this->species_list_vol_emis[s];
	return this->GetSpeciesIndex(species_name);
      }
  }


  //! Returns the sources at current date.
  /*!
    \return The sources at current date.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 4>&
  Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetSource_i()
  {
    if (source_splitting)
      return Source_i;
    else
      return this->VolumeEmission_i;
  }


  //! Returns the sources at next date.
  /*!
    \return The sources at next date.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 4>&
  Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetSource_f()
  {
    if (source_splitting)
      return Source_f;
    else
      return this->VolumeEmission_f;
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Moves model input-data to the current date.
  /*! This method prepares the model for a time integration from the current
    date. It reads input data to related to chemistry be read before InitStep
    and Forward.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitAllData()
  {
    /*** Photolysis rates ***/

    if (this->option_manage["photolysis_rate"])
      InitPhotolysis(this->current_date, PhotolysisRate_f);

    /*** Attenuation ***/

    if (this->option_manage["attenuation"])
      this->InitData("meteo", "Attenuation", FileAttenuation_i,
		     FileAttenuation_f, this->current_date, Attenuation_f);
    else
      Attenuation_f.Fill(1.);

    /*** Forced concentrations ***/

    if (this->option_manage["forced_concentration"])
      for (int i = 0; i < Ns_forced; i++)
	this->InitData("forced_concentration", species_list_forced[i],
		       FileForcedConcentration_i, FileForcedConcentration_f,
		       this->current_date, i, ForcedConcentration_f);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRY_CXX
#endif
