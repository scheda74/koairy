// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Vivien Mallet, Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFCHEMISTRY_CXX


#include "GaussianPuffChemistry.hxx"


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
  GaussianPuffChemistry<T, ClassChemistry>
  ::GaussianPuffChemistry(string config_file):
    GaussianPuffTransport<T>(config_file)
  {
  }


  //! Destructor.
  template<class T, class ClassChemistry>
  GaussianPuffChemistry<T, ClassChemistry>::~GaussianPuffChemistry()
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
  void GaussianPuffChemistry<T, ClassChemistry>::ReadConfiguration()
  {
    GaussianPuffTransport<T>::ReadConfiguration();

    /*** Gaussian options ***/

    this->config.SetSection("[gaussian]");

    this->config.PeekValue("With_chemistry",
                           this->option_process["with_chemistry"]);

    if (this->option_process["with_chemistry"])
      {
        this->config.PeekValue("With_photolysis",
                               this->option_process["with_photolysis"]);
        this->config.PeekValue("With_forced_concentration",
                               this->option_process
                               ["with_forced_concentration"]);
        this->config.PeekValue("With_puff_interaction",
                               this->option_process
                               ["with_puff_interaction"]);
        this->config.PeekValue("With_chemical_subcycling",
                               this->option_process
                               ["with_chemical_subcycling"]);
        if (this->option_process["with_chemical_subcycling"])
          this->config.PeekValue("Ncycle", "> 0", Ncycle);
        else
          Ncycle = 1;
      }

    /*** Output options ***/

    this->config.SetSection("[output]");

    if (this->option_process["with_chemistry"])
      {
       	this->config.PeekValue("With_output_plume_mass",
	       this->option_process["with_output_plume_mass"]);
	if (this->option_process["with_output_plume_mass"])
	  this->config.PeekValue("File_mass", file_mass);
	this->config.PeekValue("With_output_plume_concentration",
	       this->option_process["with_output_plume_concentration"]);
	this->config.PeekValue("With_output_plume_coordinate",
	       this->option_process["with_output_plume_coordinate"]);
	this->config.PeekValue("With_output_plume_number",
	       this->option_process["with_output_plume_number"]);
      }
  }


  //! Allocates memory.
  /*! Allocates the grids and the concentration Data for aerosol species.
   */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::Allocate()
  {
    GaussianPuffTransport<T>::Allocate();
  }


  //! Model initialization.
  /*! It reads the configuration.
   */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::Init()
  {
    GaussianPuffTransport<T>::Init();
  }


  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::InitStep()
  {
    GaussianPuffTransport<T>::InitStep();
  }


  //! Puffs initialization.
  /*!
    For each source, it creates the corresponding puff (in case of a puff
    source) or series of puffs (for a continuous source). The continuous
    source is discretized with the model time step.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::InitPuff()
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
	T source_water = 0.;
	T volume_prev = 0.;
	string source_id;


        // Getting source plume rise parameters.
        if (this->PointEmissionManager->HasPlumeRise(emission))
          this->PointEmissionManager->
            GetPlumeRiseParam(velocity, temperature, diameter, emission);
        source_id = this->PointEmissionManager->GetEmissionSourceId(emission);
	source_water = this->PointEmissionManager->GetEmissionSourceWater(emission);
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

            // Creating one puff for each location with all species.
            Npoint = point_emission.extent(0);
            for (int index = 0; index < Npoint; index++)
              {
                Puff<T>* puff
                  = new PuffChemistry<T>(time_puff, velocity,
                                         temperature, diameter,
                                         quantity,
                                         point_emission(index, 0),
                                         point_emission(index, 1),
                                         point_emission(index, 2),
										 source_water,
										 volume_prev,
                                         this->species_list,
                                         photolysis_reaction_list,
                                         source_id);
                puff->InitPuff();
                this->PuffList.push_back(puff);
                this->Npuff++;
              }
          }
      }
    this->puff_index = 0;
    this->current_puff = this->PuffList.begin();
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////

  //! Returns whether there is chemistry used or not.
  /*!  \return true if option_chemistry is set to 'yes', false otherwise.
   */
  template<class T, class ClassChemistry>
  bool GaussianPuffChemistry<T, ClassChemistry>::WithChemistry()
  {
    return this->option_process["with_chemistry"];
  }


  //! Returns whether there is photolysis used or not.
  /*!  \return true if option_process["with_photolysis"] is set to 'yes',
    false otherwise.
  */
  template<class T, class ClassChemistry>
  bool GaussianPuffChemistry<T, ClassChemistry>::WithPhotolysis()
  {
    return this->option_process["with_photolysis"];
  }


  ////////////////////////////////////
  // BACKGROUND METEOROLOGICAL DATA //
  ///////////////////////////////////


  //! Initializes meteorological conditions.
  /*! It sets the meteorological data from a configuration file.  The
    situation is described in a dedicated section "[situation]" in which one
    finds the following entries:
    <ul>
    <li> Temperature: the temperature of ambient air (Celsius degrees).
    <li> Wind_angle: the angle between positive x-axis and wind,
    counterclockwise (degrees).
    <li> Wind: the wind velocity (m/s).
    <li> Stability: stability class in [A, F].
    </ul>
    \param meteo ConfigStream instance through which all entries (Temperature,
    Wind_angle, ...) may be read to set the meteorological situation.
    \param show_meteo indicates whether the meteorological data is to be
    displayed on screen.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {
    GaussianPuffTransport<T>::InitMeteo(meteo, show_meteo);

    /*** For chemistry ***/

    if (this->option_process["with_chemistry"])
      {
        meteo.PeekValue("Attenuation", background_attenuation_);
        meteo.PeekValue("Pressure", background_pressure_);
        meteo.PeekValue("Specific_humidity", background_specific_humidity_);

        InitChemistry(meteo);
      }
  }


  /*! Initializes photolysis parameters.
    \param Nr The number of photolysis reactions.
    \param photolysis_list  List of species with photolysis reactions.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::InitPhotolysis(int Nr, vector<string> photolysis_list)
  {
    if (this->option_process["with_chemistry"]
        && this->option_process["with_photolysis"])
      {
        Nr_photolysis = Nr;
        photolysis_reaction_list = photolysis_list;
        photolysis_rate_.resize(Nr_photolysis);
        background_concentration_.resize(this->Ns);

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
        Chemistry_.Init(*this);
#endif
      }
  }


  /*! Initializes background concentrations and photolysis parameters.
    \param meteo ConfigStream instance through which background concentrations
    and photolysis rates may be read.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::InitChemistry(ConfigStream& meteo)
  {
    // Longitude and latitude (for chemistry only).
    this->config.SetSection("[domain]");
    this->config.PeekValue("Longitude", background_longitude);
    this->config.PeekValue("Latitude", background_latitude);

    vector<string>::iterator iter;

    // Reads the list of species with photolysis.
    if (this->option_process["with_photolysis"])
      {
        ConfigStream species_stream(this->file_species);
        species_stream.SetSection("[photolysis]");
        while (!species_stream.IsEmpty())
          photolysis_reaction_list.push_back(species_stream.GetElement());
        Nr_photolysis = int(photolysis_reaction_list.size());
        photolysis_rate_.resize(Nr_photolysis);

        meteo.Find("Photolysis_rate");
        vector<string> photolysis =  split(meteo.GetLine());
        for (int i = 0; i < Nr_photolysis; i++)
          {
            iter = find(photolysis.begin(), photolysis.end(),
                        photolysis_reaction_list[i]);
            if (iter++ == photolysis.end() || iter == photolysis.end())
              throw string("Unable to find photolysis rate for")
                + string(" species \"") + photolysis_reaction_list[i] + "\".";
            photolysis_rate_(i) = to_num<T>(*iter);
          }
      }

    background_concentration_.resize(this->Ns);
    background_concentration_ = 0.;

    meteo.Find("Background_concentration");
    vector<string> concentration =  split(meteo.GetLine());
    for (iter = concentration.begin(); iter != concentration.end();
         iter++)
      for (int i = 0; i < this->Ns; i++)
        if (*iter == this->species_list[i])
          {
            iter++;
            background_concentration_(i) = to_num<T>(*iter);
          }

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
    Chemistry_.Init(*this);
#endif
  }


  //! Sets the current meteorological data to puff data, if available.
  /*! It sets the current meteorological data to puff data, if available,
    and to background meteorological data otherwise.
    \param puff the puff.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetCurrentMeteo(Puff<T>* puff)
  {
    GaussianPuffTransport<T>::SetCurrentMeteo(puff);

    if (this->option_process["with_chemistry"])
      {
        if (puff->HasMeteo())
          puff->GetAdditionalMeteo(attenuation_, pressure_,
                                   specific_humidity_);
        else
          {
            attenuation_ = background_attenuation_;
            pressure_ = background_pressure_;
            specific_humidity_ = background_specific_humidity_;
            this->longitude_ = background_longitude;
            this->latitude_ = background_latitude;
          }
      }
  }


  //! Sets the current meteorological data to puff data, if available.
  /*! It sets the current meteorological data to puff data, if available,
    and to background meteorological data otherwise.
    \param puff the puff.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetPuffCurrentMeteo(T interaction_coefficient, Puff<T>* puff)
  {
    GaussianPuffTransport<T>::SetCurrentMeteo(puff);
      
    if (this->option_process["with_chemistry"])
      {
	if (puff->HasMeteo())
	  {
	    T water_puff = 0.;
	    puff->GetAdditionalMeteo(attenuation_, pressure_,
				     specific_humidity_tmp);
	    water_puff = puff->GetEmittedWater();
	    water_puff *= interaction_coefficient;

	    specific_humidity_ = specific_humidity_tmp + (water_puff * 287. * this->temperature_ 
	    			   / pressure_);

	  }
	else
	  {
	    attenuation_ = background_attenuation_;
	    pressure_ = background_pressure_;
	    specific_humidity_ = background_specific_humidity_;
	    this->longitude_ = background_longitude;
	    this->latitude_ = background_latitude;
	  }
      }
  }





  ////////////////////////////////////////
  // ACCESS METHODS FOR PUFF ATTRIBUTES //
  ///////////////////////////////////////


  //! Sets values of additional meteorological data for a given puff.
  /*! It sets the additionnal meteorological data (needed for chemistry)
    for a given puff.
    \param index: puff index.
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetPuffAdditionalMeteo(int index,
                           T attenuation, T pressure,
                           T specific_humidity)
  {
    this->SetCurrentPuff(index);
    (*this->current_puff)->SetAdditionalMeteo(attenuation, pressure,
                                              specific_humidity);
  }


  //! Sets values of photolysis rates for a given puff.
  /*! It sets the photolysis rates for a given puff.
    \param index: puff index.
    \param photolysis_rate: array containing the photolysis rates for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetPuffPhotolysisRate(int index,
                          Array<T, 1> photolysis_rate)
  {
    this->SetCurrentPuff(index);

    if ((int) photolysis_rate.size() != Nr_photolysis)
      throw string("The number of photolysis reactions is set to ")
        + to_str<int>(Nr_photolysis) + " but the array is of size "
        + to_str<int>(photolysis_rate.size());

    for (int r = 0; r < Nr_photolysis; r++)
      (*this->current_puff)->SetPhotolysisRate(photolysis_rate(r), r);
  }


  //! Sets values of species background concentrations for a given puff.
  /*! It sets the background concentrations (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetPuffBackgroundConcentration(int index, Array<T, 1> concentration)
  {
    this->SetCurrentPuff(index);

    if ((int) concentration.size() != this->Ns)
      throw string("The number of species is set to ")
        + to_str<int>(this->Ns) + " but the array is of size "
        + to_str<int>(concentration.size());

    for (int s = 0; s < this->Ns; s++)
      (*this->current_puff)->SetBackgroundConcentration(concentration(s), s);
  }

  //! Gets values of species background concentrations for a given puff.
  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::GetPuffBackgroundConcentration(int index, int s)
  {
    this->SetCurrentPuff(index);
    
    return (*this->current_puff)->GetBackgroundConcentration(s);
  }

  //! Sets species quantity for a given puff.
  /*! It sets the species quantity (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::SetPuffQuantity(int index, Array<T, 1> quantity)
  {
    this->SetCurrentPuff(index);
    for(int s = 0; s < this->Ns; s++)
      {
	(*this->current_puff)->SetQuantity(quantity(s), s);
      }
  }


  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::GetPuffQuantity(int index, int s)
  {
    this->SetCurrentPuff(index);
    return (*this->current_puff)->GetQuantity(s);
   }




  ///////////////////////////
  // COMPUTATIONAL METHODS //
  //////////////////////////


  //!  Performs the chemistry.
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::Chemistry()
  {
    // Definition of arrays.
    Array<T, 2> quantity_list(this->Npuff, this->Ns);
    Array<T, 2> interaction_coefficient(this->Npuff, this->Npuff);
    Array<T, 1> concentration_list(this->Ns);
    Array<T, 1> overlap_volume(this->Ns);
    Array<T, 1> background_concentration(this->Ns);
    Array<T, 1> photolysis_rate(Nr_photolysis);
    Array<T, 1> source(this->Ns);
    Data<T, 1> species_quantity(this->Ns);

    // Initialization of arrays.
    source = 0.;
    quantity_list = 0.;
    interaction_coefficient = 0.;
    species_quantity.SetZero();

    // Other variables.
    int s, r, alpha, beta, i;
    typename list<Puff<T>* >::iterator iter, iter2;
    Puff<T>* puff;

    // For interactions.
    T tmp, tmp2;
    Array<vector<int>, 1 > PuffInteractionList(this->Npuff);
    vector<int> PuffList_tmp;
    int Ninteraction;

    // Getting quantities, computing interaction coefficients for all puffs.
    alpha = 0;
    for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
      {
        // List of quantities for all puffs.
        for (s = 0; s < this->Ns; s++)
          quantity_list(alpha, s) = (*iter)->GetQuantity(s);

        // Interaction coefficient for all puff pairs.
        if (this->option_process["with_puff_interaction"])
          {
            beta = 0;
            PuffList_tmp.clear();
            for (iter2 = this->PuffList.begin();
                 iter2 != this->PuffList.end(); iter2++)
              {
                tmp = ComputeInteractionCoefficient(*iter, *iter2);
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
            ComputeInteractionCoefficient(*iter, *iter);
        alpha++;
      }

    // Chemistry for all puffs.
    alpha = 0;
    for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
      {
        puff = *iter;
        // SetCurrentMeteo(puff);

	//VR

	SetPuffCurrentMeteo(interaction_coefficient(alpha, alpha), puff);
	//

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
              background_concentration(s) = background_concentration_(s);

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
              if (concentration_list(s) != 0.)
                overlap_volume(s) = quantity_list(alpha, s)
                  / concentration_list(s);
              else
                overlap_volume(s) = 1. /
                  (interaction_coefficient(alpha, alpha) * tmp2);
            else
              overlap_volume(s) = 1. / interaction_coefficient(alpha, alpha);

            // Total concentrations (puff and background).
            concentration_list(s) += background_concentration(s);
          }
        for (r = 0; r < Nr_photolysis; r++)
          if (puff->HasMeteo())
            photolysis_rate(r) = puff->GetPhotolysisRate(r);
          else
            photolysis_rate(r) = photolysis_rate_(r);

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
        // Chemistry for background concentrations only.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           this->temperature_, pressure_, source,
                           photolysis_rate,
                           T(this->next_date.
                             GetSecondsFrom(this->current_date)),
                           attenuation_, specific_humidity_,
                           this->temperature_, pressure_, source,
                           photolysis_rate, this->longitude_,
                           this->latitude_, background_concentration);

        // Chemistry for total concentrations.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           this->temperature_, pressure_, source,
                           photolysis_rate,
                           T(this->next_date.
                             GetSecondsFrom(this->current_date)),
                           attenuation_, specific_humidity_,
                           this->temperature_, pressure_, source,
                           photolysis_rate, this->longitude_,
                           this->latitude_, concentration_list);
#endif
        // New puff quantities.
        for (s = 0; s < this->Ns; s++)
          {
            T quantity = 0.;
            // Substracting background concentrations.
            concentration_list(s) -= background_concentration(s);
            quantity = concentration_list(s) * overlap_volume(s);

            // Setting puff quantity.
            puff->SetQuantity(quantity, s);
            species_quantity(s) += quantity;

            // Setting puff background concentration.
            puff->SetBackgroundConcentration(background_concentration(s), s);
            background_concentration_(s) = background_concentration(s);
          }
        alpha++;
      }
    if (this->option_process["with_output_plume_mass"])
      FormatBinary<float>().Append(species_quantity, file_mass);
  }


  //!  Performs the chemistry.
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>
  ::Chemistry(vector<list<int> > PuffCellList,
              vector<T> PuffCellVolume,
              Array<vector<T>, 1 >& PuffCellConcentration)
  {
    int Ncell = PuffCellVolume.size();
    int Ntot = this->Npuff + Ncell;

    // Definition of arrays.
    Array<T, 2> quantity_list(Ntot, this->Ns);
    Array<T, 2> interaction_coefficient(Ntot, Ntot);
    Array<T, 1> concentration_list(this->Ns);
    Array<T, 1> overlap_volume(this->Ns);
    Array<T, 1> background_concentration(this->Ns);
    Array<T, 1> photolysis_rate(Nr_photolysis);
    Array<T, 1> source(this->Ns);
    Data<T, 1> species_quantity(this->Ns);

    // Initialization of arrays.
    source = 0.;
    quantity_list = 0.;
    interaction_coefficient = 0.;
    species_quantity.SetZero();

    // Other variables.
    int s, r, alpha, beta, i;
    typename list<Puff<T>* >::iterator iter, iter2;
    Puff<T>* puff;

    // For interactions.
    T tmp, tmp2;
    Array<vector<int>, 1 > PuffInteractionList(Ntot);
    vector<int> PuffList_tmp;
    int Ninteraction;

    // Getting quantities, computing interaction coefficients for all puffs.
    alpha = 0;
    for (iter = this->PuffList.begin(); iter != this->PuffList.end(); iter++)
      {
        // List of quantities for all puffs.
        for (s = 0; s < this->Ns; s++)
          quantity_list(alpha, s) = (*iter)->GetQuantity(s);


        // Interaction coefficient for all puff pairs.
        if (this->option_process["with_puff_interaction"])
          {
            beta = 0;
            PuffList_tmp.clear();
            for (iter2 = this->PuffList.begin();
                 iter2 != this->PuffList.end(); iter2++)
              {
                tmp = ComputeInteractionCoefficient(*iter, *iter2);
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
            ComputeInteractionCoefficient(*iter, *iter);
        alpha++;
      }

    // Getting quantities, computing interaction coefficients for cells.
    typename list<int>::iterator iter3;
    list<int> pufflist_tmp;
    int index;
    for (i = 0; i < Ncell; i++)
      {
        index = this->Npuff + i;
        pufflist_tmp.clear();
        pufflist_tmp = PuffCellList[i];

        // Background quantities.
        for (s = 0; s < this->Ns; s++)
          quantity_list(index, s) = PuffCellVolume[i]
            * PuffCellConcentration(s)[i];

        // Interaction coefficient between puff and cell.
        for (alpha = 0; alpha < this->Npuff; alpha++)
          {
            iter3 = find(pufflist_tmp.begin(), pufflist_tmp.end(), alpha);
            if (iter3 != pufflist_tmp.end())
              {
                interaction_coefficient(alpha, index) = 1. / PuffCellVolume[i];
                interaction_coefficient(index, alpha) = 1. / PuffCellVolume[i];
                PuffInteractionList(alpha).push_back(index);
                PuffInteractionList(index).push_back(alpha);
              }
          }
        // Interaction coefficient between two indentical cells.
        interaction_coefficient(index, index) = 1. / PuffCellVolume[i];
        PuffInteractionList(index).push_back(index);

      }

    // Chemistry for all puffs and cells.
    int puff_index;
    T quantity;
    for (alpha = 0; alpha < Ntot; alpha++)
      {
        // Looking for the appropriate puff to take the meteorological data.
        puff_index = 0;
        if (alpha < this->Npuff)
          puff_index = alpha;
        else
          puff_index = PuffInteractionList(alpha)[0];

        this->SetCurrentPuff(puff_index);
        puff = *this->current_puff;
        SetCurrentMeteo(puff);
        overlap_volume = 0.;
        PuffList_tmp = PuffInteractionList(alpha);
        Ninteraction = PuffList_tmp.size();

        // Computing overlap concentrations and overlap volumes.
        for (s = 0; s < this->Ns; s++)
          {
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

            if (concentration_list(s) < 0.)
              concentration_list(s) = 0.;

            // Overlap volume.
            if (this->option_process["with_puff_interaction"])
              if (concentration_list(s) != 0.)
                overlap_volume(s) = quantity_list(alpha, s)
                  / concentration_list(s);
              else
                overlap_volume(s) = 1. /
                  (interaction_coefficient(alpha, alpha) * tmp2);
            else
              overlap_volume(s) = 1. / interaction_coefficient(alpha, alpha);
          }
        for (r = 0; r < Nr_photolysis; r++)
          if (puff->HasMeteo())
            photolysis_rate(r) = puff->GetPhotolysisRate(r);
          else
            photolysis_rate(r) = photolysis_rate_(r);

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
        // Chemistry for total concentrations.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                                 attenuation_, specific_humidity_,
                                 this->temperature_, pressure_, source,
                                 photolysis_rate,
                                 T(this->next_date.GetNumberOfSeconds()),
                                 attenuation_, specific_humidity_,
                                 this->temperature_, pressure_, source,
                                 photolysis_rate, this->longitude_,
                                 this->latitude_, concentration_list);
#endif
        // New puff quantities.
        if (alpha < this->Npuff)
          {
            for (s = 0; s < this->Ns; s++)
              {
                quantity = concentration_list(s) * overlap_volume(s);
                species_quantity(s) += quantity;
                puff->SetQuantity(quantity, s);
              }
          }
        // New cell quantities.
        else
          {
            i = alpha - this->Npuff;

            // Getting background concentrations (without puffs).
            for (s = 0; s < this->Ns; s++)
              background_concentration(s) = PuffCellConcentration(s)[i];

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
            // Chemistry for background concentrations.
            Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                               attenuation_, specific_humidity_,
                               this->temperature_, pressure_, source,
                               photolysis_rate,
                               T(this->next_date.GetNumberOfSeconds()),
                               attenuation_, specific_humidity_,
                               this->temperature_, pressure_, source,
                               photolysis_rate, this->longitude_,
                               this->latitude_, background_concentration);
#endif
            // Getting the concentration perturbation due to puffs.
            for (s = 0; s < this->Ns; s++)
              {
                quantity = concentration_list(s) * overlap_volume(s);
                PuffCellConcentration(s)[i] =
                  quantity * interaction_coefficient(alpha, alpha)
                  - background_concentration(s);
              }
          }
      }
    if (this->option_process["with_output_plume_mass"])
      FormatBinary<float>().Append(species_quantity, file_mass);
  }


  //! Performs one step forward.
  template<class T, class ClassChemistry>
  void GaussianPuffChemistry<T, ClassChemistry>::Forward()
  {
    this->Advection();
    this->Diffusion();

    int s, i, j, k;
    this->Concentration.SetZero();

    if (this->option_process["with_chemistry"])
      {
        Chemistry();
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              for (i = 0; i < this->Nx; i++)
                this->Concentration(s, k, j, i)
                  = background_concentration_(s);
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
      }
    this->AddTime(this->Delta_t);
    this->step++;
  }


  //! Computes concentration at a given point.
  /*!
    \return The concentration at the point.
  */
  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::GetConcentration(int species, T z, T y, T x)
  {
    T concentration = 0.;

    if (this->option_process["with_chemistry"])
      concentration = background_concentration_(species);

    this->SubtractTime(this->Delta_t);
    for (typename list<Puff<T>* >::iterator iter = this->PuffList.begin();
         iter != this->PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        if (this->current_time >= puff->GetReleaseTime() &&
            (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == species))
          {
            SetCurrentMeteo(puff);
            concentration += this->ComputePuffConcentration(puff, species, x, y, z);
          }
      }
    this->AddTime(this->Delta_t);
    return concentration;
  }


  //! Computes the mean concentration over a given volume.
  /*!
    \return The concentration over the given volume.
  */
  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::GetIntegratedConcentration(int species, T z, T y, T x,
                               T lz, T ly, T lx)
  {
    T concentration = 0.;

    if (this->option_process["with_chemistry"])
      concentration = background_concentration_(species);

    this->SubtractTime(this->Delta_t);

    for (typename list<Puff<T>* >::iterator iter = this->PuffList.begin();
         iter != this->PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        if (this->current_time >= puff->GetReleaseTime() &&
            (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == species))
          {
            SetCurrentMeteo(puff);
            concentration += this->ComputePuffIntegral(puff, species, x, y, z,
                                                       lx, ly, lz);
          }
      }
    this->AddTime(this->Delta_t);
    return concentration;
  }


  //!  Computes the interaction coefficient between two puffs.
  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::ComputeInteractionCoefficient(Puff<T>* puff_alpha,
                                  Puff<T>* puff_beta)
  {
    // 1. / (2. * sqrt(pi))
    const T coeff_pi = 0.28209479177387814;

    // Puff alpha.
    T sigmax_alpha, sigmay_alpha, sigmaz_alpha, x_alpha, y_alpha, z_alpha;
    T initial_sigma_y_2_alpha, initial_sigma_z_2_alpha, initial_sigma_x_2_alpha;
    sigmax_alpha = puff_alpha->GetSigma_x();
    sigmay_alpha = puff_alpha->GetSigma_y();
    sigmaz_alpha = puff_alpha->GetSigma_z();

    if (this->option_process["with_ADMS_dispersion_formula"])
      {
	initial_sigma_x_2_alpha = 0.;
	initial_sigma_y_2_alpha = 0.;
	initial_sigma_z_2_alpha = 0.;
      }
    else
      {
	initial_sigma_x_2_alpha = puff_alpha->GetInitialSigma_y();
	initial_sigma_y_2_alpha = puff_alpha->GetInitialSigma_y();
	initial_sigma_z_2_alpha = puff_alpha->GetInitialSigma_z();
      }

    sigmax_alpha = sqrt(sigmax_alpha * sigmax_alpha + initial_sigma_x_2_alpha);
    sigmay_alpha = sqrt(sigmay_alpha * sigmay_alpha + initial_sigma_y_2_alpha);
    sigmaz_alpha = sqrt(sigmaz_alpha * sigmaz_alpha + initial_sigma_z_2_alpha);

    x_alpha = puff_alpha->GetX();
    y_alpha = puff_alpha->GetY();
    z_alpha = puff_alpha->GetZ();

    // Puff beta.
    T sigmax_beta, sigmay_beta, sigmaz_beta, x_beta, y_beta, z_beta;
    T initial_sigma_y_2_beta, initial_sigma_z_2_beta, initial_sigma_x_2_beta;
    sigmax_beta = puff_beta->GetSigma_x();
    sigmay_beta = puff_beta->GetSigma_y();
    sigmaz_beta = puff_beta->GetSigma_z();

    if (this->option_process["with_ADMS_dispersion_formula"])
      {
	initial_sigma_x_2_beta = 0.;
	initial_sigma_y_2_beta = 0.;
	initial_sigma_z_2_beta = 0.;
      }
    else
      {
	initial_sigma_x_2_beta = puff_beta->GetInitialSigma_y();
	initial_sigma_y_2_beta = puff_beta->GetInitialSigma_y();
	initial_sigma_z_2_beta = puff_beta->GetInitialSigma_z();
      }

    sigmax_beta = sqrt(sigmax_beta * sigmax_beta + initial_sigma_x_2_beta);
    sigmay_beta = sqrt(sigmay_beta * sigmay_beta + initial_sigma_y_2_beta);
    sigmaz_beta = sqrt(sigmaz_beta * sigmaz_beta + initial_sigma_z_2_beta);

    x_beta = puff_beta->GetX();
    y_beta = puff_beta->GetY();
    z_beta = puff_beta->GetZ();

    if (abs(x_alpha - x_beta) > 2. * (sigmax_alpha + sigmax_beta)
        || abs(y_alpha - y_beta) > 2. * (sigmay_alpha + sigmay_beta)
        || abs(z_alpha - z_beta) > 2. * (sigmaz_alpha + sigmaz_beta)
        || sigmax_alpha == 0. || sigmax_beta == 0.
        || sigmay_alpha == 0. || sigmay_beta == 0.
        || sigmaz_alpha == 0. || sigmaz_beta == 0.)
      return 0.;
    else
      {
        T c_alpha_beta, coeff, alpha, beta, delta;
        T interaction_coefficient = 1.;

        // Direction x.
        alpha = 1. / (2. * sigmax_alpha * sigmax_alpha);
        beta = 1. / (2. * sigmax_beta * sigmax_beta);
        delta = x_beta - x_alpha;

        c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
        coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
        interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);

        // Direction y.
        alpha = 1. / (2. * sigmay_alpha * sigmay_alpha);
        beta = 1. / (2. * sigmay_beta * sigmay_beta);
        delta = y_beta - y_alpha;

        c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
        coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
        interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);

        // Direction z.
        alpha = 1. / (2. * sigmaz_alpha * sigmaz_alpha);
        beta = 1. / (2. * sigmaz_beta * sigmaz_beta);
        delta = z_beta - z_alpha;

        c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
        coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
        interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);

        return interaction_coefficient;
      }
  }

  //!  Computes the interaction coefficient between two puffs.
  template<class T, class ClassChemistry>
  T GaussianPuffChemistry<T, ClassChemistry>
  ::ComputePuffOverlap(int alpha,
		       int beta)
  {
    Puff<T>* puff_alpha;
    this->SetCurrentPuff(alpha);
    puff_alpha = (*this->current_puff);
    Puff<T>* puff_beta;
    this->SetCurrentPuff(beta);
    puff_beta = (*this->current_puff);

    // 1. / (2. * sqrt(pi))
    const T coeff_pi = 0.28209479177387814;
    // Puff alpha.
    T sigmax_alpha, sigmay_alpha, sigmaz_alpha, x_alpha, y_alpha, z_alpha;
    // T initial_sigma_y_2_alpha, initial_sigma_z_2_alpha, initial_sigma_x_2_alpha;
    sigmax_alpha = puff_alpha->GetSigma_x();
    sigmay_alpha = puff_alpha->GetSigma_y();
    sigmaz_alpha = puff_alpha->GetSigma_z();

    x_alpha = puff_alpha->GetX();
    y_alpha = puff_alpha->GetY();
    z_alpha = puff_alpha->GetZ();

    // Puff beta.
    T sigmax_beta, sigmay_beta, sigmaz_beta, x_beta, y_beta, z_beta;
    sigmax_beta = puff_beta->GetSigma_x();
    sigmay_beta = puff_beta->GetSigma_y();
    sigmaz_beta = puff_beta->GetSigma_z();

    x_beta = puff_beta->GetX();
    y_beta = puff_beta->GetY();
    z_beta = puff_beta->GetZ();

    if (abs(x_alpha - x_beta) > 2. * (sigmax_alpha + sigmax_beta)
	|| abs(y_alpha - y_beta) > 2. * (sigmay_alpha + sigmay_beta)
	|| abs(z_alpha - z_beta) > 2. * (sigmaz_alpha + sigmaz_beta)
	|| sigmax_alpha == 0. || sigmax_beta == 0.
	|| sigmay_alpha == 0. || sigmay_beta == 0.
	|| sigmaz_alpha == 0. || sigmaz_beta == 0.)
      return 0.;
    else
      {
	T c_alpha_beta, coeff, alpha, beta, delta;
	T interaction_coefficient = 1.;
    
	// Direction x.
	alpha = 1. / (2. * sigmax_alpha *sigmax_alpha);
	beta = 1. / (2. * sigmax_beta *sigmax_beta);
	delta = x_beta - x_alpha;
    
	c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
	coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
	interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);
    
	// Direction y.
	alpha = 1. / (2. * sigmay_alpha *sigmay_alpha);
	beta = 1. / (2. * sigmay_beta *sigmay_beta);
	delta = y_beta - y_alpha;
    
	c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
	coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
	interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);

	// Direction z.
	alpha = 1. / (2. * sigmaz_alpha *sigmaz_alpha);
	beta = 1. / (2. * sigmaz_beta *sigmaz_beta);
	delta = z_beta - z_alpha;
    
	c_alpha_beta = alpha * beta * delta * delta / (alpha + beta);
	coeff =  sqrt((4. * alpha * beta) / (alpha + beta));
	interaction_coefficient *= coeff_pi * coeff * exp(- c_alpha_beta);

	return interaction_coefficient;
      }
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFCHEMISTRY_CXX
#endif
