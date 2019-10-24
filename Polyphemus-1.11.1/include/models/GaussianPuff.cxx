// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Vivien Mallet, Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_CXX


#include "BriggsFormula.hxx"
#include "GaussianPuff.hxx"



namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds the model. Nothing else is performed.
  */
  template<class T>
  GaussianPuff<T>::GaussianPuff():
    BaseModel<T>(), Npuff(0), cos_angle_(0.), sin_angle_(0.)
  {
  }


  //! Main constructor.
  /*!
    \param config_file configuration filename.
  */
  template<class T>
  GaussianPuff<T>::GaussianPuff(string config_file):
    BaseModel<T>(config_file), Npuff(0), cos_angle_(0.), sin_angle_(0.)
  {
  }


  //! Destructor.
  template<class T>
  GaussianPuff<T>::~GaussianPuff()
  {
    ClearPuffList();
    delete this->PointEmissionManager;
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T>
  void GaussianPuff<T>::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    /*** Domain ***/

    this->config.SetSection("[domain]");

    // Land category ("rural" or "urban").
    string land_category;
    this->config.PeekValue("Land_category", "rural | urban", land_category);
    option_rural = land_category == "rural";

    // Time (Night or Day).
    string day_night;
    this->config.PeekValue("Time", "day | night", day_night);
    option_day = day_night == "day";

    /*** Gaussian options ***/

    this->config.SetSection("[gaussian]");

    this->config.PeekValue("With_radioactive_decay",
                           option_radioactive_decay);
    this->config.PeekValue("With_biological_decay",
                           option_biological_decay);
    this->config.PeekValue("With_scavenging",
                           option_scavenging);
    this->config.PeekValue("With_dry_deposition",
                           option_dry_deposition);
    this->config.PeekValue("With_plume_rise",
                           option_plume_rise);
    this->config.PeekValue("With_increasing_sigma",
                           option_increasing_sigma);
    this->config.PeekValue("With_chemistry", option_chemistry);

    if (option_chemistry)
      {
        bool option_subcycle;
        this->config.PeekValue("With_puff_interaction", option_interaction);
        this->config.PeekValue("With_chemical_subcycling", option_subcycle);
        if (option_subcycle)
          this->config.PeekValue("Ncycle", "> 0", Ncycle);
        else
          Ncycle = 1;
      }

    // Option useful only if there is chemistry.
    with_background_concentration = false;

    // Parameterization to compute standard deviations ("Briggs", "Doury").
    string sigma_formula, sigma_formula_above;
    this->config.PeekValue("Sigma_parameterization",
                           "Briggs | Doury | similarity_theory",
                           sigma_formula);
    option_briggs = sigma_formula == "Briggs";
    option_doury = sigma_formula == "Doury";

    // Parameterization for standard deviation above BL.
    this->config.PeekValue("Above_BL", "Gillani | none", sigma_formula_above);
    option_gillani = sigma_formula_above == "Gillani";

    // Alternative formulae for similarity theory (for elevated sources).
    if (sigma_formula == "similarity_theory")
      this->config.PeekValue("With_HPDM", option_hpdm);

    // Plume rise parameterization.
    if (option_plume_rise)
      {
        this->config.PeekValue("With_plume_rise_breakup", option_breakup);
        string plume_rise_formula;
        this->config.PeekValue("Plume_rise_parameterization",
                               "HPDM | Holland | Concawe",
                               plume_rise_formula);
        option_holland = plume_rise_formula == "Holland";
        option_concawe = plume_rise_formula == "Concawe";
      }

    // Gets the deposition model.
    if (option_dry_deposition)
      {
        this->config.SetSection("[deposition]");
        this->config.PeekValue("Deposition_model", "Chamberlain | Overcamp",
                               deposition_model);
      }

    if (option_chemistry)
      {
        this->config.SetSection("[output]");
        this->config.PeekValue("With_output_plume_mass", option_mass);
        if (option_mass)
          this->config.PeekValue("File_mass", file_mass);
      }

    /*** Species data ***/

    ConfigStream species_stream(this->file_species);

    // Reads the species half-lives in the species file.
    if (option_radioactive_decay)
      {
        half_life_time.resize(this->Ns);

        species_stream.SetSection("[half_life]");
        for (int i = 0; i < this->Ns; i++)
          species_stream.PeekValue(this->species_list[i], half_life_time(i));

        // Conversion from days to seconds.
        T factor = 24. * 3600.;
        for (int i = 0; i < this->Ns; i++)
          half_life_time(i) *= factor;
      }

    // Reads the species biological half-lives in the species file.
    if (option_biological_decay)
      {
        biological_half_life_time.resize(this->Ns);

        T day_value, night_value;
        for (int i = 0; i < this->Ns; i++)
          {
            species_stream.SetSection("[half_life_time]");
            species_stream.Find(this->species_list[i]);
            species_stream.GetNumber(day_value);
            species_stream.GetNumber(night_value);
            if (option_day)
              biological_half_life_time(i) = day_value;
            else
              biological_half_life_time(i) = night_value;
          }
      }
  }


  //! Allocates memory.
  /*! Allocates the grids and the concentration Data for aerosol species.
   */
  template<class T>
  void GaussianPuff<T>::Allocate()
  {
    BaseModel<T>::Allocate();
  }


  //! Model initialization.
  /*! It reads the configuration.
   */
  template<class T>
  void GaussianPuff<T>::Init()
  {
    BaseModel<T>::Init();
  }


  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T>
  void GaussianPuff<T>::InitStep()
  {
    BaseModel<T>::InitStep();
    InitPuff();
  }


  //! Puff source initialization.
  /*!
    Puff source initialization. It initializes the sources
    and creates the corresponding puffs.
  */
  template<class T>
  void GaussianPuff<T>::InitSource()
  {
    // File containing puffs data.
    string file_puff;
    this->config.SetSection("[gaussian]");
    this->config.PeekValue("File_puff", file_puff);
    this->config.PeekValue("Delta_t_puff", Delta_t_puff);

    InitSource(file_puff, Delta_t_puff);
  }


  //! Source initialization.
  /*! It sets all sources from a text file. Each source is described in a
    dedicated section "[source]" containing the following entries:
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
    \param delta_t_puff time step between two puffs.
  */
  template<class T>
  void GaussianPuff<T>::InitSource(string file_puff, T delta_t_puff)
  {
    // List of emissions.
    this->PointEmissionManager = new BasePointEmission<T>();
    BasePointEmission<T>& emissionManager = *this->PointEmissionManager;

    emissionManager.Init(file_puff, this->species_list);
    Delta_t_puff = delta_t_puff;

    int Nemis = emissionManager.GetNumberEmission();
    for (int emission = 0; emission < Nemis; emission++)
      if (emissionManager.GetEmissionType(emission) == "continuous_line")
        {
          const vector<int>& emitted_species_index =
            emissionManager.GetEmittedSpeciesIndex(emission);

          list<Array<T, 1> > coordinate_list;
          emissionManager.GetEmissionCoordinates(coordinate_list, emission);

          int Ns_emitted = emitted_species_index.size();
          vector<T> rate(Ns_emitted);
          int id_section = 0;
          typedef typename list<Array<T, 1> >::const_iterator T_Iter;
          for (T_Iter iter = coordinate_list.begin();
               iter != coordinate_list.end(); iter++, id_section++)
            {
              T width = emissionManager.GetWidth(emission, id_section);
              T velocity, temperature, diameter;
              emissionManager.GetPlumeRiseParam(velocity, temperature,
                                                diameter, emission);
              const Array<T, 1>& coordinate = *iter;
              const Date& date_beg = emissionManager.GetDateBeg(emission);
              const Date& date_end = emissionManager.GetDateEnd(emission);

              int N = 20;
              int N1 = max(floor(width), 1);
              int N2 = floor(N / N1);
              for (int i = 0; i < int(rate.size()); i++)
                rate[i] = emissionManager.GetRate(emission, i, id_section)
                  / T(N1 * N2);

              T x1 = coordinate(0);
              T y1 = coordinate(1);
              T z1 = coordinate(2);
              T x2 = coordinate(3);
              T y2 = coordinate(4);
              T z2 = coordinate(5);
              T norm = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
              T cos_theta = abs(y1 - y2) / norm;
              T sin_theta = abs(x1 - x2) / norm;

              for (int j = 0; j < N1; j++)
                {
                  T x1_step = x1, x2_step = x2, y1_step = y1, y2_step = y2;
                  if (N1 > 1)
                    {
                      T distance = width * (j / double(N1 - 1) - 0.5);
                      T x_distance = cos_theta * distance;
                      T y_distance = sin_theta * distance;
                      x1_step += x_distance;
                      x2_step += x_distance;
                      y1_step += y_distance;
                      y2_step += y_distance;
                    }

                  for (int k = 0; k < N2; k++)
                    {
                      T step_fraction = (k / double(N2 - 1));
                      T x = x2_step + (x1_step - x2_step) * step_fraction;
                      T y = y2_step + (y1_step - y2_step) * step_fraction;
                      T z = z2 + (z1 - z2) * step_fraction;
                      emissionManager.
                        AddContinuousEmission(x, y, z, rate, date_beg,
                                              date_end, emitted_species_index,
                                              diameter, velocity,
                                              temperature);
                    }
                }
            }
        }
  }


  //! Gets the sources coordinates.
  /*!
    It gets the coordinates (abscissa, ordinate, height) of all point sources.
    \param source_coordinates Array to contain all sources coordinates.
    \note For each source i, source_coordinates(i, 0) is the source abscissa,
    source_coordinates(i, 1) is the ordinate and source_coordinates(i, 2) is
    the source height.
  */
  template<class T>
  void GaussianPuff<T>::GetSourceCoordinates(Array<T, 2>& source_coordinates)
  {
    int Nemis = this->PointEmissionManager->GetNumberEmission();
    source_coordinates.resize(Nemis, 3);
    for (int emission = 0; emission < Nemis; emission++)
      {
        T abscissa, ordinate, height;
        this->PointEmissionManager->GetEmissionCoordinates(abscissa, ordinate,
                                                           height, emission);
        source_coordinates(emission, 0) = abscissa;
        source_coordinates(emission, 1) = ordinate;
        source_coordinates(emission, 2) = height;
      }
  }


  //! Sets the sources coordinates.
  /*!
    It sets the coordinates (abscissa, ordinate, height) of all point sources.
    \param source_coordinates Array that contains all new sources coordinates.
    \note For each source i, source_coordinates(i, 0) is the source abscissa,
    source_coordinates(i, 1) is the ordinate and source_coordinates(i, 2) is
    the source height.
  */  template<class T>
  void GaussianPuff<T>::SetSourceCoordinates(Array<T, 2> source_coordinates)
  {
    int Nemis = this->PointEmissionManager->GetNumberEmission();
    if (Nemis != source_coordinates.extent(0))
      throw string("The number of coordinates and sources must be equal.");
    for (int emission = 0; emission < Nemis; emission++)
      this->PointEmissionManager->
        SetEmissionCoordinates(source_coordinates(emission, 0),
                               source_coordinates(emission, 1),
                               source_coordinates(emission, 2), emission);
  }


  //! Puffs initialization.
  /*!
    For each source, it creates the corresponding puff (in case of a puff
    source) or series of puffs (for a continuous source). The continuous
    source is discretized with the model time step.
  */
  template<class T>
  void GaussianPuff<T>::InitPuff()
  {
    int Nemis = this->PointEmissionManager->GetNumberEmission();
    for (int emission = 0; emission < Nemis; emission++)
      {
        if (this->PointEmissionManager->GetEmissionType(emission)
            != "continuous_line")
          {
            vector<int> emitted_species_index =
              this->PointEmissionManager->GetEmittedSpeciesIndex(emission);
            int Ns_emitted = emitted_species_index.size();
            int Npoint;
            T time_puff;
            Array<T, 2> point_emission;
            T velocity = 0.;
            T temperature = 0.;
            T diameter = 0.;
            string source_id;

            // Getting source plume rise parameters.
            if (this->PointEmissionManager->HasPlumeRise(emission))
              this->PointEmissionManager->
                GetPlumeRiseParam(velocity, temperature, diameter, emission);

            source_id = this->PointEmissionManager->GetEmissionSourceId(emission);
            time_puff = this->current_date.GetSecondsFrom(this->Date_min);
            int time_step = int(time_puff / this->Delta_t);
            int Nt_puff = max(int(Delta_t_puff / this->Delta_t), 1);
            Date puff_next_date = this->current_date;
            puff_next_date.AddSeconds(this->Delta_t * Nt_puff);
            // If the source is emitting at time step t.
            if (this->PointEmissionManager->
                IsEmitting(this->current_date, this->next_date, emission)
                && time_step % Nt_puff == 0)
              {
                if (!option_chemistry)
                  for (int species = 0; species < Ns_emitted; species++)
                    {
                      // Getting source emission parameters for the species.
                      int s = emitted_species_index[species];
                      this->PointEmissionManager->
                        GetEmission(this->current_date, puff_next_date,
                                    species, emission, point_emission);
                      Npoint = point_emission.extent(0);
                      for (int index = 0; index < Npoint; index++)
                        {
                          Puff<T>* puff
                            = new Puff<T>(time_puff, velocity,
                                          temperature, diameter,
                                          point_emission(index, 3),
                                          point_emission(index, 0),
                                          point_emission(index, 1),
                                          point_emission(index, 2),
                                          s, source_id);
                          puff->InitPuff();
                          PuffList.push_back(puff);
                          Npuff++;
                        }
                    }
                else
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
                                                 this->species_list,
                                                 photolysis_reaction_list,
                                                 source_id);
                        puff->InitPuff();
                        PuffList.push_back(puff);
                        Npuff++;
                      }
                  }
              }
          }
      }
    puff_index = 0;
    current_puff = PuffList.begin();
  }


  //! Initializes puffs position and resets time.
  /*! It reinitializes puffs position, date, time.
   */
  template<class T>
  void GaussianPuff<T>::InitPosition()
  {

    ClearPuffList();

    // Initializes time, step.
    this->current_time = 0.;
    this->SetDate(this->Date_min);
    this->step = 0;
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Sets the simulation beginning date.
  /*!
    \param date_min the simulation beginning date.
  */
  template<class T>
  void GaussianPuff<T>::SetDateMin(Date date_min)
  {
    this->Date_min = date_min;
  }


  //! Sets the time step.
  /*!
    \param delta_t time step in seconds.
  */
  template<class T>
  void GaussianPuff<T>::SetTimeStep(T delta_t)
  {
    this->Delta_t = delta_t;
  }


  //! Sets the number of time steps.
  /*!
    \param Nt the number of time steps.
  */
  template<class T>
  void GaussianPuff<T>::SetNt(int Nt)
  {
    this->Nt = Nt;
  }


  //! Returns whether there is similarity theory used or not.
  /*!  \return true if sigma parameterization is 'similarity_theory', false
    otherwise.
  */
  template<class T>
  bool GaussianPuff<T>::WithSimilarity()
  {
    if (option_briggs || option_doury)
      return 0;
    else
      return 1;
  }


  //! Returns whether there is deposition used or not.
  /*!  \return true if option_dry_deposition is set to 'yes', false otherwise.
   */
  template<class T>
  bool GaussianPuff<T>::WithDeposition()
  {
    return option_dry_deposition;
  }


  //! Returns whether there is scavenging used or not.
  /*!  \return true if option_scavenging is set to 'yes', false otherwise.
   */
  template<class T>
  bool GaussianPuff<T>::WithScavenging()
  {
    return option_scavenging;
  }


  //! Returns whether there is chemistry used or not.
  /*!  \return true if option_chemistry is set to 'yes', false otherwise.
   */
  template<class T>
  bool GaussianPuff<T>::WithChemistry()
  {
    return option_chemistry;
  }


  //////////////////////////
  // PUFF LIST MANAGEMENT //
  /////////////////////////


  //! Returns list of puffs.
  /*!
    \return the list of puffs.
  */
  template<class T>
  const list<Puff<T>* >& GaussianPuff<T>::GetPuffList() const
  {
    return PuffList;
  }


  //! Returns the number of puffs.
  /*!
    \return the number of puffs.
  */
  template<class T>
  int GaussianPuff<T>::GetInputSourceCount() const
  {
    return PuffList.size();
  }


  //! Returns the number of puffs.
  /*!
    Returns the number of discretized source count, which equals
    to the number of input puff sources since no discretization
    is required in this particular model.
    \return the number of puffs.
  */
  template<class T>
  int GaussianPuff<T>::GetGaussianSourceCount() const
  {
    return this->GetInputSourceCount();
  }


  //! Sets the list of puffs.
  /*!
    \param pufflist the list of puffs.
  */
  template<class T>
  void GaussianPuff<T>::SetPuffList(list<Puff<T>* > pufflist)
  {
    ClearPuffList();
    PuffList = pufflist;
    Npuff = pufflist.size();
  }


  //! Sets the current puff index and iterator to a given puff.
  /*! It sets current puff index to the index given as parameter, and
    current puff points to the corresponding puff in PuffList. It returns
    the current puff.
    \param index puff index.
  */
  template<class T>
  void GaussianPuff<T>::SetCurrentPuff(int index)
  {
    if (puff_index < index)
      while (puff_index != index)
        {
          if (current_puff != PuffList.end())
            {
              puff_index ++;
              current_puff ++;
            }
          else
            throw string("Puff index: ") + to_str<int>(index)
              + " out of range.";
        }
    else if (puff_index > index)
      {
        puff_index = 0;
        current_puff = PuffList.begin();
        while (puff_index != index)
          {
            if (current_puff != PuffList.end())
              {
                puff_index ++;
                current_puff ++;
              }
            else
              throw string("Puff index: ") + to_str<int>(index)
                + " out of range.";
          }
      }
  }


  //! Returns the index corresponding to a given puff.
  /*! It returns the index corresponding to a given puff in the PuffList,
    without modifying the current puff or puff index.
    \param puff the puff.
  */
  template<class T>
  int GaussianPuff<T>
  ::GetPuffIndex(typename list<Puff<T>* >::iterator puff)
  {
    if (current_puff == puff)
      return puff_index;
    else
      {
        int index = 0;
        typename list<Puff<T>* >::iterator iter;
        iter = PuffList.begin();
        while (iter != puff)
          if (iter != PuffList.end())
            {
              index ++;
              iter ++;
            }
          else
            throw string("Puff index: ") + to_str<int>(index)
              + " out of range.";
        return index;
      }
  }


  //! Erases a given puff from the puff list.
  /*!
    \param index index of the puff to be erased.
  */
  template<class T>
  void GaussianPuff<T>::ErasePuff(int index)
  {
    SetCurrentPuff(index);
    typename list<Puff<T>* >::iterator next_puff;
    next_puff = current_puff;
    next_puff++;
    delete *current_puff;
    PuffList.erase(current_puff);
    puff_index++;
    current_puff = next_puff;
    Npuff = Npuff - 1;
  }


  //! Clears the puff list.
  /*!
   */
  template<class T>
  void GaussianPuff<T>::ClearPuffList()
  {
    int Npuff_tot = Npuff;
    for (int i = 0; i < Npuff_tot; i++)
      ErasePuff(0);
    PuffList.clear();
    Npuff = 0;
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
  template<class T>
  void GaussianPuff<T>::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {
    const T pi = 3.14159265358979323846264;

    /*** Main meteorological variables ***/

    meteo.PeekValue("Temperature", background_temperature_);
    meteo.PeekValue("Wind_angle", background_wind_angle_);
    meteo.PeekValue("Wind", background_wind_);

    background_temperature_ = background_temperature_ + 273.15;
    T wind_angle_rad = background_wind_angle_ / 180. * pi;
    background_cos_angle_ = cos(wind_angle_rad);
    background_sin_angle_ = sin(wind_angle_rad);
    meteo.PeekValue("Boundary_height", "positive",
                    background_boundary_height_);

    if (option_briggs || option_doury)
      {
        if (option_plume_rise)
          if (option_breakup && !option_holland)
            {
              meteo.PeekValue("Friction_velocity", friction_velocity_);
              meteo.PeekValue("Convective_velocity", convective_velocity_);
            }
        meteo.PeekValue("Stability", "A | B | C | D | E | F",
                        background_stability_class_);
        background_stability_class_ = upper_case(background_stability_class_);
        if (background_stability_class_ == "A")
          background_stability_ = 0;
        else if (background_stability_class_ == "B")
          background_stability_ = 1;
        else if (background_stability_class_ == "C")
          background_stability_ = 2;
        else if (background_stability_class_ == "D")
          background_stability_ = 3;
        else if (background_stability_class_ == "E")
          background_stability_ = 4;
        else if (background_stability_class_ == "F")
          background_stability_ = 5;
      }
    else
      {
        meteo.PeekValue("Friction_velocity", background_friction_velocity_);
        meteo.PeekValue("Convective_velocity",
                        background_convective_velocity_);
        meteo.PeekValue("LMO", background_lmo_);
        meteo.PeekValue("Coriolis", background_coriolis_);

        // Matching between LMO and stability classes.
        ComputeStabilityClass(background_lmo_,
                              background_stability_,
                              background_stability_class_);
      }

    if (show_meteo)
      cout << "\t" << background_temperature_
           << "\t\t" << background_wind_angle_
           << "\t\t" << background_wind_
           << "\t\t" << background_stability_class_ << endl;

    /*** Scavenging coefficients ***/

    if (option_scavenging)
      InitScavenging(meteo);

    /*** Deposition velocities ***/

    if (option_dry_deposition)
      InitDeposition(meteo);

    /*** For chemistry ***/
    if (option_chemistry)
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
  template<class T>
  void GaussianPuff<T>::InitPhotolysis(int Nr, vector<string> photolysis_list)
  {
    Nr_photolysis = Nr;
    photolysis_reaction_list = photolysis_list;
    photolysis_rate_.resize(Nr_photolysis);
    background_concentration_.resize(this->Ns);
    if (option_chemistry)
      Chemistry_.Init(*this);
  }


  /*! Initializes background concentrations and photolysis parameters.
    \param meteo ConfigStream instance through which background concentrations
    and photolysis rates may be read.
  */
  template<class T>
  void GaussianPuff<T>::InitChemistry(ConfigStream& meteo)
  {
    // Longitude and latitude (for chemistry only).
    this->config.SetSection("[domain]");
    this->config.PeekValue("Longitude", background_longitude);
    this->config.PeekValue("Latitude", background_latitude);

    // Reads the list of species with photolysis.
    ConfigStream species_stream(this->file_species);
    species_stream.SetSection("[photolysis]");
    while (!species_stream.IsEmpty())
      photolysis_reaction_list.push_back(species_stream.GetElement());
    Nr_photolysis = int(photolysis_reaction_list.size());
    photolysis_rate_.resize(Nr_photolysis);

    meteo.Find("Photolysis_rate");
    vector<string> photolysis =  split(meteo.GetLine());
    vector<string>::iterator iter;
    for (int i = 0; i < Nr_photolysis; i++)
      {
        iter = find(photolysis.begin(), photolysis.end(),
                    photolysis_reaction_list[i]);
        if (iter++ == photolysis.end() || iter == photolysis.end())
          throw string("Unable to find photolysis rate for")
            + string(" species \"") + photolysis_reaction_list[i] + "\".";
        photolysis_rate_(i) = to_num<T>(*iter);
      }

    background_concentration_.resize(this->Ns);
    background_concentration_ = 0.;

    // Background concentrations must be saved.
    with_background_concentration = true;

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
    Chemistry_.Init(*this);
  }


  /*! Initializes scavenging coefficients.
    \param meteo ConfigStream instance through which scavenging coefficients
    may be read.
  */
  template<class T>
  void GaussianPuff<T>::InitScavenging(ConfigStream& meteo)
  {
    scavenging_coefficient.resize(this->Ns);

    meteo.Find("Scavenging_coefficient");
    vector<string> scav_coefficient =  split(meteo.GetLine());
    vector<string>::iterator iter;
    for (int i = 0; i < this->Ns; i++)
      {
        iter = find(scav_coefficient.begin(), scav_coefficient.end(),
                    this->species_list[i]);
        if (iter++ == scav_coefficient.end() ||
            iter == scav_coefficient.end())
          throw string("Unable to find scavenging coefficient for")
            + string(" species \"") + this->species_list[i] + "\".";
        scavenging_coefficient(i) = to_num<T>(*iter);
      }
  }


  /*! Initializes deposition velocities.
    \param meteo ConfigStream instance through which deposition velocities
    may be read.
  */
  template<class T>
  void GaussianPuff<T>::InitDeposition(ConfigStream& meteo)
  {
    deposition_velocity.resize(this->Ns);

    meteo.Find("Deposition_velocity");
    vector<string> dep_velocity = split(meteo.GetLine());
    vector<string>::iterator iter;
    for (int  i = 0; i < this->Ns; i++)
      {
        iter = find(dep_velocity.begin(), dep_velocity.end(),
                    this->species_list[i]);
        if (iter++ == dep_velocity.end() ||
            iter == dep_velocity.end())
          throw string("Unable to find deposition velocity for")
            + string(" species \"") + this->species_list[i] + "\".";
        deposition_velocity(i) = to_num<T>(*iter);
      }
  }


  //! Sets values of meteorological parameters.
  /*! It sets the meteorological data.
    \param temperature: the temperature of ambient air (Kelvin degrees).
    \param wind_angle: the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind: the wind velocity (m/s).
    \param stability: stability class in [A, F].
  */
  template<class T>
  void GaussianPuff<T>::SetBackgroundMeteo(T temperature, T wind_angle,
                                           T wind, string stability_class,
                                           T longitude, T latitude,
                                           T boundary_height)
  {
    background_temperature_ = temperature;
    background_wind_angle_ = wind_angle;
    background_cos_angle_ = cos(background_wind_angle_);
    background_sin_angle_ = sin(background_wind_angle_);
    background_wind_ = wind;
    background_longitude = longitude;
    background_latitude = latitude;
    option_day = IsDay(longitude, latitude, this->current_date);

    background_stability_class_ = upper_case(stability_class);
    if (background_stability_class_ == "A")
      background_stability_ = 0;
    else if (background_stability_class_ == "B")
      background_stability_ = 1;
    else if (background_stability_class_ == "C")
      background_stability_ = 2;
    else if (background_stability_class_ == "D")
      background_stability_ = 3;
    else if (background_stability_class_ == "E")
      background_stability_ = 4;
    else if (background_stability_class_ == "F")
      background_stability_ = 5;
    else
      throw string("Stability class should be in [A, F], but \"")
        + background_stability_class_ + "\" was provided.";
    background_boundary_height_ = boundary_height;
  }


  //! Sets values of meteorological parameters.
  /*! It sets the meteorological data.
    \param temperature: the temperature of ambient air (Kelvin degrees).
    \param wind_angle: the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind: the wind velocity (m/s).
    \param friction_velocity: friction velocity (m/s).
    \param convective_velocity: convective velocity (m/s).
    \param boundary_height: boundary layer height (m).
    \param lmo: Monin Obukhov length (m).
    \param coriolis: Coriolis parameter.
  */
  template<class T>
  void GaussianPuff<T>::SetBackgroundMeteo(T temperature, T wind_angle,
                                           T wind, T friction_velocity,
                                           T convective_velocity,
                                           T boundary_height, T lmo,
                                           T coriolis, T longitude,
                                           T latitude)
  {
    background_temperature_ = temperature;
    background_wind_angle_ = wind_angle;
    background_cos_angle_ = cos(background_wind_angle_);
    background_sin_angle_ = sin(background_wind_angle_);
    background_wind_ = wind;
    background_friction_velocity_ = friction_velocity;
    background_convective_velocity_ = convective_velocity;
    background_boundary_height_ = boundary_height;
    background_lmo_ = lmo;
    background_coriolis_ = coriolis;
    background_longitude = longitude;
    background_latitude = latitude;
    option_day = IsDay(longitude, latitude, this->current_date);
  }


  //! Sets values of additional meteorological data.
  /*! It sets the additionnal meteorological data (needed for chemistry).
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void GaussianPuff<T>::SetBackgroundAdditionalMeteo(T attenuation,
                                                     T pressure,
                                                     T specific_humidity)
  {
    background_attenuation_ = attenuation;
    background_pressure_ = pressure;
    background_specific_humidity_ = specific_humidity;
  }


  //! Sets values of photolysis rates.
  /*! It sets the photolysis rates.
    \param photolysis_rate: array containing the photolysis rates for all
    concerned species.
  */
  template<class T>
  void GaussianPuff<T>
  ::SetBackgroundPhotolysisRate(Array<T, 1> photolysis)
  {
    if (photolysis_rate_.size() != Nr_photolysis)
      throw string("The number of photolysis reactions is set to ")
        + to_str<int>(Nr_photolysis) + " but the array is of size "
        + to_str<int>(photolysis_rate_.size());

    photolysis_rate_ = photolysis;
  }


  //! Sets values of species background concentrations.
  /*! It sets the background concentrations (homogeneous).
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T>
  void GaussianPuff<T>
  ::SetBackgroundConcentration(Array<T, 1> concentration)
  {
    if (concentration.size() != this->Ns)
      throw string("The number of species is set to ")
        + to_str<int>(this->Ns) + " but the array is of size "
        + to_str<int>(concentration.size());

    background_concentration_ = concentration;
  }


  //! Sets the current meteorological data to puff data, if available.
  /*! It sets the current meteorological data to puff data, if available,
    and to background meteorological data otherwise.
    \param puff the puff.
  */
  template<class T>
  void GaussianPuff<T>::SetCurrentMeteo(Puff<T>* puff)
  {
    if (option_doury || option_briggs)
      {
        if (puff->HasMeteo())
          {
            puff->GetMeteo(temperature_, wind_, cos_angle_, sin_angle_,
                           stability_, boundary_height_,
                           longitude_, latitude_, isday_, rural_);
          }
        else
          {
            temperature_ = background_temperature_;
            cos_angle_ = background_cos_angle_;
            sin_angle_ = background_sin_angle_;
            wind_ = background_wind_;
            boundary_height_ = background_boundary_height_;
            stability_ = background_stability_;
            isday_ = option_day;
            rural_ = option_rural;
          }
      }
    else
      {
        if (puff->HasMeteo())
          {
            puff->GetMeteo(temperature_, wind_, cos_angle_, sin_angle_,
                           friction_velocity_, convective_velocity_,
                           boundary_height_, lmo_, coriolis_,
                           longitude_, latitude_, isday_, rural_);
          }
        else
          {
            temperature_ = background_temperature_;
            cos_angle_ = background_cos_angle_;
            sin_angle_ = background_sin_angle_;
            wind_ = background_wind_;
            friction_velocity_ = background_friction_velocity_;
            convective_velocity_ = background_convective_velocity_;
            boundary_height_ = background_boundary_height_;
            lmo_ = background_lmo_;
            coriolis_ = background_coriolis_;
            isday_ = option_day;
            rural_ = option_rural;
          }

        // Matching between LMO and stability classes.
        ComputeStabilityClass(lmo_, stability_, stability_class_);
      }
    if (option_chemistry)
      if (puff->HasMeteo())
        puff->GetAdditionalMeteo(attenuation_, pressure_,
                                 specific_humidity_);
      else
        {
          attenuation_ = background_attenuation_;
          pressure_ = background_pressure_;
          specific_humidity_ = background_specific_humidity_;
          longitude_ = background_longitude;
          latitude_ = background_latitude;
        }

    if (stability_ != 4 && stability_ != 5)
      inversion_height_ = boundary_height_;
    else
      inversion_height_ = 0.;
  }


  //! Returns true if the given option is activated.
  /*!
    \param option the option name.
    \return True if the option is activated.
  */
  template<class T>
  bool GaussianPuff<T>::IsOption(string option_name)
  {
    if (option_name == "gillani")
      return option_gillani;
    if (option_name == "hpdm")
      return option_hpdm;
    if (option_name == "rural")
      return option_rural;
    if (option_name == "day")
      return option_day;
    if (option_name == "scavenging")
      return option_scavenging;
    if (option_name == "dry_deposition")
      return option_dry_deposition;
    if (option_name == "plume_rise")
      return option_plume_rise;
    if (option_name == "plume_rise_breakup")
      return option_breakup;
    throw string("Unknown Boolean option name: \"") + option_name + "\".";
  }


  //! Gets the value of model parameter \a name.
  /*!
    \param name name of the model parameter.
    \return The value of the model parameter.
  */
  template<class T>
  T GaussianPuff<T>::GetModelParameter(string name)
  {
    if (name == "temperature")
      return temperature_;
    if (name == "wind_angle")
      return wind_angle_;
    if (name == "wind")
      return wind_;
    if (name == "inversion_height")
      return inversion_height_;
    if (name == "friction_velocity")
      return friction_velocity_;
    if (name == "convective_velocity")
      return convective_velocity_;
    if (name == "boundary_height")
      return boundary_height_;
    if (name == "lmo")
      return lmo_;
    if (name == "coriolis")
      return coriolis_;
    if (name == "cos_angle")
      return cos_angle_;
    if (name == "sin_angle")
      return sin_angle_;
    if (name == "Nt")
      return this->Nt;
    throw string("Unknown numerical value: \"") + name + "\".";
  }


  //! Sets model parameter \a name to \a value.
  /*!
    \param name name of the model parameter.
    \param value the new value for the model parameter.
  */
  template<class T>
  void GaussianPuff<T>::SetModelParameter(string name, T value)
  {
    if (name == "temperature")
      temperature_ = value;
    else if (name == "wind_angle")
      wind_angle_ = value;
    else if (name == "wind")
      wind_ = value;
    else if (name == "inversion_height")
      inversion_height_ = value;
    else if (name == "friction_velocity")
      friction_velocity_ = value;
    else if (name == "convective_velocity")
      convective_velocity_ = value;
    else if (name == "boundary_height")
      boundary_height_ = value;
    else if (name == "lmo")
      lmo_ = value;
    else if (name == "coriolis")
      coriolis_ = value;
    else if (name == "cos_angle")
      cos_angle_ = value;
    else if (name == "sin_angle")
      sin_angle_ = value;
    else if (name == "Nt")
      this->Nt = value;
    else
      throw string("Unknown model parameter name: \"") + name + "\".";
  }


  ////////////////////////////////////////
  // ACCESS METHODS FOR PUFF ATTRIBUTES //
  ///////////////////////////////////////


  //! Sets values of meteorological parameters for a given puff.
  /*! It sets the meteorological data for a given puff.
    \param index: puff index.
    \param temperature: the temperature of ambient air (Kelvin degrees).
    \param wind_angle: the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind: the wind velocity (m/s).
    \param stability: stability class in [A, F].
    \param longitude the longitude coordinate.
    \param latitude the latitude coordinate.
    \param isday: boolean equal to true if it is daytime, false otherwise.
    \param rural true in rural environment.
    \param inversion_height: inversion height (meters). Default value: 0.
  */
  template<class T>
  void GaussianPuff<T>::SetGaussianMeteo(int index,
                                         T temperature, T wind_angle, T wind,
                                         string stability_class,
                                         T longitude, T latitude,
                                         bool isday, bool rural,
                                         T boundary_height)
  {
    SetCurrentPuff(index);
    (*current_puff)->SetMeteo(temperature, wind_angle, wind, stability_class,
                              longitude, latitude, isday, rural,
                              boundary_height);
  }


  //! Sets values of meteorological parameters for a given puff.
  /*! It sets the meteorological data for a given puff.
    \param index: puff index.
    \param temperature: the temperature of ambient air (Kelvin degrees).
    \param wind_angle: the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind: the wind velocity (m/s).
    \param friction_velocity: friction velocity (m/s).
    \param convective_velocity: convective velocity (m/s).
    \param boundary_height: boundary layer height (m).
    \param lmo: Monin Obukhov length (m).
    \param coriolis: Coriolis parameter.
    \param longitude the longitude coordinate.
    \param latitude the latitude coordinate.
    \param isday: boolean equal to true if it is daytime, false otherwise.
    \param rural true in rural environment.
  */
  template<class T>
  void GaussianPuff<T>
  ::SetGaussianMeteo(int index, T temperature, T wind_angle, T wind,
                     T friction_velocity,
                     T convective_velocity, T boundary_height, T lmo,
                     T coriolis, T longitude, T latitude, bool isday, bool rural)
  {
    SetCurrentPuff(index);
    (*current_puff)->SetMeteo(temperature, wind_angle, wind,
                              friction_velocity, convective_velocity,
                              boundary_height, lmo, coriolis,
                              longitude, latitude, isday, rural);
  }


  //! Sets values of additional meteorological data for a given puff.
  /*! It sets the additionnal meteorological data (needed for chemistry)
    for a given puff.
    \param index: puff index.
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void GaussianPuff<T>::SetPuffAdditionalMeteo(int index,
                                               T attenuation, T pressure,
                                               T specific_humidity)
  {
    SetCurrentPuff(index);
    (*current_puff)->SetAdditionalMeteo(attenuation, pressure,
                                        specific_humidity);
  }


  //! Sets values of photolysis rates for a given puff.
  /*! It sets the photolysis rates for a given puff.
    \param index: puff index.
    \param photolysis_rate: array containing the photolysis rates for all
    concerned species.
  */
  template<class T>
  void GaussianPuff<T>::SetGaussianPhotolysisRate(int index,
                                                  const Array<T, 1>& photolysis_rate)
  {
    SetCurrentPuff(index);

    if ((int) photolysis_rate.size() != Nr_photolysis)
      throw string("The number of photolysis reactions is set to ")
        + to_str<int>(Nr_photolysis) + " but the array is of size "
        + to_str<int>(photolysis_rate.size());

    for (int r = 0; r < Nr_photolysis; r++)
      (*current_puff)->SetPhotolysisRate(photolysis_rate(r), r);
  }


  //! Sets values of deposition velocities for a given puff.
  /*! It sets the deposition velocities for a given puff.
    \param index: puff index.
    \param deposition_velocity: array containing the deposition velocities
    for all species (set to 0. if there is no deposition).
  */
  template<class T>
  void GaussianPuff<T>::SetGaussianDepositionVelocity(int index,
                                                      const Array<T, 1>& deposition_velocity)
  {
    SetCurrentPuff(index);

    if ((int) deposition_velocity.size() != this->Ns)
      throw string("The number of species with deposition is set to ")
        + to_str<int>(this->Ns) + " but the array is of size "
        + to_str<int>(deposition_velocity.size());

    for (int s = 0; s < this->Ns; s++)
      (*current_puff)->SetDepositionVelocity(deposition_velocity(s), s);
  }


  //! Sets values of scavenging coefficients for a given puff.
  /*! It sets the scavenging coefficients for a given puff.
    \param index: puff index.
    \param scavenging_coefficient: array containing the scavenging
    coefficients for all species (set to 0. if there is no scavenging).
  */
  template<class T>
  void GaussianPuff<T>::SetGaussianScavengingCoefficient(int index,
                                                         const Array<T, 1>& scavenging_coefficient)
  {
    SetCurrentPuff(index);

    if ((int) scavenging_coefficient.size() != this->Ns)
      throw string("The number of species with deposition is set to ")
        + to_str<int>(this->Ns) + " but the array is of size "
        + to_str<int>(scavenging_coefficient.size());

    for (int s = 0; s < this->Ns; s++)
      (*current_puff)->SetScavengingCoefficient(scavenging_coefficient(s), s);
  }


  //! Sets values of species background concentrations for a given puff.
  /*! It sets the background concentrations (homogeneous) for a given puff.
    \param concentration: array containing the concentrations for all
    concerned species.
  */
  template<class T>
  void GaussianPuff<T>
  ::SetGaussianBackgroundConcentration(int index,
                                       const Array<T, 1>& concentration)
  {
    SetCurrentPuff(index);

    if ((int) concentration.size() != this->Ns)
      throw string("The number of species is set to ")
        + to_str<int>(this->Ns) + " but the array is of size "
        + to_str<int>(concentration.size());

    for (int s = 0; s < this->Ns; s++)
      (*current_puff)->SetBackgroundConcentration(concentration(s), s);
  }


  //! Gets the new position of the puff center for a given puff index.
  /*!
    \param index puff index.
    \param x abscissa (meters).
    \param y ordinate (meters).
    \param z height (meters).
    \param distance distance (meters).
    \param time time since release (s).
  */  template<class T>
  void GaussianPuff<T>::GetGaussianPosition(int index, T& x, T& y, T& z,
                                            T& distance, T& time)
  {
    SetCurrentPuff(index);
    Puff<T>* puff = *current_puff;
    x = puff->GetX();
    y = puff->GetY();
    z = puff->GetZ();
    distance = puff->GetDistance();
    time = puff->GetPuffTime();
  }


  //! Returns the standard deviations for a given puff index.
  /*!
    \param index puff index.
    \param sigma_x standard deviation in the downwind direction (m).
    \param sigma_y standard deviation in the crosswind direction (m).
    \param sigma_z vertical standard deviation (m).
  */
  template<class T>
  void GaussianPuff<T>::GetGaussianSigma(int index, T& sigma_x,
                                         T& sigma_y, T& sigma_z)
  {
    SetCurrentPuff(index);
    Puff<T>* puff = *current_puff;
    sigma_x = puff->GetSigma_x();
    sigma_y = puff->GetSigma_y();
    sigma_z = puff->GetSigma_z();
  }


  //! Returns the release time of a given puff.
  /*!
    \return The puff release time (s).
  */
  template<class T>
  T GaussianPuff<T>::GetPuffReleaseTime(int index)
  {
    SetCurrentPuff(index);
    return (*current_puff)->GetReleaseTime();
  }


  //! Returns the quantity of species s in a given puff.
  /*!
    \return The quantity of species s in a given puff.
  */
  template<class T>
  T GaussianPuff<T>::GetPuffQuantity(int index, int s)
  {
    SetCurrentPuff(index);
    if ((*current_puff)->GetNs() == this->Ns
        || (*current_puff)->GetSpeciesIndex() == s)
      return (*current_puff)->GetQuantity(s);
    else
      return 0.;
  }


  ///////////////////////////
  // COMPUTATIONAL METHODS //
  //////////////////////////


  //! Computes the new position of the puff center after advection.
  template<class T>
  void GaussianPuff<T>::Advection()
  {
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        SetCurrentMeteo(*iter);
        Advection(*iter);
      }
  }


  //!  Computes the horizontal and vertical standard deviations.
  template<class T>
  void GaussianPuff<T>::Diffusion()
  {
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        SetCurrentMeteo(*iter);
        Diffusion(*iter);
      }
  }


  //! Computes the puff loss factors.
  template<class T>
  void GaussianPuff<T>::ComputeLossFactor()
  {
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        SetCurrentMeteo(*iter);
        ComputeLossFactor(*iter);
      }
  }


  //! Performs the puff chemistry.
  template<class T>
  void GaussianPuff<T>::ComputeChemistry()
  {
    // Definition of arrays.
    Array<T, 2> quantity_list(Npuff, this->Ns);
    Array<T, 2> interaction_coefficient(Npuff, Npuff);
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
    concentration_list = 0.;

    // Other variables.
    int s, r, alpha, beta, i;
    typename list<Puff<T>* >::iterator iter, iter2;
    Puff<T>* puff;

    // For interactions.
    T tmp, tmp2;
    Array<vector<int>, 1 > PuffInteractionList(Npuff);
    vector<int> PuffList_tmp;
    int Ninteraction;

    // Getting quantities, computing interaction coefficients for all puffs.
    alpha = 0;
    for (iter = PuffList.begin(); iter != PuffList.end(); iter++)
      {
        // List of quantities for all puffs.
        for (s = 0; s < this->Ns; s++)
          quantity_list(alpha, s) = (*iter)->GetQuantity(s);

        // Interaction coefficient for all puff pairs.
        if (option_interaction)
          {
            beta = 0;
            PuffList_tmp.clear();
            for (iter2 = PuffList.begin(); iter2 != PuffList.end(); iter2++)
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
    for (iter = PuffList.begin(); iter != PuffList.end(); iter++)
      {
        puff = *iter;
        SetCurrentMeteo(puff);
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
            if (option_interaction)
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
            if (option_interaction)
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

        // Chemistry for background concentrations only.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate,
                           T(this->next_date.
                             GetSecondsFrom(this->current_date)),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate, longitude_,
                           latitude_, background_concentration);

        // Chemistry for total concentrations.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate,
                           T(this->next_date.
                             GetSecondsFrom(this->current_date)),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate, longitude_,
                           latitude_, concentration_list);

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
    if (option_mass)
      FormatBinary<float>().Append(species_quantity, file_mass);
  }


  //! Performs the puff chemistry.
  template<class T>
  void GaussianPuff<T>::ComputeChemistry(const vector<list<int> >& CellPuff,
                                         const vector<T>& CellVolume,
                                         Array<vector<T>, 1>& CellConcentration)
  {
    int Ncell = CellVolume.size();
    int Ntot = Npuff + Ncell;

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
    for (iter = PuffList.begin(); iter != PuffList.end(); iter++)
      {
        // List of quantities for all puffs.
        for (s = 0; s < this->Ns; s++)
          quantity_list(alpha, s) = (*iter)->GetQuantity(s);


        // Interaction coefficient for all puff pairs.
        if (option_interaction)
          {
            beta = 0;
            PuffList_tmp.clear();
            for (iter2 = PuffList.begin(); iter2 != PuffList.end(); iter2++)
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
        index = Npuff + i;
        pufflist_tmp.clear();
        pufflist_tmp = CellPuff[i];

        // Background quantities.
        for (s = 0; s < this->Ns; s++)
          quantity_list(index, s) = CellVolume[i]
            * CellConcentration(s)[i];

        // Interaction coefficient between puff and cell.
        for (alpha = 0; alpha < Npuff; alpha++)
          {
            iter3 = find(pufflist_tmp.begin(), pufflist_tmp.end(), alpha);
            if (iter3 != pufflist_tmp.end())
              {
                interaction_coefficient(alpha, index) = 1. / CellVolume[i];
                interaction_coefficient(index, alpha) = 1. / CellVolume[i];
                PuffInteractionList(alpha).push_back(index);
                PuffInteractionList(index).push_back(alpha);
              }
          }
        // Interaction coefficient between two indentical cells.
        interaction_coefficient(index, index) = 1. / CellVolume[i];
        PuffInteractionList(index).push_back(index);

      }

    // Chemistry for all puffs and cells.
    int puff_index;
    T quantity;
    for (alpha = 0; alpha < Ntot; alpha++)
      {
        // Looking for the appropriate puff to take the meteorological data.
        puff_index = 0;
        if (alpha < Npuff)
          puff_index = alpha;
        else
          puff_index = PuffInteractionList(alpha)[0];

        SetCurrentPuff(puff_index);
        puff = *current_puff;
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
            if (option_interaction)
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
            if (option_interaction)
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

        // Chemistry for total concentrations.
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate,
                           T(this->next_date.GetNumberOfSeconds()),
                           attenuation_, specific_humidity_,
                           temperature_, pressure_, source,
                           photolysis_rate, longitude_,
                           latitude_, concentration_list);

        // New puff quantities.
        if (alpha < Npuff)
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
            i = alpha - Npuff;

            // Getting background concentrations (without puffs).
            for (s = 0; s < this->Ns; s++)
              background_concentration(s) = CellConcentration(s)[i];

            // Chemistry for background concentrations.
            Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
                               attenuation_, specific_humidity_,
                               temperature_, pressure_, source,
                               photolysis_rate,
                               T(this->next_date.GetNumberOfSeconds()),
                               attenuation_, specific_humidity_,
                               temperature_, pressure_, source,
                               photolysis_rate, longitude_,
                               latitude_, background_concentration);

            // Getting the concentration perturbation due to puffs.
            for (s = 0; s < this->Ns; s++)
              {
                quantity = concentration_list(s) * overlap_volume(s);
                CellConcentration(s)[i] =
                  quantity * interaction_coefficient(alpha, alpha)
                  - background_concentration(s);
              }
          }
      }
    if (option_mass)
      FormatBinary<float>().Append(species_quantity, file_mass);
  }


  //! Performs one step forward.
  template<class T>
  void GaussianPuff<T>::Forward()
  {
    Advection();
    Diffusion();

    int s, i, j, k;
    this->Concentration.SetZero();

    if (option_chemistry && with_background_concentration)
      {
        ComputeChemistry();
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              for (i = 0; i < this->Nx; i++)
                this->Concentration(s, k, j, i)
                  = background_concentration_(s);
      }
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        SetCurrentMeteo(puff);
        for (s = 0; s < this->Ns; s++)
          if (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == s)
            {
              // Computing loss factors.
              ComputeLossFactor(puff, s);
              // Loop on all points to compute concentration.
              for (k = 0; k < this->Nz; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration(s, k, j, i) +=
                      ComputePuffConcentration(puff, s, this->GridX4D(i),
                                               this->GridY4D(j),
                                               this->GridZ4D(k));
            }
      }
    this->AddTime(this->Delta_t);
    this->step++;
  }

  //! Computes plume rise for all newly emitted puffs.
  template<class T>
  void GaussianPuff<T>::ComputePlumeRise()
  {
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        T time = puff->GetPuffTime();
        // Newly emitted puffs only.
        if (time == 0.)
          {
            SetCurrentMeteo(puff);
            bool is_volume_source = puff->IsVolumeSource();
            if (option_plume_rise && !is_volume_source)
              ComputePuffEffectiveHeight(puff);
            else
              {
                T z_c = puff->GetZ();
                if (inversion_height_ > 0. && z_c >= inversion_height_)
                  {
                    puff->SetHeight(0.);
                    puff->SetHeightAboveBL(z_c);
                  }
              }
          }
      }
  }


  //! Computes plume rise for a given source.
  template<class T>
  void GaussianPuff<T>::ComputePuffEffectiveHeight(Puff<T>* puff)
  {
    // Computing the effective source height.
    T source_height = puff->GetZ();
    T plume_rise, velocity, temperature, diameter;
    velocity = puff->GetSourceVelocity();
    temperature = puff->GetSourceTemperature();
    diameter = puff->GetSourceDiameter();

    if (option_holland)
      plume_rise = ComputeHollandPlumeRise(temperature_, wind_,
                                           velocity, temperature,
                                           diameter);
    else if (option_concawe)
      plume_rise = ComputeConcawePlumeRise(temperature_, wind_,
                                           inversion_height_,
                                           convective_velocity_,
                                           friction_velocity_,
                                           stability_, velocity,
                                           temperature, diameter,
                                           source_height, option_breakup);
    else
      plume_rise = ComputeHPDMPlumeRise(temperature_, wind_,
                                        inversion_height_,
                                        convective_velocity_,
                                        friction_velocity_,
                                        stability_, velocity,
                                        temperature, diameter,
                                        source_height, option_breakup);

    // Computing the additional plume spread due to plume rise.
    T initial_sigma_y_2 = puff->GetInitialSigma_y();
    T initial_sigma_z_2 = puff->GetInitialSigma_z();
    T sigma_y_pr = plume_rise / 3.5;
    T sigma_z_pr;
    if (option_gillani)
      {
        T dTdz;
        ComputePotentialTemperatureGradient(stability_, dTdz);
        sigma_z_pr = 15. * exp(- 117. * dTdz);
      }
    else
      sigma_z_pr = plume_rise / 2.;
    puff->SetInitialSigma_y(initial_sigma_y_2 + sigma_y_pr * sigma_y_pr);
    puff->SetInitialSigma_z(initial_sigma_z_2 + sigma_z_pr * sigma_z_pr);

    // Computing effective height.
    T effective_height = 0.;
    T effective_height_above = 0.;
    T penetration_factor = 0.;

    if (inversion_height_ <= 0.)
      {
        // There is no inversion.
        effective_height = source_height + plume_rise;
        penetration_factor = 0.;
      }
    else if (source_height + 0.5 * plume_rise >= inversion_height_)
      {
        // Plume is above inversion layer.
        penetration_factor = 1.;
        effective_height = 0.;
        effective_height_above = source_height + plume_rise;
      }
    else if (source_height + 1.5 * plume_rise >= inversion_height_)
      {
        // Partial penetration in the inversion layer.
        penetration_factor
          = 1.5 - (inversion_height_ - source_height) / plume_rise;
        plume_rise = (0.62 + 0.38 * penetration_factor) *
          (inversion_height_ - source_height);
        effective_height = source_height + plume_rise;
        effective_height_above = source_height +
          (1.38 + 0.38 * penetration_factor)
          * (inversion_height_ - source_height);
      }
    else
      {
        // Plume is inside the boundary layer.
        effective_height = source_height + plume_rise;
        penetration_factor = 0.;
      }
    puff->SetPenetrationFactor(penetration_factor);
    puff->SetHeight(effective_height);
    puff->SetHeightAboveBL(effective_height_above);
  }


  //! Computes concentration at a given point.
  /*!
    \return The concentration at the point.
  */
  template<class T>
  T GaussianPuff<T>::GetConcentration(int species, T z, T y, T x)
  {
    T concentration = 0.;

    if (option_chemistry && with_background_concentration)
      concentration = background_concentration_(species);

    this->SubtractTime(this->Delta_t);
    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        if (this->current_time >= puff->GetReleaseTime() &&
            (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == species))
          {
            SetCurrentMeteo(puff);
            concentration += ComputePuffConcentration(puff, species, x, y, z);
          }
      }
    this->AddTime(this->Delta_t);
    return concentration;
  }


  //! Computes the mean concentration over a given volume.
  /*!
    \return The concentration over the given volume.
  */
  template<class T>
  T GaussianPuff<T>::GetIntegratedConcentration(int species, T z, T y, T x,
                                                T lz, T ly, T lx)
  {
    T concentration = 0.;

    if (option_chemistry && with_background_concentration)
      concentration = background_concentration_(species);

    this->SubtractTime(this->Delta_t);

    for (typename list<Puff<T>* >::iterator iter = PuffList.begin();
         iter != PuffList.end(); iter++)
      {
        Puff<T>* puff = *iter;
        if (this->current_time >= puff->GetReleaseTime() &&
            (puff->GetNs() == this->Ns || puff->GetSpeciesIndex() == species))
          {
            SetCurrentMeteo(puff);
            concentration += ComputePuffIntegral(puff, species, x, y, z,
                                                 lx, ly, lz);
          }
      }
    this->AddTime(this->Delta_t);
    return concentration;
  }


  //! Performs advection for a given puff.
  /*!
    \param index the puff index.
  */
  template<class T>
  void GaussianPuff<T>::Advection(int index)
  {
    SetCurrentPuff(index);
    SetCurrentMeteo(*current_puff);
    Advection(*current_puff);
  }


  //! Performs advection for a given puff.
  /*!
    \param puff the puff.
  */
  template<class T>
  void GaussianPuff<T>::Advection(Puff<T>* puff)
  {
    // Puff position.
    T x_c, y_c, z_c, distance, time;
    x_c = puff->GetX();
    y_c = puff->GetY();
    z_c = puff->GetZ();

    distance = puff->GetDistance();
    time = puff->GetPuffTime();

    // Advection.
    x_c += wind_ * cos_angle_ * this->Delta_t;
    y_c += wind_ * sin_angle_ * this->Delta_t;
    distance += wind_ * this->Delta_t;
    time += this->Delta_t;
    puff->SetPuffPosition(x_c, y_c, z_c, distance, time);
  }


  //! Performs diffusion for a given puff.
  /*!
    \param index the puff index.
  */
  template<class T>
  void GaussianPuff<T>::Diffusion(int index)
  {
    SetCurrentPuff(index);
    SetCurrentMeteo(*current_puff);
    Diffusion(*current_puff);
  }


  //! Performs diffusion for a given puff.
  /*!
    \param puff the puff.
  */
  template<class T>
  void GaussianPuff<T>::Diffusion(Puff<T>* puff)
  {
    T z, sigma_x, sigma_y, sigma_z;
    T distance_x = puff->GetDistance();
    T time = puff->GetPuffTime();
    T z_c = puff->GetZ();
    T z_above = puff->GetHeightAboveBL();
    if (z_c != 0.)
      z = z_c;
    else
      z = z_above;

    // Former standard deviations.
    sigma_x = puff->GetSigma_x();
    sigma_y = puff->GetSigma_y();
    sigma_z = puff->GetSigma_z();

    if (!option_increasing_sigma || sigma_y == 0.)
      ComputeSigma(distance_x, time, z, sigma_x, sigma_y, sigma_z);
    else
      {
        // Former distance / time with current meteo.
        T previous_distance, previous_time;
        ComputeTime(z, sigma_y, previous_distance, previous_time);

        // New distance / time to compute standard deviations.
        T next_distance, next_time;
        next_distance = previous_distance + this->Delta_t * wind_;
        next_time = previous_time + this->Delta_t;
        ComputeSigma(next_distance, next_time, z, sigma_x, sigma_y, sigma_z);
      }
    puff->SetSigma_x(sigma_x);
    puff->SetSigma_y(sigma_y);
    ComputeSigma(distance_x, time, z, sigma_x, sigma_y, sigma_z);
    puff->SetSigma_z(sigma_z);
  }


  //! Computes the species concentration due to the puff at a given point.
  /*!
    \param index puff index.
    \param x abscissa of the point where concentration is computed (m).
    \param y ordinate of the point where concentration is computed (m).
    \param z height of the point where concentration is computed (m).
    \return The concentration at (\a x, \a y, \a z).
  */
  template<class T>
  T GaussianPuff<T>::ComputePuffConcentration(int index, int s, T x, T y, T z)
  {
    SetCurrentPuff(index);
    SetCurrentMeteo(*current_puff);
    if ((*current_puff)->GetNs() == this->Ns
        || (*current_puff)->GetSpeciesIndex() == s)
      return ComputePuffConcentration(*current_puff, s, x, y, z);
    else
      return 0.;
  }


  //! Computes the species concentration due to the puff at a given point.
  /*!
    \param puff the puff.
    \param x abscissa of the point where concentration is computed (m).
    \param y ordinate of the point where concentration is computed (m).
    \param z height of the point where concentration is computed (m).
    \return The concentration at (\a x, \a y, \a z).
  */
  template<class T>
  T GaussianPuff<T>::ComputePuffConcentration(Puff<T>* puff, int s,
                                              T x, T y, T z)
  {

    // Puff quantity.
    T quantity = puff->GetQuantity(s);

    // Puff position.
    T z_c = puff->GetZ();
    T z_above = puff->GetHeightAboveBL();
    T penetration_factor = puff->GetPenetrationFactor();
    T overcamp_factor = puff->GetReflectionFactor(s);
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
    T dx = (x - x_c) * cos_angle_ + (y - y_c) * sin_angle_;
    // Crosswind distance from puff center.
    T dy = (x_c - x) * sin_angle_ + (y - y_c) * cos_angle_;

    // Downwind contribution.
    T Fx = exp(-dx * dx  / (2. * sigma_x * sigma_x));

    // Crosswind contribution.
    T Fy = exp(-dy * dy  / (2. * sigma_y * sigma_y));

    // Ground reflection: plume is below the inversion layer.
    T FzG = 0.;
    T sigma_z_2 = sigma_z * sigma_z;
    if (inversion_height_ == 0. || z_c != 0.)
      if (z_c - sigma_z <= 0.)
        {
          T dzg = z + z_c;
          FzG = overcamp_factor * exp(-dzg * dzg / (2. * sigma_z_2));
        }
    T Fz = 0.;
    T concentration = 0.;

    // If there is an inversion layer.
    if (inversion_height_ > 0.)
      {
        T FzI = 0.;
        // Part of the puff that is below inversion layer.
        if (z_c != 0.)
          {
            // Vertical contribution.
            T dzr = z - z_c;
            Fz = exp(-dzr * dzr / (2. * sigma_z_2));

            // Reflection occurs only when puff touches the inversion layer.
            if (z_c + sigma_z >= inversion_height_)
              {
                T dzi = z + z_c - 2. * inversion_height_;
                T dzig = z - z_c + 2. * inversion_height_;
                T dzgi = z - z_c - 2. * inversion_height_;

                FzI += exp(-dzi * dzi / (2. * sigma_z_2))
                  + exp(-dzgi * dzgi / (2. * sigma_z_2))
                  + exp(-dzig * dzig / (2. * sigma_z_2));
              }

            // Plume only impacts below the inversion layer.
            if (z <= inversion_height_)
              {
                // Far field.
                if (sigma_z > 1.5 * inversion_height_)
                  concentration +=  quantity * Fx * Fy
                    / (2 * pi * sigma_x * sigma_y * inversion_height_);
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

            if (option_gillani)
              {
                initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
                sigma_z = sqrt(initial_sigma_z_2 * (1. + 2.3 * sqrt(time)));
                sigma_z_2 = sigma_z * sigma_z;
              }
            // Reflection on the inversion layer.
            T dzi = z + z_above - 2. * inversion_height_;
            FzI += exp(-dzi * dzi / (2. * sigma_z_2));
            // Plume only impacts above the inversion layer.
            if (z >= inversion_height_)
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

        if (option_gillani && z_c >= boundary_height_)
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


  //! Computes the integral of a puff over a given volume.
  template<class T>
  T GaussianPuff<T>::ComputePuffIntegral(int index, int s,
                                         T x, T y, T z,
                                         T lx, T ly, T lz)
  {
    SetCurrentPuff(index);
    SetCurrentMeteo(*current_puff);
    if ((*current_puff)->GetNs() == this->Ns
        || (*current_puff)->GetSpeciesIndex() == s)
      return ComputePuffIntegral(*current_puff, s, x, y, z, lx, ly, lz);
    else
      return 0.;
  }


  //! Computes the integral of a puff over a given volume.
  template<class T>
  T GaussianPuff<T>::ComputePuffIntegral(Puff<T>* puff, int s,
                                         T x, T y, T z,
                                         T lx, T ly, T lz)
  {
    // sqrt(2.)
    const T sqrt_2 = 1.41421356237;

    // Puff quantity.
    T quantity = puff->GetQuantity(s);
    T overcamp_factor = puff->GetReflectionFactor(s);
    T concentration = 0.;

    // Puff center coordinates.
    T z_c = puff->GetZ();
    T y_ = puff->GetY();
    T x_ = puff->GetX();

    // Conversion.
    T x_c = x_ * cos_angle_ + y_ * sin_angle_;
    T y_c = y_ * cos_angle_ - x_ * sin_angle_;

    // Standard deviations.
    T sigma_x = puff->GetSigma_x();
    T sigma_y = puff->GetSigma_y();
    T sigma_z = puff->GetSigma_z();
    T initial_sigma_y_2 = puff->GetInitialSigma_y();
    T initial_sigma_z_2 = puff->GetInitialSigma_z();
    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
    sigma_z = sqrt(sigma_z * sigma_z + initial_sigma_z_2);

    // Integration bounds.
    T x1, x2, y1, y2, z1, z2;
    x1 = x * cos_angle_ + y * sin_angle_ - lx / 2.;
    x2 = x * cos_angle_ + y * sin_angle_ + lx / 2.;

    y1 = y * cos_angle_ - x * sin_angle_ - ly / 2.;
    y2 = y * cos_angle_ - x * sin_angle_ + ly / 2.;

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


  //!  Computes the horizontal and vertical standard deviations.
  template<class T>
  void GaussianPuff<T>::ComputeSigma(T distance, T transfer_time,
                                     T z, T& sigma_x, T& sigma_y,
                                     T& sigma_z, bool diff)
  {
    // Horizontal diffusion parameter.
    if (option_briggs)
      {
        if (rural_)
          sigma_y =
            ComputeRuralPlumeHorizontalSigma(distance, stability_);
        else
          sigma_y =
            ComputeUrbanPlumeHorizontalSigma(distance, stability_);
        sigma_x = sigma_y;
      }
    else if (option_doury)
      {
        if (isday_)
          sigma_y =
            ComputeDouryNormalPlumeHorizontalSigma(transfer_time);
        else if (wind_ > 3.)
          sigma_y =
            ComputeDouryNormalPlumeHorizontalSigma(transfer_time);
        else
          sigma_y =
            ComputeDouryLowPlumeHorizontalSigma(transfer_time);
        sigma_x = sigma_y;
      }
    else
      {
        sigma_x = ComputeDownwindSigma(transfer_time, z,
                                       friction_velocity_,
                                       convective_velocity_,
                                       boundary_height_,
                                       lmo_, coriolis_,
                                       option_hpdm, stability_);
        sigma_y = ComputeCrosswindSigma(transfer_time, z,
                                        friction_velocity_,
                                        convective_velocity_,
                                        boundary_height_,
                                        lmo_, coriolis_,
                                        option_hpdm, stability_);
      }
    // Vertical diffusion parameter.
    if (option_briggs)
      if (rural_)
        if (!diff)
          sigma_z =
            ComputeRuralPlumeVerticalSigma(distance, stability_);
        else
          sigma_z =
            DifferentiateRuralPlumeVerticalSigma(distance, stability_);
      else if (!diff)
        sigma_z =
          ComputeUrbanPlumeVerticalSigma(distance, stability_);
      else
        sigma_z =
          DifferentiateUrbanPlumeVerticalSigma(distance, stability_);
    else if (option_doury)
      {
        if (isday_)
          if (!diff)
            sigma_z =
              ComputeDouryNormalPlumeVerticalSigma(transfer_time);
          else
            sigma_z =
              DifferentiateDouryNormalPlumeVerticalSigma(transfer_time);
        else if (wind_ > 3.)
          if (!diff)
            sigma_z =
              ComputeDouryNormalPlumeVerticalSigma(transfer_time);
          else
            sigma_z =
              DifferentiateDouryNormalPlumeVerticalSigma(transfer_time);
        else if (!diff)
          sigma_z =
            ComputeDouryLowPlumeVerticalSigma(transfer_time);
        else
          sigma_z =
            DifferentiateDouryLowPlumeVerticalSigma(transfer_time);
      }
    else if (!diff)
      sigma_z = ComputeVerticalSigma(transfer_time, z,
                                     friction_velocity_,
                                     convective_velocity_,
                                     boundary_height_, lmo_, coriolis_,
                                     option_hpdm, temperature_, stability_);
    else
      sigma_z = DifferentiateVerticalSigma(transfer_time, z, wind_,
                                           friction_velocity_,
                                           convective_velocity_,
                                           boundary_height_,
                                           lmo_, coriolis_, option_hpdm,
                                           temperature_, stability_);
  }


  //!  Computes the time/distance corresponding to given standard deviations.
  template<class T>
  void GaussianPuff<T>::ComputeTime(T z, T sigma_y,
                                    T& distance, T& transfer_time)
  {
    // Horizontal diffusion parameter.
    if (option_briggs)
      {
        if (rural_)
          distance =
            ComputeRuralPlumeHorizontalDistance(sigma_y, stability_);
        else
          distance =
            ComputeUrbanPlumeHorizontalDistance(sigma_y, stability_);
        transfer_time = distance / wind_;
      }
    else if (option_doury)
      {
        transfer_time = ComputeDouryNormalPlumeHorizontalTime(sigma_y);
        distance = transfer_time * wind_;
      }
    else
      {
        transfer_time = ComputeHorizontalTime(sigma_y, z,
                                              friction_velocity_,
                                              convective_velocity_,
                                              boundary_height_, lmo_,
                                              coriolis_,
                                              option_hpdm, stability_);
        distance = transfer_time * wind_;
      }
  }


  //! Computes the loss factor for a given puff, and all species.
  /*! Computes the loss factor for a given puff, and all species.
    \param index puff index.
  */
  template<class T>
  void GaussianPuff<T>::ComputeLossFactor(int index)
  {
    SetCurrentPuff(index);
    for (int s = 0; s < this->Ns; s++)
      if ((*current_puff)->GetNs() == this->Ns
          || (*current_puff)->GetSpeciesIndex() == s)
        ComputeLossFactor(*current_puff, s);
  }



  //! Computes the loss factor for a given puff, and all species.
  /*! Computes the loss factor for a given puff, and all species.
    \param puff the puff.
  */
  template<class T>
  void GaussianPuff<T>::ComputeLossFactor(Puff<T>* puff)
  {
    for (int s = 0; s < this->Ns; s++)
      if (puff->GetNs() == this->Ns
          || puff->GetSpeciesIndex() == s)
        ComputeLossFactor(puff, s);
  }


  //! Computes the loss factor.
  /*! Computes the loss factor.
    \param puff the puff.
    \param s species index.
  */  template<class T>
  void GaussianPuff<T>::ComputeLossFactor(Puff<T>* puff, int s)
  {
    T loss_factor = 1.;
    T overcamp_factor = 1.;
    T distance = puff->GetDistance();
    T transfer_time = puff->GetPuffTime();
    T quantity = puff->GetQuantity(s);
    T z = puff->GetZ();
    T rad, bio, scav, dep;
    if (option_radioactive_decay)
      rad = half_life_time(s);
    else
      rad = 0.;
    if (option_biological_decay)
      bio = biological_half_life_time(s);
    else
      bio = 0.;
    if (puff->HasMeteo())
      {
        if (option_scavenging)
          scav = puff->GetScavengingCoefficient(s);
        else
          scav = 0.;
        if (option_dry_deposition)
          dep = puff->GetDepositionVelocity(s);
        else
          dep = 0.;
      }
    else
      {
        if (option_scavenging)
          scav = scavenging_coefficient(s);
        else
          scav = 0.;
        if (option_dry_deposition)
          dep = deposition_velocity(s);
        else
          dep = 0.;
      }
    ComputeLossFactor(distance, transfer_time, z, rad, bio, scav, dep,
                      loss_factor, overcamp_factor);
    puff->SetQuantity(loss_factor * quantity, s);
    puff->SetReflectionFactor(overcamp_factor, s);
  }


  //!  Computes the loss factor.
  template<class T>
  void GaussianPuff<T>::ComputeLossFactor(T distance, T transfer_time,
                                          T z, T rad, T bio, T scav, T dep,
                                          T& loss_factor, T& overcamp_factor)
  {
    T radioactive_factor = 1.;
    T biological_factor = 1.;
    T scavenging_factor = 1.;
    T chamberlain_factor = 1.;
    overcamp_factor = 1.;

    const T pi = 3.14159265358979323846264;
    // log(2.)
    const T log2 = 0.69314718055994529;

    // Radioactive decay.
    if (option_radioactive_decay && rad != 0.)
      radioactive_factor = exp(-log2 * this->Delta_t / rad);
    // Biological decay.
    if (option_biological_decay && bio != 0.)
      biological_factor = exp(-log2 * this->Delta_t / bio);

    // Scavenging.
    if (option_scavenging)
      scavenging_factor = exp(- scav * this->Delta_t);

    // Dry deposition.
    if (option_dry_deposition)
      {
        if (deposition_model == "Chamberlain")
          {
            T z2 = 0.5 * z * z;
            T dsigma_x, dsigma_y, dsigma_z;
            T fdry = 0;
            ComputeSigma(distance, transfer_time, z,
                         dsigma_x, dsigma_y, dsigma_z);
            fdry = exp(-z2 / (dsigma_z * dsigma_z)) / dsigma_z;
            chamberlain_factor = exp(-dep * sqrt(2. / pi)
                                     * fdry * this->Delta_t);
          }
        else
          {
            T dsigma_z, sigma_x, sigma_y, sigma_z;
            bool diffz = 1;
            ComputeSigma(distance, transfer_time, z,
                         sigma_x, sigma_y, sigma_z);
            ComputeSigma(distance, transfer_time,
                         z, sigma_x, sigma_y, dsigma_z, diffz);
            overcamp_factor = 1. - 2. * dep
              / (dep + wind_ * z * dsigma_z / sigma_z);
          }
      }
    loss_factor = radioactive_factor * scavenging_factor *
      biological_factor * chamberlain_factor;
  }


  //!  Computes the interaction coefficient between two puffs.
  template<class T>
  T GaussianPuff<T>::ComputeInteractionCoefficient(Puff<T>* puff_alpha,
                                                   Puff<T>* puff_beta)
  {
    // 1. / (2. * sqrt(pi))
    const T coeff_pi = 0.28209479177387814;

    // Puff alpha.
    T sigmax_alpha, sigmay_alpha, sigmaz_alpha, x_alpha, y_alpha, z_alpha;
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

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_CXX
#endif
