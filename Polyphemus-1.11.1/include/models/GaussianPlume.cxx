// Copyright (C) 2006-2012, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Irène Korsakissok,
// Régis Briant
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

// This file is part of a Gaussian plume model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_CXX


#include "GaussianPlume.hxx"
#include "BriggsFormula.hxx"


namespace Polyphemus
{

  template<class T>
  const T GaussianPlume<T>::pi = 3.14159265358979323846264;

  template<class T>
  const T GaussianPlume<T>::log_2 = log(2.);

  template<class T>
  const T GaussianPlume<T>::cos_89 = cos(89. * pi / 180.);

  template<class T>
  const T GaussianPlume<T>::sin_89 = sin(89. * pi / 180.);


  /////////////////
  // CONSTRUCTOR //
  /////////////////


  //! Main constructor.
  /*!
    \param config_file Configuration file.
  */
  template<class T>
  GaussianPlume<T>::GaussianPlume(string config_file):
    BaseModel<T>(config_file), Nsource(0), cos_angle_(0.), sin_angle_(0.),
    option_cartesian(true), option_infinite_plume(true),
    computed_source_(-1), current_plume(0), previous_plume(0)
  {
  }


  //! Destructor.
  template<class T>
  GaussianPlume<T>::~GaussianPlume()
  {
    DeleteSource();
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Gets a Boolean option.
  /*!
    \param option the option name.
    \return The option value.
  */
  template<class T>
  bool GaussianPlume<T>::IsOption(string option)
  {
    if (option == "gillani")
      return option_gillani;
    if (option == "hpdm")
      return option_hpdm;
    if (option == "rural")
      return option_rural;
    if (option == "day")
      return option_day;
    if (option == "scavenging")
      return option_scavenging;
    if (option == "dry_deposition")
      return option_dry_deposition;
    if (option == "plume_rise")
      return option_plume_rise;
    if (option == "plume_rise_breakup")
      return option_breakup;
    if (option == "infinite_plume")
      return option_infinite_plume;
    throw string("Unknown Boolean option: \"") + option + "\".";
  }


  //! Sets a Boolean option.
  /*!
    \param option the option name.
    \param value the option value.
  */
  template<class T>
  void GaussianPlume<T>::SetOption(string option, bool value)
  {
    if (option == "gillani")
      option_gillani = value;
    else if (option == "hpdm")
      option_hpdm = value;
    else if (option == "rural")
      option_rural = value;
    else if (option == "day")
      option_day = value;
    else if (option == "scavenging")
      option_scavenging = value;
    else if (option == "dry_deposition")
      option_dry_deposition = value;
    else if (option == "plume_rise")
      option_plume_rise = value;
    else if (option == "plume_rise_breakup")
      option_breakup = value;
    else if (option == "cartesian")
      option_cartesian = value;
    else if (option == "infinite_plume")
      option_infinite_plume = value;
    else
      throw string("Unknown Boolean option: \"") + option + "\".";
  }


  //! Gets a string option.
  /*!
    \param option the option name.
    \return The option value.
  */
  template<class T>
  string GaussianPlume<T>::GetStringOption(string option)
  {
    if (option == "parameterization_std")
      if (option_briggs)
        return "Briggs";
      else if (option_doury)
        return "Doury";
      else
        return "similarity_theory";
    if (option == "deposition_model")
      return deposition_model;
    if (option == "stability_class")
      return stability_class_;
    throw string("Unknown string option: \"") + option + "\".";
  }


  //! Sets a string option.
  /*!
    \param option the option name.
    \param value on exit, the option value.
  */
  template<class T>
  void GaussianPlume<T>::SetStringOption(string option, string value)
  {
    if (option == "parameterization_std")
      if (value == "Briggs")
        {
          option_briggs = true;
          option_doury = false;
        }
      else if (value == "Doury")
        {
          option_briggs = false;
          option_doury = true;
        }
      else if (value == "similarity_theory")
        {
          option_briggs = false;
          option_doury = false;
        }
      else
        throw string("Unknown parameterization for standard deviation: \"")
          + value + "\".";
    else if (option == "deposition_model")
      {
        deposition_model = value;
        if (deposition_model != "Chamberlain"
            && deposition_model != "Overcamp")
          throw string("The deposition model must be \"Chamberlain\"")
            + " or \"Overcamp\".";
      }
    else if (option == "stability_class")
      {
        stability_class_ = value;
        if (stability_class_ == "A")
          stability_ = 0;
        else if (stability_class_ == "B")
          stability_ = 1;
        else if (stability_class_ == "C")
          stability_ = 2;
        else if (stability_class_ == "D")
          stability_ = 3;
        else if (stability_class_ == "E")
          stability_ = 4;
        else if (stability_class_ == "F")
          stability_ = 5;
        else
          throw string("Unknown stability class: \"")
            + stability_class_ + "\".";
      }
    else
      throw string("Unknown string option: \"") + option + "\".";
  }


  //! Gets the value of model parameter \a name.
  /*!
    \param name name of the model parameter.
    \return The value of the model parameter.
  */
  template<class T>
  T GaussianPlume<T>::GetModelParameter(string name)
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
    if (name == "stability")
      return stability_;
    if (name == "rural")
      return option_rural;
    throw string("Unknown numerical value: \"") + name + "\".";
  }


  //! Sets model parameter \a name to \a value.
  /*!
    \param name name of the model parameter.
    \param value the new value for the model parameter.
  */
  template<class T>
  void GaussianPlume<T>::SetModelParameter(string name, T value)
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
      throw string("Unknown numerical value: \"") + name + "\".";
  }


  //! Gets source data.
  /*!
    \param name the name of the data.
    \param value on exit, the values for all sources.
  */
  template<class T>
  void GaussianPlume<T>::GetSourceData(string name, Array<T, 1>& value)
  {
    value.resize(Nsource);

    source_iterator iter;
    int i = 0;
    if (name == "rate")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetRate();
    else if (name == "velocity")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetVelocity();
    else if (name == "temperature")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetTemperature();
    else if (name == "z")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetZ();
    else if (name == "y")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetY();
    else if (name == "x")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetX();
    else if (name == "diameter")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        value(i) = (*iter)->GetDiameter();
    else
      throw string("Unknown source data: \"") + name + "\".";
  }


  //! Sets source data.
  /*!
    \param name the name of the data to be set.
    \param value The values.
  */
  template<class T>
  void GaussianPlume<T>::SetSourceData(string name, const Array<T, 1>& value)
  {
    source_iterator iter;
    int i = 0;
    if (name == "rate")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetRate(value(i));
    else if (name == "velocity")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetVelocity(value(i));
    else if (name == "temperature")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetTemperature(value(i));
    else if (name == "z")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetZ(value(i));
    else if (name == "y")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetY(value(i));
    else if (name == "x")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetX(value(i));
    else if (name == "diameter")
      for (iter = SourceList.begin(); iter != SourceList.end(); iter++, i++)
        (*iter)->SetDiameter(value(i));
    else
      throw string("Unknown source data: \"") + name + "\".";
  }


  //! Gets species data.
  /*!
    \param name the name of the data.
    \param value on exit, the values.
  */
  template<class T>
  void GaussianPlume<T>::GetSpeciesData(string name, Array<T, 1>& value)
  {
    if (name == "half_life_time")
      value = half_life_time;
    else if (name == "biological_half_life_time")
      value = biological_half_life_time;
    else if (name == "scavenging_coefficient")
      value = scavenging_coefficient;
    else if (name == "deposition_velocity")
      value = deposition_velocity;
    else
      throw string("Unknown species data: \"") + name + "\".";
  }


  //! Sets species data.
  /*!
    \param name the name of the data to be set.
    \param value The values.
  */
  template<class T>
  void GaussianPlume<T>::SetSpeciesData(string name, const Array<T, 1>& value)
  {
    if (name == "half_life_time")
      half_life_time = value;
    else if (name == "biological_half_life_time")
      biological_half_life_time = value;
    else if (name == "scavenging_coefficient")
      scavenging_coefficient = value;
    else if (name == "deposition_velocity")
      deposition_velocity = value;
    else
      throw string("Unknown species data: \"") + name + "\".";
  }


  //! Returns whether there is deposition used or not.
  /*! \return true if option_dry_deposition is set to 'yes', false otherwise.
   */
  template<class T>
  bool GaussianPlume<T>::WithDeposition()
  {
    return option_dry_deposition;
  }


  //! Returns whether there is scavenging used or not.
  /*! \return true if option_scavenging is set to 'yes', false otherwise.
   */
  template<class T>
  bool GaussianPlume<T>::WithScavenging()
  {
    return option_scavenging;
  }


  //! Returns whether there is similarity theory used or not.
  /*! \return true if sigma parameterization is 'similarity_theory', false
    otherwise.
  */
  template<class T>
  bool GaussianPlume<T>::WithSimilarity()
  {
    return !option_briggs && !option_doury;
  }


  //! Returns whether the NO2 chemistry is activated or not.
  /*! \return true if the NO2 chemistry is activated, false otherwise.
   */
  template<class T>
  bool GaussianPlume<T>::WithChemistry()
  {
    return option_NO2_chemistry || option_OH_chemistry;
  }


  //! Get concentration of a given species at a given point.
  /*!
    \param species species index.
    \param x abscissa (m).
    \param y ordinate (m).
    \param z height (m).
    \return The value of concentration at the given point (mass/m^3).
  */
  template<class T>
  T GaussianPlume<T>::GetConcentration(int species, T z, T y, T x)
  {
    int Ncoordinate = coordinate_list_.extent(0);
    for (int i = 0; i < Ncoordinate; ++i)
      if (coordinate_list_(i, 0) == z
          && coordinate_list_(i, 1) == y
          && coordinate_list_(i, 2) == x)
        return Concentration_list_(i, species);

    throw "Cannot get the concentration at the point (" + to_str(z) + ", "
      + to_str(y) + ", " + to_str(x) + ") because it is not in the list "
      "of points to be computed.";
  }


  //! Return concentration of the whole domain.
  /*!
    \return The concentration matrix for the whole domain (mass/m^3).
  */
  template<class T>
  Data<T, 4>&  GaussianPlume<T>::GetConcentration()
  {
    return this->Concentration;
  }


  //! Sets the current plume source.
  /*! \param index the plume source index to be set as current.
   */
  template<class T>
  void GaussianPlume<T>::SetCurrentPlume(int index)
  {
    CheckSourceIndex(index);
    current_plume = SourceList[index];
  }


  //! Sets meteorological parameters for a given plume source.
  /*!
    \param index plume source index.
    \param temperature the temperature of ambient air (Kelvin degrees).
    \param wind_angle the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind the wind velocity (m/s).
    \param stability_class stability class in [A, F].
    \param longitude the longitude coordinate.
    \param latitude the latitude coordinate.
    \param isday true if it is daytime.
    \param rural true in rural environment.
    \param boundary_height boundary layer height (meters).
  */
  template<class T>
  void GaussianPlume<T>::SetGaussianMeteo(int index,
                                          T temperature, T wind_angle, T wind,
                                          string stability_class,
                                          T longitude, T latitude,
                                          bool isday, bool rural,
                                          T boundary_height)
  {
    temperature_ = temperature;
    wind_angle_ = wind_angle;
    wind_ = wind;
    boundary_height_ = boundary_height;
    stability_class_ = upper_case(stability_class);
    if (stability_class_ == "A")
      stability_ = 0;
    else if (stability_class_ == "B")
      stability_ = 1;
    else if (stability_class_ == "C")
      stability_ = 2;
    else if (stability_class_ == "D")
      stability_ = 3;
    else if (stability_class_ == "E")
      stability_ = 4;
    else if (stability_class_ == "F")
      stability_ = 5;
    else
      throw string("Stability class should be in [A, F], but \"")
        + stability_class_ + "\" was provided.";
    option_day = isday;
    option_rural = rural;
  }

  //! Sets the meteorological data for a given plume source.
  /*!
    \param index plume source index
    \param temperature the temperature of ambient air (Kelvin degrees).
    \param wind_angle the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind the wind velocity (m/s).
    \param friction_velocity friction velocity (m/s).
    \param convective_velocity convective velocity (m/s).
    \param boundary_height boundary layer height (m).
    \param lmo Monin Obukhov length (m).
    \param coriolis Coriolis parameter.
    \param isday boolean equal to true if it is daytime, false otherwise.
    \param rural true in rural environment.
  */
  template<class T>
  void GaussianPlume<T>
  ::SetGaussianMeteo(int index, T temperature, T wind_angle, T wind,
                     T friction_velocity,
                     T convective_velocity, T boundary_height, T lmo,
                     T coriolis, T longitude, T latitude, bool isday,
                     bool rural)
  {
    temperature_ = temperature;
    wind_angle_ = wind_angle;
    friction_velocity_ = friction_velocity;
    convective_velocity_ = convective_velocity;
    coriolis_ = coriolis;
    option_day = isday;
    option_rural = rural;
    wind_ = wind;
    boundary_height_ = boundary_height;
    lmo_ = lmo;
  }


  //! Sets values of deposition velocities for a given plume.
  /*! It sets the deposition velocities for a given plume.
    \param[in] index plume source index.
    \param[in] deposition_velocity_in array containing the deposition velocities
    for all species (set to 0. if there is no deposition).
  */
  template<class T>
  void GaussianPlume<T>::SetGaussianDepositionVelocity(int index,
                                                       const Array<T, 1>& deposition_velocity_in)
  {
    SetCurrentPlume(index);
    deposition_velocity = deposition_velocity_in;
  }


  //! Sets values of scavenging coefficients for a given plume.
  /*! It sets the scavenging coefficients for a given plume.
    \param[in] index plume source index.
    \param[in] scavenging_coefficient_in array containing the scavenging
    coefficients for all species (set to 0. if there is no scavenging).
  */
  template<class T>
  void GaussianPlume<T>::
  SetGaussianScavengingCoefficient(int index,
                                   const Array<T, 1>& scavenging_coefficient_in)
  {
    SetCurrentPlume(index);
    scavenging_coefficient = scavenging_coefficient_in;
  }


  //! Gets the position of the plume center for a given plume source index.
  /*!
    \param index plume source index.
    \param x abscissa (meters).
    \param y ordinate (meters).
    \param z height (meters).
    \param distance distance (meters).
    \param time time since release (s).
  */  template<class T>
  void GaussianPlume<T>::GetGaussianPosition(int index, T& x, T& y, T& z,
                                             T& distance, T& time)
  {
    SetCurrentPlume(index);
    if (current_plume->GetEmissionType() == "continuous")
      {
        x = current_plume->GetX();
        y = current_plume->GetY();
        z = current_plume->GetZ();
      }
    else if (current_plume->GetEmissionType() == "continuous_line")
      {
        x = 0.5 * (current_plume->GetX() + current_plume->GetX2());
        y = 0.5 * (current_plume->GetY() + current_plume->GetY2());
        z = 0.5 * (current_plume->GetZ() + current_plume->GetZ2());
      }

    // Distance and time are not used.
    distance = 0.;
    time = 0.;
  }


  //! Returns the standard deviations for a given plume index.
  /*!
    \param index plume source index.
    \param sigma_x standard deviation in the downwind direction (m).
    \param sigma_y standard deviation in the crosswind direction (m).
    \param sigma_z vertical standard deviation (m).
  */
  template<class T>
  void GaussianPlume<T>::GetGaussianSigma(int index, T& sigma_x,
                                          T& sigma_y, T& sigma_z)
  {
    SetCurrentPlume(index);
    sigma_y = current_plume->GetSigma_y();
    sigma_z = current_plume->GetSigma_z();
  }


  //! Returns the emission type of a given source index.
  /*!
    \param index plume source index.
    \return the emission type of a given source index.
  */
  template<class T>
  string GaussianPlume<T>::GetEmissionType(int index)
  {
    SetCurrentPlume(index);
    return current_plume->GetEmissionType();
  }


  //! Returns the parameter named \a name of the source of index \a index.
  /*!
    \param name name of the source parameter.
    \param index index of the source plume.
    \return The source parameter value.
  */
  template<class T>
  T GaussianPlume<T>::GetSourceParameter(string name, int index)
  {
    SetCurrentPlume(index);
    if (name == "X")
      return current_plume->GetX();
    if (name == "Y")
      return current_plume->GetY();
    if (name == "Z")
      return current_plume->GetZ();
    if (name == "X2")
      return current_plume->GetX2();
    if (name == "Y2")
      return current_plume->GetY2();
    if (name == "Z2")
      return current_plume->GetZ2();
    if (name == "rate")
      return current_plume->GetRate();
    if (name == "width")
      return current_plume->GetWidth();
    if (name == "VehicleVelocity")
      return current_plume->GetVehicleVelocity();
    if (name == "Area")
      return current_plume->GetArea();
    if (name == "Density")
      return current_plume->GetDensity();
    throw string("Unknown source parameter: \"") + name + "\".";
  }


  //! Returns the species coordinate of a given source index.
  /*!
    \param index plume source index.
    \return the species coordinates of a given source index.
  */
  template<class T>
  int GaussianPlume<T>::GetSpecies(int index)
  {
    SetCurrentPlume(index);
    int species = current_plume->GetSpeciesIndex();
    return species;
  }


  //! Sets the simulation beginning date.
  /*!
    \param date_min the simulation beginning date.
  */
  template<class T>
  void GaussianPlume<T>::SetDateMin(Date date_min)
  {
    this->Date_min = date_min;
  }


  //! Sets the simulation time step.
  /*!
    \param delta_t time step in seconds.
  */
  template<class T>
  void GaussianPlume<T>::SetTimeStep(T delta_t)
  {
    this->Delta_t = delta_t;
  }


  //! Returns the number of plume source.
  /*!
    \return the number of plume source.
  */
  template<class T>
  int GaussianPlume<T>::GetInputSourceCount() const
  {
    return Nsource;
  }


  //! Gets the time step.
  /*!
    \param delta_t time step in seconds.
  */
  template<class T>
  T GaussianPlume<T>::GetTimeStep()
  {
    return this->Delta_t;
  }


  //! Returns the number of sources, after discretization.
  /*!
    \return the number of sources, after discretization.
  */
  template<class T>
  int GaussianPlume<T>::GetGaussianSourceCount() const
  {
    return SourceList.size();
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
  void GaussianPlume<T>::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

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
    this->config.PeekValue("File_source", file_source);
    this->config.PeekValue("File_correction", correction_coefficient_path_);
    this->config.PeekValue("With_NO2_chemistry", option_NO2_chemistry);
    this->config.PeekValue("With_OH_chemistry", option_OH_chemistry);
    this->config.PeekValue("With_high_width_accuracy",
                           option_high_width_precision);

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
    if (sigma_formula == "similarity_theory" || read_all)
      this->config.PeekValue("With_HPDM", option_hpdm);

    // Plume rise parameterization.
    if (option_plume_rise || read_all)
      {
        this->config.PeekValue("With_plume_rise_breakup", option_breakup);
        string plume_rise_formula;
        this->config.PeekValue("Plume_rise_parameterization",
                               "HPDM | Holland | Concawe",
                               plume_rise_formula);
        option_holland = plume_rise_formula == "Holland";
        option_concawe = plume_rise_formula == "Concawe";
      }

    this->config.PeekValue("Compute_domain", option_compute_domain);
    this->config.PeekValue("Compute_list", option_compute_list);
    this->config.PeekValue("Npmax", Np_max);
    this->config.PeekValue("Discretization_step", discretization_step);

    // Read coordinate list.
    if (option_compute_list)
      {
        string coordinates_computed_points_file;
        this->config.PeekValue("Coordinates_computed_points",
                               coordinates_computed_points_file);

        list<vector<string> > param_list;
        vector<string> param;
        ConfigStream coordinate_stream(coordinates_computed_points_file);
        while (!coordinate_stream.IsEmpty())
          param_list.push_back(split(coordinate_stream.GetLine()));

        int N = int(param_list.size());
        coordinate_list_.resize(N, 3);

        for (int i = 0; i < N; i++)
          {
            param = param_list.back();
            param_list.pop_back();
            coordinate_list_(i, 0) = to_num<T>(param[0]);

            if (option_cartesian)
              {
                coordinate_list_(i, 1) = to_num<T>(param[1]);
                coordinate_list_(i, 2) = to_num<T>(param[2]);
              }
            else
              LatLonToCartesian(to_num<T>(param[2]), to_num<T>(param[1]),
                                coordinate_list_(i, 2), coordinate_list_(i, 1));
          }
        Concentration_list_.Resize(N, this->Ns);
      }

    // Gets the deposition model.
    if (option_dry_deposition || read_all)
      {
        this->config.SetSection("[deposition]");
        this->config.PeekValue("Deposition_model", "Chamberlain | Overcamp",
                               deposition_model);

        if (deposition_model == "Chamberlain")
          this->config.PeekValue("Nchamberlain", "> 0", Nchamberlain);
      }
    else
      Nchamberlain = 0;

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
  /*! Allocates the grids and the concentration Data for gaseous species.
   */
  template<class T>
  void GaussianPlume<T>::Allocate()
  {
    BaseModel<T>::Allocate();
  }


  //! Sources initialization.
  /*! It sets all sources from a text file. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Rate: the rate (mass unit per second),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees),
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species_name: species name.
    </ul>
    \param source_file file that describes the sources.
  */
  template<class T>
  void GaussianPlume<T>::InitSource(string source_file)
  {
    // List of emissions.
    this->PointEmissionManager = new BasePointEmission<T>();
    BasePointEmission<T>& emission = *this->PointEmissionManager;

    emission.Init(source_file, this->species_list);

    T rate, velocity, temperature, diameter, x1, y1, z1;
    string species_name;

    int Nemis = emission.GetNumberEmission();
    for (int i = 0; i < Nemis; i++)
      {
        const vector<int>& emitted_species_index =
          emission.GetEmittedSpeciesIndex(i);
        int Ns_emitted = emitted_species_index.size();
        int s;

        // Loop on all emitted species.
        for (int species = 0; species < Ns_emitted; species++)
          {
            if (emission.GetEmissionType(i) == "continuous")
              {
                s = emitted_species_index[species];
                rate = emission.GetRate(i, species);
                emission.GetPlumeRiseParam(velocity, temperature,
                                           diameter, i);
                emission.GetEmissionCoordinates(x1, y1, z1, i);
                PlumeSource<T>* ContinuousSource =
                  new PlumeSource<T>(rate, velocity, temperature, diameter,
                                     x1, y1, z1, s, i);

                SourceList.push_back(ContinuousSource);
                Nsource++;
              }
            else if (emission.GetEmissionType(i) == "continuous_line")
              {
                s = emitted_species_index[species];
                emission.GetPlumeRiseParam(velocity, temperature,
                                           diameter, i);
                list<Array<T, 1> > coordinate_list;
                Array<T, 1> coordinate(6);
                emission.GetEmissionCoordinates(coordinate_list, i);
                int id_section = 0;
                for (typename list<Array<T, 1> >::const_iterator iter
                       = coordinate_list.begin();
                     iter != coordinate_list.end(); iter++)
                  {
                    rate = emission.GetRate(i, species, id_section);
                    T width = emission.GetWidth(i, id_section);
                    coordinate = *iter;

                    T VehicleVelocity = emission.GetVehicleVelocity(i, id_section);
                    T Area = emission.GetArea(i, id_section);
                    T Density = emission.GetDensity(i, id_section);

                    PlumeSource<T>* LineSource =
                      new PlumeLineSource<T>(rate, coordinate(0),
                                             coordinate(1), coordinate(2),
                                             coordinate(3), coordinate(4),
                                             coordinate(5), width,
                                             VehicleVelocity, Area, Density,
                                             s, i, id_section);
                    SourceList.push_back(LineSource);
                    Nsource++;
                    id_section++;
                  }
              }
            else
              throw string("Source type must be \"continuous\" or ")
                + "\"continuous_line\" but is "
                + emission.GetEmissionType(i);
          }
      }

    if (!SourceList.empty())
      current_plume = *SourceList.begin();
    Nsource_temp = Nsource;
  }


   //! Sources removal.
  template<class T>
  void GaussianPlume<T>::DeleteSource()
  {
    EmptyPointSourceList();

    for (int i = 0; i < (int) SourceList.size(); ++i)
        delete SourceList[i];
    SourceList.clear();
    Nsource = 0;
    current_plume = 0;

    delete this->PointEmissionManager;
  }


  //! Model initialization.
  /*! It reads the configuration.
    \param read_all_input Should all input data and options be read in the
    configuration files? If so, even if unnecessary data will be read,
    e.g. the Monin-Obukhov length while the Doury parameterization is
    used. This enables to switch et a later time to other parameterizations
    (e.g., to "similarity_theory").
  */
  template<class T>
  void GaussianPlume<T>::Init(bool read_all_input)
  {
    read_all = read_all_input;

    BaseModel<T>::Init();

    /*** Initializes the sources ***/

    InitSource(file_source);
  }


  //! Initializes meteorological conditions.
  /*! It sets the meteorological data from a configuration file.  The
    situation is described in a dedicated section "[situation]" in which one
    finds the following entries:
    <ul>
    <li> Temperature: the temperature of ambient air (Celsius degrees).
    <li> Wind_angle: the angle between positive x-axis and wind,
    counterclockwise (degrees).
    <li> Wind: the wind velocity (m/s).
    <li> Inversion_height: the inversion height (m).
    <li> Stability: stability class in [A, F].
    </ul>
    \param meteo ConfigStream instance through which all entries (Temperature,
    Wind_angle, ...) may be read to set the meteorological situation.
    \param show_meteo indicates whether the meteorological data is to be
    displayed on screen.
  */
  template<class T>
  void GaussianPlume<T>::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {
    meteo.PeekValue("Temperature", temperature_);
    meteo.PeekValue("Wind_angle", wind_angle_);
    meteo.PeekValue("Wind", wind_);
    meteo.PeekValue("Boundary_height", "positive", boundary_height_);
    if (option_briggs || option_doury || read_all)
      {
        if (option_plume_rise && !read_all)
          if (option_breakup && !option_holland)
            {
              meteo.PeekValue("Friction_velocity", friction_velocity_);
              meteo.PeekValue("Convective_velocity", convective_velocity_);
            }
        meteo.PeekValue("Stability", "A | B | C | D | E | F",
                        stability_class_);
        stability_class_ = upper_case(stability_class_);
        if (stability_class_ == "A")
          stability_ = 0;
        else if (stability_class_ == "B")
          stability_ = 1;
        else if (stability_class_ == "C")
          stability_ = 2;
        else if (stability_class_ == "D")
          stability_ = 3;
        else if (stability_class_ == "E")
          stability_ = 4;
        else if (stability_class_ == "F")
          stability_ = 5;
      }
    if (!option_briggs && !option_doury || read_all)
      {
        meteo.PeekValue("Friction_velocity", friction_velocity_);
        meteo.PeekValue("Convective_velocity", convective_velocity_);
        meteo.PeekValue("LMO", lmo_);
        meteo.PeekValue("Coriolis", coriolis_);
      }

    if (show_meteo)
      cout << "\t" << temperature_
           << "\t\t" << wind_angle_
           << "\t\t" << wind_
           << "\t\t" << stability_class_ << endl;

    temperature_ = temperature_ + 273.15;

    /*** Scavenging coefficients ***/

    if (option_scavenging)
      InitScavenging(meteo);

    /*** Deposition velocities ***/

    if (option_dry_deposition)
      InitDeposition(meteo);

    if (option_NO2_chemistry)
      {
        // *** Hard coded value set to Paris, France.
        // *** It should be editable in configuration files...
        T latitude = 48.8566;
        T longitude = 2.35194;

        if (option_compute_domain)
          {
            Array<T, 1> K(2);
            k2_domain.resize(this->Ny, this->Nx);
            T zenith_angle = ComputeZenithAngle(this->GetCurrentDate()
                                                .GetNumberOfSeconds(),
                                                latitude, longitude);

            // Computes the kinetic constant of the reaction of NO with O3
            // in ppb^-1.seconds^-1.
            k1_ = exp(-27.29454887930734 - 1310. / temperature_);
            k1_ = k1_ * 2.46 * pow(10, 10);

            k2_domain = ComputeNO2PhotolysisRate(zenith_angle);
          }

        if (option_compute_list)
          {
            Array<T, 1> K(2);
            k2_list.resize(this->Concentration_list_.GetArray().extent(0));
            T zenith_angle = ComputeZenithAngle(this->GetCurrentDate().
                                                GetNumberOfSeconds(),
                                                latitude, longitude);

            // Computes the kinetic constant of the reaction of NO with O3
            // in ppb^-1.seconds^-1.
            k1_ = exp(-27.29454887930734 - 1310. / temperature_);
            k1_ = k1_ * 2.46 * pow(10, 10);

            k2_list = ComputeNO2PhotolysisRate(zenith_angle);
          }
      }

    if (option_OH_chemistry)
      {
        meteo.PeekValue("Relative_humidity", RelativeHumidity_);
        meteo.PeekValue("Pressure", Pressure_);
      }
  }


  /*! Initializes scavenging coefficients.
    \param meteo ConfigStream instance through which scavenging coefficients
    may be read.
  */
  template<class T>
  void GaussianPlume<T>::InitScavenging(ConfigStream& meteo)
  {

    scavenging_coefficient.resize(this->Ns);
    meteo.Find("Scavenging_coefficient");
    vector<string> scav_coefficient = split(meteo.GetLine());
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
  void GaussianPlume<T>::InitDeposition(ConfigStream& meteo)
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


  //////////////////
  // COMPUTATIONS //
  //////////////////


  //! Computations performed before any concentration can be computed.
  /*!  After this method is called, one may call 'GetConcentration' or
    'Compute'. But never before!
  */
  template<class T>
  void GaussianPlume<T>::InitCompute()
  {
    T wind_angle_rad = wind_angle_ / 180. * pi;
    cos_angle_ = cos(wind_angle_rad);
    sin_angle_ = sin(wind_angle_rad);

    // If some sources were added during the previous meteo situation.
    if (Nsource_temp != Nsource)
      {
        this->EraseDiscretizedSource();
        this->RestoreRate();
      }

    if (!option_briggs && !option_doury)
      // Matching between LMO and stability classes.
      ComputeStabilityClass(lmo_, stability_, stability_class_);
    if (stability_ != 4 && stability_ != 5)
      inversion_height_ = boundary_height_;
    else
      inversion_height_ = 0.;

    this->Discretization();

    if (option_plume_rise)
      ComputePlumeRise();

    InitCorrectionCoefficients();
  }


  //! Computes plume rise.
  template<class T>
  void GaussianPlume<T>::ComputePlumeRise()
  {
    // Loop over all sources.
    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); iter++)
      if ((*iter)->GetEmissionType() == "continuous")
        ComputeSourcePlumeRise(**iter);
  }


  //! Computes plume rise for a given source.
  template<class T>
  void GaussianPlume<T>::ComputeSourcePlumeRise(PlumeSource<T>& source)
  {
    // Computing the effective source height.
    T source_height = source.GetZ();
    T plume_rise;
    if (option_holland)
      plume_rise = ComputeHollandPlumeRise(temperature_, wind_,
                                           source.GetVelocity(),
                                           source.GetTemperature(),
                                           source.GetDiameter());
    else if (option_concawe)
      plume_rise = ComputeConcawePlumeRise(temperature_, wind_,
                                           inversion_height_,
                                           convective_velocity_,
                                           friction_velocity_,
                                           stability_,
                                           source.GetVelocity(),
                                           source.GetTemperature(),
                                           source.GetDiameter(),
                                           source_height, option_breakup);
    else
      plume_rise = ComputeHPDMPlumeRise(temperature_, wind_,
                                        inversion_height_,
                                        convective_velocity_,
                                        friction_velocity_,
                                        stability_,
                                        source.GetVelocity(),
                                        source.GetTemperature(),
                                        source.GetDiameter(),
                                        source_height, option_breakup);

    // Computing the additional plume spread due to plume rise.
    T initial_sigma_y_2 = source.GetSigma_y_2();
    T initial_sigma_z_2 = source.GetSigma_z_2();
    T sigma_y_pr = plume_rise / 3.5;
    T sigma_z_pr;
    if (option_gillani)
      {
        T dTdz;
        ComputePotentialTemperatureGradient(stability_, dTdz);
        dTdz = dTdz - 0.986;
        sigma_z_pr = 15. * exp(- 117. * dTdz);
      }
    else
      sigma_z_pr = plume_rise / 2.;
    source.SetSigma_y_2(initial_sigma_y_2 + sigma_y_pr * sigma_y_pr);
    source.SetSigma_z_2(initial_sigma_z_2 + sigma_z_pr * sigma_z_pr);

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
    source.SetPenetrationFactor(penetration_factor);
    source.SetHeight(effective_height);
    source.SetHeightAboveBL(effective_height_above);
  }


  template<class T>
  void GaussianPlume<T>::ComputeLossFactor(T distance, T transfer_time,
                                           T z, T rad, T bio, T scav, T dep,
                                           T& loss_factor, T& overcamp_factor)
  {
    T radioactive_factor = 1.;
    T biological_factor = 1.;
    T scavenging_factor = 1.;
    T chamberlain_factor = 1.;
    overcamp_factor = 1.;

    // Radioactive decay.
    if (option_radioactive_decay && rad != 0.)
      radioactive_factor = exp(-log_2 * transfer_time / rad);

    // Biological decay.
    if (option_biological_decay && bio != 0.)
      biological_factor = exp(-log_2 * transfer_time / bio);

    // Scavenging.
    if (option_scavenging)
      scavenging_factor = exp(- scav * transfer_time);

    // Dry deposition.
    if (option_dry_deposition)
      if (deposition_model == "Chamberlain")
        {
          T delta_x = distance / T(Nchamberlain);
          T z2 = 0.5 * z * z;
          T dsigma_y, dsigma_z;
          T x(0.);
          T exp_factor = - dep / wind_ * sqrt(2. / pi);

          T fdry = 0;
          for (int i = 1; i < Nchamberlain + 1; i++)
            {
              x += delta_x;
              T time = x / wind_;
              ComputeSigma(x, time, z, dsigma_y, dsigma_z);
              fdry += exp(-z2 / (dsigma_z * dsigma_z)) / dsigma_z;
            }
          chamberlain_factor = exp(exp_factor * fdry * delta_x);
        }
      else
        {
          T dsigma_z, sigma_y, sigma_z;
          bool diffz = 1;
          ComputeSigma(distance, transfer_time, z, sigma_y, sigma_z);
          ComputeSigma(distance, transfer_time, z, sigma_y, dsigma_z, diffz);
          overcamp_factor = 1. - 2. * dep /
            (dep + wind_ * z * dsigma_z / sigma_z);
        }
    loss_factor = radioactive_factor * scavenging_factor *
      biological_factor * chamberlain_factor;
  }


  //!  Computes the loss factor.
  template<class T>
  void GaussianPlume<T>::ComputeLossFactor(T distance, T z, int species,
                                           T& loss_factor, T& overcamp_factor)
  {
    T rad, bio, scav, dep;
    T transfer_time = distance / wind_;
    if (option_radioactive_decay)
      rad = half_life_time(species);
    else
      rad = 0.;
    if (option_biological_decay)
      bio = biological_half_life_time(species);
    else
      bio = 0.;
    if (option_scavenging)
      scav = scavenging_coefficient(species);
    else
      scav = 0.;
    if (option_dry_deposition)
      dep = deposition_velocity(species);
    else
      dep = 0.;

    ComputeLossFactor(distance, transfer_time, z, rad, bio, scav, dep,
                      loss_factor, overcamp_factor);
  }


  //! Computes concentration of a given species at a given point.
  /*!
    \param species species index.
    \param x abscissa (m).
    \param y ordinate (m).
    \param z height (m).
    \return The value of concentration at the given point (mass/m^3).
  */
  template<class T>
  T GaussianPlume<T>::ComputeGaussianConcentration(int species, T z, T y, T x)
  {
    // Distances downwind and crosswind from the source (m).
    T distance_x, distance_x1, distance_x2, distance_y, distance_y1,
      distance_y2, time;
    distance_x = 0;

    // Horizontal and vertical standard deviations (m).
    T sigma_y, sigma_z, initial_sigma_y_2, initial_sigma_z_2;

    // Effective height of release (m), inside and above BL.
    T effective_height, effective_height_above;

    // Fraction of the plume that penetrates above BL.
    T penetration_factor;

    // Rate (g/s).
    T rate = 0.;

    // Output concentration.
    T concentration = 0.;

    // Loss factors.
    T loss_factor, overcamp_factor;

    // Loop over all sources.
    int begin = 0, end = SourceList.size();
    if (computed_source_ != -1)
      {
        begin = computed_source_;
        end = begin + 1;
      }

    for (int i = begin; i < end; ++i)
      {
        PlumeSource<T>& source = *SourceList[i];
        if (source.GetSpeciesIndex() != species)
          continue;

        // Minimum volume.
        T minimum_volume = pi * source.GetDiameter()
          * source.GetDiameter() * source.GetVelocity();

        if (source.GetEmissionType() == "continuous")
          {
            // Downwind distance from source.
            distance_x = (x - source.GetX()) * cos_angle_
              + (y - source.GetY()) * sin_angle_;

            // Distance from downwind axis.
            distance_y = (source.GetX() - x) * sin_angle_
              + (y - source.GetY()) * cos_angle_;
            distance_y = abs(distance_y);

            // Effective height.
            effective_height = source.GetHeight();
            penetration_factor = source.GetPenetrationFactor();
            effective_height_above = source.GetHeightAboveBL();

            // If there is anything to compute.
            if (IsPlume(distance_x))
              {
                // Initial plume spread.
                initial_sigma_y_2 = source.GetSigma_y_2();
                initial_sigma_z_2 = source.GetSigma_z_2();

                if (effective_height != 0.)
                  {
                    // Computing standard deviations.
                    time = distance_x / wind_;
                    ComputeSigma(distance_x, time, effective_height,
                                 sigma_y, sigma_z);

                    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
                    if (option_gillani && effective_height
                        >= boundary_height_)
                      {
                        initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
                        sigma_z = sqrt(initial_sigma_z_2
                                       * (1. + 2.3 * sqrt(time)));
                      }
                    else
                      sigma_z = sqrt(sigma_z * sigma_z
                                     + initial_sigma_z_2);

                    // Computing loss factor.
                    ComputeLossFactor(distance_x, effective_height,
                                      species, loss_factor,
                                      overcamp_factor);

                    rate = source.GetRate() * (1. - penetration_factor);
                    concentration += loss_factor
                      * ComputePlumeConcentration(wind_, inversion_height_
                                                  , effective_height, rate
                                                  , distance_y, z,
                                                  sigma_y, sigma_z,
                                                  overcamp_factor,
                                                  minimum_volume);
                  }
                // Fraction of the plume above the boundary layer.
                if (effective_height_above != 0.)
                  {
                    // Computing standard deviations.
                    time = distance_x / wind_;
                    ComputeSigma(distance_x, time, boundary_height_,
                                 sigma_y, sigma_z);
                    sigma_y = sqrt(sigma_y * sigma_y + initial_sigma_y_2);
                    if (option_gillani)
                      {
                        initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
                        sigma_z = sqrt(initial_sigma_z_2
                                       * (1. + 2.3 * sqrt(time)));
                      }
                    else
                      sigma_z = sqrt(sigma_z * sigma_z +
                                     initial_sigma_z_2);

                    // Computing loss factor.
                    ComputeLossFactor(distance_x, effective_height_above,
                                      species, loss_factor,
                                      overcamp_factor);

                    rate = source.GetRate() * penetration_factor;
                    concentration += loss_factor
                      * ComputePlumeConcentration(wind_, inversion_height_
                                                  , effective_height_above
                                                  , rate, distance_y, z,
                                                  sigma_y, sigma_z,
                                                  overcamp_factor,
                                                  minimum_volume);
                  }
              }
          }
        else if (source.GetEmissionType() == "continuous_line")
          {
            T width = source.GetWidth();

            if (width == 0.)
              concentration +=
                LineSourceConcentration(species, z, y, x, source);
            else
              {
                T range = width * 0.5;
                concentration +=
                  RombergLineSourceConcentration(-range, range,
                                                 species, z, y, x, source);
              }
          }
      }
    return concentration;
  }

  //!  Computes the horizontal and vertical standard deviations.
  template<class T>
  void GaussianPlume<T>::ComputeSigma(T distance, T transfer_time,
                                      T z, T& sigma_y,
                                      T& sigma_z, bool diff)
  {
    // Horizontal diffusion parameter.
    if (option_briggs)
      {
        if (option_rural)
          sigma_y =
            ComputeRuralPlumeHorizontalSigma(distance, stability_);
        else
          sigma_y =
            ComputeUrbanPlumeHorizontalSigma(distance, stability_);
      }
    else if (option_doury)
      {
        if (option_day)
          sigma_y =
            ComputeDouryNormalPlumeHorizontalSigma(transfer_time);
        else if (wind_ > 3.)
          sigma_y =
            ComputeDouryNormalPlumeHorizontalSigma(transfer_time);
        else
          sigma_y =
            ComputeDouryLowPlumeHorizontalSigma(transfer_time);
      }
    else
      sigma_y = ComputeCrosswindSigma(transfer_time, z,
                                      friction_velocity_,
                                      convective_velocity_,
                                      boundary_height_,
                                      lmo_, coriolis_, option_hpdm,
                                      stability_);
    // Vertical diffusion parameter.
    if (option_briggs)
      if (option_rural)
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
        if (option_day)
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


  //! Computes concentrations in the whole domain.
  template<class T>
  void GaussianPlume<T>::Compute()
  {
    this->InitCompute();

    // Computes concentrations in the whole domain.
    if (option_compute_domain)
      {
        int i, j, k, species;

        for (species = 0; species < this->Ns; species++)
          for (k = 0; k < this->Nz; k++)
            for (j = 0; j < this->Ny; j++)
              for (i = 0; i < this->Nx; i++)
                this->Concentration(species, k, j, i)
                  = ComputeGaussianConcentration(species, this->GridZ4D(k),
                                                 this->GridY4D(j),
                                                 this->GridX4D(i));

        if (option_NO2_chemistry)
          ComputeNO2Chemistry(this->Concentration);

        if (option_OH_chemistry)
          ComputeOHChemistry(this->Concentration);

        if (option_NO2_chemistry || option_OH_chemistry)
          NullifyNegativeConcentration(this->Concentration.GetArray());
      }

    // Computes concentrations in the point list.
    if (option_compute_list)
      {
        for (int i = 0; i < coordinate_list_.extent(0); i++)
          for (int species = 0; species < this->Ns; species++)
            Concentration_list_(i, species)
              = ComputeGaussianConcentration(species, coordinate_list_(i, 0),
                                             coordinate_list_(i, 1),
                                             coordinate_list_(i, 2));

        if (option_NO2_chemistry)
          ComputeNO2Chemistry(this->Concentration_list_);

        if (option_OH_chemistry)
          ComputeOHChemistry(this->Concentration_list_);

        if (option_NO2_chemistry || option_OH_chemistry)
          NullifyNegativeConcentration(Concentration_list_.GetArray());
      }
  }


  //! Sets the source to be computed. By default, all sources are computed.
  /*!
    \param source_index index of the source to be computed. -1 means all
    sources are computed. Default is -1.
  */
  template<class T>
  void GaussianPlume<T>::RestrictComputationToSource(int source_index)
  {
    if (source_index != -1)
      CheckSourceIndex(source_index);
    computed_source_ = source_index;
  }


  //! Throws an exception if \a source_index is out of range.
  template<class T>
  void GaussianPlume<T>::CheckSourceIndex(int source_index)
  {
    if (source_index < 0 || source_index >= int(SourceList.size()))
      throw "In 'GaussianPlume<T>': Source index out of range.";
  }


  //! Restore old rates.
  template<class T>
  void GaussianPlume<T>::RestoreRate()
  {
    T combination_factor, rate;
    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); iter++)
      {
        if ((*iter)->GetEmissionType() == "continuous_line")
          {
            combination_factor = (*iter)->GetCombinationFactor();
            rate = (*iter)->GetRate();
            if (combination_factor != 1)
              {
                (*iter)->SetRate(rate / (1 - combination_factor));
                (*iter)->SetCombinationFactor(1);
              }
          }
      }
  }


  //! Discretize a line source.
  template<class T>
  void GaussianPlume<T>::Discretization()
  {
    // Cosine and sine of the angle between initial axis and source axis.
    T cos_gamma, sin_gamma, cos_theta, sin_theta;

    Nsource_temp = Nsource;

    int begin = 0, end = SourceList.size();
    if (computed_source_ != -1)
      {
        begin = computed_source_;
        end = begin + 1;
      }

    for (int i = begin; i < end; ++i)
      {
        PlumeSource<T>& source = *SourceList[i];

        // If we have a line source then it is discretized.
        if (source.GetEmissionType() == "continuous_line")
          {
            // Coordinates of the line source in the initial system of
            // coordinate.
            T x, y, z, x1, y1, z1, x2, y2, z2, rate = 0., velocity,
              temperature, diameter, length, delta_x, delta_y, delta_z, width,
              VehicleVelocity, Area, Density;
            int s  = source.GetSpeciesIndex();
            int Np;

            x1 = source.GetX();
            x2 = source.GetX2();
            y1 = source.GetY();
            y2 = source.GetY2();
            z1 = source.GetZ();
            z2 = source.GetZ2();

            // Gets the source width.
            width = source.GetWidth();

            // Get the source Velocity.
            VehicleVelocity = source.GetVehicleVelocity();

            // Get the source Area.
            Area = source.GetArea();

            // Get the source Density.
            Density = source.GetDensity();

            // The line source length must not be 0.
            length = sqrt((x2 - x1) * (x2 - x1) +
                          (y2 - y1) * (y2 - y1) +
                          (z2 - z1) * (z2 - z1));

            // If both extremities of the source are at the same height.
            if (z1 == z2)
              {
                if (length != 0)
                  {
                    cos_gamma = abs(y2 - y1) / length;
                    sin_gamma = abs(x2 - x1) / length;
                  }
                else
                  throw string("There is a line source that has its")
                    + " length equal to 0.";

                // If gamma is between 0 and -90 degrees.
                if (x2 >= x1 && y2 > y1)
                  sin_gamma = -sin_gamma;

                // Else if gamma is between -90 and -180 degrees.
                else if (x2 >= x1 && y2 <= y1)
                  {
                    cos_gamma = -cos_gamma;
                    sin_gamma = -sin_gamma;
                  }
                // Else if gamma is between 90 and 180 degrees.
                else if (x2 < x1 && y2 <= y1)
                  cos_gamma = -cos_gamma;
              }
            else
              throw "The line source must have the same height"
                " at both extremities.";

            // theta = angle_ - gamma: wind angle in the source system
            // of coordinates.
            T cos_theta, sin_theta, tan_theta = 0.;
            cos_theta = cos_angle_ * cos_gamma + sin_angle_ * sin_gamma;
            sin_theta = sin_angle_ * cos_gamma - cos_angle_ * sin_gamma;
            if (cos_theta != 0.)
              tan_theta = sin_theta / cos_theta;

            // If theta is equal to -90 or 90 degree.
            if (abs(cos_theta) < 0.00000001)
              {
                // then theta = -89 degree or 89 degree.
                cos_theta = cos_89;
                if (sin_theta > 0.)
                  sin_theta = sin_89;
                else
                  sin_theta = -sin_89;
              }

            // If the wind is parallel to the road (plus or minus 10 degrees).
            if (abs(sin_theta) > sin(80. * pi / 180.))
              {
                // Coefficient for the combination between point sources and
                // line sources
                T combination_factor = 0.;
                T angle;

                // Compute the value of the wind angle between 0 and 360
                // degrees.
                if (sin_theta > 0.)
                  angle = acos(cos_theta) * 180. / pi;
                else
                  angle = 360. - acos(cos_theta) * 180. / pi;

                // Below, 'combination_factor' varies linearly between 0 and 1
                // when the wind becomes parallel to the road.
                // If combination_factor = 0, only line sources are used.
                // If combination_factor = 1, only point sources are used.

                // Angle between 80 and 90 degrees.
                if (sin_theta >= sin(80. * pi / 180.) && cos_theta >= 0.)
                  combination_factor = angle / 10. - 8.;
                // Angle between 90 and 100 degrees.
                else if (sin_theta >= sin(100. * pi / 180.) && cos_theta < 0.)
                  combination_factor = - angle / 10. + 10.;
                // Angle between 260 and 270 degrees.
                else if (sin_theta <= sin(260. * pi / 180.) && cos_theta <= 0.)
                  combination_factor = angle / 10. - 26.;
                // Angle between 270 and 280 degrees.
                else if (sin_theta <= sin(280. * pi / 180.) && cos_theta > 0.)
                  combination_factor = - angle / 10. + 28.;
                else
                  throw "The line source should not be discretized.";

                Np = min(Np_max, (int(length) + 1) / discretization_step);
                delta_x = (x2 - x1) / (Np - 1);
                delta_y = (y2 - y1) / (Np - 1);
                delta_z = (z2 - z1) / (Np - 1);
                T range = width * 0.5;
                for (T step = -range; step <= range; step++)
                  {
                    T x1_step = x1 + step * cos_gamma;
                    T y1_step = y1 + step * sin_gamma;

                    for (int j = 0; j < Np; j++)
                      {
                        // Incrementation of coordinates.
                        x = x1_step + float(j) * delta_x;
                        y = y1_step + float(j) * delta_y;
                        z = z1 + float(j) * delta_z;

                        // Rate is multiplied by the combination_factor.
                        rate = source.GetRate() * length
                          * combination_factor / Np;
                        velocity = source.GetVelocity();
                        temperature = source.GetTemperature();
                        diameter = 2 * length / Np;

                        PlumeSource<T>* ContinuousSource =
                          new PlumeSource<T>(rate, velocity, temperature,
                                             diameter, x, y, z, s,
                                             Nsource_temp + 1);
                        T sigma_z = PlumeLineSource<T>::vertical_spread(z);
                        ContinuousSource->SetSigma_z_2(pow(sigma_z, 2));
                        SourceList.push_back(ContinuousSource);
                        Nsource_temp++;
                      }
                  }
                // The rate for the line source is multiplied by the
                // combination_factor.
                source.SetRate((1. - combination_factor)
                               * source.GetRate());
                source.SetCombinationFactor(combination_factor);
              }
          }
      }
  }


  //!  Computes the correction factor for the downwind correction.
  template<class T>
  T GaussianPlume<T>::ComputeCorrection(T cos_theta,
                                        T sin_theta, T tan_theta, T y1_source,
                                        T y2_source, T y_source, T  x_source,
                                        T distance_x, T effective_height,
                                        T sigma_y)
  {
    T correction = 0;

    if (ReceptorDownwindPlume(x_source, y_source, y1_source,
                              y2_source, cos_theta, sin_theta))
      {
        // Scale to fit the estimation.
        distance_x = 10 * distance_x;

        T tan_theta = sin_theta / cos_theta;
        // Wind angle in degrees.
        int angle = floor(abs(atan(tan_theta) * 180 / pi) + 0.5);

        if (angle >= threshold_downwind_correction_1_)
          {
            // angle between threshold_downwind_correction_1_
            // and threshold_downwind_correction_2_ included.
            if (angle <= threshold_downwind_correction_2_)
              {
                int id = angle - threshold_downwind_correction_1_;

                correction = (max_1_[id] * exp(-0.5 * pow((mu_1_[id] -
                                                           distance_x)
                                                          / sigma_1_[id], 2)) +
                              max_2_[id] * exp(-0.5 * pow((mu_2_[id] -
                                                           distance_x)
                                                          / sigma_2_[id], 2)));
              }
            // Else angle is between threshold_downwind_correction_2_ excluded
            // and 90 degrees included.
            else
              {
                int id = angle - threshold_downwind_correction_2_ - 1;
                if (angle == 90)
                  // An angle value of 90 degrees is a special since
                  // the range [0,90] is partitioned with open ended
                  // intervals [a, b[
                  id -= 1;

                if (distance_x < mu_[id])
                  correction = (a_[id] + d_[id]) *
                    (1 - exp(- b_[id] * distance_x)) - d_[id];
                else
                  correction =  beta_[id]
                    * exp(-c_[id] * (distance_x - mu_[id]));
              }
          }
      }
    return correction;
  }


  //! Return true if the receptor is downwind of the source.
  template<class T>
  bool GaussianPlume<T>::ReceptorDownwindPlume(T x_source, T y_source,
                                               T y1_source, T y2_source,
                                               T cos_theta, T sin_theta)
  {
    bool res = 0;
    T tan_theta = sin_theta / cos_theta;//

    // Theta in the first quadrant.
    if (cos_theta > 0 && sin_theta >= 0)
      {
        if (tan_theta * x_source + y1_source < y_source &&
            tan_theta * x_source  + y2_source > y_source && x_source >= 0)
          res = 1;
      }
    // Theta in the second quadrant.
    else if (cos_theta < 0 && sin_theta > 0)
      {
        if (tan_theta * x_source + y1_source < y_source &&
            tan_theta * x_source  + y2_source > y_source && x_source <= 0)
          res = 1;
      }
    // Theta in the third quadrant.
    else if (cos_theta < 0 && sin_theta < 0)
      {
        if (tan_theta * x_source + y1_source < y_source &&
            tan_theta * x_source  + y2_source > y_source && x_source <= 0)
          res = 1;
      }
    // Theta in the fourth quadrant.
    else if (cos_theta > 0 && sin_theta < 0)
      {
        if (tan_theta * x_source + y1_source < y_source &&
            tan_theta * x_source  + y2_source > y_source && x_source >= 0)
          res = 1;
      }
    return res;
  }


  //!  Computes the extremities correction factors.
  template<class T>
  T GaussianPlume<T>::ComputeExtremitiesCorrection(T rate, T loss_factor,
                                                   T wind, T inversion_height_,
                                                   T effective_height,
                                                   T distance_x,
                                                   T distance_y,
                                                   T z, T initial_sigma_z_2,
                                                   T overcamp_factor,
                                                   T minimum_volume,
                                                   T cos_theta, int species)
  {
    T extremity = 0, sigma_y, sigma_z,
      time, initial_sigma_y, velocity = 10.1,
      section = 1.25664e-05, initial_sigma_y_2;

    // Corrective coeficients based on estimations.
    T Lambda = 0;
    for (int i = 0; i < int(lambda_.size()); i++)
      Lambda = Lambda + lambda_[i] * pow(cos_theta, i);
    Lambda = exp(Lambda);

    rate = rate / Lambda;

    initial_sigma_y_2 = section / pi;

    // Computing standard deviation sigma_y, with
    // distance_x.
    time = distance_x / wind;
    ComputeSigma(distance_x, time, effective_height,
                 sigma_y, sigma_z);
    sigma_y = sqrt(sigma_y * sigma_y +
                   initial_sigma_y_2);
    if (option_gillani && effective_height
        >= boundary_height_)
      {
        initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
        sigma_z = sqrt(initial_sigma_z_2
                       * (1. + 2.3 * sqrt(time)));
      }
    else
      sigma_z = sqrt(sigma_z * sigma_z
                     + initial_sigma_z_2);

    if (effective_height != 0. && IsPlume(distance_x))
      {
        minimum_volume = section  * velocity * 4;

        // Computing loss factor.
        ComputeLossFactor(distance_x, effective_height,
                          species, loss_factor,
                          overcamp_factor);
        extremity = loss_factor
          * ComputePlumeConcentration(wind, inversion_height_,
                                      effective_height, rate, distance_y, z,
                                      sigma_y, sigma_z, overcamp_factor,
                                      minimum_volume);

      }
    return extremity;
  }


  //!  Read correction coefficients from the correction coefficients file.
  template<class T>
  vector<T> GaussianPlume<T>::ReadCoefficients(string section, int length)
  {
    vector<T> coef;
    vector<string> coefficients_line;
    ConfigStreams coefficients(correction_coefficient_path_);

    if (option_rural)
      coefficients.Find("[rural]");
    else
      coefficients.Find("[urban]");
    coefficients.Find("[" + stability_class_ + "]");
    coefficients.Find(section);
    coefficients_line = split(coefficients.GetLine());
    for (int i = 0; i < length; i++)
      coef.push_back(to_num<T>(coefficients_line[i]));
    return coef;
  }


  /*! \brief Initialization of the correction coefficients for the line source
    correction. */
  template<class T>
  void GaussianPlume<T>::InitCorrectionCoefficients()
  {
    // Correction coefficients are read only if there is a line source.
    bool haveLineSource = false;
    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); iter++)
      if ((*iter)->GetEmissionType() == "continuous_line")
        {
          haveLineSource = true;
          break;
        }
    if (!haveLineSource)
      return;

    int id_1, id_2;
    // Different range of angles for the correction regarding of
    // the stability class.
    threshold_downwind_correction_1_ = 0;
    threshold_downwind_correction_2_ = 0;

    if (option_rural)
      {
        if (stability_ == 0)
          {
            threshold_downwind_correction_1_ = 19;
            threshold_downwind_correction_2_ = 53;
            threshold_extremities_correction_ = 20;
          }
        else if (stability_ == 1)
          {
            threshold_downwind_correction_1_ = 23;
            threshold_downwind_correction_2_ = 60;
            threshold_extremities_correction_ = 30;
          }
        else if (stability_ == 2)
          {
            threshold_downwind_correction_1_ = 30;
            threshold_downwind_correction_2_ = 67;
            threshold_extremities_correction_ = 30;
          }
        else if (stability_ == 3)
          {
            threshold_downwind_correction_1_ = 40;
            threshold_downwind_correction_2_ = 74;
            threshold_extremities_correction_ = 45;
          }
        else if (stability_ == 4)
          {
            threshold_downwind_correction_1_ = 50;
            threshold_downwind_correction_2_ = 75;
            threshold_extremities_correction_ = 55;
          }
        else if (stability_ == 5)
          {
            threshold_downwind_correction_1_ = 71;
            threshold_downwind_correction_2_ = 81;
            threshold_extremities_correction_ = 70;
          }
      }
    else
      {
        if (stability_ == 0 || stability_ == 1)
          {
            threshold_downwind_correction_1_ = 5;
            threshold_downwind_correction_2_ = 43;
            threshold_extremities_correction_ = 10;
          }
        else if (stability_ == 2)
          {
            threshold_downwind_correction_1_ = 9;
            threshold_downwind_correction_2_ = 54;
            threshold_extremities_correction_ = 10;
          }
        else if (stability_ == 3)
          {
            threshold_downwind_correction_1_ = 13;
            threshold_downwind_correction_2_ = 62;
            threshold_extremities_correction_ = 20;
          }
        else if (stability_ == 4 || stability_ == 5)
          {
            threshold_downwind_correction_1_ = 20;
            threshold_downwind_correction_2_ = 69;
            threshold_extremities_correction_ = 20;
          }
      }
    id_1 = threshold_downwind_correction_2_ -
      threshold_downwind_correction_1_ + 1;
    id_2 = 90 - threshold_downwind_correction_2_ - 1;

    // Correction coefficients for downwind receptors. Angles between
    // 40 degrees and 75 degrees.
    max_1_ = ReadCoefficients("[max_1]", id_1);
    mu_1_ = ReadCoefficients("[mu_1]", id_1);
    sigma_1_ = ReadCoefficients("[sigma_1]", id_1);
    max_2_ = ReadCoefficients("[max_2]", id_1);
    mu_2_ = ReadCoefficients("[mu_2]", id_1);
    sigma_2_ = ReadCoefficients("[sigma_2]", id_1);

    // Correction coefficients for downwind receptors. Angles between
    // 75 degrees and 90 degrees.
    mu_ = ReadCoefficients("[mu_]", id_2);
    a_ = ReadCoefficients("[a_]", id_2);
    b_ = ReadCoefficients("[b_]", id_2);
    d_ = ReadCoefficients("[d_]", id_2);
    beta_ = ReadCoefficients("[beta_]", id_2);
    c_ = ReadCoefficients("[b_]", id_2);

    // Coefficients for the correction at extremities of the source.
    lambda_ = ReadCoefficients("[lambda]", 6);
  }


  //! Computes the NO2 chemistry for the whole domain.
  template<class T>
  void GaussianPlume<T>::ComputeNO2Chemistry(Data<T, 4>& concentration)
  {
    T O3_i, NO_i, NO2_i, tmp, k2_k1;
    int nz = concentration.GetLength(1);
    int ny = concentration.GetLength(2);
    int nx = concentration.GetLength(3);
    int Id_O3 = BaseModel<T>::GetSpeciesIndex("O3");
    int Id_NO2 = BaseModel<T>::GetSpeciesIndex("NO2");
    int Id_NO = BaseModel<T>::GetSpeciesIndex("NO");

    for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
          {
            k2_k1 = k2_domain(j, i) / k1_;
            // Initialization of concentrations and convertion in ppb.
            O3_i = concentration(Id_O3, k, j, i) / (0.0409 * 48.);
            NO2_i = concentration(Id_NO2, k, j, i) / (0.0409 * 46.1);
            NO_i = concentration(Id_NO, k, j, i) / (0.0409 * 30.1);

            tmp = (-k2_k1 + O3_i - NO_i +
                   sqrt((k2_k1 - O3_i + NO_i) * (k2_k1 - O3_i + NO_i) +
                        4 * k2_k1 * (O3_i + NO2_i))) / 2.;

            // Final concentration with a convertion in mug.m^(-3).
            concentration(Id_O3, k, j, i) = tmp * 0.0409 * 48.;
            concentration(Id_NO2, k, j, i) = (NO2_i + O3_i - tmp) *
              0.0409 * 46.1;
            concentration(Id_NO, k, j, i) = (NO_i - O3_i + tmp) *
              0.0409 * 30.1;
          }
  }


  //! Computes the NO2 chemistry for a coordinate_list.
  template<class T>
  void GaussianPlume<T>::ComputeNO2Chemistry(Data<T, 2>& concentration)
  {
    T O3_i, NO_i, NO2_i, tmp, k2_k1;
    int Id_O3 = BaseModel<T>::GetSpeciesIndex("O3");
    int Id_NO2 = BaseModel<T>::GetSpeciesIndex("NO2");
    int Id_NO = BaseModel<T>::GetSpeciesIndex("NO");

    int Ncoordinate = coordinate_list_.extent(0);
    for (int i = 0; i < Ncoordinate; ++i)
      {
        k2_k1 = k2_list(i) / k1_;

        // Initialization of concentrations and convertion in ppb.
        O3_i = concentration(i, Id_O3) / (0.0409 * 48.);
        NO2_i = concentration(i, Id_NO2) / (0.0409 * 46.1);
        NO_i = concentration(i, Id_NO) / (0.0409 * 30.1);

        tmp = (-k2_k1 + O3_i - NO_i +
               sqrt((k2_k1 - O3_i + NO_i) * (k2_k1 - O3_i + NO_i) +
                    4 * k2_k1 * (O3_i + NO2_i))) / 2.;

        // Final concentration with a convertion in mug.m^(-3).
        concentration(i, Id_O3) = tmp * 0.0409 * 48.;
        concentration(i, Id_NO2) = (NO2_i + O3_i - tmp) * 0.0409 * 46.1;
        concentration(i, Id_NO) = (NO_i - O3_i + tmp) * 0.0409 * 30.1;
      }
  }


  //! Computes the NO2 photolysis rate and the kinetic constant of the
  // reaction of NO with O3.
  template<class T>
  T GaussianPlume<T>::ComputeNO2PhotolysisRate(T zenith_angle, T Attenuation)
  {
    T K = 0.;
    // Compute the NO2 photolysis rate.
    if (option_day)
      {
        if (zenith_angle >= 0. && zenith_angle < 10.)
          K = 0.9310260000000001e-2 + zenith_angle
            * (zenith_angle
               * (-0.7822279432831311e-6 + zenith_angle
                  * -0.1302720567168795e-7));
        else if (zenith_angle < 20.)
          K = 0.9219010000000000e-2 + (zenith_angle - 10.)
            * (-0.1955272056716901e-4 + (zenith_angle - 10.)
               * (-0.1173044113433769e-5 + (zenith_angle - 10.)
                  * 0.3771617015067078e-8));
        else if (zenith_angle < 30.)
          K = 0.8909950000000000e-2 + (zenith_angle - 20.)
            * (-0.4188211773132428e-4 + (zenith_angle - 20.)
               * (-0.1059895602981758e-5 + (zenith_angle - 20.)
                  * (-0.5859262388581815e-8)));
        else if (zenith_angle < 40.)
          K = 0.8379279999999999e-2 + (zenith_angle - 30.)
            * (-0.6483780850753392e-4 + (zenith_angle - 30.)
               * (-0.1235673474639213e-5 + (zenith_angle - 30.)
                  * (-0.7024567460738029e-8)));
        else if (zenith_angle < 50.)
          K = 0.7600310000000000e-2 + (zenith_angle - 40.)
            * (-0.9165864823853972e-4 + (zenith_angle - 40.)
               * (-0.1446410498461367e-5 + (zenith_angle - 40.)
                  * (-0.9202467768466835e-8)));
        else if (zenith_angle < 60.)
          K = 0.6529880000000000e-2 + (zenith_angle - 50.)
            * (-0.1233475985383066e-3 + (zenith_angle - 50.)
               * (-0.1722484531515342e-5 + (zenith_angle - 50.)
                  * (-0.1612556146540100e-7)));
        else if (zenith_angle < 70.)
          K = 0.5108030000000000e-2 + (zenith_angle - 60.)
            * (-0.1626349576082332e-3 + (zenith_angle - 60.)
               * (-0.2206251375477548e-5 + (zenith_angle - 60.)
                  * 0.3226471363007382e-7));
        else if (zenith_angle < 78.)
          K = 0.3293320000000000e-2 + (zenith_angle - 70.)
            * (-0.1970805710287543e-3 + (zenith_angle - 70.)
               * (-0.1238309966574737e-5 + (zenith_angle - 70.)
                  * (0.2027078243961372e-6)));
        else if (zenith_angle < 86.)
          K = 0.1741210000000000e-2 + (zenith_angle - 78.)
            * (-0.1779736282099126e-3 + (zenith_angle - 78.)
               * (0.3626677818932555e-5 + (zenith_angle - 78.)
                  * (-0.7448311471194499e-7)));
        else if (zenith_angle < 90.)
          K = 0.5113930000000000e-3 + (zenith_angle - 86.)
            * (-0.1342475411316713e-3 + (zenith_angle - 86.)
               * (0.1839083065842406e-5 + (zenith_angle - 86.)
                  * (0.2490309929270573e-5)));
        else
          K = 0.1632080000000000e-3;
        if (Attenuation < 0.99999)
          K *= Attenuation;
      }
    return K;
  }


  //! Computes the zenith angle.
  template<class T>
  T GaussianPlume<T>::ComputeZenithAngle(int nb_sec, T latitude, T longitude)
  {
    const T nb_sec_in_a_day = 60 * 60 * 24;
    T cos_zenith;

    // Computes the real solar time.
    int nb_day_in_year = this->Date_min.GetNumberOfDays();
    T nb_day = nb_sec / nb_sec_in_a_day;
    T t = 2. * pi * (1.1 + nb_day) / nb_day_in_year;
    T corheu = ((0.000075 + 0.001868 * cos(t) - 0.032077 * sin(t))
                - 0.014615 * cos(2. * t) - 0.040849 * sin(2. * t)) * 12. / pi;
    T clock_hour = (nb_day - floor(nb_day)) * 24.;
    T solar_hour = clock_hour + longitude / 15 + corheu
      - 24 * floor(clock_hour + longitude / 15. + corheu / 24.) - 12.;

    // Computes the solar declination.
    T solar_height = solar_hour * pi / 12.;
    T declin = 0.006918 - 0.399912 * cos(t) + 0.070257 * sin(t)
      - 0.006758 * cos(2. * t) + 0.000907 * sin(2. * t)
      - 0.002697 * cos(3. * t) + 0.001480 * sin(3. * t);

    cos_zenith = sin(declin) * sin(latitude * pi / 180.) +
      cos(declin) * cos(latitude * pi / 180.) * cos(solar_height);

    return acos(abs(cos_zenith)) * 180. / pi;
  }


  //! Compute Line source term for a given line source.
  template<class T>
  T GaussianPlume<T>::LineSourceConcentration(int species, T z, T y, T x,
                                              PlumeSource<T>& iter, T step)
  {
    // Distances downwind and crosswind from the source (m).
    T distance_x, distance_x1, distance_x2, distance_y1,
      distance_y2, time;

    // Horizontal and vertical standard deviations (m).
    T initial_sigma_y_2, initial_sigma_z_2;

    // Effective height of release (m), inside and above BL.
    T effective_height;

    // Fraction of the plume that penetrates above BL.
    T penetration_factor;

    // Rate (g/s/m).
    T rate;

    // Output concentration.
    T concentration = 0.;

    // Loss factors.
    T loss_factor, overcamp_factor;

    // Coordinates of the line source in the initial system of
    // coordinate.
    T x1, x2, y1, y2, z1, z2;

    // Coordinates in the source system of coordinate.
    T x_source, y_source, y1_source, y2_source, x1_source,
      x2_source;

    // Cosine and sine of the angle between initial axis and
    // source axis.
    T cos_gamma, sin_gamma;

    // Horizontal and vertical standard deviations (m).
    T sigma_y, sigma_z, sigma_y1 = 0., sigma_y2 = 0., sigma_z1 = 0.,
      sigma_z2 = 0.;

    x1 = iter.GetX();
    x2 = iter.GetX2();
    y1 = iter.GetY();
    y2 = iter.GetY2();
    z1 = iter.GetZ();
    z2 = iter.GetZ2();

    T length = sqrt((x2 - x1) * (x2 - x1) +
                    (y2 - y1) * (y2 - y1) +
                    (z2 - z1) * (z2 - z1));

    if (length == 0.)
      throw "There is a line source that has its length equal to 0.";

    if (z1 != z2)
      throw "The line source must have the same height at both extremities.";

    // First rotation of axis, to be in the source system of
    // coordinate.
    cos_gamma = abs(y2 - y1) / length;
    sin_gamma = abs(x2 - x1) / length;

    // If gamma is between 0 and -90 degrees.
    if (x2 >= x1 && y2 > y1)
      sin_gamma = -sin_gamma;
    // Else if gamma is between -90 and -180 degrees.
    else if (x2 >= x1 && y2 <= y1)
      {
        cos_gamma = -cos_gamma;
        sin_gamma = -sin_gamma;
      }
    // Else if gamma is between 90 and 180 degrees.
    else if (x2 < x1 && y2 <= y1)
      cos_gamma = -cos_gamma;

    // Rotation and translation of axis into the source system
    // of coordinates.
    x1_source = 0.;
    y1_source = 0.;
    x2_source = (x2 - x1) * cos_gamma + (y2 - y1) * sin_gamma;
    y2_source = (x1 - x2) * sin_gamma + (y2 - y1) * cos_gamma;
    x_source = (x - x1) * cos_gamma + (y - y1) * sin_gamma + step;
    y_source = (x1 - x) * sin_gamma + (y - y1) * cos_gamma;

    // theta = wind angle in the source system of coordinates
    //       = angle_ - gamma
    T cos_theta = cos_angle_ * cos_gamma + sin_angle_ * sin_gamma;
    T sin_theta = sin_angle_ * cos_gamma - cos_angle_ * sin_gamma;
    T tan_theta = 0.;
    if (cos_theta != 0.)
      tan_theta = sin_theta / cos_theta;

    // If theta equals +/- 90 degrees, then theta is set to +/- 89 degrees.
    if (abs(cos_theta) < 0.00000001)
      {
        cos_theta = cos_89;
        if (sin_theta > 0)
          sin_theta = sin_89;
        else
          sin_theta = -sin_89;
      }

    // Distance from the point to the first extremity of the
    // source.
    distance_x1 = (x_source - x1_source) * cos_theta +
                  (y_source - y1_source) * sin_theta;
    distance_y1 = (x1_source - x_source) * sin_theta
                  + (y_source - y1_source) * cos_theta;
    distance_y1 = abs(distance_y1);

    // Distance from the point to the second extremity of the
    // source.
    distance_x2 = (x_source - x2_source) * cos_theta +
                  (y_source - y2_source) * sin_theta;
    distance_y2 = (x2_source - x_source) * sin_theta
                  + (y_source - y2_source) * cos_theta;
    distance_y2 = abs(distance_y2);

    // Downwind distance from source.
    distance_x = x_source / cos_theta;

    // Effective height.
    effective_height = iter.GetHeight();
    penetration_factor = iter.GetPenetrationFactor();

    // If there is anything to compute.
    if (IsPlume(distance_x))
      {
        // Initial plume spread.
        initial_sigma_y_2 = iter.GetSigma_y_2();
        initial_sigma_z_2 = iter.GetSigma_z_2();

        if (effective_height != 0.)
          {
            if (distance_x1 >= 0.)
              {
                // Computing standard deviation sigma_y1, with
                // distance_x1.
                time = distance_x1 / wind_;
                ComputeSigma(distance_x1, time, effective_height,
                             sigma_y1, sigma_z1);
                sigma_y1 = sqrt(sigma_y1 * sigma_y1 +
                                initial_sigma_y_2);
              }

            if (distance_x2 >= 0.)
              {
                // Computing standard deviation sigma_y2, with
                // distance_x2.
                time = distance_x2 / wind_;
                ComputeSigma(distance_x2, time, effective_height,
                             sigma_y2, sigma_z2);
                sigma_y2 = sqrt(sigma_y2 * sigma_y2 +
                                initial_sigma_y_2);
              }

            // Computing standard deviations sigma_y and sigma_z,
            // with distance_x.
            time = distance_x / wind_;
            ComputeSigma(distance_x, time, effective_height,
                         sigma_y, sigma_z);

            if (option_gillani && effective_height >= boundary_height_)
              {
                initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
                sigma_z = sqrt(initial_sigma_z_2 * (1. + 2.3 * sqrt(time)));
              }
            else
              sigma_z = sqrt(sigma_z * sigma_z + initial_sigma_z_2);

            // Computing loss factor.
            ComputeLossFactor(distance_x, effective_height, species,
                              loss_factor, overcamp_factor);

            // Minimum volume.
            T minimum_volume = sqrt(2 * pi) * initial_sigma_z_2 * wind_;

            rate = iter.GetRate() * (1. - penetration_factor);
            concentration += loss_factor
              * ComputePlumeLineConcentration(wind_, inversion_height_,
                                              effective_height, rate,
                                              x_source, y_source, z,
                                              sigma_y1, sigma_y2, sigma_z,
                                              overcamp_factor,
                                              y1_source, y2_source,
                                              cos_theta, sin_theta,
                                              minimum_volume);

            if (option_briggs)
              {
                // Downwind correction.
                concentration /=
                  ComputeCorrection(cos_theta, sin_theta,
                                    tan_theta, y1_source, y2_source,
                                    y_source, x_source, distance_x,
                                    iter.GetHeight(), iter.GetSigma_y_2()) + 1;
                // Extremities correction.
                if (abs(cos_theta) <=
                    cos(threshold_extremities_correction_ * pi / 180.))
                  {
                    T extremity_correction_1 =
                      ComputeExtremitiesCorrection(rate, loss_factor, wind_,
                                                   inversion_height_,
                                                   effective_height,
                                                   distance_x1, distance_y1,
                                                   z, iter.GetSigma_z_2(),
                                                   overcamp_factor,
                                                   minimum_volume,
                                                   abs(cos_theta), species);
                    if (concentration < extremity_correction_1)
                      extremity_correction_1 = concentration;

                    T extremity_correction_2 =
                      ComputeExtremitiesCorrection(rate, loss_factor, wind_,
                                                   inversion_height_,
                                                   effective_height,
                                                   distance_x2, distance_y2,
                                                   z, iter.GetSigma_z_2(),
                                                   overcamp_factor,
                                                   minimum_volume,
                                                   abs(cos_theta), species);
                    if (concentration < extremity_correction_2)
                      extremity_correction_2 = concentration;

                    concentration = abs(concentration + extremity_correction_2
                                        - extremity_correction_1);
                  }
              }
          }

      }
    return concentration;
  }


  //!  Computes a concentration Romberg integration over the line source width.
  template<class T>
  T GaussianPlume<T>::RombergLineSourceConcentration(T a, T b, int species,
                                                     T z, T y, T x, PlumeSource<T>& iter)
  {
    // Concentration to be computed with a Romberg integration.
    T concentration = 0.;

    // Romberg integration parameters.
    T del, sum, tnm, x_, dif, dift, ho, hp, error = 0., eps = 0.000001;
    int k, ns, jmax = 20;
    if (option_high_width_precision)
      k = 5;
    else
      k = 4;
    vector<T> c(k + jmax), d(k + jmax), s(jmax + k), h(jmax + k);

    h[1] = 1.;
    for (int j = 1; j != jmax + 1; j++)
      {
        if (j == 1)
          s[j] = 0.5 * (b - a)
            * (LineSourceConcentration(species, z, y, x, iter, a)
               + LineSourceConcentration(species, z, y, x, iter, b));
        else
          {
            tnm  = pow(2, j - 2);
            del = (b - a) / tnm;
            x_ = a + 0.5 * del;
            sum = 0.;
            for (int l = 1; l < tnm + 1; l++)
              {
                sum += LineSourceConcentration(species, z, y, x, iter, x_);
                x_ += del;
              }
            s[j] = 0.5 * (s[j] + (b - a) * sum / tnm);
          }
        if (j >= k)
          {
            ns = 1;
            dif = abs(h[j - k + 2]);
            for (int i = 0; i < k; i++)
              {
                dift = abs(h[j - k + 1 + i]);
                if (dift < dif)
                  {
                    ns = i;
                    dif = dift;
                  }
                c[i] = s[j - k + 1 + i];
                d[i] = c[i];
              }
            concentration = s[j - k + 1 + ns];
            ns--;
            for (int m = 1; m < k - 1; m++)
              {
                for (int i = 1; i <= k - m; i++)
                  {
                    ho = h[j - k + 1 + i];
                    hp = h[j - k + 1 + i + m];
                    if (ho - hp == 0.)
                      throw "From "
                        "'GaussianPlume<T>::RombergLineSourceConcentration': "
                        "Division by zero in Romberg integration.";
                    d[i] = hp * (c[i + 1] - d[i]) / (ho - hp);
                    c[i] = ho * (c[i + 1] - d[i]) / (ho - hp);
                  }
                if (2 * ns < k - m)
                  error = c[ns + 1];
                else
                  {
                    error = d[ns];
                    ns++;
                  }
                concentration += error;
              }

            if (abs(error) <= eps * abs(concentration))
              return concentration;
          }
        s[j + 1] = s[j];
        h[j + 1] = 0.25 * h[j];
      }
    throw "From 'GaussianPlume<T>::RombergLineSourceConcentration': "
      "Too many steps in Romberg integration.";
  }


  //! Returns true if plume is active at the given distance from its source.
  /*!
    \param[in] distance_from_source distance from the plume source.
    \returns true if the plume is active at this distance, false otherwise.
  */
  template<class T>
  bool GaussianPlume<T>::IsPlume(T distance_from_source)
  {
    return distance_from_source > 0.
      && (option_infinite_plume
          || distance_from_source < this->Delta_t * wind_);
  }


  //! Sources initialization.
  /*! It run the methods InitSource.
   */
  template<class T>
  void GaussianPlume<T>::InitSource()
  {
    InitSource(file_source);
  }


  //! Converts cartesian coordinates to longitude/latitude.
  /*!
    \param x abscissa (meters).
    \param y ordinate (meters).
    \param lon longitude (degrees).
    \param lat latitude (degrees).
  */
  template<class T>
  void GaussianPlume<T>::CartesianToLatLon(T x, T y, T& lon, T& lat)
  {
    const T earth_radius = 6371229.;
    lat = (y / earth_radius) * (180. / pi);
    lon = x / (earth_radius * cos(lat * pi / 180.)) * (180. / pi);
  }


  //! Converts longitude/latitude to cartesian coordinates.
  /*!
    \param lon longitude (degrees).
    \param lat latitude (degrees).
    \param x abscissa (meters).
    \param y ordinate (meters).
  */
  template<class T>
  void  GaussianPlume<T>::LatLonToCartesian(T lon, T lat, T& x, T& y)
  {
    const T earth_radius = 6371229.;
    x = earth_radius * cos(lat * pi / 180.) * (lon * pi / 180.);
    y = earth_radius * (lat * pi / 180.);
  }


  //! Retrieves the sources coordinates.
  /*!
    Retrieves the coordinates (abscissa, ordinate, height) of all point
    sources.
    \param[out] source_coordinates array to contain all sources coordinates.
    \note For each source i, source_coordinates(i, 0) is the source abscissa,
    source_coordinates(i, 1) is the ordinate and source_coordinates(i, 2) is
    the source height.
  */
  template<class T>
  void GaussianPlume<T>::GetSourceCoordinates(Array<T, 2>& coordinates)
  {
    int Ncoordinate = SourceList.size();

    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); iter++)
      if ((*iter)->GetEmissionType() == "continuous_line")
        ++Ncoordinate;
    coordinates.resize(Ncoordinate, 3);

    int current_coordinate = 0;
    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); ++iter)
      {
        const PlumeSource<T>& source = **iter;

        coordinates(current_coordinate, 0) = source.GetX();
        coordinates(current_coordinate, 1) = source.GetY();
        coordinates(current_coordinate, 2) = source.GetZ();
        ++current_coordinate;

        if (source.GetEmissionType() == "continuous_line")
          {
            coordinates(current_coordinate, 0) = source.GetX2();
            coordinates(current_coordinate, 1) = source.GetY2();
            coordinates(current_coordinate, 2) = source.GetZ2();
            ++current_coordinate;
          }
      }
  }


  //! Sets the sources coordinates.
  /*!
    Sets the coordinates (abscissa, ordinate, height) of all point sources.
    \param[in] source_coordinates array containing the new set of sources
    coordinates.
    \note For each source i, source_coordinates(i, 0) is the source abscissa,
    source_coordinates(i, 1) is the ordinate and source_coordinates(i, 2) is
    the source height.
  */
  template<class T>
  void GaussianPlume<T>::SetSourceCoordinates(const Array<T, 2>& coordinates)
  {
    int current_coordinates = 0;
    for (source_iterator iter = SourceList.begin();
         iter != SourceList.end(); iter++)
      {
        PlumeSource<T>& source = **iter;

        source.SetX(coordinates(current_coordinates, 0));
        source.SetY(coordinates(current_coordinates, 1));
        source.SetZ(coordinates(current_coordinates, 2));
        ++current_coordinates;
        if (source.GetEmissionType() == "continuous_line")
          {
            source.SetX2(coordinates(current_coordinates, 0));
            source.SetY2(coordinates(current_coordinates, 1));
            source.SetZ2(coordinates(current_coordinates, 2));
            ++current_coordinates;
          }
      }
  }


  //! Prepare for the concentration computation with the plume in grid model.
  /*! \warning This initialization method must be called before
    'ComputeGaussianConcentration' and 'Compute'.
  */
  template<class T>
  void GaussianPlume<T>::InitStep()
  {
    T wind_angle_rad = wind_angle_ / 180. * pi;
    cos_angle_ = cos(wind_angle_);
    sin_angle_ = sin(wind_angle_);

    if (int(SourceList.size()) != Nsource)
      throw string("The number of sources must be ") + to_str<int>(Nsource)
        + " but is " + to_str<int>(SourceList.size()) + ".";

    if (!option_briggs && !option_doury)
      // Matching between LMO and stability classes.
      ComputeStabilityClass(lmo_, stability_, stability_class_);
    if (stability_ != 4 && stability_ != 5)
      inversion_height_ = boundary_height_;
    else
      inversion_height_ = 0.;

    this->Discretization();

    if (option_plume_rise)
      ComputeSourcePlumeRise(*current_plume);
  }


  //! Add a line source to the source list.
  template<class T>
  void GaussianPlume<T>::AddLineSource(int index_, T x_1, T y_1, T z_1, T x_2,
                                       T y_2, T z_2, T rate_, T width_,
                                       T VehicleVelocity_, T Area_,
                                       T Density_, int s_)
  {
    SetCurrentPlume(index_);
    PlumeSource<T>* LineSource =
      new PlumeLineSource<T>(rate_, x_1, y_1, z_1, x_2, y_2, z_2, width_, s_,
                             VehicleVelocity_, Area_, Density_,
                             current_plume->GetIdSource(),
                             current_plume->GetIdSection());
    SourceList.push_back(LineSource);
    ++Nsource;
    Nsource_temp = Nsource;
  }


  //! Erases discretized sources from the source list.
  template<class T>
  void GaussianPlume<T>::EraseDiscretizedSource()
  {
    for (int i = Nsource_temp; i != Nsource; i--)
      {
        PlumeSource<T>* source_to_erase = SourceList.back();
        if (current_plume == source_to_erase)
          current_plume = 0;
        delete source_to_erase;
        SourceList.pop_back();
      }
    if (int(SourceList.size()) != Nsource)
      throw string("The number of sources must be ") + to_str<int>(Nsource)
        + " but is " + to_str<int>(SourceList.size()) + ".";
    Nsource_temp = Nsource;
  }


  //! Erases a list of line sources from the source list.
  template<class T>
  void GaussianPlume<T>::EraseSource(vector<int> index_list)
  {
    std::sort(index_list.begin(), index_list.end(), std::greater<int>());
    for (vector<int>::const_iterator it = index_list.begin();
         it != index_list.end(); ++it)
      {
        int i = *it;
        PlumeSource<T>* plume_to_erase = SourceList[i];
        if (current_plume == plume_to_erase)
          {
            if (i + 1 < int(SourceList.size()))
              current_plume = SourceList[i + 1];
            else
              current_plume = 0;
          }
        SourceList.erase(SourceList.begin() + i);
        delete plume_to_erase;
        Nsource--;
        Nsource_temp = Nsource;
      }
  }


  //! Returns standard deviations sigma_z.
  template<class T>
  T GaussianPlume<T>::ComputeSigmaZ(int index_, T x, T y)
  {
    SetCurrentPlume(index_);
    T sigma_z, sigma_y;
    T effective_height = current_plume->GetHeight();
    T initial_sigma_z_2 = current_plume->GetSigma_z_2();
    T distance_x = 0.;
    if (current_plume->GetEmissionType() == "continuous_line")
      {
        T x1 = current_plume->GetX();
        T x2 = current_plume->GetX2();
        T y1 = current_plume->GetY();
        T y2 = current_plume->GetY2();
        T z1 = current_plume->GetZ();
        T z2 = current_plume->GetZ2();

        // First rotation of axis, to be in the source system of
        // coordinate.

        // The line source length must not be 0.
        T length = sqrt(pow(x2 - x1, 2) +
                        pow(y2 - y1, 2) +
                        pow(z2 - z1, 2));
        if (length == 0)
          throw string("There is a line source that has its")
            + " length equal to 0.";

        T cos_gamma, sin_gamma;
        cos_gamma = abs(y2 - y1) / length;
        sin_gamma = abs(x2 - x1) / length;

        // If gamma is between 0 and -90 degrees.
        if (x2 >= x1 && y2 > y1)
          sin_gamma = -sin_gamma;
        // Else if gamma is between -90 and -180 degrees.
        else if (x2 >= x1 && y2 <= y1)
          {
            cos_gamma = -cos_gamma;
            sin_gamma = -sin_gamma;
          }
        // Else if gamma is between 90 and 180 degrees.
        else if (x2 < x1 && y2 <= y1)
          cos_gamma = -cos_gamma;

        T x_source = (x - x1) * cos_gamma + (y - y1) * sin_gamma;

        // theta = angle_ - gamma
        //       = wind angle in the source system of coordinates.
        T cos_theta = cos_angle_ * cos_gamma + sin_angle_ * sin_gamma;

        // If theta is too close to +/-90°, then theta is set to +/-89°.
        if (abs(cos_theta) < 0.00000001)
          cos_theta = cos_89;

        // Downwind distance from source.
        distance_x = x_source / cos_theta;
      }
    else if (current_plume->GetEmissionType() == "continuous")
      distance_x = (x - current_plume->GetX()) * cos_angle_
        + (y - current_plume->GetY()) * sin_angle_;

    T time = distance_x / wind_;
    ComputeSigma(distance_x, time, effective_height, sigma_y, sigma_z);

    if (option_gillani && effective_height >= boundary_height_)
      {
        initial_sigma_z_2 = max(initial_sigma_z_2, 9.);
        sigma_z = sqrt(initial_sigma_z_2 * (1. + 2.3 * sqrt(time)));
      }
    else
      sigma_z = sqrt(pow(sigma_z, 2) + initial_sigma_z_2);

    return sigma_z;
  }


  // Compute Chemistry.
  template<class T>
  void GaussianPlume<T>::Chemistry_Plume(const Date& current_date,
                                         const Array<T, 2>& latitude,
                                         const Array<T, 2>& longitude,
                                         Data<T, 4> BackgroundConcentration,
                                         Data<T, 4>& Concentration_out)
  {
    if (!WithChemistry())
      return;

    int Ns = BackgroundConcentration.GetLength(0);
    int Nz = BackgroundConcentration.GetLength(1);
    int Ny = BackgroundConcentration.GetLength(2);
    int Nx = BackgroundConcentration.GetLength(3);
    Data<T, 4> TotalConcentration(Ns, Nz, Ny, Nx);

    // NO2 Chemistry.
    if (option_NO2_chemistry)
      {
        // Compute Kinetic constant of the reaction of NO with O3 and NO2
        // photolysis rate.
        k2_domain.resize(Ny, Nx);
        int nb_sec = current_date.GetNumberOfSeconds();
        for (int j = 0; j < Ny; j++)
          for (int i = 0; i < Nx; i++)
            {
              T zenith_angle = ComputeZenithAngle(nb_sec, latitude(j, i),
                                                  longitude(j, i));
              k2_domain(j, i) = ComputeNO2PhotolysisRate(zenith_angle,
                                                         Attenuation_grid(0, j, i));
            }

        // Compute the kinetic constant of the reaction of NO with O3
        // in ppb^-1.seconds^-1.
        k1_ = exp(-27.29454887930734 - 1310. / temperature_);
        k1_ = k1_ * 2.46e10;

        // NO2 Chemistry.
        TotalConcentration.GetArray() =
          Concentration_out.GetArray() + BackgroundConcentration.GetArray();

        // Chemistry for total concentrations.
        ComputeNO2Chemistry(TotalConcentration);

        // Chemistry for background concentrations only.
        ComputeNO2Chemistry(BackgroundConcentration);

        // Writes the new concentrations.
        Concentration_out.GetArray() = TotalConcentration.GetArray()
          - BackgroundConcentration.GetArray();
        NullifyNegativeConcentration(Concentration_out.GetArray(), "NO");
        NullifyNegativeConcentration(Concentration_out.GetArray(), "NO2");
      }

    // OH Chemistry.
    if (option_OH_chemistry)
      {
        TotalConcentration.GetArray() = Concentration_out.GetArray()
          + BackgroundConcentration.GetArray();

        // Chemistry for total concentrations.
        ComputeOHChemistry(TotalConcentration, &latitude, &longitude);

        // Chemistry for background concentrations only.
        ComputeOHChemistry(BackgroundConcentration, &latitude,
                           &longitude);

        // Writes the new concentrations.
        Concentration_out.GetArray() = TotalConcentration.GetArray() -
          - BackgroundConcentration.GetArray();
        NullifyNegativeConcentration(Concentration_out.GetArray());
      }
  }


  // Compute Chemistry.
  template<class T>
  void GaussianPlume<T>::
  Chemistry_Plume_list(const Date& current_date,
                       Data<T, 2> BackgroundConcentration)
  {
    if (!WithChemistry())
      return;

    int Ns = BackgroundConcentration.GetLength(1);
    Data<T, 2> TotalConcentration(this->coordinate_list_.extent(0), Ns);

    // NO2 Chemistry.
    if (option_NO2_chemistry)
      {
        // Compute Kinetic constant of the reaction of NO with O3 and
        // NO2 photolysis rate.
        k2_list.resize(this->coordinate_list_.extent(0));
        T zenith_angle, latitude, longitude;
        int nb_sec = current_date.GetNumberOfSeconds();
        for (int i = 0; i < this->coordinate_list_.extent(0); i++)
          {
            CartesianToLatLon(coordinate_list_(i, 2), coordinate_list_(i, 1)
                              , longitude, latitude);
            zenith_angle = ComputeZenithAngle(nb_sec, latitude, longitude);
            k2_list(i) = ComputeNO2PhotolysisRate(zenith_angle,
                                                  Attenuation_list(i));
          }

        // Compute the kinetic constant of the reaction of NO with O3
        // in ppb^-1.seconds^-1.
        k1_ = exp(-27.29454887930734 - 1310. / temperature_);
        k1_ = k1_ * 2.46 * pow(10, 10);

        // NO2 chemistry.
        TotalConcentration.GetArray() = Concentration_list_.GetArray()
          + BackgroundConcentration.GetArray();

        // Chemistry for total concentrations.
        ComputeNO2Chemistry(TotalConcentration);

        // Chemistry for background concentrations only.
        ComputeNO2Chemistry(BackgroundConcentration);

        Concentration_list_.GetArray() = TotalConcentration.GetArray()
          - BackgroundConcentration.GetArray();
        NullifyNegativeConcentration(Concentration_list_.GetArray(), "NO");
        NullifyNegativeConcentration(Concentration_list_.GetArray(), "NO2");
      }

    // OH Chemistry.
    if (option_OH_chemistry)
      {
        TotalConcentration.GetArray() = Concentration_list_.GetArray()
          + BackgroundConcentration.GetArray();

        // Chemistry for total concentrations.
        ComputeOHChemistry(TotalConcentration, true);

        // Chemistry for background concentrations only.
        ComputeOHChemistry(BackgroundConcentration, true);

        Concentration_list_.GetArray() = TotalConcentration.GetArray()
          - BackgroundConcentration.GetArray();
        NullifyNegativeConcentration(Concentration_list_.GetArray());
      }
  }

  //! Reset to zero any negative concentrations of the given species.
  /*!
    \param[in, out] concentrationData concentrations to be reset to zero on
    negative concentrations.
    \param[in] species_name name of the species to consider.
  */
  template<class T>
  template<int Dim>
  void GaussianPlume<T>::NullifyNegativeConcentration(
                                                      Array<T, Dim>& concentration,
                                                      const string& species_name) const
  {
    int species = BaseModel<T>::GetSpeciesIndex(species_name);
    Array<T, Dim> concentration_subarray =
      concentration(Range(species, species));
    NullifyNegativeConcentration(concentration);
  }


  //! Reset to zero any negative concentrations.
  /*!
    \param[in, out] concentrationData concentrations to be reset to zero on
    negative concentrations.
  */
  template<class T>
  template<int Dim>
  void GaussianPlume<T>::NullifyNegativeConcentration(
                                                      Array<T, Dim>& concentration) const
  {
    for (typename Array<T, Dim>::iterator it = concentration.begin();
         it != concentration.end(); ++it)
      {
        T& concentration = *it;
        if (concentration < 0.)
          concentration = 0.;
      }
  }


  //! Return the coordinate list of point to be computed.
  /*!
   */
  template<class T>
  Array<T, 2> GaussianPlume<T>::GetCoordinateList()
  {
    return this->coordinate_list_;
  }


  //! Sets the coordinate list to zero.
  /*!
   */
  template<class T>
  void GaussianPlume<T>::SetZeroConcentrationList()
  {
    this->Concentration_list_.SetZero();
  }


  //! This method is not implemented and throws an exception on call.
  /*! It is not possible to compute the mean concentration over a given
   *  volume in the case of the "continuous plume in grid" model.
   */
  template<class T>
  T GaussianPlume<T>::GetIntegratedConcentration(int species, T z, T y, T x,
                                                 T lz, T ly, T lx)
  {
    throw "The method GaussianPlume<T>::GetIntegratedConcentration is not "
      "implemented. Do not use the spatial averaged option with "
      "\"continuous plume in grid\" model";
  }


  //! Computes concentrations in the whole domain.
  template<class T>
  void GaussianPlume<T>::Compute_list(int index_)
  {
    SetCurrentPlume(index_);
    int species = current_plume->GetSpeciesIndex();
    for (int i = 0; i < coordinate_list_.extent(0); i++)
      {
        T point_concentration =
          ComputeGaussianConcentration(species, coordinate_list_(i, 0),
                                       coordinate_list_(i, 1),
                                       coordinate_list_(i, 2));

        this->Concentration_list_(i, species) += point_concentration;
      }
  }


  //! Fill the Point source list that needs to be computed.
  template<class T>
  void GaussianPlume<T>::FillPointSourceList(int begin, int end)
  {
    if (begin == end)
      return;
    if (begin > end)
      throw string("Invalid index range [") + to_str(begin) + ", "
        + to_str(end) + "]: the first index must be inferior to the second.";

    CheckSourceIndex(begin);
    if (end > begin + 1)
      CheckSourceIndex(end - 1);

    for (int i = begin; i < end; ++i)
      PointSourceList.push_back(SourceList[i]);
  }


  //! Fill the Point source list that needs to be computed.
  template<class T>
  void GaussianPlume<T>::FillPointSourceList(int index)
  {
    CheckSourceIndex(index);
    PointSourceList.push_back(SourceList[index]);
  }


  //! Empty the Point source list that needs to be computed.
  template<class T>
  void GaussianPlume<T>::EmptyPointSourceList()
  {
    PointSourceList.clear();
  }


  //! Save current_plume iterator
  template<class T>
  void GaussianPlume<T>::SaveCurrentPlume()
  {
    previous_plume = current_plume;
  }


  //! Restore current_plume iterator
  template<class T>
  void GaussianPlume<T>::RestoreCurrentPlume()
  {
    current_plume = previous_plume;
  }


  template<class T>
  void GaussianPlume<T>::SplitSourceList(int rank, int Nproc)
  {
    vector<PlumeSource<T>* > rankSourceList;
    int N = SourceList.size();
    for (int i = rank % Nproc; i < N; i += Nproc)
      rankSourceList.push_back(SourceList[i]);
    SourceList.swap(rankSourceList);
    Nsource = Nsource_temp = SourceList.size();
  }


  template<class T>
  Array<T, 2>& GaussianPlume<T>::ConcentrationList()
  {
    return Concentration_list_.GetArray();
  }


  template<class T>
  const Array<T, 2>& GaussianPlume<T>::ConcentrationList() const
  {
    return Concentration_list_.GetArray();
  }


  template<class T>
  void GaussianPlume<T>::MultiplySourcesRate(T coefficient)
  {
    for (typename vector<PlumeSource<T>* >::const_iterator iter =
           SourceList.begin(); iter != SourceList.end(); ++iter)
      (*iter)->SetRate((*iter)->GetRate() * coefficient);
  }

  // Compute Chemistry for OH.
  template<class T>
  void GaussianPlume<T>::ComputeOHChemistry(Data<T, 4>& concentration,
                                            const Array<T, 2>* latitude_matrix,
                                            const Array<T, 2>* longitude_matrix)
  {
    T time = 0., H2O_i = 0.;
    T Na = 6.02213e23;
    int nz = concentration.GetLength(1);
    int ny = concentration.GetLength(2);
    int nx = concentration.GetLength(3);
    int Id_BUTA = BaseModel<T>::GetSpeciesIndex("BUTA");
    int Id_NO2 = BaseModel<T>::GetSpeciesIndex("NO2");
    int Id_HCHO = BaseModel<T>::GetSpeciesIndex("HCHO");
    int Id_O3 = BaseModel<T>::GetSpeciesIndex("O3");

    T plume_length[][6] = {{ 53,  53,  67,  96, 175, 175 },
                           { 67, 107, 161, 240, 409, 744 }
    };

    if (longitude_matrix == 0)
      {
        time = (200. + 0.5 *
                plume_length[int(this->GetModelParameter("rural"))]
                [int(this->GetModelParameter("stability"))]) / wind_;

        // In ppmv;
        T t = 1 - 373.15 / temperature_;
        T Ps = 1e2 * 101325. * exp(13.3185 * t - 1.976 * pow(t, 2) - 0.6445 *
                                   pow(t, 3) - 0.1299 * pow(t, 4));
        H2O_i = 1e6 * RelativeHumidity_ * Ps / Pressure_;

        // Conversion in µg.m^-3.
        H2O_i = H2O_i * 40.9 * 18.;
      }

    // Reaction rates in cm^3.molecule^-1.s^-1.
    T k_but = 6.78e-11; // at T = 295
    T k_HCHO = 9.0e-12; // at T = 298

    // In ppm^-1.min^-1. (RIVAD)
    T k_O1D = 4.45e10;
    T k_H2O = 3.4e5;
    T k_no2 = 1.4e4;

    // Conversion in µg^-1.m^3.s^-1.
    k_but = (k_but * Na / 54.) * 1e-12;
    k_HCHO = (k_HCHO * Na / 30.03) * 1e-12;
    k_O1D = k_O1D / (60. * 40.9 * 16.);
    k_H2O = k_H2O / (60. * 40.9 * 18.);
    k_no2 = k_no2 / (60. * 40.9 * 46.1);

    T O1D, OH, BUTA, NO2, HCHO, O3;
    T BUTA_i, NO2_i, HCHO_i, O3_i;
    T k_O3 = 0., latitude, longitude, zenith_angle;
    T nb_sec = this->current_date.GetNumberOfSeconds();
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        {
          if (longitude_matrix == 0)
            {
              CartesianToLatLon(this->GridX4D(i), this->GridY4D(j),
                                longitude, latitude);

              zenith_angle = ComputeZenithAngle(nb_sec, latitude, longitude);

              // In s^-1.
              k_O3 = ComputeO3PhotolysisRate(zenith_angle);
            }
          else
            {
              latitude = (*latitude_matrix)(j, i);
              longitude = (*longitude_matrix)(j, i);

              zenith_angle = ComputeZenithAngle(nb_sec, latitude, longitude);
            }

          for (int k = 0; k < nz; k++)
            {
              if (longitude_matrix == 0)
                {
                  // In s^-1.
                  k_O3 = ComputeO3PhotolysisRate(zenith_angle,
                                                 Attenuation_grid(k, j, i));

                  time = (200. + 0.5 * plume_length[int(Rural_grid(k, j, i))]
                          [int(Stability_grid(k, j, i))]) / Wind_grid(k, j, i);

                  // Compute third body.
                  // Conversion = Avogadro*1d-6/Perfect gas constant.
                  // PRESS in Pascal, SUMM in molecules/cm3, TEMP in Kelvin
                  T SumM = Pressure_grid(k, j, i) * 7.243 * pow(10, 16)
                    / Temperature_grid(k, j, i);

                  // Number of water molecules computed from the massic
                  // fraction (absolute humidity) in mollecules.cm^-3
                  H2O_i = 29. * SumM * SpecificHumidity_grid(k, j, i)
                    / (18. + 11. * SpecificHumidity_grid(k, j, i));

                  // Conversion in µg.m^-3.
                  H2O_i = H2O_i * 40.9 * 18. / 2.46 * pow(10, -13);
                }

              BUTA_i = concentration(Id_BUTA, k, j, i);
              NO2_i = concentration(Id_NO2, k, j, i);
              HCHO_i = concentration(Id_HCHO, k, j, i);
              O3_i = concentration(Id_O3, k, j, i);

              if (BUTA_i + NO2_i + HCHO_i != 0.)
                {
                  // Compute O1D and OH
                  O1D = k_O3 * O3_i / (k_O1D + k_H2O * H2O_i);
                  OH = 2. * k_H2O * O1D * H2O_i / (k_but * BUTA_i +
                                                   k_no2 * NO2_i +
                                                   k_HCHO * HCHO_i);

                  // Update NO2, HCHO, BUTA and O3 concentrations.
                  NO2 = NO2_i * exp(-k_no2 * OH * time);
                  HCHO = HCHO_i * exp(-k_HCHO * OH * time);
                  BUTA = BUTA_i * exp(-k_but * OH * time);

                  concentration(Id_BUTA, k, j, i) = BUTA;
                  concentration(Id_NO2, k, j, i) = NO2;
                  concentration(Id_HCHO, k, j, i) = HCHO;
                }
            }
        }
  }


  // Compute Chemistry for OH.
  template<class T>
  void GaussianPlume<T>::ComputeOHChemistry(Data<T, 2>& concentration,
                                            bool option_plume_in_grid)
  {
    T time = 0., H2O_i = 0.;
    T Na = 6.02213e23;
    int Id_BUTA = BaseModel<T>::GetSpeciesIndex("BUTA");
    int Id_NO2 = BaseModel<T>::GetSpeciesIndex("NO2");
    int Id_HCHO = BaseModel<T>::GetSpeciesIndex("HCHO");
    int Id_O3 = BaseModel<T>::GetSpeciesIndex("O3");

    T plume_length[][6] = {{ 53,  53,  67,  96, 175, 175 },
                           { 67, 107, 161, 240, 409, 744 }
    };

    if (!option_plume_in_grid)
      {
        time = (200. + 0.5 *
                plume_length[int(this->GetModelParameter("rural"))]
                [int(this->GetModelParameter("stability"))]) / wind_;

        // In ppmv;
        T t = 1 - 373.15 / temperature_;
        T Ps = 1e2 * 101325. * exp(13.3185 * t - 1.976 * pow(t, 2)
                                   - 0.6445 * pow(t, 3) - 0.1299 * pow(t, 4)) * pow(10, 2);
        H2O_i = 1e6 * RelativeHumidity_ * Ps / Pressure_;

        // Conversion in µg.m^-3.
        H2O_i = H2O_i * 40.9 * 18.;
      }

    // Reaction rates in cm^3.molecule^-1.s^-1.
    T k_but = 6.78e-11; // at T = 295
    T k_HCHO = 9.0e-12; // at T = 298

    // In ppm^-1.min^-1. (RIVAD)
    T k_O1D = 4.45e10;
    T k_H2O = 3.4e5;
    T k_no2 = 1.4e4;

    // Conversion in µg^-1.m^3.s^-1.
    k_but = (k_but * Na / 54.) * 1e-12;
    k_HCHO = (k_HCHO * Na / 30.03) * 1e-12;
    k_O1D = k_O1D / (60. * 40.9 * 16.);
    k_H2O = k_H2O / (60. * 40.9 * 18.);
    k_no2 = k_no2 / (60. * 40.9 * 46.1);

    T O1D, OH, BUTA, NO2, HCHO, O3;
    T BUTA_i, NO2_i, HCHO_i, O3_i;
    T k_O3 = 0., latitude, longitude, zenith_angle;
    T nb_sec = this->current_date.GetNumberOfSeconds();

    for (int i = 0; i < coordinate_list_.extent(0); i++)
      {
        if (option_cartesian)
          CartesianToLatLon(coordinate_list_(i, 2), coordinate_list_(i, 1),
                            longitude, latitude);
        else
          {
            latitude = coordinate_list_(i, 1);
            longitude = coordinate_list_(i, 2);
          }

        zenith_angle = ComputeZenithAngle(nb_sec, latitude, longitude);

        if (option_plume_in_grid)
          {
            // In s^-1.
            k_O3 = ComputeO3PhotolysisRate(zenith_angle, Attenuation_list(i));

            time = (200. + 0.5 * plume_length[int(Rural_list(i))]
                    [int(Stability_list(i))]) / Wind_list(i);

            // Compute third body.
            // Conversion = Avogadro*1d-6/Perfect gas constant.
            // PRESS in Pascal, SUMM in molecules/cm3, TEMP in Kelvin
            T SumM = Pressure_list(i) * 7.243 * pow(10, 16)
              / Temperature_list(i);

            // Number of water molecules computed from the massic fraction
            // (absolute humidity) H2O_i in mollecules.cm^-3
            H2O_i = 29. * SumM * SpecificHumidity_list(i)
              / (18. + 11. * SpecificHumidity_list(i));

            // Conversion in µg.m^-3.
            H2O_i = H2O_i * 40.9 * 18. / 2.46 * pow(10, -13);
          }
        else
          {
            // In s^-1.
            k_O3 = ComputeO3PhotolysisRate(zenith_angle);
          }

        BUTA_i = concentration(i, Id_BUTA);
        NO2_i = concentration(i, Id_NO2);
        HCHO_i = concentration(i, Id_HCHO);
        O3_i = concentration(i, Id_O3);

        if (BUTA_i + NO2_i + HCHO_i != 0.)
          {
            // Compute O1D and OH
            O1D = k_O3 * O3_i / (k_O1D + k_H2O * H2O_i);
            OH = 2. * k_H2O * O1D * H2O_i / (k_but * BUTA_i + k_no2 * NO2_i
                                             + k_HCHO * HCHO_i);

            // Update NO2, HCHO, BUTA and O3 concentrations.
            NO2 = NO2_i * exp(-k_no2 * OH * time);
            HCHO = HCHO_i * exp(-k_HCHO * OH * time);
            BUTA = BUTA_i * exp(-k_but * OH * time);

            concentration(i, Id_BUTA) = BUTA;
            concentration(i, Id_NO2) = NO2;
            concentration(i, Id_HCHO) = HCHO;
          }
      }
  }


  //! Compute the O3 photolysis rate.
  template<class T>
  T GaussianPlume<T>::ComputeO3PhotolysisRate(T zenith_angle,
                                              T Attenuation)
  {
    T O3_photolysis_rate;

    // Compute the O3 photolysis rate.
    if (!option_day)
      O3_photolysis_rate = 0.;
    else
      {
        if (zenith_angle < 10.)
          {
            O3_photolysis_rate = -0.8776629099833464e-10;
            O3_photolysis_rate = -0.1165033709001661e-7
              + zenith_angle * O3_photolysis_rate;
            O3_photolysis_rate = zenith_angle * O3_photolysis_rate;
            O3_photolysis_rate = 0.3523480000000000e-4
              + zenith_angle * O3_photolysis_rate;
          }
        else if (zenith_angle  < 20.)
          {
            O3_photolysis_rate = 0.1474988729949909e-9;
            O3_photolysis_rate = -0.1428332581996665e-7
              + (zenith_angle - 0.10e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.2593366290998327e-6
              + (zenith_angle - 0.10e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.3398200000000000e-4
              + (zenith_angle - 0.10e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 30.)
          {
            O3_photolysis_rate = 0.1300707990183827e-9;
            O3_photolysis_rate = -0.9858359630116941e-8
              + (zenith_angle - 0.20e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.5007534836006686e-6
              + (zenith_angle - 0.20e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.3010780000000000e-4
              + (zenith_angle - 0.20e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 40.)
          {
            O3_photolysis_rate = 0.1988179309314732e-9;
            O3_photolysis_rate = -0.5956235659565481e-8
              + (zenith_angle - 0.30e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.6588994364974921e-6
              + (zenith_angle - 0.30e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.2424450000000000e-4
              + (zenith_angle - 0.30e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 50.)
          {
            O3_photolysis_rate = 0.2219574772557277e-9;
            O3_photolysis_rate = 0.8302268378835231e-11
              + (zenith_angle - 0.40e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.7183787704093613e-6
              + (zenith_angle - 0.40e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.1725870000000000e-4
              + (zenith_angle - 0.40e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 60.)
          {
            O3_photolysis_rate = 0.1913521600455895e-9;
            O3_photolysis_rate = 0.6667026586050136e-8
              + (zenith_angle - 0.50e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.6516254818650674e-6
              + (zenith_angle - 0.50e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.1029770000000000e-4
              + (zenith_angle - 0.50e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 70.)
          {
            O3_photolysis_rate = 0.1602388256152816e-10;
            O3_photolysis_rate = 0.1240759138741867e-7
              + (zenith_angle - 0.60e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.4608793021303539e-6
              + (zenith_angle - 0.60e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.4639500000000000e-5
              + (zenith_angle - 0.60e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 78.)
          {
            O3_photolysis_rate = -0.3089359890783960e-9;
            O3_photolysis_rate = 0.1288830786428400e-7
              + (zenith_angle - 0.70e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.2079203096133002e-6
              + (zenith_angle - 0.70e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.1287490000000000e-5
              + (zenith_angle - 0.70e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 86.)
          {
            O3_photolysis_rate = -0.2034952628632162e-9;
            O3_photolysis_rate = 0.5473844126395715e-8
              + (zenith_angle - 0.78e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.6102309368797090e-7
              + (zenith_angle - 0.78e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.2908040000000000e-6
              + (zenith_angle - 0.78e2) * O3_photolysis_rate;
          }
        else if (zenith_angle  < 90.)
          {
            O3_photolysis_rate = 0.1623544916128788e-9;
            O3_photolysis_rate = 0.5899578175158973e-9
              + (zenith_angle - 0.86e2) * O3_photolysis_rate;
            O3_photolysis_rate = -0.1251267813581064e-7
              + (zenith_angle - 0.86e2) * O3_photolysis_rate;
            O3_photolysis_rate = 0.4875570000000000e-7
              + (zenith_angle - 0.86e2) * O3_photolysis_rate;
          }
        else
          {
            O3_photolysis_rate = 0.1853500000000000e-7;
          }
        if (Attenuation < 0.99999)
          O3_photolysis_rate *= Attenuation;
      }
    return O3_photolysis_rate;
  }


  //! Sets meteorological parameters for the chemistry.
  template<class T>
  void GaussianPlume<T>::
  SetChemistryMeteorologicalParameters(Array<T, 3> SpecificHumidity_grid_,
                                       Array<T, 3> Temperature_grid_,
                                       Array<T, 3> Pressure_grid_,
                                       Array<int, 3> Stability_grid_,
                                       Array<int, 3> Rural_grid_,
                                       Array<T, 3> Wind_grid_,
                                       Array<T, 3> Attenuation_grid_)
  {
    SpecificHumidity_grid.resize(SpecificHumidity_grid_.shape());
    Temperature_grid.resize(Temperature_grid_.shape());
    Pressure_grid.resize(Pressure_grid_.shape());
    Stability_grid.resize(Stability_grid_.shape());
    Rural_grid.resize(Rural_grid_.shape());
    Wind_grid.resize(Wind_grid_.shape());
    Attenuation_grid.resize(Attenuation_grid_.shape());

    SpecificHumidity_grid = SpecificHumidity_grid_;
    Temperature_grid = Temperature_grid_;
    Pressure_grid = Pressure_grid_;
    Stability_grid = Stability_grid_;
    Rural_grid = Rural_grid_;
    Wind_grid = Wind_grid_;
    Attenuation_grid = Attenuation_grid_;
  }


  //! Sets meteorological parameters for the chemistry.
  template<class T>
  void GaussianPlume<T>::
  SetChemistryMeteorologicalParameters(Array<T, 1> SpecificHumidity_list_,
                                       Array<T, 1> Temperature_list_,
                                       Array<T, 1> Pressure_list_,
                                       Array<int, 1> Stability_list_,
                                       Array<int, 1> Rural_list_,
                                       Array<T, 1> Wind_list_,
                                       Array<T, 1> Attenuation_list_)
  {
    SpecificHumidity_list.resize(SpecificHumidity_list_.shape());
    Temperature_list.resize(Temperature_list_.shape());
    Pressure_list.resize(Pressure_list_.shape());
    Stability_list.resize(Stability_list_.shape());
    Rural_list.resize(Rural_list_.shape());
    Wind_list.resize(Wind_list_.shape());
    Attenuation_list.resize(Attenuation_list_.shape());

    SpecificHumidity_list = SpecificHumidity_list_;
    Temperature_list = Temperature_list_;
    Pressure_list = Pressure_list_;
    Stability_list = Stability_list_;
    Rural_list = Rural_list_;
    Wind_list = Wind_list_;
    Attenuation_list = Attenuation_list_;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_CXX
#endif
