// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Ir?e Korsakissok, Vivien Mallet, Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFTRANSPORT_CXX


#include "PuffTransport.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////



  //! Main constructor.
  /*!
    \param time_puff time at which the puff is released (seconds).
    \param quantity total mass released by the source (mass unit).
    \param source_abscissa abscissa (meters).
    \param source_ordinate ordinate (meters).
    \param source_height height (meters).
    \param species_index index associated to the species.
  */
  template<class T>
  Puff<T>::Puff(T release_time, T velocity, T temperature, T diameter,
		T quantity, T abscissa, T ordinate, T height, T source_water,
		T volume_prev,
		int species_index, string source_id):
    release_time_(release_time), quantity_(quantity),
    source_velocity_(velocity), source_temperature_(temperature),
    source_diameter_(diameter), source_abscissa_(abscissa),
    source_ordinate_(ordinate), source_height_(height), 
    source_water_(source_water), volume_prev_(volume_prev),
    source_id_(source_id),
    z_above_(0.), penetration_factor_(0.),
    sigma_x_(0.), sigma_y_(0.), sigma_z_(0.),
    sigma_y_2_(diameter * diameter / 4.), sigma_z_2_(0.),
    Ns(1), species_index_(species_index),
    scavenging_coefficient_(0.), deposition_velocity_(0.),
    reflection_factor_(1.), has_meteo(0)
  {
  }

  //! Alternative constructor to model a puff over the volume source.
  /*!
    \param time_puff time at which the puff is released (seconds).
    \param quantity total mass released by the source (mass unit).
    \param source_abscissa abscissa (meters).
    \param source_ordinate ordinate (meters).
    \param source_height height (meters).
    \param species_index index associated to the species.
    \param width volume source width (meters).
    \param length volume source length (meters).
    \param height volume source height (meters).
  */
  template<class T>
  Puff<T>::Puff(T release_time, T velocity, T temperature, T width,
                T length, T quantity,
		T abscissa, T ordinate, T height, T source_water,
		T volume_prev,
		int species_index, string source_id):
    release_time_(release_time), quantity_(quantity),
    source_velocity_(velocity), source_temperature_(temperature),
    source_diameter_(0.), source_abscissa_(abscissa),
    source_ordinate_(ordinate), source_height_(height), source_water_(source_water),
    volume_prev_(volume_prev), source_id_(source_id),
    z_above_(0.), penetration_factor_(0.),
    sigma_x_(0.), sigma_y_(0.), sigma_z_(0.),
    sigma_y_2_((width * width + length * length) / 18.49),
    sigma_z_2_(height * height / 18.49),
    Ns(1), species_index_(species_index),
    scavenging_coefficient_(0.), deposition_velocity_(0.),
    reflection_factor_(1.), has_meteo(0)
  {
  }

  //! Alternative constructor to model a puff over the volume source.
  /*!
    \param time_puff time at which the puff is released (seconds).
    \param quantity total mass released by the source (mass unit).
    \param source_abscissa abscissa (meters).
    \param source_ordinate ordinate (meters).
    \param source_height height (meters).
    \param species_index index associated to the species.
    \param width volume source width (meters).
    \param length volume source length (meters).
    \param height volume source height (meters).
    \warning source_height should not be zero.
  */
  template<class T>
  Puff<T>::Puff(T release_time, T velocity, T temperature, T diameter, T width,
                T length, 
		T quantity, T abscissa, T ordinate, T height, T source_water,
		T volume_prev,
		int species_index, bool is_volume_source, string source_id):
    release_time_(release_time), quantity_(quantity),
    source_velocity_(velocity), source_temperature_(temperature),
    source_diameter_(diameter), source_abscissa_(abscissa),
    source_ordinate_(ordinate), source_water_(source_water),
    volume_prev_(volume_prev), source_id_(source_id), 
    z_above_(0.), penetration_factor_(0.),
    sigma_x_(0.), sigma_y_(0.), sigma_z_(0.),
    Ns(1), species_index_(species_index),
    scavenging_coefficient_(0.), deposition_velocity_(0.),
    reflection_factor_(1.), has_meteo(0), is_volume_source_(is_volume_source)
  {
    // cout << "velocity in Puff.cxx: " << velocity << endl;
    // cout << "temperature in Puff.cxx: " << temperature << endl;
    if (is_volume_source)
      {
        sigma_y_2_ = (width * width + length * length) / 18.49;
        sigma_x_2_ = sigma_y_2_;
        sigma_z_2_ = height * height / 18.49;
        source_height_ = height / 2.;
      }
    else
      {
        sigma_y_2_ = diameter * diameter / 4.;
        sigma_x_2_ = 0.;
        sigma_z_2_ = 0.;
        source_height_ = height;
      }
  }

  //! Legacy constructor for GaussianPuff_aer.
  /*!
    \param time_puff time at which the puff is released (seconds).
    \param quantity total mass released by the source (mass unit).
    \param source_abscissa abscissa (meters).
    \param source_ordinate ordinate (meters).
    \param source_height height (meters).
    \param species_index index associated to the species.
  */
  template<class T>
  Puff<T>::Puff(T release_time, T velocity, T temperature, T diameter,
                T quantity, T abscissa, T ordinate, T height,
                int species_index, string source_id):
    release_time_(release_time), quantity_(quantity),
    source_velocity_(velocity), source_temperature_(temperature),
    source_diameter_(diameter), source_abscissa_(abscissa),
    source_ordinate_(ordinate), source_height_(height),
    source_id_(source_id),
    z_above_(0.), penetration_factor_(0.),
    sigma_x_(0.), sigma_y_(0.), sigma_z_(0.),
    sigma_y_2_(diameter * diameter / 4.), sigma_z_2_(0.),
    Ns(1), species_index_(species_index),
    scavenging_coefficient_(0.), deposition_velocity_(0.),
    reflection_factor_(1.), has_meteo(0)
  {
  }

  //! Destructor.
  template<class T>
  Puff<T>::~Puff()
  {
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Initializes the puff.
  /*!
    Sets the puff center coordinates equal to those of the release point.
  */
  template<class T>
  void Puff<T>::InitPuff()
  {
    x_ = source_abscissa_;
    y_ = source_ordinate_;
    z_ = source_height_;
    distance_ = 0.;
    puff_time_ = 0.;
    has_meteo = false;
    is_volume_source_ = false;
    sigma_x_ = 0.;
    sigma_y_ = 0.;
    sigma_z_ = 0.;
  }


  ////////////////////////////////////////
  // ACCESS METHODS FOR PUFF ATTRIBUTES //
  ///////////////////////////////////////


  //! Returns the time at which the puff is released.
  /*!
    \return The time at which the puff is released (seconds).
  */
  template<class T>
  inline T Puff<T>::GetReleaseTime() const
  {
    return release_time_;
  }


  //! Returns the puff source id
  template<class T>
  inline string Puff<T>::GetSourceId() const
  {
    return source_id_;
  }

    //! Returns the puff source id
  template<class T>
  inline T Puff<T>::GetPuffVolume() const
  {
    return volume_prev_;
  }

  //! Returns the puff emitted water quantity.
  /*!
    \return The water quantity (mass unit).
  */
  template<class T>
  inline T Puff<T>::GetEmittedWater() const
  {
    return source_water_;
  }

  //! Returns the efflux speed.
  /*!
    \return The efflux speed.
  */
  template<class T>
  inline T Puff<T>::GetSourceVelocity() const
  {
    return source_velocity_;
  }


  //! Returns the temperature of emissions.
  /*!
    \return The temperature of emissions.
  */
  template<class T>
  inline T Puff<T>::GetSourceTemperature() const
  {
    return source_temperature_;
  }


  //! Returns the source diameter.
  /*!
    \return The source diameter.
  */
  template<class T>
  inline T Puff<T>::GetSourceDiameter() const
  {
    return source_diameter_;
  }

    //! Returns the source height.
  /*!
    \return The source height.
  */
  template<class T>
  inline T Puff<T>::GetSourceHeight() const
  {
    return source_height_;
  }


  //! Returns the time since puff release.
  /*!
    \return The time since puff was released (seconds).
  */
  template<class T>
  T Puff<T>::GetPuffTime() const
  {
    return puff_time_;
  }


  //! Returns the species index.
  /*!
    \return The species index.
  */
  template<class T>
  int Puff<T>::GetSpeciesIndex() const
  {
    return species_index_;
  }


  //! Returns the number of species.
  /*!
    \return The number of species.
  */
  template<class T>
  int Puff<T>::GetNs()
  {
    return 1;
  }

  template<class T>
  int Puff<T>::GetNs_aer()
  {
    return 1;
  }


  //! Returns the puff distance from release point.
  /*!
    \return The puff distance from release point (meters).
  */
  template<class T>
  T Puff<T>::GetDistance() const
  {
    return distance_;
  }


  //! Returns the puff quantity.
  /*!
    \return The puff quantity (mass unit).
  */
  template<class T>
  inline T Puff<T>::GetQuantity(int s) const
  {
    if (this->species_index_ == s)
      return quantity_;
    else
      throw string("Species ") + to_str<int>(s) + " not in current puff.";
  }

  template<class T>
  inline T Puff<T>::GetQuantity(int s, int b) const
  {
    if (this->species_index_ == s)
      return quantity_;
    else
      throw string("Aerosol species ") + to_str<int>(s) + " not in current puff.";
  }

  template<class T>
  inline T Puff<T>::GetNumberQuantity(int b) const
  {
    if (this->species_index_ == b)
      return quantity_;
    else
      throw string("Aerosol species ") + to_str<int>(b) + " not in current puff.";
  }



  //! Returns the puff current abscissa.
  /*!
    \return The puff current abscissa (meters).
  */
  template<class T>
  inline T Puff<T>::GetX() const
  {
    return x_;
  }


  //! Returns the puff current ordinate.
  /*!
    \return The puff current ordinate (meters).
  */
  template<class T>
  inline T Puff<T>::GetY() const
  {
    return y_;
  }


  //! Returns the puff current altitude.
  /*!
    \return The puff current altitude (meters).
  */
  template<class T>
  inline T Puff<T>::GetZ() const
  {
    return z_;
  }


  //! Returns the effective height of the plume part above BL.
  /*!
    \return The plume effective height.
  */
  template<class T>
  inline T Puff<T>::GetHeightAboveBL() const
  {
    return z_above_;
  }


  //! Returns the total height of the puff.
  /*!
    \return The puff height.
  */
  template<class T>
  inline T Puff<T>::GetHeight() const
  {
    if (z_ != 0.)
      return z_;
    else
      return z_above_;
  }


  //! Gets the fraction of the puff above the boundary layer.
  /*!
    \return the penetration factor.
  */
  template<class T>
  T Puff<T>::GetPenetrationFactor() const
  {
    return penetration_factor_;
  }


  //! Gets the fraction of the puff reflected by the ground.
  /*!
    \return the reflection factor.
  */
  template<class T>
  T Puff<T>::GetReflectionFactor(int s) const
  {
    if (this->species_index_ == s)
      return reflection_factor_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }

  //! Gets the fraction of the puff reflected by the ground.
  /*!
    \return the reflection factor.
  */
  template<class T>
  T Puff<T>::GetReflectionFactor(int s, int b) const
  {
    if (this->species_index_ == s)
      return reflection_factor_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }


  //! Sets the puff height.
  /*!
    \param z the puff height.
  */
  template<class T>
  void Puff<T>::SetHeight(T z)
  {
    z_ = z;
  }

  //! Sets the puff volume.
  /*!
    \param volume_prev the volume of the puff at previous timestep.
  */
  template<class T>
  void Puff<T>::SetPuffVolume(T volume_prev)
  {
    volume_prev_ = volume_prev;
  }



  //! Sets the puff height above the BL.
  /*!
    \param z the puff height above the BL.
  */
  template<class T>
  void Puff<T>::SetHeightAboveBL(T z_above)
  {
    z_above_ = z_above;
  }


  //! Sets the fraction of the puff above the boundary layer.
  /*!
    \param P the penetration factor.
  */
  template<class T>
  void Puff<T>::SetPenetrationFactor(T P)
  {
    penetration_factor_ = P;
  }


  //! Sets the fraction of the puff reflected by the ground.
  /*!
    \param alpha reflection factor.
  */
  template<class T>
  void Puff<T>::SetReflectionFactor(T alpha, int s)
  {
    if (this->species_index_ == s)
      reflection_factor_ = alpha;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }

  //! Sets the fraction of the puff reflected by the ground.
  /*!
    \param alpha reflection factor.
  */
  template<class T>
  void Puff<T>::SetReflectionFactor(T alpha, int s, int b)
  {
    if (this->species_index_ == s)
      reflection_factor_ = alpha;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }


  //! Returns the horizontal downwind standard deviation.
  /*!
    \return Horizontal standard deviation downwind in meters.
  */
  template<class T>
  inline T Puff<T>::GetSigma_x() const
  {
    return sigma_x_;
  }


  //! Returns the horizontal crosswind standard deviation.
  /*!
    \return Horizontal standard deviation crosswind in meters.
  */
  template<class T>
  inline T Puff<T>::GetSigma_y() const
  {
    return sigma_y_;
  }


  //! Returns the vertical standard deviation.
  /*!
    \return Vertical standard deviation in meters.
  */
  template<class T>
  inline T Puff<T>::GetSigma_z() const
  {
    return sigma_z_;
  }


  //! Returns the square of the puff initial horizontal spread.
  /*!
    \return The square of the puff initial horizontal spread.
  */
  template<class T>
  inline T Puff<T>::GetInitialSigma_y() const
  {
    return sigma_y_2_;
  }

  template<class T>
  inline T Puff<T>::GetInitialSigma_x() const
  {
    return sigma_x_2_;
  }

  //! Returns the square of the puff initial vertical spread.
  /*!
    \return The square of the puff initial vertical spread.
  */
  template<class T>
  inline T Puff<T>::GetInitialSigma_z() const
  {
    return sigma_z_2_;
  }



  //! Sets the puff quantity.
  template<class T>
  void Puff<T>::SetQuantity(T quantity, int s)
  {

    if (this->species_index_ == s)
      quantity_ = quantity;
    else
      throw string("Species ") + to_str<int>(s) + " not in current puff.";
  }

  //! Sets the puff aerosol quantity.
  template<class T>
  void Puff<T>::SetQuantity(T quantity, int s, int b)
  {

    throw string("\"Puff<T>::SetQuantity(T&, int&, int&)\"")
      + " is not defined.";
  
  }

  //! Sets the puff aerosol quantity.
  template<class T>
  void Puff<T>::SetNumberQuantity(T quantity, int b)
  {

    throw string("\"Puff<T>::SetQuantity(T&, int&, int&)\"")
      + " is not defined.";
  
  }

  //! Sets the horizontal standard deviation in the downwind direction.
  /*!  \param sigma_x horizontal standard deviation in the downwind direction
    (meters).
  */
  template<class T>
  void Puff<T>::SetSigma_x(T sigma_x)
  {
    if (sigma_x > sigma_x_)
      sigma_x_ = sigma_x;
  }


  //! Sets the standard deviation in the crosswind direction.
  /*!  \param sigma_y horizontal standard deviation in the crosswind direction
    (meters).
  */
  template<class T>
  void Puff<T>::SetSigma_y(T sigma_y)
  {
    if (sigma_y > sigma_y_)
      sigma_y_ = sigma_y;
  }


  //! Sets the vertical standard deviation.
  /*!
    \param sigma_z vertical standard deviation (meters).
  */
  template<class T>
  void Puff<T>::SetSigma_z(T sigma_z)
  {
    if (sigma_z > sigma_z_)
      sigma_z_ = sigma_z;
  }

  //! Sets the puff initial horizontal spread.
  /*!
    \param sigma_y_2 The puff initial horizontal spread.
  */
  template<class T>
  void Puff<T>::SetInitialSigma_y(T sigma_y_2)
  {
    sigma_y_2_ = sigma_y_2;
  }

  //! Sets the puff initial horizontal spread.
  /*!
    \param sigma_x_2 The puff initial horizontal spread.
  */
  template<class T>
  void Puff<T>::SetInitialSigma_x(T sigma_x_2)
  {
    sigma_x_2_ = sigma_x_2;
  }


  //! Sets the puff initial vertical spread.
  /*!
    \param sigma_z_2 The puff initial vertical spread.
  */
  template<class T>
  void Puff<T>::SetInitialSigma_z(T sigma_z_2)
  {
    sigma_z_2_ = sigma_z_2;
  }


  //! Returns true if meteo has been set, false otherwise.
  /*!
    \return true if meteo has been set, false otherwise.
  */
  template<class T>
  bool Puff<T>::HasMeteo() const
  {
    return has_meteo;
  }


  //! Gets values of meteorological parameters.
  /*! It gets the meteorological data.
    \param temperature: the temperature of ambient air (Kelvin degrees).
    \param wind_angle: the angle between positive x-axis and wind,
    counterclockwise (radians).
    \param wind: the wind velocity (m/s).
    \param stability: stability class in [0, 5].
  */
  template<class T>
  void Puff<T>::GetMeteo(T& temperature, T& wind,
                         T& cos_angle, T& sin_angle,
                         int& stability, T& boundary_height,
                         T& longitude, T& latitude,
                         bool& isday, bool& rural)
  {
    temperature = temperature_;
    cos_angle = cos_angle_;
    sin_angle = sin_angle_;
    wind = wind_;
    boundary_height = boundary_height_;
    stability = stability_;
    longitude = longitude_;
    latitude = latitude_;
    isday = isday_;
    rural = rural_;
  }


  //! Gets values of meteorological parameters.
  /*! It gets the meteorological data.
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
  void Puff<T>::GetMeteo(T& temperature, T& wind,
                         T& cos_angle, T& sin_angle,
                         T& friction_velocity, T& convective_velocity,
                         T& boundary_height, T& lmo, T& coriolis,
                         T& longitude, T& latitude,
                         bool& isday, bool& rural)
  {
    temperature = temperature_;
    cos_angle = cos_angle_;
    sin_angle = sin_angle_;
    wind = wind_;
    friction_velocity = friction_velocity_;
    convective_velocity = convective_velocity_;
    boundary_height = boundary_height_;
    lmo = lmo_;
    coriolis = coriolis_;
    longitude = longitude_;
    latitude = latitude_;
    isday = isday_;
    rural = rural_;
  }


  //! Gets values of additional meteorological data.
  /*! It gets the additionnal meteorological data (needed for chemistry).
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void Puff<T>::GetAdditionalMeteo(T& attenuation, T& pressure,
                                   T& specific_humidity) const
  {
    throw string("\"Puff<T>::GetAdditionalMeteo(T&, T&, T&)\"")
      + " is not defined.";
  }

  template<class T>
  void Puff<T>::GetAdditionalMeteo(T& liquid_water_content) const
  {
    throw string("\"Puff<T>::GetAdditionalMeteo(T&)\"")
      + " is not defined.";
  }


  //! Gets scavenging coefficient.
  template<class T>
  T Puff<T>::GetScavengingCoefficient(int s) const
  {
    if (this->species_index_ == s)
      return scavenging_coefficient_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }

  //! Gets scavenging coefficient.
  template<class T>
  T Puff<T>::GetScavengingCoefficient(int s, int b) const
  {
    if (this->species_index_ == s)
      return scavenging_coefficient_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }


  //! Gets deposition velocity.
  template<class T>
  T Puff<T>::GetDepositionVelocity(int s) const
  {
    if (this->species_index_ == s)
      return deposition_velocity_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }

  //! Gets deposition velocity.
  template<class T>
  T Puff<T>::GetDepositionVelocity(int s, int b) const
  {
    if (this->species_index_ == s)
      return deposition_velocity_;
    else
      throw string("Species index ") + to_str(s) + " not in puff.";
  }


  //! Gets the species photolysis rate for a given species.
  template<class T>
  T Puff<T>::GetPhotolysisRate(int s) const
  {
    throw string("\"Puff<T>::GetPhotolysisRate(int)\"")
      + " is not defined.";
  }


  //! Gets the species background concentration for a given species.
  template<class T>
  T Puff<T>::GetBackgroundConcentration(int s) const
  {
    throw string("\"Puff<T>::GetBackgroundConcentration(int)\"")
      + " is not defined.";
  }

  //! Gets the species background concentration for a given species.
  template<class T>
  T Puff<T>::GetPreviousBackgroundConcentration(int s) const
  {
    throw string("\"Puff<T>::GetBackgroundConcentration(int)\"")
      + " is not defined.";
  }

  //! Gets the aerosol species background concentration for a given species.
  template<class T>
  T Puff<T>::GetBackgroundConcentration(int s, int b) const
  {
    throw string("\"Puff<T>::GetBackgroundConcentration(int, int)\"")
      + " is not defined.";
  }

  //! Gets the aerosol species background concentration for a given species.
  template<class T>
  T Puff<T>::GetPreviousBackgroundConcentration(int s, int b) const
  {
    throw string("\"Puff<T>::GetBackgroundConcentration(int, int)\"")
      + " is not defined.";
  }

  //! Gets the aerosol species background concentration for a given species.
  template<class T>
  T Puff<T>::GetBackgroundNumberConcentration(int b) const
  {
    throw string("\"Puff<T>::GetBackgroundNumberConcentration(int, int)\"")
      + " is not defined.";
  }

  //! Gets the aerosol species background concentration for a given species.
  template<class T>
  T Puff<T>::GetPreviousBackgroundNumberConcentration(int b) const
  {
    throw string("\"Puff<T>::GetBackgroundNumberConcentration(int, int)\"")
      + " is not defined.";
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
  void Puff<T>::SetMeteo(T temperature, T wind_angle, T wind,
                         string stability_class, T longitude, T latitude,
                         bool isday, bool rural, T boundary_height)
  {
    temperature_ = temperature;
    wind_angle_ = wind_angle;
    cos_angle_ = cos(wind_angle_);
    sin_angle_ = sin(wind_angle_);
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
    longitude_ = longitude;
    latitude_ = latitude;
    has_meteo = 1;
    isday_ = isday;
    rural_ = rural;
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
  void Puff<T>::SetMeteo(T temperature, T wind_angle, T wind,
                         T friction_velocity, T convective_velocity,
                         T boundary_height, T lmo, T coriolis,
                         T longitude, T latitude, bool isday, bool rural)
  {
    temperature_ = temperature;
    wind_angle_ = wind_angle;
    cos_angle_ = cos(wind_angle_);
    sin_angle_ = sin(wind_angle_);
    wind_ = wind;
    boundary_height_ = boundary_height;
    friction_velocity_ = friction_velocity;
    convective_velocity_ = convective_velocity;
    boundary_height_ = boundary_height;
    lmo_ = lmo;
    coriolis_ = coriolis;
    longitude_ = longitude;
    latitude_ = latitude;
    has_meteo = 1;
    isday_ = isday;
    rural_ = rural;
  }

  //! Sets values of additional meteorological data.
  /*! It sets the additionnal meteorological data (needed for chemistry).
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void Puff<T>::SetAdditionalMeteo(T attenuation, T pressure,
                                   T specific_humidity)
  {
    throw string("\"Puff<T>::SetAdditionalMeteo(T, T, T)\"")
      + " is not defined.";
  }

  template<class T>
  void Puff<T>::SetAdditionalMeteo(T liquid_water_content)
  {
    throw string("\"Puff<T>::SetAdditionalMeteo(T)\"")
      + " is not defined.";
  }


  //! Sets the species scavenging coefficient.
  /*!
    \param lambda the species scavenging coefficient ( s^(-1) ).
  */
  template<class T>
  void Puff<T>::SetScavengingCoefficient(T lambda, int s)
  {
    if (this->species_index_ == s)
      scavenging_coefficient_ = lambda;
  }


  // //! Returns the puff LWC.
  // /*!
  //   \return the LWC.
  // */
  // template<class T>
  // T Puff<T>::GetLiquidWaterContent() const
  // {
  //     throw string("\"Puff<T>::GetLiquidWaterContent(T)\"")
  //     + " is not defined.";;
  // }

  //! Returns the puff specific humidity.
  /*!
    \return the specific humidity.
  */
  template<class T>
  T Puff<T>::GetSpecificHumidity() const
  {
      throw string("\"Puff<T>::GetSpecificHumidity(T)\"")
      + " is not defined.";;
  }


  //! Sets the species deposition velocity.
  /*!
    \param lambda the species deposition velocity (m/s).
  */
  template<class T>
  void Puff<T>::SetDepositionVelocity(T vd, int s)
  {
    if (this->species_index_ == s)
      deposition_velocity_ = vd;
  }


  //! Sets the species photolysis rate for a given species.
  /*!
    \param rate the species photolysis rate.
  */
  template<class T>
  void Puff<T>::SetPhotolysisRate(T rate, int s)
  {
    throw string("\"Puff<T>::SetPhotolysisRate(T, int)\"")
      + " is not defined.";
  }


  //! Sets the species background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void Puff<T>::SetBackgroundConcentration(T concentration, int s)
  {
    throw string("\"Puff<T>::SetBackgroundConcentration(T, int)\"")
      + " is not defined.";
  }

  //! Sets the species background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void Puff<T>::SetPreviousBackgroundConcentration(T concentration, int s)
  {
    throw string("\"Puff<T>::SetBackgroundConcentration(T, int)\"")
      + " is not defined.";
  }

  //! Sets the aerosol species background concentration for a given species.
  /*!
    \param concentration the aerosol species background concentration.
  */
  template<class T>
  void Puff<T>::SetBackgroundConcentration(T concentration, int s, int b)
  {
    throw string("\"Puff<T>::SetBackgroundConcentration(T, int, int)\"")
      + " is not defined.";
  }

  //! Sets the aerosol species background concentration for a given species.
  /*!
    \param concentration the aerosol species background concentration.
  */
  template<class T>
  void Puff<T>::SetPreviousBackgroundConcentration(T concentration, int s, int b)
  {
    throw string("\"Puff<T>::SetBackgroundConcentration(T, int, int)\"")
      + " is not defined.";
  }

  template<class T>
  void Puff<T>::SetBackgroundNumberConcentration(T concentration, int b)
  {
    throw string("\"Puff<T>::SetBackgroundNumberConcentration(T, int, int)\"")
      + " is not defined.";
  }

  template<class T>
  void Puff<T>::SetPreviousBackgroundNumberConcentration(T concentration, int b)
  {
    throw string("\"Puff<T>::SetBackgroundNumberConcentration(T, int, int)\"")
      + " is not defined.";
  }

  //! Sets the new position of the puff center.
  /*!
    \param x new abscissa (meters).
    \param y new ordinate (meters).
    \param z new height (meters).
    \param distance new distance (meters).
    \param time new time since release (s).
  */  template<class T>
  void Puff<T>::SetPuffPosition(T x, T y, T z, T distance, T time)
  {
    x_ = x;
    y_ = y;
    z_ = z;
    distance_ = distance;
    puff_time_ = time;
  }

  //! Returns true if the puff is from a volume source, false otherwise.
  /*!
    \return true if the puff is from a volume source, false otherwise.
  */
  template<class T>
  bool Puff<T>::IsVolumeSource() const
  {
    return is_volume_source_;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFTRANSPORT_CXX
#endif
