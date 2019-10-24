// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Ir�ne Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_CXX


#include "ContinuousEmission.hxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  ContinuousEmission<T>::ContinuousEmission()
  {
  }


  //! Destructor.
  template<class T>
  ContinuousEmission<T>::~ContinuousEmission()
  {
    // Nothing.
  }


  //! Continuous emission initialization.
  /*! It reads the configuration. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species: list of species names,
    <li> Date beg: the release date,
    <li> Date end: the release ending date,
    <li> Rate: the list of release rates (mass unit per seconds),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees).
    </ul>
    \param species_list the list of species names.
    \param config the ConfigStream instance containing the sources parameters.
  */
  template<class T>
  void ContinuousEmission<T>::Init(ConfigStream& config,
				   vector<string> species_list)
  {
    config.PeekValue("Abscissa", x);
    config.PeekValue("Ordinate", y);
    config.PeekValue("Altitude", z);
    config.PeekValue("Velocity", this->velocity_);
    config.PeekValue("Temperature", this->temperature_);
    config.PeekValue("IsVolumeSource", this->is_volume_source_);
    if (this->is_volume_source_)
      {
        config.PeekValue("Width", "positive", this->width_);
        config.PeekValue("Length", "positive", this->length_);
        this->diameter_ = sqrt(this->width_ * this->width_ + this->length_ * this->length_);
      }
    else
      config.PeekValue("Diameter", "positive", this->diameter_);

    this->temperature_ += 273.15;

    date_end = config.PeekValue("Date_end");
    PointEmissionUnit<T>::Init(config, species_list);
    config.Find("Rate");
    vector<string> tmp = split(config.GetLine());
    int Nrate = tmp.size();
    if (Nrate != this->Ns_emis)
      throw string("Error in \"ContinuousEmission<T>::Init()\": ")
        + " there must be one rate per emitted species.";
    for (int i = 0; i < Nrate; i++)
      rate.push_back(to_num<T>(tmp[i]));
  }


  //! Continuous emission initialization.
  /*! It reads the configuration. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species: list of species names,
    <li> Date beg: the release date,
    <li> Date end: the release ending date,
    <li> Rate: the list of release rates (mass unit per seconds),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees),
    </ul>
  */
  template<class T>
  void ContinuousEmission<T>::Init(T abscissa_, T ordinate_, T height_,
                                   const vector<T>&  rate_, Date date_beg_,
                                   Date date_end_, const vector<int>& species,
                                   T diameter, T velocity, T temperature)
  {
    this->type = "continuous";
    this->emitted_species_index = species;
    this->Ns_emis = this->emitted_species_index.size();
    this->date_beg = date_beg_;
    this->velocity_ = velocity;
    this->temperature_ = temperature;
    this->diameter_ = diameter;
    x = abscissa_;
    y = ordinate_;
    z = height_;
    rate = rate_;
    date_end = date_end_;
  }


  //! Gets the emission coordinates.
  /*!
    Gets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
  */
  template<class T>
  void ContinuousEmission<T>::GetCoordinates(T& abscissa, T& ordinate,
                                             T& height) const
  {
    abscissa = x;
    ordinate = y;
    height = z;
  }


  //! Sets the emission coordinates.
  /*!
    Sets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
  */
  template<class T>
  void ContinuousEmission<T>::SetCoordinates(T abscissa, T ordinate, T height)
  {
    x = abscissa;
    y = ordinate;
    z = height;
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool ContinuousEmission<T>::HasEnded(Date date)
  {
    return date >= date_end;
  }


  //! Returns the ending date of the emission.
  /*!
    \return The ending date of the emission.
  */
  template<class T>
  Date ContinuousEmission<T>::GetEndDate() const
  {
    return date_end;
  }


  //! Gets the species rate.
  /*!
    \return The species rate.
  */
  template<class T>
  T ContinuousEmission<T>::GetRate(int species) const
  {
    return rate[species];
  }


  //! Multiplies all rates by some factor.
  /*!
    \param[in] factor the factor by which to multiply all emission rates.
  */
  template<class T>
  void ContinuousEmission<T>::MultiplyRate(T factor)
  {
    for (int i = 0; i < int(rate.size()); i++)
      rate[i] *= factor;
  }


  //! Multiplies all rates of a given species by some factor.
  /*!
    \param[in] species index of the species.
    \param[in] factor the factor by which to multiply all emission rates.
  */
  template<class T>
  void ContinuousEmission<T>::MultiplyRate(int species, T factor)
  {
    rate[species] *= factor;
  }


  //! Shifts the emission times by some value.
  /*! The initial and final release times are both shifted.
    \param[in] shift the time shift in seconds (a positive value will postpone
    the release).
  */
  template<class T>
  void ContinuousEmission<T>::ApplyTimeShift(T shift)
  {
    this->date_beg.AddSeconds(shift);
    date_end.AddSeconds(shift);
  }


  //! Shifts the emission altitudes by some value.
  /*!
    \param[in] shift the altitude shift in meters.
  */
  template<class T>
  void ContinuousEmission<T>::ApplyAltitudeShift(T shift)
  {
    z += shift;
  }


  //! Gets emitted quantity for one species during a given time interval.
  /*!
    \param date_beg beginning date of the current time interval
    \param date_end ending date of the current time interval
    \param s species index
    \param point_emission a 2D array of size (1, 4) containing the source
    coordinates (x, y, z) and the species quantity.
  */
  template<class T>
  void ContinuousEmission<T>::GetEmission(Date current_date,
                                          Date next_date, int s,
                                          Array<T, 2>& point_emission)
  {
    point_emission.resize(1, 4);
    point_emission(0, 0) = x;
    point_emission(0, 1) = y;
    point_emission(0, 2) = z;

    T delta_t = next_date.GetSecondsFrom(current_date);
    T quantity = delta_t
      - max(0., this->date_beg.GetSecondsFrom(current_date))
      - max(0., next_date.GetSecondsFrom(date_end));
    quantity *= rate[s];
    point_emission(0, 3) = quantity;
  }


  //! Gets the timespan of the emission.
  /*!
    \return The emission timespan.
  */
  template<class T>
  T ContinuousEmission<T>::GetEmissionTimespan(const Date& current_date,
                                               const Date& next_date) const
  {
    return (next_date.GetSecondsFrom(current_date)
            - max(0., this->date_beg.GetSecondsFrom(current_date))
            - max(0., next_date.GetSecondsFrom(date_end)));
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_CXX
#endif
