// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFEMISSION_CXX


#include "PuffEmission.hxx"
#include "PointEmissionUnit.cxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  PuffEmission<T>::PuffEmission()
  {
  }


  //! Model initialization.
  /*! It reads the configuration. Each source is described in a
    dedicated section "[source]" in which one finds the following entries:
    <ul>
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species: list of species names,
    <li> Date beg: the release date,
    <li> Quantity: the list of released quantitys (mass unit),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees),
    </ul>
    \param config The ConfigStream instance containing the sources parameters.
  */
  template<class T>
  void PuffEmission<T>::Init(ConfigStream& config,
                             vector<string> species_list)
  {
    config.PeekValue("Abscissa", x);
    config.PeekValue("Ordinate", y);
    config.PeekValue("Altitude", z);
    config.PeekValue("Velocity", this->velocity_);
    config.PeekValue("Temperature", this->temperature_);
    config.PeekValue("Diameter", "positive", this->diameter_);
    this->temperature_ += 273.15;

    PointEmissionUnit<T>::Init(config, species_list);
    date_end = this->date_beg;
    config.Find("Quantity");
    vector<string> tmp = split(config.GetLine());
    int Nquantity = tmp.size();
    if (Nquantity != this->Ns_emis)
      throw string("Error in \"PuffEmission<T>::Init()\": ")
        + " there must be one quantity per emitted species.";
    for (int i = 0; i < Nquantity; i++)
      quantity.push_back(to_num<T>(tmp[i]));
  }


  //! Gets the emission coordinates.
  /*!
    Gets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
  */
  template<class T>
  void PuffEmission<T>::GetCoordinates(T& abscissa, T& ordinate,
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
  void PuffEmission<T>::SetCoordinates(T abscissa, T ordinate, T height)
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
  bool PuffEmission<T>::HasEnded(Date date)
  {
    return date > date_end;
  }


  //! Gets emitted quantities for one species.
  /*!
    \return The emitted quantities for one species.
  */
  template<class T>
  void PuffEmission<T>::GetEmission(Date current_date,
                                    Date next_date, int s,
                                    Array<T, 2>& point_emission)
  {
    GetEmission(s, point_emission);
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
  void PuffEmission<T>::GetEmission(int s, Array<T, 2>& point_emission)
  {
    point_emission.resize(1, 4);
    point_emission(0, 0) = x;
    point_emission(0, 1) = y;
    point_emission(0, 2) = z;
    point_emission(0, 3) = quantity[s];
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFEMISSION_CXX
#endif
