// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMELINESOURCE_CXX


#include "PlumeLineSource.hxx"


namespace Polyphemus
{


  //! Computes the initial vertical spread of the line source (m).
  /*!
    \param z height of the source (m).
    \return The initial vertical spread of the line source.
  */
  template<class T>
  T  PlumeLineSource<T>::vertical_spread(T z)
  {
    return 0.7 * z;
  }


  //! Computes the initial vertical spread of the line source (m).
  /*!
    \param z height of the source (m).
    \param width width of the source (m).
    \param area area of the line source (m²).
    \param density density of the line source (veh/m).
    \param vehicle_velocity vehicle velocity of the line source (m/s).
    \return The initial vertical spread of the line source.
  */
  template<class T>
  T  PlumeLineSource<T>::vertical_spread(T z, T width, T area, T density,
                                         T vehicle_velocity)
  {
    if (width > 0.)
      return 0.3 * sqrt(density * pow(vehicle_velocity, 2) * area / width);
    else
      return vertical_spread(z);
  }


  //! Main constructor.
  /*!
    \param rate source emission rate.
    \param x abscissa of the source first extremities.
    \param y ordinate of the source first extremities.
    \param z height of the source first extremities.
    \param x2 abscissa of the source second extemities.
    \param y2 ordinate of the source second extremities.
    \param z2 height of the source second extremities.
    \param width width of the source.
    \param species_index index associated to the species.
    \param id1 index of the source.
    \param id2 index of the section.
  */
  template<class T>
  PlumeLineSource<T>::PlumeLineSource(T rate, T x, T y, T z, T x2, T y2, T z2,
                                      T width, T VehicleVelocity, T Area,
                                      T Density, int species_index,
                                      int id_source, int id_section):
    PlumeSource<T>(rate, 0., 0., 0., x, y, z, species_index, id_source),
    x2_(x2), y2_(y2), z2_(z2), combination_factor_(1), width_(width),
    VehicleVelocity_(VehicleVelocity), Area_(Area), Density_(Density),
    id_section_(id_section)
  {
    if (z != z2)
      throw "Both extremities of each road sections must be at the same "
        "height.";
    this->sigma_z_2_ = pow(vertical_spread(z, width, Area, Density,
                                           VehicleVelocity), 2);
    this->type_ = "continuous_line";
  }


  //! Destructor.
  template<class T>
  PlumeLineSource<T>::~PlumeLineSource()
  {
    // Nothing.
  }


  //! Returns the source second extremity along x.
  /*!
    \return The source second extremity abscissa.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetX2() const
  {
    return x2_;
  }


  //! Sets the source second extremity along x.
  /*!
    \param x the source second extremity abscissa.
  */
  template<class T>
  inline void PlumeLineSource<T>::SetX2(T x2)
  {
    x2_ = x2;
  }


  //! Returns the source second extremity along y.
  /*!
    \return The source second extremity ordinate.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetY2() const
  {
    return y2_;
  }


  //! Sets the source second extremity along y.
  /*!
    \param x the source second extremity ordinate.
  */
  template<class T>
  inline void PlumeLineSource<T>::SetY2(T y2)
  {
    y2_ = y2;
  }


  //! Returns the source second extremity along z.
  /*!
    \return The source second extremity height.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetZ2() const
  {
    return z2_;
  }


  //! Sets the source second extremity along z.
  /*!
    \param x the source second extremity height.
  */
  template<class T>
  inline void PlumeLineSource<T>::SetZ2(T z2)
  {
    z2_ = z2;
  }


  //! Returns the combination factor.
  /*!
    \return The combination factor.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetCombinationFactor() const
  {
    return combination_factor_;
  }


  //! Sets the combination factor.
  /*!
    \param combination_factor the combination factor.
  */
  template<class T>
  void PlumeLineSource<T>::SetCombinationFactor(T combination_factor)
  {
    combination_factor_ = combination_factor;
  }


  //! Returns width of the line source.
  /*!
    \return The width of the line source.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetWidth() const
  {
    return width_;
  }


  //! Returns the source section identifier.
  /*!
    \return the source section identifier
  */
  template<class T>
  inline int PlumeLineSource<T>::GetIdSection() const
  {
    return id_section_;
  }


  //! Returns VehicleVelocity of the line source (m/s).
  /*!
    \return The VehicleVelocity of the line source.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetVehicleVelocity() const
  {
    return VehicleVelocity_;
  }

  //! Returns Area of the line source (m²).
  /*!
    \return The Area of the line source.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetArea() const
  {
    return Area_;
  }

  //! Returns density of the line source (veh/m).
  /*!
    \return The density of the line source.
  */
  template<class T>
  inline T PlumeLineSource<T>::GetDensity() const
  {
    return Density_;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PLUMELINESOURCE_CXX
#endif
