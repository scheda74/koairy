// Copyright (C) 2006-2012, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMESOURCE_CXX


#include "PlumeSource.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*!
    \param rate source emission rate.
    \param velocity efflux speed of gases.
    \param temperature temperature of emissions.
    \param diameter source diameter.
    \param x abscissa.
    \param y ordinate.
    \param z height.
    \param species_index index associated to the species.
    \param id_source index of the source.
  */
  template<class T>
  PlumeSource<T>::PlumeSource(T rate, T velocity, T temperature, T diameter,
                              T x, T y, T z, int species_index, int id_source):
    rate_(rate), velocity_(velocity), temperature_(temperature),
    diameter_(diameter), x_(x), y_(y), z_(z), effective_height_(z_),
    effective_height_above_(0.), penetration_factor_(0.),
    sigma_y_2_(diameter * diameter / 4.), sigma_z_2_(0.),
    species_index_(species_index), type_("continuous"), id_source_(id_source)
  {
  }


  //! Destructor.
  template<class T>
  PlumeSource<T>::~PlumeSource()
  {
  }


  //! Returns the source rate.
  /*!
    \return The source rate.
  */
  template<class T>
  inline T PlumeSource<T>::GetRate() const
  {
    return rate_;
  }


  //! Sets the source rate.
  /*!
    \param rate The source rate.
  */
  template<class T>
  void PlumeSource<T>::SetRate(T rate)
  {
    rate_ = rate;
  }


  //! Sets the fraction of the plume above the boundary layer.
  /*!
    \param P the penetration factor.
  */
  template<class T>
  void PlumeSource<T>::SetPenetrationFactor(T P)
  {
    penetration_factor_ = P;
  }


  //! Gets the fraction of the plume above the boundary layer.
  /*!
    \return the penetration factor.
  */
  template<class T>
  T PlumeSource<T>::GetPenetrationFactor() const
  {
    return penetration_factor_;
  }


  //! Returns the efflux speed.
  /*!
    \return The efflux speed.
  */
  template<class T>
  inline T PlumeSource<T>::GetVelocity() const
  {
    return velocity_;
  }


  //! Sets the efflux speed.
  /*!
    \param velocity the new efflux speed.
  */
  template<class T>
  inline void PlumeSource<T>::SetVelocity(T velocity)
  {
    velocity_ = velocity;
  }


  //! Returns the temperature of emissions.
  /*!
    \return The temperature of emissions.
  */
  template<class T>
  inline T PlumeSource<T>::GetTemperature() const
  {
    return temperature_;
  }


  //! Sets the temperature of emissions.
  /*!
    \param temperature the new temperature of emissions.
  */
  template<class T>
  inline void PlumeSource<T>::SetTemperature(T temperature)
  {
    temperature_ = temperature;
  }


  //! Returns the source diameter.
  /*!
    \return The source diameter.
  */
  template<class T>
  inline T PlumeSource<T>::GetDiameter() const
  {
    return diameter_;
  }


  //! Sets the source diameter.
  /*!
    \param diameter the new source diameter.
  */
  template<class T>
  inline void PlumeSource<T>::SetDiameter(T diameter)
  {
    diameter_ = diameter;
  }


  //! Returns the source position along x.
  /*!
    \return The source abscissa.
  */
  template<class T>
  inline T PlumeSource<T>::GetX() const
  {
    return x_;
  }


  //! Sets the source position along x.
  /*!
    \param x the new source abscissa.
  */
  template<class T>
  inline void PlumeSource<T>::SetX(T x)
  {
    x_ = x;
  }


  //! Returns the source position along y.
  /*!
    \return The source ordinate.
  */
  template<class T>
  inline T PlumeSource<T>::GetY() const
  {
    return y_;
  }


  //! Sets the source position along y.
  /*!
    \param the new source ordinate.
  */
  template<class T>
  inline void PlumeSource<T>::SetY(T y)
  {
    y_ = y;
  }


  //! Returns the source height.
  /*!
    \return The source height.
  */
  template<class T>
  inline T PlumeSource<T>::GetZ() const
  {
    return z_;
  }


  //! Sets the source height.
  /*!
    \param z the new source height.
  */
  template<class T>
  inline void PlumeSource<T>::SetZ(T z)
  {
    z_ = z;
  }


  //! Returns the plume effective height (with plume rise).
  /*!
    \return The plume effective height.
  */
  template<class T>
  inline T PlumeSource<T>::GetHeight() const
  {
    return effective_height_;
  }


  //! Returns the effective height of the plume part above BL.
  /*!
    \return The plume effective height.
  */
  template<class T>
  inline T PlumeSource<T>::GetHeightAboveBL() const
  {
    return effective_height_above_;
  }


  //! Sets the source height.
  /*!
    \param rate The source height.
  */
  template<class T>
  void PlumeSource<T>::SetHeight(T height)
  {
    effective_height_ = height;
  }


  //! Sets the plume height above the BL.
  /*!
    \param effective_height_above the plume height above the BL.
  */
  template<class T>
  void PlumeSource<T>::SetHeightAboveBL(T effective_height_above)
  {
    effective_height_above_ = effective_height_above;
  }


  //! Returns the plume initial horizontal spread.
  /*!
    \return The plume initial horizontal spread.
  */
  template<class T>
  inline T PlumeSource<T>::GetSigma_y() const
  {
    return sqrt(sigma_y_2_);
  }


  //! Returns the plume initial vertical spread.
  /*!
    \return The plume initial vertical spread.
  */
  template<class T>
  inline T PlumeSource<T>::GetSigma_z() const
  {
    return sqrt(sigma_z_2_);
  }


  //! Sets the square of the plume initial horizontal spread.
  /*!
    \param sigma_y_2 The square of the plume initial horizontal spread.
  */
  template<class T>
  void PlumeSource<T>::SetSigma_y_2(T sigma_y_2)
  {
    sigma_y_2_ = sigma_y_2;
  }


  //! Sets the square of the plume initial vertical spread.
  /*!
    \param sigma_z_2 The square of the plume initial vertical spread.
  */
  template<class T>
  void PlumeSource<T>::SetSigma_z_2(T sigma_z_2)
  {
    sigma_z_2_ = sigma_z_2;
  }


  //! Returns the square of the plume initial horizontal spread.
  /*!
    \return The square of the plume initial horizontal spread.
  */
  template<class T>
  inline T PlumeSource<T>::GetSigma_y_2() const
  {
    return sigma_y_2_;
  }


  //! Returns the square of the plume initial vertical spread.
  /*!
    \return The square of the plume initial vertical spread.
  */
  template<class T>
  inline T PlumeSource<T>::GetSigma_z_2() const
  {
    return sigma_z_2_;
  }


  //! Returns the species index.
  /*!
    \return The species index.
  */
  template<class T>
  inline int PlumeSource<T>::GetSpeciesIndex() const
  {
    return species_index_;
  }


  //! Returns the emission type.
  /*!
    \return The emission type.
  */
  template<class T>
  inline string PlumeSource<T>::GetEmissionType() const
  {
    return type_;
  }


  //! Returns the source second extremity along x.
  /*!
    \return The source second extremity abscissa.
  */
  template<class T>
  inline T PlumeSource<T>::GetX2() const
  {
    throw "'PlumeSource<T>::GetX2() const' is undefined.";
  }


  //! Returns the source second extremity along y.
  /*!
    \return The source second extremity ordinate.
  */
  template<class T>
  inline T PlumeSource<T>::GetY2() const
  {
    throw "'PlumeSource<T>::GetY2() const' is undefined.";
  }


  //! Returns the source second extremity along z.
  /*!
    \return The source second extremity height.
  */
  template<class T>
  inline T PlumeSource<T>::GetZ2() const
  {
    throw "'PlumeSource<T>::GetZ2() const' is undefined.";
  }


  //! Returns the source second extremity along x.
  /*!
    \return The source second extremity abscissa.
  */
  template<class T>
  inline void PlumeSource<T>::SetX2(T x_2)
  {
    throw "'PlumeSource<T>::SetX2(T) const' is undefined.";
  }


  //! Returns the source second extremity along y.
  /*!
    \return The source second extremity ordinate.
  */
  template<class T>
  inline void PlumeSource<T>::SetY2(T y_2)
  {
    throw "'PlumeSource<T>::SetY2(T) const' is undefined.";
  }


  //! Returns the source second extremity along z.
  /*!
    \return The source second extremity height.
  */
  template<class T>
  inline void PlumeSource<T>::SetZ2(T z_2)
  {
    throw "'PlumeSource<T>::SetZ2(T) const' is undefined.";
  }


  //! Returns the combination factor.
  /*!
    \return The combination factor.
  */
  template<class T>
  inline T PlumeSource<T>::GetCombinationFactor() const
  {
    throw "'PlumeSource<T>::GetCombinationFactor() const' is undefined.";
  }


  //! Sets the combination factor.
  /*!
    \param combination_factor the combination factor.
  */
  template<class T>
  void PlumeSource<T>::SetCombinationFactor(T combination_factor)
  {
    throw "'PlumeSource<T>::SetCombinationFactor() const' is undefined.";
  }


  //! Returns the line source width.
  /*!
    \return The line source width.
  */
  template<class T>
  inline T PlumeSource<T>::GetWidth() const
  {
    throw "'PlumeSource<T>::GetWidth()' is undefined.";
  }


  //! Returns the source identifier.
  /*!
    \return the source identifier
  */
  template<class T>
  inline int PlumeSource<T>::GetIdSource() const
  {
    return id_source_;
  }


  //! Returns the source section identifier.
  /*!
    \return the source section identifier
  */
  template<class T>
  int PlumeSource<T>::GetIdSection() const
  {
    throw "'PlumeSource<T>::GetIdSection() const' is undefined.";
  }


  //! Returns the line source VehicleVelocity.
  /*!
    \return The line source Vehiclevelocity.
  */
  template<class T>
  inline T PlumeSource<T>::GetVehicleVelocity() const
  {
    throw "'PlumeSource<T>::GetVehicleVelocity() const' is undefined.";
  }


  //! Returns the line source Area.
  /*!
    \return The line source Area.
  */
  template<class T>
  inline T PlumeSource<T>::GetArea() const
  {
    throw "'PlumeSource<T>::GetArea() const' is undefined.";
  }


  //! Returns the line source density.
  /*!
    \return The line source density.
  */
  template<class T>
  inline T PlumeSource<T>::GetDensity() const
  {
    throw "'PlumeSource<T>::GetDensity() const' is undefined.";
  }
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PLUMESOURCE_CXX
#endif
