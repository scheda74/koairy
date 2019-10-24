// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFCHEMISTRY_CXX


#include "PuffChemistry.hxx"


namespace Polyphemus
{


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
  PuffChemistry<T>::PuffChemistry(T release_time, T velocity,  T temperature,
				  T diameter, vector<T> quantity, T abscissa,
				  T ordinate, T height, T source_water, T volume_prev,
				  vector<string> species_list,
				  vector<string> photolysis_reaction_list, string source_id):
    Puff<T>(release_time, velocity,  temperature, diameter, 
	    0., abscissa, ordinate, height, source_water, volume_prev, 0, source_id),
    quantity_list(quantity), species_list(species_list),
    photolysis_reaction_list(photolysis_reaction_list)
  {
    this->Ns = species_list.size();
    Nr_photolysis = photolysis_reaction_list.size();
  }

  //! Alternative constructor
  template<class T>
  PuffChemistry<T>::PuffChemistry(T release_time, T velocity,  T temperature,
				  T diameter, T width, T length,
                                  vector<T> quantity, T abscissa,
				  T ordinate, T height, T source_water,
				  T volume_prev,
				  vector<string> species_list,
				  vector<string> photolysis_reaction_list,
                                  bool is_volume_source, string source_id):
    Puff<T>(release_time, velocity,  temperature, diameter, width, length,
	    0., abscissa, ordinate, height, source_water, volume_prev, 0, is_volume_source, source_id),
    quantity_list(quantity), species_list(species_list),
    photolysis_reaction_list(photolysis_reaction_list)
  {
    this->Ns = species_list.size();
    Nr_photolysis = photolysis_reaction_list.size();
  }

  //! Legacy constructor for GaussianPuff.
  /*!
    \param time_puff time at which the puff is released (seconds).
    \param quantity total mass released by the source (mass unit).
    \param source_abscissa abscissa (meters).
    \param source_ordinate ordinate (meters).
    \param source_height height (meters).
    \param species_index index associated to the species.
  */
  template<class T>
  PuffChemistry<T>::PuffChemistry(T release_time, T velocity,  T temperature,
                                  T diameter, vector<T> quantity, T abscissa,
                                  T ordinate, T height,
                                  vector<string> species_list,
                                  vector<string> photolysis_reaction_list,
                                  string source_id):
    Puff<T>(release_time, velocity,  temperature, diameter,
            0., abscissa, ordinate, height, 0, source_id),
    quantity_list(quantity), species_list(species_list),
    photolysis_reaction_list(photolysis_reaction_list)
  {
    this->Ns = species_list.size();
    Nr_photolysis = photolysis_reaction_list.size();
  }

  //! Initializes the puff.
  template<class T>
  void PuffChemistry<T>::InitPuff()
  {
    Puff<T>::InitPuff();
    scavenging_coefficient.resize(this->Ns);
    deposition_velocity.resize(this->Ns);
    reflection_factor.resize(this->Ns);
    photolysis_rate.resize(Nr_photolysis);
    background_concentration.resize(this->Ns);
    previous_background_concentration.resize(this->Ns);
    for (int i = 0; i < this->Ns; i++)
      {
        scavenging_coefficient(i) = 0.;
        deposition_velocity(i) = 0.;
        reflection_factor(i) = 1.;
      }
    for (int r = 0; r < Nr_photolysis; r++)
      photolysis_rate(r) = 0.;
  }


  //! Returns the number of species.
  /*!
    \return The number of species.
  */
  template<class T>
  inline int PuffChemistry<T>::GetNs()
  {
    return this->Ns;
  }


  //! Returns the puff quantity for a given species.
  /*!
    \return The species quantity (mass unit).
  */
  template<class T>
  inline T PuffChemistry<T>::GetQuantity(int s) const
  {
    return quantity_list[s];
  }


  //! Gets values of additional meteorological data.
  /*! It gets the additionnal meteorological data (needed for chemistry).
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void PuffChemistry<T>::GetAdditionalMeteo(T& attenuation, T& pressure,
                                            T& specific_humidity) const
  {
    attenuation = attenuation_;
    pressure = pressure_;
    specific_humidity = specific_humidity_;
  }

  //! Returns the puff specific humidity.
  /*!
    \return the specific humidity.
  */
  template<class T>
  T PuffChemistry<T>::GetSpecificHumidity() const
  {
    return specific_humidity_;
  }


  //! Gets scavenging coefficient for a given species.
  template<class T>
  inline T PuffChemistry<T>::GetScavengingCoefficient(int s) const
  {
    return scavenging_coefficient(s);
  }


  //! Gets deposition velocity for a given species.
  template<class T>
  inline T PuffChemistry<T>::GetDepositionVelocity(int s) const
  {
    return deposition_velocity(s);
  }


  //! Gets the fraction of the puff reflected by the ground.
  /*!
    \return the reflection factor.
  */
  template<class T>
  T PuffChemistry<T>::GetReflectionFactor(int s) const
  {
    return reflection_factor(s);
  }


  //! Gets photolysis for a given species.
  template<class T>
  inline T PuffChemistry<T>::GetPhotolysisRate(int s) const
  {
    return photolysis_rate(s);
  }


  //! Gets background concentration for a given species.
  template<class T>
  inline T PuffChemistry<T>::GetBackgroundConcentration(int s) const
  {
    return background_concentration(s);
  }

  //! Gets background concentration for a given species.
  template<class T>
  inline T PuffChemistry<T>::GetPreviousBackgroundConcentration(int s) const
  {
    return previous_background_concentration(s);
  }


  //! Sets the quantity for a given species.
  /*!
    \param quantity the species quantity (mass unit).
  */
  template<class T>
  void PuffChemistry<T>::SetQuantity(T quantity, int s)
  {
    quantity_list[s] = quantity;
  }


  //! Sets values of additional meteorological data.
  /*! It sets the additionnal meteorological data (needed for chemistry).
    \param attenuation: the attenuation.
    \param pressure: the pressure (Pa).
    \param specific_humidity: the specific humidity.
  */
  template<class T>
  void PuffChemistry<T>::SetAdditionalMeteo(T attenuation, T pressure,
                                            T specific_humidity)
  {
    attenuation_ = attenuation;
    pressure_ = pressure;
    specific_humidity_ = specific_humidity;
  }


  //! Sets the species scavenging coefficient for a given species.
  /*!
    \param lambda the species scavenging coefficient ( s^(-1) ).
  */
  template<class T>
  void PuffChemistry<T>::SetScavengingCoefficient(T lambda, int s)
  {
    scavenging_coefficient(s) = lambda;
  }


  //! Sets the species deposition velocity for a given species.
  /*!
    \param lambda the species deposition velocity (m/s).
  */
  template<class T>
  void PuffChemistry<T>::SetDepositionVelocity(T vd, int s)
  {
    deposition_velocity(s) = vd;
  }


  //! Sets the fraction of the puff reflected by the ground.
  /*!
    \param alpha reflection factor.
  */
  template<class T>
  void PuffChemistry<T>::SetReflectionFactor(T alpha, int s)
  {
    reflection_factor(s) = alpha;
  }


  //! Sets the species photolysis rate for a given species.
  /*!
    \param rate the species photolysis rate.
  */
  template<class T>
  void PuffChemistry<T>::SetPhotolysisRate(T rate, int s)
  {
    photolysis_rate(s) = rate;
  }


  //! Sets background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void PuffChemistry<T>::SetBackgroundConcentration(T concentration, int s)
  {
    background_concentration(s) = concentration;
  }

  //! Sets background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void PuffChemistry<T>::SetPreviousBackgroundConcentration(T concentration, int s)
  {
    previous_background_concentration(s) = concentration;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFCHEMISTRY_CXX
#endif
