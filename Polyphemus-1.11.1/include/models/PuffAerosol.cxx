// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFAEROSOL_CXX


#include "PuffAerosol.hxx"


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
  PuffAerosol<T>::PuffAerosol(T release_time, T velocity,  T temperature,
                              T diameter, T width, T length,
                              vector<T> quantity, T abscissa,
                              T ordinate, T height, T source_water,
							  T volume_prev,
                              vector<string> species_list,
                              vector<string> photolysis_reaction_list,
                              vector<string> species_list_aer,
                              vector<string> bin_list,
                              Array<T, 2> quantity_aer, bool is_volume_source, string source_id):
    PuffChemistry<T>(release_time, velocity,  temperature, diameter,
                     width, length,
                     quantity, abscissa, ordinate, height, source_water, volume_prev, species_list,
                     photolysis_reaction_list, is_volume_source, source_id)
  {
    this->Ns = species_list.size();
    this->Nr_photolysis = photolysis_reaction_list.size();
    Nbin_aer = bin_list.size() - 1;
    Ns_aer = species_list_aer.size();
    quantity_list_aer.resize(Ns_aer, Nbin_aer);
    quantity_list_aer = quantity_aer;
  }

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
  PuffAerosol<T>::PuffAerosol(T release_time, T velocity,  T temperature,
                              T diameter, T width, T length,
                              vector<T> quantity, T abscissa,
			      T ordinate, T height, T source_water,
			      T volume_prev,
			      vector<string> species_list,
                              vector<string> photolysis_reaction_list,
                              vector<string> species_list_aer,
                              vector<string> bin_list,
                              Array<T, 2> quantity_aer, Array<T, 1> quantity_aer_number,
			      bool is_volume_source, string source_id):
    PuffChemistry<T>(release_time, velocity,  temperature, diameter, 
                     width, length, 
                     quantity, abscissa, ordinate, height, source_water,
		     volume_prev, species_list, 
                     photolysis_reaction_list, is_volume_source, source_id)
  {
    this->Ns = species_list.size();
    this->Nr_photolysis = photolysis_reaction_list.size();
    Nbin_aer = bin_list.size() - 1;
    Ns_aer = species_list_aer.size();
    quantity_list_aer.resize(Ns_aer, Nbin_aer);
    quantity_list_number.resize(Nbin_aer);
    quantity_list_aer = quantity_aer;
    quantity_list_number = quantity_aer_number;
  }
  


  //! Initializes the puff.
  template<class T>
  void PuffAerosol<T>::InitPuff()
  {
    PuffChemistry<T>::InitPuff();
    reflection_factor_aer.resize(this->Ns_aer, this->Nbin_aer);
    background_concentration_aer.resize(this->Ns_aer, this->Nbin_aer);
    background_concentration_aer = 0.;
    background_concentration_number.resize(this->Nbin_aer);
    background_concentration_number = 0.;
    previous_background_concentration_aer.resize(this->Ns_aer, this->Nbin_aer);
    previous_background_concentration_aer = 0.;
    previous_background_concentration_number.resize(this->Nbin_aer);
    previous_background_concentration_number = 0.;
    
    reflection_factor_aer = 1.;
  }


  //! Returns the number of species.
  /*!
    \return The number of species.
  */
  template<class T>
  inline int PuffAerosol<T>::GetNs_aer()
  {
    return this->Ns_aer;
  }

  //! Returns the puff quantity for a given species.
  /*!
    \return The species quantity (mass unit).
  */
  template<class T>
  inline T PuffAerosol<T>::GetQuantity(int s, int b) const
  {
    return quantity_list_aer(s, b);
  }

  //! Returns the puff quantity for a given aerosol bin.
  /*!
    \return The bin quantity (number unit).
  */
  template<class T>
  inline T PuffAerosol<T>::GetNumberQuantity(int b) const
  {
    return quantity_list_number(b);
  }


  //! Gets values of additional meteorological data.
  /*! It gets the additionnal meteorological data (needed for chemistry).
    \param liquid_water_content:
  */
  template<class T>
  void PuffAerosol<T>::GetAdditionalMeteo(T& liquid_water_content) const
  {
    liquid_water_content = liquid_water_content_;
  }


  //! Gets scavenging coefficient for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetScavengingCoefficient(int s, int b) const
  {
    return scavenging_coefficient_aer(s, b);
  }


  //! Gets deposition velocity for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetDepositionVelocity(int s, int b) const
  {
    return deposition_velocity_aer(s, b);
  }


  // //! Gets the fraction of the puff reflected by the ground.
  // /*!
  //   \return the reflection factor.
  // */
  template<class T>
  T PuffAerosol<T>::GetReflectionFactor(int s, int b) const
  {
    return reflection_factor_aer(s, b);
  }


  //! Gets background concentration for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetBackgroundConcentration(int s, int b) const
  {
    return background_concentration_aer(s, b);
  }

 //! Gets background concentration for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetPreviousBackgroundConcentration(int s, int b) const
  {
    return previous_background_concentration_aer(s, b);
  }

  //! Gets background number concentration for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetBackgroundNumberConcentration(int b) const
  {
    return background_concentration_number(b);
  }

  //! Gets background number concentration for a given species.
  template<class T>
  inline T PuffAerosol<T>::GetPreviousBackgroundNumberConcentration(int b) const
  {
    return previous_background_concentration_number(b);
  }



  //! Sets the quantity for a given aerosol species.
  /*!
    \param quantity the species quantity (mass unit).
  */
  template<class T>
  void PuffAerosol<T>::SetQuantity(T quantity, int s, int b)
  {
    quantity_list_aer(s, b) = quantity;
  }

  //! Sets the quantity for a given aerosol bin.
  /*!
    \param quantity the bin quantity (number unit).
  */
  template<class T>
  void PuffAerosol<T>::SetNumberQuantity(T quantity, int b)
  {
    quantity_list_number(b) = quantity;
  }


  //! Sets values of additional meteorological data.
  /*! It sets the additionnal meteorological data (needed for chemistry).
    \param liquid_water_content: the liquid water content (?).
  */
  template<class T>
  void PuffAerosol<T>::SetAdditionalMeteo(T liquid_water_content)
  {
    liquid_water_content_ = liquid_water_content;
  }


  //! Sets the fraction of the puff reflected by the ground.
  /*!
    \param alpha reflection factor.
  */
  template<class T>
  void PuffAerosol<T>::SetReflectionFactor(T alpha, int s, int b)
  {
    reflection_factor_aer(s, b) = alpha;
  }


  //! Sets background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void PuffAerosol<T>::SetBackgroundConcentration(T concentration, int s, int b)
  {
    background_concentration_aer(s, b) = concentration;
  }

  //! Sets background concentration for a given species.
  /*!
    \param concentration the species background concentration.
  */
  template<class T>
  void PuffAerosol<T>::SetPreviousBackgroundConcentration(T concentration, int s, int b)
  {
    previous_background_concentration_aer(s, b) = concentration;
  }

  //! Sets background number concentration for a given species.
  /*!
    \param concentration the species background number concentration.
  */
  template<class T>
  void PuffAerosol<T>::SetBackgroundNumberConcentration(T concentration, int b)
  {
    background_concentration_number(b) = concentration;
  }

  //! Sets background number concentration for a given species.
  /*!
    \param concentration the species background number concentration.
  */
  template<class T>
  void PuffAerosol<T>::SetPreviousBackgroundNumberConcentration(T concentration, int b)
  {
    previous_background_concentration_number(b) = concentration;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFAEROSOL_CXX
#endif
