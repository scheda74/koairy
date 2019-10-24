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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFCHEMISTRY_HXX



#include "PuffTransport.cxx"
#include <list>

namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////

  using namespace std;


  //////////////
  // GAUSSIAN //
  //////////////


  //! Class that stores all data about a puff with chemistry.
  template<class T>
  class PuffChemistry: public Puff<T>
  {

  protected:

    //! List of species quantities.
    vector<T> quantity_list;

    //! List of species.
    vector<string> species_list;

    //! Number of photolysis reactions.
    int Nr_photolysis;

    //! List of species with photolysis reactions.
    vector<string> photolysis_reaction_list;

    //! Attenuation.
    T attenuation_;

    //! Pressure.
    T pressure_;

    //! Specific humidity.
    T specific_humidity_;

    //! List of species scavenging coefficients ( s^(-1) ).
    Array<T, 1> scavenging_coefficient;

    //! List of species deposition velocities (m/s).
    Array<T, 1> deposition_velocity;

    //! List of species reflection factors.
    Array<T, 1> reflection_factor;

    //! List of photolysis rates.
    Array<T, 1> photolysis_rate;

    //! List of background species concentrations.
    Array<T, 1> background_concentration;
    Array<T, 1> previous_background_concentration;

  public:

    /*** Constructor ***/

    PuffChemistry(T time_puff, T velocity, T temperature, T diameter,
		  vector<T> quantity, T abscissa, T ordinate, T height,
		  T source_water,T volume_prev,
		  vector<string> species,
		  vector<string> photolysis_reaction_list, string source_id);

    PuffChemistry(T time_puff, T velocity, T temperature, T diameter, T width, T length,
		  vector<T> quantity, T abscissa, T ordinate, T height,
		  T source_water, T volume_prev,
		  vector<string> species,
		  vector<string> photolysis_reaction_list,
                  bool is_volume_source, string source_id);

    PuffChemistry(T time_puff, T velocity, T temperature, T diameter,
                  vector<T> quantity, T abscissa, T ordinate, T height,
                  vector<string> species,
                  vector<string> photolysis_reaction_list,
                  string source_id);


    /*** Methods ***/

    void InitPuff();

    int GetNs();
    T GetQuantity(int s) const;
    T GetSpecificHumidity() const;
    void GetAdditionalMeteo(T& attenuation, T& pressure,
                            T& specific_humidity) const;
    T GetScavengingCoefficient(int s) const;
    T GetDepositionVelocity(int s) const;
    T GetReflectionFactor(int s) const;
    T GetPhotolysisRate(int s) const;
    T GetBackgroundConcentration(int s) const;
    T GetPreviousBackgroundConcentration(int s) const;

    void SetQuantity(T quantity, int s);
    void SetAdditionalMeteo(T attenuation, T pressure, T specific_humidity);
    void SetScavengingCoefficient(T lambda, int s);
    void SetDepositionVelocity(T vd, int s);
    void SetReflectionFactor(T alpha, int s);
    void SetPhotolysisRate(T rate, int s);
    void SetBackgroundConcentration(T concentration, int s);
    void SetPreviousBackgroundConcentration(T concentration, int s);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFCHEMISTRY_HXX
#endif
