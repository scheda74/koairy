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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFAEROSOL_HXX



#include "PuffChemistry.cxx"

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
  class PuffAerosol: public PuffChemistry<T>
  {

  protected:

    //! List of aerosol species quantities.
    Array<T, 2> quantity_list_aer;
    Array<T, 1> quantity_list_number;

    T liquid_water_content_;
    T ph_;

    //! List of species scavenging coefficients ( s^(-1) ).
    Array<T, 2> scavenging_coefficient_aer;

    //! List of species deposition velocities (m/s).
    Array<T, 2> deposition_velocity_aer;

    //! List of species reflection factors.
    Array<T, 2> reflection_factor_aer;

    // //! List of background species concentrations.
    Array<T, 2> background_concentration_aer;
    Array<T, 2> previous_background_concentration_aer;

    // //! List of background species number concentrations.
    Array<T, 1> background_concentration_number;
    Array<T, 1> previous_background_concentration_number;


    int Nbin_aer;
    int Ns_aer;

  public:

    /*** Constructor ***/

    PuffAerosol(T time_puff, T velocity, T temperature, T diameter,
		vector<T> quantity, T abscissa, T ordinate, T height,
		T source_water, T volume_prev,
		vector<string> species,
                vector<string> photolysis_reaction_list,
                vector<string> species_list_aer,
                vector<string> bin_list,
                Array<T, 2> quantity_aer, string source_id);

    PuffAerosol(T time_puff, T velocity, T temperature, T diameter,
                T width, T length,
		vector<T> quantity, T abscissa, T ordinate, T height,
		T source_water, T volume_prev,
		vector<string> species,
                vector<string> photolysis_reaction_list,
                vector<string> species_list_aer,
                vector<string> bin_list,
                Array<T, 2> quantity_aer, bool is_volume_source, string source_id);
  
    /*** Constructor for number computation ***/
    PuffAerosol(T time_puff, T velocity, T temperature, T diameter,
		vector<T> quantity, T abscissa, T ordinate, T height,
		T source_water, T volume_prev,
		vector<string> species,
                vector<string> photolysis_reaction_list,
                vector<string> species_list_aer,
                vector<string> bin_list,
                Array<T, 2> quantity_aer, 
                Array<T, 1> quantity_aer_number, 
		string source_id);

    PuffAerosol(T time_puff, T velocity, T temperature, T diameter,
                T width, T length,
		vector<T> quantity, T abscissa, T ordinate, T height,
		T source_water, T volume_prev,
		vector<string> species,
                vector<string> photolysis_reaction_list,
                vector<string> species_list_aer,
                vector<string> bin_list,
                Array<T, 2> quantity_aer, 
                Array<T, 1> quantity_aer_number, 
		bool is_volume_source, string source_id);

    /*** Methods ***/

    void InitPuff();

    int GetNs_aer();
    T GetQuantity(int s, int b) const;
    T GetNumberQuantity(int b) const;
    void GetAdditionalMeteo(T& liquid_water_content) const;
    T GetScavengingCoefficient(int s, int b) const;
    T GetDepositionVelocity(int s, int b) const;
    T GetReflectionFactor(int s, int b) const;
    T GetBackgroundConcentration(int s, int b) const;
    T GetPreviousBackgroundConcentration(int s, int b) const;
    T GetBackgroundNumberConcentration(int b) const;
    T GetPreviousBackgroundNumberConcentration(int b) const;

    void SetQuantity(T quantity, int s, int b);
    void SetNumberQuantity(T quantity, int b);
    void SetAdditionalMeteo(T liquid_water_content);
    void SetScavengingCoefficient(T lambda, int s, int b);
    void SetDepositionVelocity(T vd, int s, int b);
    void SetReflectionFactor(T alpha, int s, int b);
    void SetBackgroundConcentration(T concentration, int s, int b);
    void SetPreviousBackgroundConcentration(T concentration, int s, int b);
    void SetBackgroundNumberConcentration(T concentration, int b);
    void SetPreviousBackgroundNumberConcentration(T concentration, int b);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFAEROSOL_HXX
#endif
