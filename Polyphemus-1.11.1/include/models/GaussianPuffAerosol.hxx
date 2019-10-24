// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFAEROSOL_HXX


#include "PuffAerosol.cxx"
#include "GaussianPuffChemistry.cxx"

namespace Polyphemus
{


  using namespace std;


  /////////////////////////////
  // GAUSSIAN PUFF CHEMISTRY //
  /////////////////////////////


  //! Gaussian puff model with gaseous chemistry.
  template<class T, class ClassChemistry>
  class GaussianPuffAerosol: public GaussianPuffChemistry<T, ClassChemistry>
  {

  protected:

    /********* Loss processes *********/

    //! List of species scavenging coefficients ( s^(-1) ).
    Array<T, 2> scavenging_coefficient_aer;

    //! List of species deposition velocities (m/s).
    Array<T, 2> deposition_velocity_aer;

    /********* Chemistry *********/

    //! Directory to save the puff mass.
    string directory_puff_mass;
    T delta_t_output;

    //! Directory to save the puff concentration.
    string directory_puff_concentration;

    //! Directory to save the puff background concentration.
    string directory_puff_background_concentration;


    //! Directory to save the puff coordinate
    string directory_puff_coordinate;

    //! Directory to save the puff number concentration
    string directory_puff_number;

    //! Directory to save the puff number concentration
    string directory_puff_number_concentration;

    /*** Species ***/

    //! List of aerosol bins.
    vector<string> bin_list;
    vector<float> Mass_Density_aer;

   //! List of aerosol fraction.
    vector<string> fraction_list;
    vector<string> fraction_list_2;
    //! Bins bounds (in m).
    Array<T, 1> Fractionbound_aer;
    Array<T, 1> Fractionbound_aer_2;
    
    //! Number of aerosol fraction sections.
    int Nfraction_aer;
    int Nfraction_aer_2;

    //! List of photolysis rates.
    Array<T, 1> photolysis_rate_;

    //! List of background concentrations.
    Array<T, 2> background_concentration_aer_;

    //! List of Number background concentrations.
    Array<T, 1> background_concentration_number_;

    // /********* Background Meteorological data ************/

    //! Liquid water content.
    T background_liquid_water_content_;
    T background_ph_;


    // /********* Current Meteorological data ************/

    T relative_humidity_;
    T liquid_water_content_;
    T ph_;

  public:

    using GaussianPuffChemistry<T, ClassChemistry>::Chemistry;

    /*** Constructors and destructor ***/

    //    GaussianPuffChemistry();
    GaussianPuffAerosol(string config_file);
    virtual ~GaussianPuffAerosol();

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    virtual void Init();
    void InitSource(string file_puff, T delta_t_puff);
    void InitPuff();

    /*** Background Meteorological data ***/

    virtual    void InitMeteo(ConfigStream& meteo, bool show_meteo);
    virtual    void InitPhotolysis(int Nr, vector<string> photolysis_list);
    virtual    void InitChemistry(ConfigStream& meteo);
    void InitScavenging(ConfigStream& meteo);
    void InitDeposition(ConfigStream& meteo);
    void SetPuffCurrentMeteo(T interaction_coefficient, Puff<T>* puff);
    virtual void SetCurrentMeteo(Puff<T>* puff);

    /*** Access methods for Puff attributes ***/

    virtual void SetPuffAdditionalMeteo_aer(int index, T liquid_water_content);
    void SetPuffBackgroundConcentration_aer(int index, Array<T, 2> concentration_aer);
    void SetPuffBackgroundConcentration_number(int index, Array<T, 1> concentration_number);
    T GetPuffBackgroundConcentration_aer(int index, int s, int b);
    T GetPuffBackgroundConcentration_number(int index, int b);
    void SetPuffQuantity_aer(int index, Array<T,2> quantity_aer);
    void SetPuffQuantity_number(int index, Array<T,1> quantity_number);
    T GetPuffQuantity_aer(int index, int s, int b);
    T GetPuffQuantity_number(int index, int b);

    /*** Computational Methods ***/
    virtual    void Chemistry();
    void PuffOutputSaver();
    void SavePuffQuantity();
    void SavePlumeQuantity();
    void SavePlumeNumberQuantity();
    void SavePuffNumberQuantity();
    // virtual    void Chemistry(vector<list<int> > PuffCellList,
    // 			      vector<T> PuffCellVolume,
    // 			      Array<vector<T>, 1 >& PuffCellConcentration,
    // 			      Array<vector<T>, 2 >& PuffCellConcentration_aer,
    // 			      Array<vector<T>, 1 >& PuffCellConcentration_number);
    virtual    void Forward();
    virtual    void InitWetDiameter_aer(T& RelativeHumidity_, 
                       T& Temperature_, 
                       Array<T, 1>& WetDiameter_aer_);
    T ComputeDensity(Data<T,1> Conc_aer_tmp, 
		     vector<float> Rho_species, float TotalMass, int Ns);
    void ComputePuffNumberEmission_aer(int emission, Array<T, 1>& puff_number_emissions);
    T ComputeNumberPuff(int b, Array<T, 2> concentration_list_aer);
    T ComputeNumberPuffDiameter(int b, T diameter, Array<T, 2> concentration_list_aer);
    void ComputeRedistribution(Array<T, 2>& concentration_list_aer, 
			       Array<T, 1>& concentration_list_number,
			       int bin);
    T ComputeDiameter(int b,
		      Array<T, 2> concentration_list_aer, 
		      Array<T, 1> concentration_list_number);

    void ComputePrecursorContribution_aer(Array<T, 1> quantity_puf_beta,
				       Array<T, 1>& precursors_contribution);

    void ComputePrecursorContribution(Array<T, 1> quantity_puf_beta,
				   Array<T, 1>& precursors_contribution);

    T ComputeOutputNumber(int b, 
			  Array<T, 2> concentration_list_aer,
			  Array<T, 1> concentration_list_number,
			  Array<T, 2> background_concentration_aer,
			  Array<T, 1> background_concentration_number,
			  vector<float> Rho_species);

    T ComputeInputNumber(int b, 
			 Array<T, 2> concentration_list_aer,
			 Array<T, 1> delta_concentration_list_number,
			 Array<T, 2> background_concentration_aer,
			 Array<T, 1> background_concentration_number,
			 vector<float> Rho_species,
			 int write_in);

    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;
  //! relations between aerosol species index and groups index
    Array<int, 1> aerosol_species_group_relation;
    //! List of aerosol compositions
    Data<T, 3> composition_bounds;

    Array<T, 1> wet_diameter_aer;

    /*** Domain ***/
    //! 5D grid along z.
    RegularGrid<T> GridZ5D;
    //! 5D grid along y.
    RegularGrid<T> GridY5D;
    //! 5D grid along x.
    RegularGrid<T> GridX5D;

    //! 4D grid along z.
    RegularGrid<T> GridZ4D;
    //! 5D grid along y.
    RegularGrid<T> GridY4D;
    //! 5D grid along x.
    RegularGrid<T> GridX4D;

    /*** Species ***/
    //! 5D grid for aerosol bins.
    RegularGrid<T> GridB5D_aer;
    //! 5D grid for aerosol species.
    RegularGrid<T> GridS5D_aer;
    //! 4D grid for gas species.
    RegularGrid<T> GridS4D;
    //! 4D grid for number species.
    RegularGrid<T> GridS4D_number;

    void ComputeLossFactorAerosol(Puff<T>* puff, int s, int b);
    T ComputePuffConcentrationAerosol(Puff<T>* puff, int s, int b, T x, T y, T z);
    T ComputePuffIntegral_aer(int index, int s, int b, T x, T y, T z,
                              T lx, T ly, T lz);
    T ComputePuffIntegral_aer(Puff<T>* puff, int s, int b, T x, T y, T z,
			  T lz, T ly, T lx);
  
    T ComputePuffIntegral_number(int index, int b, T x, T y, T z,
                          T lx, T ly, T lz);
    T ComputePuffIntegral_number(Puff<T>* puff, int b, T x, T y, T z,
			  T lz, T ly, T lx);

    void CombineOverlappingPuff(int alpha, int beta, T volume_alpha, T volume_beta);
    void GetMassDensity_aer(Array<T, 1>& MassDensity_aer);

    T ComputePuffAdditionalLWC(T temperature, T pressure, T specific_humidity);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFAEROSOL_HXX
#endif
