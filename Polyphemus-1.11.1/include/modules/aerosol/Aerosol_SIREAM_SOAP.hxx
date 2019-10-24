// Copyright (C) 2007-2017, ENPC - INRIA - EDF R&D
// Author(s): Edouard Debry, Youngseob Kim, Stephanie Deschamps
//
// This file is part of the Size Resolved Aerosol Model (SIREAM), a component
// of the air quality modeling system Polyphemus.
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


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SOAP_HXX

extern "C"
{
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/ipc.h>
}

#include <vector>
#include <map>
#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"
#include "parameters.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////

#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _compute_coagulation_coefficient compute_coagulation_coefficient__
#else
#define _compute_coagulation_coefficient compute_coagulation_coefficient_
#endif

#define _chem_0d chem_0d_
#define _aerosol_0d aerosol_0d_
#define _dimensions dimensions_

  extern "C"
  {
    void _chem_0d(int*, int*, int*, int*, int*, int*,
                  double*, double*, int*,
                  double*, double*, double*, double*, double*, double*,
                  double*, double*, double*, double*, double*,
                  double*, double*, double*, double*, double*, int*,
                  double*, double*, double*, int*, int*, int*, int*,
                  double*, double*, double*, double*, int*, double*,
                  int*, double*, double*, int*, int*, int*, int*, int*,
						double*, double*);

    void _aerosol_0d(int*, int*, int*, int*, int*,
                     double*, double*, double*, double*, double*, double*,
                     double*, int*, int*, int*, int*, int*, double*, double*,
                     double*, double*, double*, int*, int*, int*, double*,
                     double*, double*, double*, double*, double*, double*,
                     double*, double*, double*, double*, double*, double*,
                     double*, double*, double*, double*, double*, int*, int*,
                     int*, int*, int*, int*, int*, int*, int*,
                     double*, double*, int*, double*, double*, double*,
					 double*, double*, model_config*, vector<species>*);

    void _compute_coagulation_coefficient(int*, double*,
                                          int*, int*, int*, double*);

    void _dimensions(int*, int*, int*);

  }


  //////////////////////////
  // AEROSOL_SIREAM_SOAP ///
  //////////////////////////


  //! \brief This class is a numerical solver for the chemical mechanism
  //! RACM_SIREAM_SOAP, CB05_SIREAM_SOAP and RACM2_SIREAM_SOAP.
  /*! It uses a second-order Rosenbrock method.
   */
  template<class T>
  class Aerosol_SIREAM_SOAP: public BaseModuleParallel
  {

  protected:
    //! Number of species.
    int Ns;
    //! Number of reactions.
    int Nr;
    //! Number of photolysis reactions.
    int Nr_photolysis;
    //! Number of species with volume sources.
    int Ns_source;
    //! Number of levels with volume sources.
    int Nz_source;

    //! Number of sub-cycles.
    int Ncycle;

    //! Sorted list of species names.
    vector<string> species_list;
    vector<string> species_list_aer;

    //! Molecular weights of species.
    Array<T, 1> molecular_weight;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 1> ConversionFactor;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 2> ConversionFactorJacobian;

    /*! \brief Map between photolysis reactions names and their indices in
      photolysis reactions.
    */
    map<string, int> photolysis_reaction_name;
    //! Indices of photolysis reactions among other reactions.
    Array<int, 1> photolysis_reaction_index;
    //! Indices of species with volume sources.
    Array<int, 1> source_index;

    //! Number of sub-cycles ofr aerosol processes
    int Ncycle_aer;

    //! Number of aerosol species.
    int Ns_aer;

    //! Number of aerosol bins.
    int Nbin_aer;

    //!Indice for PBiPER number
    int jBiPER;

    //!Indice for BiPER degradation kinetics
    int kBiPER;

    //! Bin index corresponding to fixed cutting diameter.
    int cutting_bin;
    /* bins from 1 to cutting_bin included are at equilibrium
       and are dynamic from cutting_bin+1 to Nbin_aer */

    //! Indices of species with heterogeneous reactions.
    Array<int, 1> heterogeneous_reaction_index;

    //! list of aerosol options.
    int Noptions_aer;
    Array<int, 1> options_aer;

  public:
    /*** Configuration ***/

    map<string, bool> option_process_aer;

    //! Numerical solver for dynamic bin condensation (etr, ros2 or ebi).
    string dynamic_condensation_solver;

    //! Cutting diameter between equilibrium and dynamic bins.
    double fixed_cutting_diameter;

    //! Sulfate condensation computation method (equilibrium, dynamic).
    string sulfate_computation;

    //! SOA computation method (equilibrium, dynamic).
    string soa_computation;

    /*! Redistribution method of lagrangian bins (number-conserving,
      interpolation).
    */
    string redistribution_method;

	//! Normalisation for mass concervation with redistribution : euler-number, hemen, euler-coupled
	string with_normalisation;
	//! Conserving mass tolerance
	double conserving_mass_tolerance;

    //! If nucleation, which nucleation model (binary, ternary).
    string nucleation_model;
    double K_nucl_factor;
    double P_nucl_factor;

    //! Aerosol wet diameter estimation (gerber or isoropia).
    string wet_diameter_estimation;

    //! Thermodynamic_model (constant, ideal or unifac).
    string Thermodynamic_model;

    //! Name of the aqueous module.
    string aqueous_module;

    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;

    //! Aerosol density (kg / m^3).
    double FixedDensity_aer;
    Array<T, 1> Density_aer;

    //! Coagulation coefficients.
    Array<int, 1> CoagulationCouple;
    Array<int, 2> CoagulationFirstIndex;
    Array<int, 2> CoagulationSecondIndex;
    Array<T, 3> CoagulationPartitionCoefficient;

    //! liquid water content threshold for clouds.
    double lwc_cloud_threshold;

    //! Options for adaptive time stepping for gaseous chemistry
    bool with_adaptive_time_step;
    int option_adaptive_time_step;
    double adaptive_time_step_tolerance, min_adaptive_time_step;

    //! Option for photolysis
    int option_photolysis_tabulation;

    //! Saturation pressure (Pascals).
    Array<T, 1> saturation_pressure;

    //! Partition coefficient (m^3.\mu g^-1).
    Array<T, 1> partition_coefficient;

    //! Vaporization enthalpy (J.mol^-1).
    Array<T, 1> vaporization_enthalpy;

    //! Accomodation coefficient, dimensionless, between 0 and 1.
    Array<T, 1> accomodation_coefficient;

    //! Surface tension (N.m^-1).
    Array<T, 1> surface_tension;

    //!  Saturation pressure (m^3.\mu g^-1).
    Array<T, 1> saturation_pressure_mass;

    //!  Saturation pressure (torr).
    Array<T, 1> saturation_pressure_torr;

    //!  Deliquescence relative humidity ().
    Array<T, 1> deliquescence_relative_humidity;

    //! Molecular weight (\mu g.mol^-1).
    Array<T, 1> molecular_weight_aer;

    //! Molecular diameter (Angstreum).
    Array<T, 1> molecular_diameter_aer;

    //! Collision factor ().
    Array<T, 1> collision_factor_aer;

    //! Mass density (\mu g.\mu m^-3).
    Array<T, 1> mass_density_aer;

    //! Gas species aerosol interacting index ().
    Array<int, 1> species_index_aerosol_interact;

    //! Gas species cloud interacting index ().
    int Ns_cloud_interact;
    vector<int> species_index_cloud_interact;

    //! Index of aerosol species managed by isorropia equilibrium.
    int Ns_isorropia;
    vector<int> isorropia_species;

    //! Index of aerosol species managed by aec equilibrium.
    int Ns_aec;
    vector<int> aec_species;

    //! Index of aerosol species managed by pankow equilibrium.
    int Ns_pankow;
    vector<int> pankow_species;

    //! Index of primary organic aerosol species.
    int Ns_poa;
    vector<int> poa_species;

    int md_species, bc_species;

    /*! \brief SOAP model parameters defined in species.h
      \param model_config class for model configuration
      \param species class for species   
     */
    //! Model configuration defined in parameter.cxx 
    model_config soap_config;
    //! Vector of the species defined in species.cxx
    vector<species> surrogate;

    //! Diffusion coefficient in the organic phase.
    double dorg;

    //! Number of aerosol layers.
    int nlayer;


    /*** Constructor ***/
    
    Aerosol_SIREAM_SOAP();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void DisplayConfiguration(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);

    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    void Forward(T current_time,
                 Data<T, 3>& Attenuation_i,
                 Data<T, 3>& SpecificHumidity_i,
                 Data<T, 3>& Temperature_i,
                 Data<T, 3>& Pressure_i,
                 Data<T, 4>& VolumeEmission_i,
                 Data<T, 4>& PhotolysisRate_i,
                 T next_time,
                 Data<T, 3>& Attenuation_f,
                 Data<T, 3>& SpecificHumidity_f,
                 Data<T, 3>& Temperature_f,
                 Data<T, 3>& Pressure_f,
                 Data<T, 4>& VolumeEmission_f,
                 Data<T, 4>& PhotolysisRate_f,
                 Array<T, 1>& Longitude,
                 Array<T, 1>& Latitude,
                 Data<T, 4>& Concentration);

    void Forward(T current_time,
                 T& attenuation,
                 T& humid,
                 T& temperature,
                 T& pressure,
                 Array<T, 1>& source,
                 Array<T, 1>& photolysis,
                 T next_time,
                 T& attenuation_f,
                 T& humid_f,
                 T& temperature_f,
                 T& pressure_f,
                 Array<T, 1>& source_f,
                 Array<T, 1>& photolysis_f,
                 T& lon, T& lat,
                 Array<T, 1>& concentration);

    void Forward(T current_time,
                 Data<T, 3>& Attenuation_i,
                 Data<T, 3>& SpecificHumidity_i,
                 Data<T, 3>& Temperature_i,
                 Data<T, 3>& Pressure_i,
                 Data<T, 4>& VolumeEmission_i,
                 Data<T, 4>& PhotolysisRate_i,
                 T next_time,
                 Data<T, 3>& Attenuation_f,
                 Data<T, 3>& SpecificHumidity_f,
                 Data<T, 3>& Temperature_f,
                 Data<T, 3>& Pressure_f,
                 Data<T, 4>& VolumeEmission_f,
                 Data<T, 4>& PhotolysisRate_f,
                 Array<T, 1>& Longitude,
                 Array<T, 1>& Latitude,
                 Data<T, 4>& Concentration,
                 Data<T, 3>& LiquidWaterContent_i,
                 Data<T, 4>& WetDiameter_aer,
                 Data<T, 5>& Concentration_aer,
                 Data<T, 3>& pH,
        		Data<T, 4>& NumberConcentration_aer);

    void Forward(T current_time,
                 T& attenuation,
                 T& humid,
                 T& temperature,
                 T& pressure,
                 Array<T, 1>& source,
                 Array<T, 1>& photolysis,
                 T next_time,
                 T& attenuation_f,
                 T& humid_f,
                 T& temperature_f,
                 T& pressure_f,
                 Array<T, 1>& source_f,
                 Array<T, 1>& photolysis_f,
                 T& lon, T& lat,
                 Array<T, 1>& concentration,
                 T& liquid_water_content,
                 Array<T, 1>& wet_diameter_aer,
                 Array<T, 2>& concentration_aer,
                 T& ph,
                 Array<T, 1>& number_concentration_aer);

    void Forward_aer(T current_time,
					 Data<T, 3>& SpecificHumidity_i,
					 Data<T, 3>& Temperature_i,
					 Data<T, 3>& Pressure_i,
					 T next_time,
					 Data<T, 4>& Concentration,
					 Data<T, 3>& LiquidWaterContent_i,
					 Data<T, 2>& Rain_i,
					 Array<T, 1>& VerticalInterface,
					 Data<T, 5>& Concentration_aer,
					 Data<T, 3>& InCloudWetDepositionFlux,
					 Data<T, 4>& InCloudWetDepositionFlux_aer,
					 Data<T, 3>& pH,
					 Data<T, 4>& NumberConcentration_aer,
					 Data<T, 3>& InCloudWetDepositionFluxNumber_aer);
	
	void Forward_aer(T current_time,
					 T& specifichumidity,
					 T& temperature,
					 T& pressure,
					 T delta_t,
					 Array<T, 1>& concentration,
					 T& liquidwatercontent,
					 T& rain,
					 Array<T, 1>& CurrentVerticalInterface,
					 Array<T, 2>& concentration_aer,
					 Array<T, 1>& incloudwetdepositionflux,
					 Array<T, 2>& incloudwetdepositionflux_aer,
                     T& ph,
                     T& lwc_avg, T& heightfog, int& ifog,
					 Array<T, 1>& number_concentration_aer,
			 Array<T, 1>& incloudwetdepositionfluxnumber_aer);


    bool IsRequired(string field);
    bool IsComputed(string field);
    void FogSettling(Array<T, 1>& LiquidWaterContent,
                     T& lwc_cloud_threshold,
                     Array<T, 1>& VerticalInterface,
                     T& lwc_avg,
                     int& nfoglay,
                     T& heightfog);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SOAP_HXX
#endif
