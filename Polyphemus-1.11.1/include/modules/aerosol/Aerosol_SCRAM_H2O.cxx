// Copyright (C) 2014, ENPC - INRIA - EDF R&D
// Author(s): Shupeng ZHU
//
// This file is part of the Size Composition Resolved Aerosol Model (SCRAM), 
// a component of the air quality modeling system Polyphemus.
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


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SCRAM_H2O_CXX

#include "Aerosol_SCRAM_H2O.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  Aerosol_SCRAM_H2O<T>::Aerosol_SCRAM_H2O()
  {
    _dimensions(&Ns,&Nr,&Nr_photolysis);
    Ncycle = 1;
    molecular_weight.resize(Ns);
    photolysis_reaction_index.resize(Nr_photolysis);
    Ncycle_aer = 1;
    Noptions_aer = 16;
	
  }


  //! Initialization of the scheme.
  /*! In case the model is incompatible with the chemical mechanism, an
    exception is thrown.
    \param Model model with the following interface:
    <ul>
    <li> GetNs()
    <li> GetNbin_aer()
    <li> GetNgroup_aer()
    <li> GetNcomposition_aer()
    <li> GetNsize_section_aer()
    <li> GetSpeciesList()
    <li> GetSpeciesFile()
    <li> GetNr_photolysis()
    <li> GetPhotolysisReactionList()
    <li> GetNs_source()
    <li> GetNz_source()
    <li> SourceGlobalIndex(int)
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SCRAM_H2O<T>::Init(ClassModel& Model)
  {
    if (Model.GetNs() != Ns)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNs())
	+ string(" species while chemistry has ") + to_str(Ns) + " species.";
    if (Model.GetNr_photolysis() != Nr_photolysis)
      throw string("Incompatibility: model manages ")
	+ to_str(Model.GetNr_photolysis())
	+ string(" photolysis reactions while chemistry has ")
	+ to_str(Nr_photolysis) + " photolysis reactions.";
    species_list = Model.GetSpeciesList();
    species_list_aer = Model.GetSpeciesList_aer();
    Ncomposition_aer=Model.GetNcomposition_aer();
    Nsize_section_aer= Model.GetNsize_section_aer();
    Nbin_aer = Model.GetNbin_aer();
    Ngroup_aer = Model.GetNgroup_aer();
    Ns_aer = Model.GetNs_aer();

    aerosol_species_group_relation.resize(Ns_aer);
    for(int i = 0; i < Ns_aer; i++)
      aerosol_species_group_relation(i)=Model.aerosol_species_group_relation(i);
    
    composition_bounds.resize(Ncomposition_aer*Ngroup_aer*2);
    int id=0;
    for(int i = 0; i < Ncomposition_aer; i++)
      for(int j = 0; j < Ngroup_aer; j++)
	for(int k = 0; k < 2; k++)
	{
	  composition_bounds(id)=Model.composition_bounds(i,j,k);
	  id++;
	}
    // Photolysis reactions.
    ConfigStream config_species(Model.GetSpeciesFile());
    int iBiPER;
    config_species.SetSection("[molecular_weight]");
    for (int i = 0; i < Ns; i++)
      {
	config_species.PeekValue(species_list[i], molecular_weight(i));
	if (species_list[i]=="BiPER")
	  {
	    iBiPER = i;
	  }
      }
    
    config_species.SetSection("[photolysis_reaction_index]");
    string species;
    for (int i = 0; i < Nr_photolysis; i++)
      {
	species = config_species.GetElement();
	if (species=="BiPER") kBiPER = i+1;
	photolysis_reaction_name[species] =
	  convert<int>(config_species.GetElement());
      }

    for (int i = 0; i < Nr_photolysis; i++)
      photolysis_reaction_index(i) =
	photolysis_reaction_name[Model.GetPhotolysisReactionList()[i]];

    map<string, double> parameter;
    map<string, double>::iterator iter;

    // Conversions.
    double Navogadro = 6.02213e23;
    ConversionFactor.resize(Ns);
    for (int i = 0; i < Ns; i++)
      ConversionFactor(i) = Navogadro * 1e-12 / molecular_weight(i);
    ConversionFactorJacobian.resize(Ns, Ns);
    for (int i = 0; i < Ns; i++)
      for (int j = 0; j < Ns; j++)
	ConversionFactorJacobian(i, j) = molecular_weight(j)
	  / molecular_weight(i);

    // Read molecular weight for aerosol species (\mu g.mol^-1).
    config_species.SetSection("[molecular_weight_aer]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    molecular_weight_aer.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  // Convert from g to \mu g.
	  molecular_weight_aer(i) = iter->second * 1.e6;
	else
	  throw string("Module: no molecular weight for ")
	    + string("aerosol species \"") + species_list_aer[i]
	    + string("\". Please provide one.");
      }
    parameter.clear();

    // Read molecular diameter (Angstreum).
    config_species.SetSection("[molecular_diameter_aer]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    molecular_diameter_aer.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  molecular_diameter_aer(i) = iter->second;
	else
	  molecular_diameter_aer(i) = 1.e10;
      }
    parameter.clear();

    // Read collision factor ().
    config_species.SetSection("[collision_factor_aer]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    collision_factor_aer.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  collision_factor_aer(i) = iter->second;
	else
	  collision_factor_aer(i) = 1.e10;
      }
    parameter.clear();

    // Read Mass density (\mu g.\mu m^-3).
    config_species.SetSection("[mass_density_aer]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    mass_density_aer.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  // Convert from g.cm^-3 to \mu g.\mu m^-3.
	  mass_density_aer(i) = iter->second * 1.e-6;
	else
	  throw string("Module: no mass density for ")
	    + string("aerosol species \"") + species_list_aer[i]
	    + string("\". Please provide one.");
      }
    parameter.clear();

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_aerosol_interact]");
    
    map<string, string> parameter2;
    map<string, string>::iterator iter2;
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	parameter2[species] = config_species.GetElement();
      }
    species_index_aerosol_interact.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter2 = parameter2.find(species_list_aer[i]);
	if(iter2 != parameter2.end())
	  {
	    int gas_index = Model.GetSpeciesIndex(iter2->second) + 1;
	    species_index_aerosol_interact(i) = gas_index;
	  }
	else
	  species_index_aerosol_interact(i) = -1;
      }
    parameter2.clear();

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_cloud_interact]");

    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int gas_index = Model.GetSpeciesIndex(species) + 1;
	species_index_cloud_interact.push_back(gas_index);
      }
    Ns_cloud_interact = int(species_index_cloud_interact.size());

    // Aerosol species managed by isorropia equilibrium.
    config_species.SetSection("[isorropia_species]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	isorropia_species.push_back(aer_index);
      }
    Ns_isorropia = int(isorropia_species.size());

    // Aerosol species managed by AEC equilibrium.
    config_species.SetSection("[aec_species]");
    jBiPER=0;
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	if (species=="PBiPER")
	  {
	    //	jBiPER = aer_index;
	  }
	aec_species.push_back(aer_index);
      }
    Ns_aec = int(aec_species.size());

    // Aerosol species managed by Pankow equilibrium.
    config_species.SetSection("[pankow_species]");

    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	pankow_species.push_back(aer_index);
      }
    Ns_pankow = int(pankow_species.size());

    // Primary organic aerosol species.
    config_species.SetSection("[poa_species]");

    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	poa_species.push_back(aer_index);
      }
    Ns_poa = int(poa_species.size());

    // Mineral dust aerosol species.
    config_species.SetSection("[mineral_dust_species]");
    species = config_species.GetElement();
    md_species = Model.GetSpeciesIndex_aer(species) + 1;

    // Black carbon aerosol species.
    config_species.SetSection("[black_carbon_species]");
    species = config_species.GetElement();
    bc_species = Model.GetSpeciesIndex_aer(species) + 1;

    // Read accomodation coefficients for aerosol species (m^3.\mu g^-1).
    config_species.SetSection("[accomodation_coefficient]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    accomodation_coefficient.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  accomodation_coefficient(i) = iter->second;
	else
	  accomodation_coefficient(i) = 0.0;
      }
    parameter.clear();

    // Read surface tension for aerosol species (m^3.\mu g^-1).
    config_species.SetSection("[surface_tension]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    surface_tension.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  surface_tension(i) = iter->second;
	else
	  surface_tension(i) = 0.0;
      }
    parameter.clear();

    // Read saturation pressure for aerosol species (Pascals).
    config_species.SetSection("[saturation_pressure]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    saturation_pressure.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  saturation_pressure(i) = iter->second;
	else
	  saturation_pressure(i) = 0.0;
      }
    parameter.clear();

    // Read partition coefficients for aerosol species (m^3.\mu g^-1).
    config_species.SetSection("[partition_coefficient]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    partition_coefficient.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  partition_coefficient(i) = iter->second;
	else
	  partition_coefficient(i) = 0.0;
      }
    parameter.clear();

    // Read vaporization enthalpy for aerosol species (J.mol^-1).
    config_species.SetSection("[vaporization_enthalpy]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    vaporization_enthalpy.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  vaporization_enthalpy(i) = iter->second;
	else
	  vaporization_enthalpy(i) = 0.0;
      }
    parameter.clear();

    // Read saturation pressure for aerosol species (m^3.\mu g^-1).
    config_species.SetSection("[saturation_pressure_mass]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    saturation_pressure_mass.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  saturation_pressure_mass(i) = iter->second;
	else
	  saturation_pressure_mass(i) = 0.0;
      }
    parameter.clear();

    // Read saturation pressure for aerosol species (torr).
    config_species.SetSection("[saturation_pressure_torr]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    saturation_pressure_torr.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  saturation_pressure_torr(i) = iter->second;
	else
	  saturation_pressure_torr(i) = 0.0;
      }
    parameter.clear();

    // Deliquescence relative humidity ().
    config_species.SetSection("[deliquescence_relative_humidity]");
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    deliquescence_relative_humidity.resize(Ns_aer);
    for (int i = 0; i < Ns_aer; i++)
      {
	iter = parameter.find(species_list_aer[i]);
	if(iter != parameter.end())
	  deliquescence_relative_humidity(i) = iter->second;
	else
	  deliquescence_relative_humidity(i) = 0.0;
      }
    parameter.clear();

    // Volumic sources.
    Ns_source = Model.GetNs_source();
    Nz_source = Model.GetNz_source();
    source_index.resize(Ns_source);
    for (int i = 0; i < Ns_source; i++)
      source_index(i) = Model.SourceGlobalIndex(i);

    // Module options.
    ConfigStream config(Model.GetConfigurationFile());
    config.SetSection("[options]");
	
    config.PeekValue("With_coagulation",
		     option_process_aer["with_coagulation"]);
    config.PeekValue("With_condensation",
		     option_process_aer["with_condensation"]);
    config.PeekValue("With_nucleation",
		     option_process_aer["with_nucleation"]);
    config.PeekValue("With_in_cloud_scavenging",
		     option_process_aer["with_in_cloud_scavenging"]);
    config.PeekValue("Collect_wet_flux_aerosol",
		     option_process_aer["collect_wet_flux_aer"]);
    config.PeekValue("Collect_wet_flux",
		     option_process_aer["collect_wet_flux"]);
    config.PeekValue("aqueous_module",aqueous_module);
    aqueous_module = lower_case(aqueous_module);
    config.PeekValue("With_heterogeneous_reactions",
		     option_process_aer["with_heterogeneous_reactions"]);
    config.PeekValue("With_kelvin_effect",
		     option_process_aer["with_kelvin_effect"]);
    config.PeekValue("With_number_concentration",
		     option_process_aer["with_number_concentration"]);
    config.PeekValue("With_volume_emission_aerosol",
		     option_process_aer["with_volume_emission_aer"]);
    config.PeekValue("With_oligomerization",
		     option_process_aer["with_oligomerization"]);
    config.PeekValue("Dynamic_condensation_solver",
		     dynamic_condensation_solver);
    dynamic_condensation_solver = lower_case(dynamic_condensation_solver);
    config.PeekValue("Fixed_cutting_diameter", fixed_cutting_diameter);
    fixed_cutting_diameter = 1.e-6 * fixed_cutting_diameter;
    config.PeekValue("Sulfate_computation", sulfate_computation);
    sulfate_computation = lower_case(sulfate_computation);
    config.PeekValue("Redistribution_method", redistribution_method);
    redistribution_method = lower_case(redistribution_method);
    config.PeekValue("With_normalisation", with_normalisation);
    with_normalisation = lower_case(with_normalisation);
    config.PeekValue("Conserving_mass_tolerance", "positive",conserving_mass_tolerance);
    config.PeekValue("Nucleation_model", nucleation_model);
    nucleation_model = lower_case(nucleation_model);
    if (nucleation_model == "unary")
      {
	config.PeekValue("Nucleation_exponent", P_nucl_factor);
	config.PeekValue("Nucleation_prefactor", K_nucl_factor);
      }
    config.PeekValue("With_fixed_density",
		     option_process_aer["with_fixed_density"]);
    config.PeekValue("Wet_diameter_estimation", wet_diameter_estimation);
    wet_diameter_estimation = lower_case(wet_diameter_estimation);
	

    config.PeekValue("With_adaptive_time_step_for_gas_chemistry",
		     with_adaptive_time_step);
    if (with_adaptive_time_step)
      {
	option_adaptive_time_step = 1;
	config.PeekValue("Adaptive_time_step_tolerance",
			 adaptive_time_step_tolerance);
	config.PeekValue("Min_adaptive_time_step", min_adaptive_time_step);
      }
    config.PeekValue("Photolysis_tabulation_option",option_photolysis);
    config.PeekValue("Organic_thermodynamic_model", Thermodynamic_model);
    Thermodynamic_model = lower_case(Thermodynamic_model);

    // Density in kg / m^3.
    config.PeekValue("Fixed_aerosol_density", FixedDensity_aer);
    Density_aer.resize(Nbin_aer);
    Density_aer = FixedDensity_aer;

    // Lwc cloud threshold.
    config.PeekValue("Lwc_cloud_threshold", lwc_cloud_threshold);

    // Reads bin bounds.
    vector<string> bin_list;
    config.SetSection("[domain]");
    config.Find("Bin_bounds");
    bin_list = split(config.GetLine());
    
    // Reads bin bounds in micrometers and converts it to meters.
    BinBound_aer.resize(Nsize_section_aer + 1);
    for (int i = 0; i < Nsize_section_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]);

    if (option_process_aer["with_coagulation"])
      {
	config.SetSection("[options]");
	config.PeekValue("Coagulation_coefficient_file", coagulation_coefficient_file);
	if(!exists(coagulation_coefficient_file))
	{
	  cout<<"Error, Coagulation Coefficient data is not found, please computed it first!"<<endl;
	  exit(1);
	}
	Read_Coagulation_Coefficient(coagulation_coefficient_file);
      }
    else
      {
	coef_size=0;
      }

    // Bin corresponding to cutting diameter (for condensation).
    cutting_bin = 0;
    while (BinBound_aer(cutting_bin) < fixed_cutting_diameter &&
	   cutting_bin < Nsize_section_aer)
      cutting_bin++;

    cutting_bin=cutting_bin*Ncomposition_aer;

    // List of aerosol options to be passed to fortran routine.
    options_aer.resize(Noptions_aer);

    options_aer(0) = option_process_aer["with_coagulation"] ? 1 : 0;
    options_aer(1) = option_process_aer["with_condensation"] ? 1 : 0;
    options_aer(2) = option_process_aer["with_nucleation"] ? 1 : 0;

    if (aqueous_module=="no")
      options_aer(3) = 0;
    else if (aqueous_module=="vsrm")
      options_aer(3) = 1;
    else if (aqueous_module=="simple")// simple
      options_aer(3) = 2;
    else
      {
	options_aer(3) = 2;
	aqueous_module = "simple";
      }

    options_aer(4) = option_process_aer["with_in_cloud_scavenging"] ? 1 : 0;
    options_aer(5) = option_process_aer["with_kelvin_effect"] ? 1 : 0;
    options_aer(6) = option_process_aer["with_fixed_density"] ? 0 : 1;
    options_aer(7) = cutting_bin;
    if (sulfate_computation == "dynamic")
      options_aer(8) = 1;
    else if (sulfate_computation == "equilibrium")
      options_aer(8) = 0;
    else
      throw string("bad string \"") + sulfate_computation +
	string("\" for sulfate computation option,\n") +
	string("possibilities are dynamic or equilibrium.");

    if (dynamic_condensation_solver == "etr")
      options_aer(9) = 1;
    else if (dynamic_condensation_solver == "ros2")
      options_aer(9) = 2;
//     else if (dynamic_condensation_solver == "ebi")
//       options_aer(9) = 3;
    else
      throw string("bad string \"") + dynamic_condensation_solver +
	string("\" for dynamic c/e solver option,\n") +
	string("possibilities are etr, ros2.");
    
    if (redistribution_method == "mass-conserving")
      options_aer(10) = 1;
    else if (redistribution_method == "interpolation")
      options_aer(10) = 2;
    else if (redistribution_method == "euler-mass")
      options_aer(10) = 3;
    else if (redistribution_method == "euler-number" and with_normalisation == "no")
      options_aer(10) = 4;
    else if (redistribution_method == "hemen" and with_normalisation == "no")
      options_aer(10) = 5;
    else if (redistribution_method == "euler-coupled" and with_normalisation == "no")
      options_aer(10) = 6;
    else if (redistribution_method == "euler-number" and with_normalisation == "yes")
      options_aer(10) = 7;
    else if (redistribution_method == "hemen" and with_normalisation == "yes")
      options_aer(10) = 8;
    else if (redistribution_method == "euler-coupled" and with_normalisation == "yes")
      options_aer(10) = 9;
    else if (redistribution_method == "moving-diameter")
      options_aer(10) = 10;

    else
      throw string("bad string \"") + redistribution_method +
	string("\" for redistribution method option,\n") +
	string("possibilities are number-conserving, interpolation,") +
	string(" euler-mass, euler-number, hemen or euler-coupled.");
    
    if (nucleation_model == "binary")
      options_aer(11) = 0;
    else if (nucleation_model == "ternary")
      options_aer(11) = 1;
    else if (nucleation_model == "unary")
      options_aer(11) = 2;
    else
      throw string("bad string \"") + nucleation_model +
	string("\" for nucleation model option,\n") +
	string("possibilities are binary or ternary.");
    
    if (wet_diameter_estimation == "isoropia"
	|| wet_diameter_estimation == "isorropia")
      options_aer(12) = 0;
    else if (wet_diameter_estimation == "gerber")
      options_aer(12) = 1;
    else
      throw string("bad string \"") + wet_diameter_estimation +
	string("\" for wet diameter estimation option,\n") +
	string("possibilities are gerber or isorropia.");

    options_aer(13) = option_process_aer["with_oligomerization"] ? 1 : 0;
	
    if (Thermodynamic_model=="unifac")
      options_aer(14)=1;
    else if (Thermodynamic_model=="constant")
      options_aer(14)=0;
    else if (Thermodynamic_model=="ideal")
      options_aer(14)=2;
    else
      throw string("bad string \"") + Thermodynamic_model +
	string("\" organic thermodynamic model option,\n") +
	string("possibilities are unifac, constant or ideal.");
    options_aer(15) = option_process_aer["with_number_concentration"] ? 1 : 0;
	
    // Index of species subject to heterogeneous reactions.
    heterogeneous_reaction_index.resize(4);
    heterogeneous_reaction_index(0) = Model.GetSpeciesIndex("HO2");
    heterogeneous_reaction_index(1) = Model.GetSpeciesIndex("NO2");
    heterogeneous_reaction_index(2) = Model.GetSpeciesIndex("NO3");
    heterogeneous_reaction_index(3) = Model.GetSpeciesIndex("N2O5");

    if (option_process_aer["with_in_cloud_scavenging"])
      {
	if (option_process_aer["collect_wet_flux_aer"])
	  {
	    if (!Model.HasField("InCloudWetDepositionFlux_aer"))
	      throw string("if field \" InCloudWetDepositionFlux_aer \" ") +
		string("does not exist in Model, in cloud wet deposition ") +
		string("fluxes cannot be collected.");
	    if (Model.D4("InCloudWetDepositionFlux_aer").GetLength(2)
		!= Model.GetNy()
		|| Model.D4("InCloudWetDepositionFlux_aer").GetLength(3)
		!= Model.GetNx())
	      throw string("\" InCloudWetDepositionFlux_aer \" has") +
		string(" wrong dimension.");
	    if (option_process_aer["with_number_concentration"])
	      {
		if(!Model.HasField("InCloudWetDepositionFluxNumber_aer"))
		  throw string("if field \" InCloudWetDepositionFluxNumber_aer \" ") +
		    string("does not exist in Model, in cloud wet deposition ") +
		    string(" number fluxes cannot be collected.");
		if (Model.D3("InCloudWetDepositionFluxNumber_aer").GetLength(1)
		    != Model.GetNy()
		    || Model.D3("InCloudWetDepositionFluxNumber_aer").GetLength(2)
		    != Model.GetNx())
		  throw string("\" InCloudWetDepositionFluxNumber_aer \" has") +
		    string(" wrong dimension.");
	      }
			
	  }
	if (option_process_aer["collect_wet_flux"])
	  {
	    if (!Model.HasField("InCloudWetDepositionFlux"))
	      throw string("if field \" InCloudWetDepositionFlux \" ") +
		string("does not exist in Model, in cloud wet deposition ") +
		string("fluxes cannot be collected.");
	    if (Model.D3("InCloudWetDepositionFlux").GetLength(1)
		!= Model.GetNy()
		|| Model.D3("InCloudWetDepositionFlux").GetLength(2)
		!= Model.GetNx())
	      throw string("\" InCloudWetDepositionFlux \" has") +
		string(" wrong dimension.");
	  }

	// Check if in-cloud scavenging takes place
	// and if an aqueous-phase module is not activated.
	// In-cloud scavenging must be accompanied with
        // an aqueous-phase module (simple or VSRM module).
     	if (aqueous_module=="no")
	  throw string("Warning: activate an aqueous module ") +
	    string("(simple or VSRM) in your configuration file, ") +
	    string("polair3d.cfg. ") +
	    string("In-cloud scavenging is calculated ") +
            string("in the aqueous module. ") +
	    string("See aerosol.f and Aerosol_SCRAM.cxx");
      }

    if (option_process_aer["with_number_concentration"] == 0
	and option_process_aer["with_fixed_density"] == 0)
      throw string("Warning: activate \"With_number_concentration\"") +
	string(" or \"With_fixed_density\" or both in your") +
	string(" configuration file.");
    if (option_process_aer["with_number_concentration"] == 0
	and option_process_aer["with_fixed_density"] == 1
	and options_aer(10) > 2)
      throw string("Warning: The euler redistributions") +
	string(" cannot be used with the option \"With_fixed_density\".");
	
	

    BaseModuleParallel::Init(Model);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Partition along axis x.
    BuildPartition_x();
#endif
    /*** Display options ***/
    config.SetSection("[display]");
    // Should the configuration be displayed on screen?
    config.PeekValue("Show_configuration",
		     option_process_aer["show_configuration"]);
    
    if (option_process_aer["show_configuration"] && GetRank() == 0)
      DisplayConfiguration(Model);

  }
  

  //! Display the configuration.
  template<class T>
  template<class ClassModel>
  void Aerosol_SCRAM_H2O<T>::DisplayConfiguration(ClassModel& Model)
  {
    species_list = Model.GetSpeciesList();
    species_list_aer = Model.GetSpeciesList_aer();
    Ncomposition_aer=Model.GetNcomposition_aer();
    Nsize_section_aer= Model.GetNsize_section_aer();
    Nbin_aer = Model.GetNbin_aer();
    Ns_aer = Model.GetNs_aer();

    ConfigStream config_species(Model.GetSpeciesFile());

    map<string, double> parameter;
    map<string, double>::iterator iter;

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_aerosol_interact]");
    
    map<string, string> parameter2;
    map<string, string>::iterator iter2;
    cout << "Interacting aerosol species : ";
    string species;
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	parameter2[species] = config_species.GetElement();
      }
    for (int i = 0; i < Ns_aer; i++)
      {
	iter2 = parameter2.find(species_list_aer[i]);
	if(iter2 != parameter2.end())
	  {
	    int gas_index = Model.GetSpeciesIndex(iter2->second) + 1;
	    species_index_aerosol_interact(i) = gas_index;
	    cout << species_list_aer[i] << "<->" << iter2->second
		 << "(" << gas_index << ") ";
	  }
      }
    parameter2.clear();
    cout << endl;

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_cloud_interact]");

    cout << "Interacting cloud species : ";
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int gas_index = Model.GetSpeciesIndex(species) + 1;
	cout << species << "(" << gas_index << ") ";
      }
    cout << endl;

    // Aerosol species managed by isorropia equilibrium.
    config_species.SetSection("[isorropia_species]");

    cout << "Isorropia species : ";
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	cout << species << "(" << aer_index << ") ";
      }
    cout << endl;

    // Aerosol species managed by AEC equilibrium.
    config_species.SetSection("[aec_species]");

    cout << "AEC species : ";
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	cout << species << "(" << aer_index << ") ";
      }
    cout << endl;

    // Aerosol species managed by Pankow equilibrium.
    config_species.SetSection("[pankow_species]");

    cout << "Pankow species : ";
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	cout << species << "(" << aer_index << ") ";
      }
    cout << endl;

    // Primary organic aerosol species.
    config_species.SetSection("[poa_species]");
    cout << "Primary organic species : ";
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int aer_index = Model.GetSpeciesIndex_aer(species) + 1;
	cout << species << "(" << aer_index << ") ";
      }
    cout << endl;
	
    // Mineral dust aerosol species.
    config_species.SetSection("[mineral_dust_species]");
    species = config_species.GetElement();
    md_species = Model.GetSpeciesIndex_aer(species) + 1;
    cout << "Mineral dust species : " << species
	 << "(" << md_species << ")" << endl;

    // Black carbon aerosol species.
    config_species.SetSection("[black_carbon_species]");
    species = config_species.GetElement();
    bc_species = Model.GetSpeciesIndex_aer(species) + 1;
    cout << "Black carbon species : " << species
	 << "(" << bc_species << ")" << endl;

    if (with_adaptive_time_step)
      cout << "Module: with_adaptive_time_step_for_gas_chemistry" << endl;

    cout << "Module: photolysis option: ";
    if (option_photolysis == 1)
      cout << "Tabulation by SPACK" << endl;
    else if (option_photolysis == 2)
      cout << "Reading the binary files obtained by FastJ" << endl;
    else
      throw string("Wrong number is given. Please put 1 or 2");

    cout << "Module: Fixed_aerosol_density = " << FixedDensity_aer << endl;
    cout << "Module: Lwc_cloud_threshold = " << lwc_cloud_threshold << endl;
    
    if (option_process_aer["with_coagulation"])
      {
	cout << "Module: reading coagulation coefficients ... ";
	cout << "done" << endl;
      }

    if (option_process_aer["with_coagulation"])
      cout << "Module: with coagulation" << endl;
    if (option_process_aer["with_condensation"])
      cout << "Module: with condensation" << endl;
    if (option_process_aer["with_nucleation"])
      cout << "Module: with nucleation" << endl;
    cout << "Module: cloud chemistry: " << aqueous_module << endl;
    if (option_process_aer["with_in_cloud_scavenging"])
      cout << "Module: with in-cloud scavenging" << endl;
    if (option_process_aer["with_heterogeneous_reactions"])
      cout << "Module: with heterogeneous reactions" << endl;
    if (option_process_aer["with_kelvin_effect"])
      cout << "Module: with kelvin effect" << endl;
    if (option_process_aer["with_fixed_density"])
      cout << "Module: with fixed density" << endl;
    if (cutting_bin == 0)
      cout << "Module: all bins at dynamic." << endl;
    else if (cutting_bin == Nbin_aer)
      cout << "Module: all bins at equilibrium." << endl;
    else
      cout << "Module: bins from 1 to " << cutting_bin <<
	"(included) at equilibrium." << endl;

    cout << "Module: dynamic condensation solver: " <<
      dynamic_condensation_solver << endl;
    
    cout << "Module: redistribution method: " <<
      redistribution_method << endl;
    
    if (option_process_aer["with_nucleation"])
      {
	cout << "Module: nucleation model: " <<
	  nucleation_model << endl;
	// deschams
	cout << "K "<< K_nucl_factor <<", P " << P_nucl_factor<< endl; 
      }
    cout << "Module: wet diameter estimation: " <<
      wet_diameter_estimation << endl;
    
    if (option_process_aer["with_oligomerization"])
      cout << "Module: with oligomerization" << endl;

    cout << "Module: Organic Thermodynamic Model :" <<
      Thermodynamic_model << endl;

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    cout << "Module: without sea salt in isorropia" << endl;
#else
    cout << "Module: with sea salt in isorropia" << endl;
#endif
  }

    
  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i")
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetSource_i()
    <li> D4("PhotolysisRate_i")
    <li> D3("Attenuation_f")
    <li> D3("SpecificHumidity_f")
    <li> D3("Temperature_f")
    <li> D3("Pressure_f")
    <li> GetSource_f()
    <li> D4("PhotolysisRate_f")
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D4("WetDiameter_aer")
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SCRAM_H2O<T>::Forward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

    Forward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
	    Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
	    Model.D3("Pressure_i"), Model.GetSource_i(),
	    Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
	    Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
	    Model.D3("Temperature_f"), Model.D3("Pressure_f"),
	    Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
	    Model.GetGridXArray1D(), Model.GetGridYArray1D(),
	    Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
	    Model.D4("WetDiameter_aer"), Model.GetConcentration_aer(),
	    Model.D3("pH"), Model.GetNumberConcentration_aer());
	
  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D2("Rain_i")
    <li> GetLayerInterface()
    <li> GetConcentration_aer()
    <li> D3("pH")
    <li> GetNumberConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SCRAM_H2O<T>::Forward_aer(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();
	
    Data<T, 3> in_cloud_wet_flux(Ns, Model.GetNy(), Model.GetNx());
    Data<T, 4> in_cloud_wet_flux_aer(Ns_aer, Nbin_aer, Model.GetNy(),
				     Model.GetNx());
    Data<T, 3> in_cloud_wet_flux_number_aer(Nbin_aer, Model.GetNy(), Model.GetNx());

    Forward_aer(T(date_i.GetNumberOfSeconds()),
		Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
		Model.D3("Pressure_i"),	T(date_f.GetSecondsFrom(date_i)),
		Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
		Model.D2("Rain_i"), Model.GetLayerInterface(),
		Model.GetConcentration_aer(), in_cloud_wet_flux,
		in_cloud_wet_flux_aer,
		Model.D3("pH"),
		Model.GetNumberConcentration_aer(),
		in_cloud_wet_flux_number_aer);
	
    if (option_process_aer["with_in_cloud_scavenging"])
      {
	if (option_process_aer["collect_wet_flux"])
	  for (int s = 0; s < Ns; s++)
	    if(Model.HasScavenging(s))
	      for (int j = 0; j < Model.GetNy(); j++)
		for (int i = 0; i < Model.GetNx(); i++)
		  Model.D3("InCloudWetDepositionFlux")
		    (Model.ScavengingIndex(s), j, i) =
		    in_cloud_wet_flux(s, j, i);
	if (option_process_aer["collect_wet_flux_aer"])
	  for (int b = 0; b < Nbin_aer; b++)
	    if(Model.HasScavenging_aer(b))
	      for (int j = 0; j < Model.GetNy(); j++)
		for (int i = 0; i < Model.GetNx(); i++)
		  {
		    for (int s = 0; s < Ns_aer; s++)
		      Model.D4("InCloudWetDepositionFlux_aer")
			(s, Model.ScavengingIndex_aer(b), j, i) =
			in_cloud_wet_flux_aer(s, b, j, i);
		    if(option_process_aer["with_number_concentration"])
		      Model.D3("InCloudWetDepositionFluxNumber_aer")
			(Model.ScavengingIndex_aer(b), j, i) =
			in_cloud_wet_flux_number_aer(b, j, i);
		  }
		  
      }
		
  }


  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param next_time time at the end of the time step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param WetDiameter_aer Aerosol wet diameters at the beginning of the time
    step.
    \param Concentration_aer aerosol concentrations.
  */
  template<class T>
  void Aerosol_SCRAM_H2O<T>::Forward(T current_time,
				      Data<T, 3>& Attenuation_i,
				      Data<T, 3>& SpecificHumidity_i,
				      Data<T, 3>& Temperature_i,
				      Data<T, 3>& Pressure_i,
				      Data<T, 4>& Source_i,
				      Data<T, 4>& PhotolysisRate_i,
				      T delta_t,
				      Data<T, 3>& Attenuation_f,
				      Data<T, 3>& SpecificHumidity_f,
				      Data<T, 3>& Temperature_f,
				      Data<T, 3>& Pressure_f,
				      Data<T, 4>& Source_f,
				      Data<T, 4>& PhotolysisRate_f,
				      Array<T, 1>& Longitude,
				      Array<T, 1>& Latitude,
				      Data<T, 4>& Concentration,
				      Data<T, 3>& LiquidWaterContent_i,
				      Data<T, 4>& WetDiameter_aer,
				      Data<T, 5>& Concentration_aer,
				      Data<T, 3>& pH,
				      Data<T, 4>& NumberConcentration_aer)
    
  {

    int icld = options_aer(3);
    int iheter = option_process_aer["with_heterogeneous_reactions"] ? 1 : 0;

    int first_index_along_x, last_index_along_x;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T,4> OrdConc;
    Array<T,5> OrdConc_aer;
    Array<T,4> OrdNumConc_aer;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    ScatterSlice_x_MPI(Concentration_aer, OrdConc_aer);
    if (option_process_aer["with_number_concentration"])
      ScatterSlice_x_MPI(NumberConcentration_aer, OrdNumConc_aer);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif
    
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
      for (int k = 0; k < Concentration.GetLength(1); k++)
	for (int j = 0; j < Concentration.GetLength(2); j++)
	  {
            Array<T, 1> Photolysis1D(Nr_photolysis);
            Array<T, 1> Photolysis1D_f(Nr_photolysis);
            Array<T, 1> Concentration1D(Ns);
            Array<T, 1> Source1D(Ns_source);
            Array<T, 1> Source1D_f(Ns_source);
            Array<T, 1> WetDiameter_aer1D(Nbin_aer);
            Array<T, 2> Concentration_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> NumberConcentration_aer1D(Nbin_aer);
            int s, b;

	    for (s = 0; s < Nr_photolysis; s++)
	      {
		Photolysis1D(s) = PhotolysisRate_i()(s, k, j, i);
		Photolysis1D_f(s) = PhotolysisRate_f()(s, k, j, i);
	      }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		Concentration1D(s) = OrdConc(i, k, j, s);
#else
		Concentration1D(s) = Concentration(s, k, j, i);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                WetDiameter_aer1D(b) = WetDiameter_aer(b, k, j, i);
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    Concentration_aer2D(s, b) = OrdConc_aer(i, b, k, j, s);
#else
                    Concentration_aer2D(s, b) = Concentration_aer(s, b, k, j, i);
#endif

                  }
              }

	    if (option_process_aer["with_number_concentration"])
	      for (b = 0; b < Nbin_aer; b++)
		{
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		  // Rank of indices must match the one defined in
		  // ScatterSlice_x_MPI.

		  if (isnan(OrdNumConc_aer(i, k, j, b)))
		    {
		      cout <<"isnan "<< b <<" "<< i <<" "<<  j<<" " <<  k;
						
		    }

					
		  NumberConcentration_aer1D(b) = OrdNumConc_aer(i, k, j, b);
#else

		  NumberConcentration_aer1D(b) = NumberConcentration_aer(b, k, j, i);
#endif
		}
	    else
	      NumberConcentration_aer1D = 0.0;
			
	    if (k < Nz_source)
	      {
		for (s = 0; s < Ns_source ; s++)
		  {
		    Source1D(s) = Source_i()(s, k, j, i);
		    Source1D_f(s) = Source_f()(s, k, j, i);
		  }
	      }
	    else
	      {
		Source1D = 0;
		Source1D_f = 0;
	      }
	    Forward(current_time, Attenuation_i(k, j, i),
		    SpecificHumidity_i(k, j, i), Temperature_i(k, j, i),
		    Pressure_i(k, j, i), Source1D, Photolysis1D,
		    delta_t, Attenuation_f(k, j, i),
		    SpecificHumidity_f(k, j, i), Temperature_f(k, j, i),
		    Pressure_f(k, j, i), Source1D_f, Photolysis1D_f,
		    Longitude(i), Latitude(j), Concentration1D, 
                    LiquidWaterContent_i(k, j, i), WetDiameter_aer1D, 
                    Concentration_aer2D, pH(k, j, i),
		    NumberConcentration_aer1D);

		    
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdConc(i, k, j, s) = Concentration1D(s);
#else
		Concentration(s, k, j, i) = Concentration1D(s);
#endif
	      }

	    for (s = 0; s < Ns_aer; s++)
              {
                for (b = 0; b < Nbin_aer; b++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    OrdConc_aer(i, b, k, j, s) = Concentration_aer2D(s, b);
#else
                    Concentration_aer(s, b, k, j, i) = Concentration_aer2D(s, b);
#endif
                  }
              }
	    if (option_process_aer["with_number_concentration"])
	      for (b = 0; b < Nbin_aer ; b++)
		{
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		  // Rank of indices must match the one defined in
		  // ScatterSlice_x_MPI.
		  OrdNumConc_aer(i, k, j, b) = NumberConcentration_aer1D(b);
#else
		  NumberConcentration_aer(b, k, j, i) = NumberConcentration_aer1D(b);
#endif
		} 


	  } // loop j
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
    GatherSlice_x_MPI(OrdConc_aer, Concentration_aer);
    if(option_process_aer["with_number_concentration"])
      GatherSlice_x_MPI(OrdNumConc_aer, NumberConcentration_aer);
#endif
  } // Forward


  
  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param next_time time at the end of the time step.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param Rain_i rain rate at the beginning of the time step.
    \param VerticalInterface height of vertical interfaces in meter.
    \param Concentration_aer aerosol concentrations.
    \param aerosol pH.
    \param NumberConcentration_aer aerosol concentrations.
  */
  template<class T>
  void Aerosol_SCRAM_H2O<T>::Forward_aer(T current_time,
					  Data<T, 3>& SpecificHumidity_i,
					  Data<T, 3>& Temperature_i,
					  Data<T, 3>& Pressure_i,
					  T delta_t,
					  Data<T, 4>& Concentration,
					  Data<T, 3>&
					  LiquidWaterContent_i,
					  Data<T, 2>& Rain_i,
					  Array<T, 1>& VerticalInterface,
					  Data<T, 5>& Concentration_aer,
					  Data<T, 3>&
					  InCloudWetDepositionFlux,
					  Data<T, 4>&
					  InCloudWetDepositionFlux_aer,
					  Data<T, 3>& pH,
					  Data<T, 4>& NumberConcentration_aer,
					  Data<T, 3>&
					  InCloudWetDepositionFluxNumber_aer)
  {
    int first_index_along_x, last_index_along_x;

    double max_c_aer=0.0;
    Array<T, 1> Index_max(5);
    Index_max=0;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T, 4> OrdConc;
    Array<T, 5> OrdConc_aer;

    // This array will store the number concentration data
    Array<T, 4> OrdNumConc_aer;
    Array<T, 3> OrdWetDepositionFluxNumber_aer;
	
    Array<T, 3> OrdWetDepositionFlux;
    Array<T, 4> OrdWetDepositionFlux_aer;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    ScatterSlice_x_MPI(Concentration_aer, OrdConc_aer);
    if (option_process_aer["with_number_concentration"])
      {
	ScatterSlice_x_MPI(NumberConcentration_aer, OrdNumConc_aer);
	ScatterSlice_x_MPI(InCloudWetDepositionFluxNumber_aer, OrdWetDepositionFluxNumber_aer);
      }
    ScatterSlice_x_MPI(InCloudWetDepositionFlux, OrdWetDepositionFlux);
    ScatterSlice_x_MPI(InCloudWetDepositionFlux_aer, OrdWetDepositionFlux_aer);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif
	
    InCloudWetDepositionFlux.SetZero();
    InCloudWetDepositionFlux_aer.SetZero();
    if (option_process_aer["with_number_concentration"])
      InCloudWetDepositionFluxNumber_aer.SetZero();
    T total_number_in=0.0;
    T total_number_out=0.0;	
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)    
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
    {
      for (int k = 0; k < Concentration.GetLength(1); k++)
	for (int j = 0; j < Concentration.GetLength(2); j++)
	  {
	    Array<T, 1> Concentration1D(Ns);
            Array<T, 2> Concentration_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> NumberConcentration_aer1D(Nbin_aer);
            Array<T, 1> InCloudWetDepositionFlux1D(Ns);
            Array<T, 2> InCloudWetDepositionFlux_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> InCloudWetDepositionFluxNumber_aer1D(Nbin_aer);
            Array<T, 1> CurrentVerticalInterface(2);
	    int s, b;
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
	    int tid = omp_get_thread_num();
	    
//	    if (tid == 0) 
		cout<<"\r"<<tid<<"=Grid point ("<<i<<" ,"<<k<<" ,"<<j<<")"<<flush;
#endif	    
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		InCloudWetDepositionFlux1D(s) = OrdWetDepositionFlux(i, j, s);
#else
		InCloudWetDepositionFlux1D(s) = InCloudWetDepositionFlux(s, j, i);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    InCloudWetDepositionFlux_aer2D(s, b) = OrdWetDepositionFlux_aer(i, b, j, s);
#else
                    InCloudWetDepositionFlux_aer2D(s, b) = InCloudWetDepositionFlux_aer(s, b, j, i);
#endif
                  }
              }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		Concentration1D(s) = OrdConc(i, k, j, s);
#else
		Concentration1D(s) = Concentration(s, k, j, i);
#endif
	      }
            double total_ms=0;
	    double total_nb=0;
            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    Concentration_aer2D(s, b) = OrdConc_aer(i, b, k, j, s);
#else
                    Concentration_aer2D(s, b) = Concentration_aer(s, b, k, j, i);
		    if(s< Ns_aer-1)
		      total_ms+=Concentration_aer2D(s, b);
#endif
                  }
              }
			
	    if (option_process_aer["with_number_concentration"])
	      {
		for (b = 0; b < Nbin_aer ; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    InCloudWetDepositionFluxNumber_aer1D(b) =
		      OrdWetDepositionFluxNumber_aer(i, j, b);
#else
                    InCloudWetDepositionFluxNumber_aer1D(b) =
		      InCloudWetDepositionFluxNumber_aer(b, j, i);
#endif
                  }
              
		for (b = 0; b < Nbin_aer; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    NumberConcentration_aer1D(b) = OrdNumConc_aer(i, k, j, b);
#else
		    NumberConcentration_aer1D(b) = NumberConcentration_aer(b, k, j, i);
		    total_nb+=NumberConcentration_aer1D(b);
		    total_number_in+= NumberConcentration_aer1D(b);		    
#endif
		  }
	      }
	    else
	      {
		InCloudWetDepositionFluxNumber_aer1D = 0.0;
		NumberConcentration_aer1D = 0.0;
	      }
			
			
            CurrentVerticalInterface(0) = VerticalInterface(k);
            CurrentVerticalInterface(1) = VerticalInterface(k+1);

            Array<T, 1> LiquidWaterContent_1D(Concentration.GetLength(1));
            int nfoglay = 0;
            T lwc_avg = 0.0;
            T heightfog = 0.0;
            int ifog = 0;
            for (int kk = 0; kk < Concentration.GetLength(1); kk++)
              LiquidWaterContent_1D(kk) = LiquidWaterContent_i(kk, j, i)
		* 1000. * Pressure_i(kk, j, i)
		/ 101325.0 * 28.97 / 0.082
		/ Temperature_i(kk, j, i);

            FogSettling(LiquidWaterContent_1D,
                        lwc_cloud_threshold, VerticalInterface,
                        lwc_avg, nfoglay, heightfog);
            if (k < nfoglay) ifog = 1;
            Forward_aer(current_time,
			SpecificHumidity_i(k, j, i), Temperature_i(k, j, i),
			Pressure_i(k, j, i),
			delta_t, Concentration1D,
			LiquidWaterContent_i(k, j, i), Rain_i(j,i),
			CurrentVerticalInterface, Concentration_aer2D,
			InCloudWetDepositionFlux1D,
			InCloudWetDepositionFlux_aer2D,
			pH(k, j, i), lwc_avg, heightfog, ifog,
			NumberConcentration_aer1D,
			InCloudWetDepositionFluxNumber_aer1D);
			
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdConc(i, k, j, s) = Concentration1D(s);
#else
		Concentration(s, k, j, i) = Concentration1D(s);
#endif
	      }
	    double total_ms_a=0;
	    double total_nb_a=0;
	    for (s = 0; s < Ns_aer; s++)
	      {
                for (b = 0; b < Nbin_aer; b++)
                  {		
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdConc_aer(i, b, k, j, s) = Concentration_aer2D(s, b);
#else
		    Concentration_aer(s, b, k, j, i) = Concentration_aer2D(s, b);

		    if(s< Ns_aer-1)
		    {
		      total_ms_a+=Concentration_aer2D(s, b);
		      if(max_c_aer<Concentration_aer2D(s, b))
		      {
			max_c_aer=Concentration_aer2D(s, b);
			Index_max(0)=s;
			Index_max(1)=b;
			Index_max(2)=k;
			Index_max(3)=j;
			Index_max(4)=i;
		      }
		    }
#endif
                  }
              }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdWetDepositionFlux(i, j, s) += InCloudWetDepositionFlux1D(s);
#else
		InCloudWetDepositionFlux(s, j, i) += InCloudWetDepositionFlux1D(s);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    OrdWetDepositionFlux_aer(i, b, j, s) +=
		      InCloudWetDepositionFlux_aer2D(s, b);
#else
                    InCloudWetDepositionFlux_aer(s, b, j, i) +=
		      InCloudWetDepositionFlux_aer2D(s, b);
						
#endif
                  }
              }
			
	    if (option_process_aer["with_number_concentration"])
	      {
		for (b = 0; b < Nbin_aer; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdWetDepositionFluxNumber_aer(i, j, b) +=
		      InCloudWetDepositionFluxNumber_aer1D(b);
#else
		    InCloudWetDepositionFluxNumber_aer(b, j, i) +=
		      InCloudWetDepositionFluxNumber_aer1D(b);
			
#endif
		  }
				
		for (b = 0; b < Nbin_aer ; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdNumConc_aer(i, k, j, b) = NumberConcentration_aer1D(b);
#else
		    NumberConcentration_aer(b, k, j, i) = NumberConcentration_aer1D(b);
		    total_nb_a+=NumberConcentration_aer1D(b);
		    total_number_out+= NumberConcentration_aer1D(b);
#endif
		  }
	      }
             if(total_nb_a*total_ms_a==0.0&&total_nb_a!=total_ms_a)
	     {
	      cout<<k<<","<<j<<","<<i<<" n0="<<total_nb<<" n1="<<
	      total_nb_a<<" m0="<<total_ms<<" m1="<<total_ms_a<<endl;
	     }
	     else
	     {
// 	       double change=abs(total_ms_a-total_ms)/total_ms*100.0;
// 	       if(change>50.0)
// 		cout<<k<<","<<j<<","<<i<<" n0="<<total_nb<<" n1="<<
// 		total_nb_a<<" m0="<<total_ms<<" ->m1="<<total_ms_a<<
// 		" = "<<change<<" %"<<endl;
	     }
	  }// loop k,j,i
    }
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
    GatherSlice_x_MPI(OrdConc_aer, Concentration_aer);
    if(option_process_aer["with_number_concentration"])
      {
	GatherSlice_x_MPI(OrdNumConc_aer, NumberConcentration_aer);
	GatherSlice_x_MPI(OrdWetDepositionFluxNumber_aer, InCloudWetDepositionFluxNumber_aer);
      }
    GatherSlice_x_MPI(OrdWetDepositionFlux, InCloudWetDepositionFlux);
    GatherSlice_x_MPI(OrdWetDepositionFlux_aer, InCloudWetDepositionFlux_aer);
#else
      cout<<"max without water="<<max_c_aer<<"("<<species_list_aer[Index_max(0)]<<","
      <<Index_max(1)<<","<<Index_max(2)<<","
      <<Index_max(3)<<","<<Index_max(4)<<")"<<endl;
      cout<<"NumberConcentration_aer";
      NumberConcentration_aer.PrintInfo();
      Array<int, 1> Index2(4);
      Index2=NumberConcentration_aer.GetMaxIndex();
      cout<<"max_index("<<Index2(0)<<","<<Index2(1)<<","<<Index2(2)<<","
      <<Index2(3)<<")"<<endl;
      cout<<"total_number_in="<<total_number_in<<endl;
      cout<<"total_number_out="<<total_number_out<<endl;      
#endif	
  } // Forward_aer






  
  template<class T>
  void Aerosol_SCRAM_H2O<T>::Forward(T current_time,
				      T& attenuation,
				      T& humid,
				      T& temperature,
				      T& pressure,
				      Array<T, 1>& source,
				      Array<T, 1>& photolysis,
				      T delta_t,
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
				      Array<T, 1>& number_concentration_aer)
  {
    int icld = options_aer(3);
    int iheter = option_process_aer["with_heterogeneous_reactions"] ? 1 : 0;
    int inum = options_aer(15);
    int idens = options_aer(6);
	
    _chem(&Ns, &Nr,
	  &Nr_photolysis, photolysis_reaction_index.data(),
	  &Ns_source, source_index.data(),
	  ConversionFactor.data(), ConversionFactorJacobian.data(),
	  &Nz_source,&lwc_cloud_threshold, molecular_weight.data(),
	  &current_time, &attenuation,
	  &humid, &temperature,
	  &pressure, source.data(),
	  photolysis.data(), &delta_t, &attenuation_f,
	  &humid_f, &temperature_f,
	  &pressure_f, source_f.data(),
	  photolysis_f.data(), &Ncycle, &lon,
	  &lat, concentration.data(),
	  &icld, &iheter, &Ns_aer, &Nbin_aer,&Ncomposition_aer,
	  &liquid_water_content,
	  BinBound_aer.data(), &FixedDensity_aer,
	  wet_diameter_aer.data(),
	  heterogeneous_reaction_index.data(),
	  concentration_aer.data(),
	  &option_adaptive_time_step, &adaptive_time_step_tolerance,
	  &min_adaptive_time_step, &option_photolysis, &jBiPER, &kBiPER,
	  &inum, &idens, number_concentration_aer.data(),
	  mass_density_aer.data());
	
  }

  template<class T>
  void Aerosol_SCRAM_H2O<T>::Forward_aer(T current_time,
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
                                          T& lwc_avg,
                                          T& heightfog,
                                          int& ifog,
					  Array<T, 1>& number_concentration_aer,
					  Array<T, 1>& incloudwetdepositionfluxnumber_aer)
  {
     __zaerosolscram_MOD_aerosol(&Ns, &Ns_isorropia, &Ns_aec, &Ns_pankow,
 	     &Ns_poa, &lwc_cloud_threshold, &current_time,
 	     &specifichumidity, &temperature,
 	     &pressure, &delta_t, concentration.data(),
 	     &Noptions_aer, options_aer.data(), &Ns_aer,
	     &Nsize_section_aer, &Nbin_aer, &Ncomposition_aer, &Ngroup_aer,
 	     &Ncycle_aer, &liquidwatercontent, &rain,
 	     BinBound_aer.data(), &FixedDensity_aer,
	     aerosol_species_group_relation.data(), composition_bounds.data(),
	     Ncoefficient.data(),index_first.data(),index_second.data(),coefficient.data(),
 	     &coef_size,CurrentVerticalInterface.data(),
 	     concentration_aer.data(),incloudwetdepositionflux.data(),
 	     incloudwetdepositionflux_aer.data(), &ph,
 	     saturation_pressure.data(), partition_coefficient.data(),
 	     vaporization_enthalpy.data(), accomodation_coefficient.data(),
 	     surface_tension.data(), saturation_pressure_mass.data(),
 	     saturation_pressure_torr.data(),
 	     deliquescence_relative_humidity.data(),
 	     molecular_weight_aer.data(), molecular_diameter_aer.data(),
 	     collision_factor_aer.data(), mass_density_aer.data(),
 	     species_index_aerosol_interact.data(),
 	     &isorropia_species.front(), &aec_species.front(),
 	     &pankow_species.front(),&poa_species.front(),&md_species,&bc_species,
 	     &Ns_cloud_interact, &species_index_cloud_interact.front(),
 	     &lwc_avg, &heightfog, &ifog, number_concentration_aer.data(),
 	     incloudwetdepositionfluxnumber_aer.data(),&conserving_mass_tolerance,
 	     &P_nucl_factor, &K_nucl_factor);
	
  }

  //! Checks whether a given field is required by this module.
  /*! Checks whether a given field must be available in the underlying model
    for this module to perform properly.
    \param field the field name.
    \return True is the field is required, false otherwise.
  */
  template<class T>
  bool Aerosol_SCRAM_H2O<T>::IsRequired(string field)
  {
    if (field == "wet_diameter_aer") return true;
    if (field == "rain") return true;
    return false;
  }


  //! Checks whether a given field is computed by this module.
  /*! Checks whether a given field is computed by this module and updated in
    the underlying model.
    \param field the field name.
    \return True is the field is computed, false otherwise.
  */
  template<class T>
  bool Aerosol_SCRAM_H2O<T>::IsComputed(string field)
  {
    if (field == "pH") return true;
    return false;
  }

  template<class T>
  void Aerosol_SCRAM_H2O<T>::FogSettling(Array<T, 1>& LiquidWaterContent,
                                          T& lwc_cloud_threshold,
                                          Array<T, 1>& VerticalInterface,
                                          T& lwc_avg,
                                          int& nfoglay,
                                          T& heightfog)
  {
    // Fog settling
    int size_lwc = LiquidWaterContent.shape()[0];  //sizeof(LiquidWaterContent);
	
    if (LiquidWaterContent(0) > lwc_cloud_threshold)
      {
	int indok = 1;
	for (int k = 0; k < size_lwc; k++)
	  {
	    if ((LiquidWaterContent(k) > lwc_cloud_threshold) && (indok == 1))
	      {
				
		heightfog = VerticalInterface(k+1);
		nfoglay = k+1;
		lwc_avg += LiquidWaterContent(k);
			
	      }
	    else
	      indok = 0;
	  }
		
	if (nfoglay > 0) 
	  lwc_avg /= nfoglay;
      }
   
  }

  //! Read coagulation coefficient repartition data 
  template<class T>
  void Aerosol_SCRAM_H2O<T>::Read_Coagulation_Coefficient(const string &input_file)
  {
    NcFile fnc(input_file.c_str(), NcFile::ReadOnly);

    if (! fnc.is_valid())
      throw string("Invalid NetCDF file pointer.");

    // Check some dimensions.
    if (Nbin_aer != int(fnc.get_dim("Nsize")->size()))
      throw string("Values of Nsize differ between NetCDF file and current object.");

    if (Nsize_section_aer != int(fnc.get_dim("Nb")->size()))
      throw string("Values of Nb differ between NetCDF file and current object.");

    if (Ncomposition_aer!= int(fnc.get_dim("Nc")->size()))
      throw string("Values of Nc differ between NetCDF file and current object.");

    int total_size=0;
    Ncoefficient.resize(Nbin_aer);
    vector<int> index_first_;
    vector<int> index_second_;
    vector<double> coefficient_;

    for (int i = 0; i < Nbin_aer; i++)
      {
        string dim_name("Ncoef_" + to_str(i));
        int Ncoef = int(fnc.get_dim(dim_name.c_str())->size());

	Ncoefficient(i)=Ncoef;
	total_size+=Ncoef;

	Array<int, 1> id1_tmp;
	Array<int, 1> id2_tmp;
	Array<double, 1> coef_tmp;

	id1_tmp.resize(Ncoef);
	id2_tmp.resize(Ncoef);
	coef_tmp.resize(Ncoef);

	NcVar *var;
	string var_name;

        var_name = "index1_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(id1_tmp.data(), Ncoef);

        var_name = "index2_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(id2_tmp.data(), Ncoef);

        var_name = "coef_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(coef_tmp.data(), Ncoef);

	for (int j = 0; j < Ncoef; j++)
	{
	  index_first_.push_back(id1_tmp(j));
	  index_second_.push_back(id2_tmp(j));
	  coefficient_.push_back(coef_tmp(j));
	}
	
      }

      index_first.resize(total_size);
      index_second.resize(total_size);
      coefficient.resize(total_size);

      for(int i=0; i<total_size; i++)
      {
	index_first(i)=index_first_[i];
	index_second(i)=index_second_[i];
	coefficient(i)=coefficient_[i];
      }
     coef_size=total_size;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SCRAM_H2O_CXX
#endif


