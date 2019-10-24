// Copyright (C) 2006-2018, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Shupeng Zhu
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

// This file is part of the Eulerian model Polair3D.

#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DAEROSOL_CXX


#include "Polair3DAerosol.hxx"
#include "fastJX.hxx"
#include "Common.cxx"

#include <string.h>
#include <blitz/array/et.h>

#define PI 3.141592653589793115997963468544185161590576171875

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Polair3DAerosol(string config_file):
    Polair3DChemistry < T, ClassAdvection, ClassDiffusion,
    ClassChemistry > (config_file)
  {

    /*** Managed data ***/

    this->option_manage["temperature"] = true;
    this->option_manage["pressure"] = true;
    this->option_manage["initial_condition_aer"] = true;
    this->option_manage["surface_temperature"] = true;
    this->option_manage["surface_pressure"] = true;
    this->option_manage["first_level_wind_module"] = true;
    this->option_manage["liquid_water_content"] = true;
    this->option_manage["snow_height"] = true;
    this->option_manage["boundary_condition_aer"] = true;
    this->option_manage["deposition_velocity_aer"] = true;
    this->option_manage["scavenging_coefficient_aer"] = true;
    this->option_manage["surface_emission_aer"] = true;
    this->option_manage["volume_emission_aer"] = true;

    /*** Pointers to 2D data ***/

    this->D2_map["SurfaceTemperature"] = &SurfaceTemperature_f;
    this->D2_map["SurfaceTemperature_i"] = &SurfaceTemperature_i;
    this->D2_map["SurfaceTemperature_f"] = &SurfaceTemperature_f;

    this->D2_map["SurfacePressure"] = &SurfacePressure_f;
    this->D2_map["SurfacePressure_i"] = &SurfacePressure_i;
    this->D2_map["SurfacePressure_f"] = &SurfacePressure_f;

    this->D2_map["FirstLevelWindModule"] = &FirstLevelWindModule_f;
    this->D2_map["FirstLevelWindModule_i"] = &FirstLevelWindModule_i;
    this->D2_map["FirstLevelWindModule_f"] = &FirstLevelWindModule_f;

    this->D2_map["SnowHeight"] = &SnowHeight_f;
    this->D2_map["SnowHeight_i"] = &SnowHeight_i;
    this->D2_map["SnowHeight_f"] = &SnowHeight_f;

    /*** Pointers to 3D data ***/

    this->D3_map["LiquidWaterContent"] = &LiquidWaterContent_i;
    this->D3_map["LiquidWaterContent_i"] = &LiquidWaterContent_i;

    this->D3_map["DepositionVelocity_aer"] = &DepositionVelocity_aer_f;
    this->D3_map["DepositionVelocity_aer_i"] = &DepositionVelocity_aer_i;
    this->D3_map["DepositionVelocity_aer_f"] = &DepositionVelocity_aer_f;

    this->D3_map["pH"] = &pH;

    this->D3_map["NumberBoundaryCondition_z_aer"]= &NumberBoundaryCondition_z_aer_i;
    this->D3_map["NumberBoundaryCondition_z_aer_i"]= &NumberBoundaryCondition_z_aer_i;

    this->D3_map["NumberSurfaceEmission_aer"] = &NumberSurfaceEmission_aer_f;
    this->D3_map["NumberSurfaceEmission_aer_i"] = &NumberSurfaceEmission_aer_i;
    this->D3_map["NumberSurfaceEmission_aer_f"] = &NumberSurfaceEmission_aer_f;

    this->D3_map["DryDepositionFluxNumber_aer"] = &DryDepositionFluxNumber_aer;
    this->D3_map["WetDepositionFluxNumber_aer"] = &WetDepositionFluxNumber_aer;
    this->D3_map["InCloudWetDepositionFluxNumber_aer"]
      = &InCloudWetDepositionFluxNumber_aer;

    this->D3_map["CloudOpticalDepth"] = &CloudOpticalDepth_i;
    this->D3_map["CloudOpticalDepth_i"] = &CloudOpticalDepth_i;
    this->D3_map["IceOpticalDepth"] = &IceOpticalDepth_i;
    this->D3_map["IceOpticalDepth_i"] = &IceOpticalDepth_i;

    /*** Pointers to 4D data ***/

    this->D4_map["BoundaryCondition_z_aer"] = &BoundaryCondition_z_aer_i;
    this->D4_map["BoundaryCondition_z_aer_i"] = &BoundaryCondition_z_aer_i;

    this->D4_map["NumberBoundaryCondition_x_aer"]= &NumberBoundaryCondition_x_aer_i;
    this->D4_map["NumberBoundaryCondition_x_aer_i"]= &NumberBoundaryCondition_x_aer_i;

    this->D4_map["NumberBoundaryCondition_y_aer"]= &NumberBoundaryCondition_y_aer_i;
    this->D4_map["NumberBoundaryCondition_y_aer_i"]= &NumberBoundaryCondition_y_aer_i;
	

    this->D4_map["ScavengingCoefficient_aer"] = &ScavengingCoefficient_aer_f;
    this->D4_map["ScavengingCoefficient_aer_i"]
      = &ScavengingCoefficient_aer_i;
    this->D4_map["ScavengingCoefficient_aer_f"]
      = &ScavengingCoefficient_aer_f;

    this->D4_map["SurfaceEmission_aer"] = &SurfaceEmission_aer_f;
    this->D4_map["SurfaceEmission_aer_i"] = &SurfaceEmission_aer_i;
    this->D4_map["SurfaceEmission_aer_f"] = &SurfaceEmission_aer_f;

    this->D4_map["NumberVolumeEmission_aer"] = &NumberVolumeEmission_aer_f;
    this->D4_map["NumberVolumeEmission_aer_i"] = &NumberVolumeEmission_aer_i;
    this->D4_map["NumberVolumeEmission_aer_f"] = &NumberVolumeEmission_aer_f;


    this->D4_map["WetDiameter_aer"] = &WetDiameter_aer;

    this->D4_map["DryDepositionFlux_aer"] = &DryDepositionFlux_aer;
    this->D4_map["WetDepositionFlux_aer"] = &WetDepositionFlux_aer;
    this->D4_map["InCloudWetDepositionFlux_aer"]
      = &InCloudWetDepositionFlux_aer;

    /*** Pointers to 5D data ***/

    this->D5_map["BoundaryCondition_y_aer"] = &BoundaryCondition_y_aer_i;
    this->D5_map["BoundaryCondition_y_aer_i"] = &BoundaryCondition_y_aer_i;

    this->D5_map["BoundaryCondition_x_aer"] = &BoundaryCondition_x_aer_i;
    this->D5_map["BoundaryCondition_x_aer_i"] = &BoundaryCondition_x_aer_i;

    this->D5_map["VolumeEmission_aer"] = &VolumeEmission_aer_f;
    this->D5_map["VolumeEmission_aer_i"] = &VolumeEmission_aer_i;
    this->D5_map["VolumeEmission_aer_f"] = &VolumeEmission_aer_f;

    this->D5_map["Source_aer"] = &Source_aer_f;
    this->D5_map["Source_aer_i"] = &Source_aer_i;
    this->D5_map["Source_aer_f"] = &Source_aer_f;

    /*** Fields bins lists ***/

    this->field_bins["DepositionVelocity_aer"] = &bin_list_dep_aer;
    this->field_bins["ScavengingCoefficient_aer"] = &bin_list_scav_aer;
  }


  //! Destructor.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::~Polair3DAerosol()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ReadConfiguration()
  {
    Polair3DChemistry < T, ClassAdvection, ClassDiffusion,
                        ClassChemistry >::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[options]");
    this->config.PeekValue("With_initial_condition_aerosol",
			   this->option_process
			   ["with_initial_condition_aer"]);
    this->config.PeekValue("With_initial_condition_number_aerosol",
			   this->option_process
			   ["with_initial_condition_number_aer"]);
    this->config.PeekValue("With_boundary_condition_aerosol",
			   this->option_process
			   ["with_boundary_condition_aer"]);
    this->config.PeekValue("With_boundary_condition_number_aerosol",
			   this->option_process
			   ["with_boundary_condition_number_aer"]);
    this->config.PeekValue("With_deposition_aerosol",
                           this->option_process["with_deposition_aer"]);
    if (this->option_process["with_deposition_aer"])
      this->config.PeekValue("Compute_deposition_aerosol",
                             this->option_process["compute_deposition_aer"]);
    else
      this->option_process["compute_deposition_aer"] = false;
    this->config.PeekValue("With_scavenging_aerosol",
                           this->option_process["with_scavenging_aer"]);
    this->config.PeekValue("With_point_emission_aerosol",
                           this->option_process["with_point_emission_aer"]);
    this->config.PeekValue("With_surface_emission_aerosol",
                           this->option_process["with_surface_emission_aer"]);
    this->config.PeekValue("With_surface_number_emission_aerosol",
			   this->option_process["with_surface_number_emission_aer"]);
    this->config.PeekValue("With_volume_emission_aerosol",
			   this->option_process["with_volume_emission_aer"]);
    this->config.PeekValue("With_volume_number_emission_aerosol",
			   this->option_process["With_volume_number_emission_aer"]);
    this->config.PeekValue("With_number_concentration",
			   this->option_process["with_number_concentration"]);
    if (this->config.Check("With_external_composition"))
      this->config.PeekValue("With_external_composition",
                             this->option_process["with_external_composition"]);
    else
      this->option_process["with_external_composition"] = false;

    if (this->option_process["with_air_density"])
      throw string("The option \"With_air_density\" can not yet be used.");
    if (this->option_process["with_number_concentration"] && 
	this->option_process["with_point_emission_aer"])
      throw string("The option with_number_concentration can not yet \"") 
	+ "\" be used with the option with_point_emission_aer.";
    if (this->option_process["with_number_concentration"] && 
	this->source_splitting)
      throw string("The option with_number_concentration can not yet \"") 
	+ "\" be used with the option source_splitting.";
    this->config.PeekValue("Collect_dry_flux_aerosol",
                           this->option_process["collect_dry_flux_aer"]);
    this->config.PeekValue("Collect_wet_flux_aerosol",
                           this->option_process["collect_wet_flux_aer"]);
    // Warning: should be set to true only if chemistry computes pH.
    this->config.PeekValue("With_pH", this->option_process["with_pH"]);

    // Liquid water content threshold above which there are clouds (in g/m3).
    this->config.PeekValue("Lwc_cloud_threshold", "positive",
			   lwc_cloud_threshold);

    this->config.PeekValue("With_fixed_density",
		     this->option_process["with_fixed_density"]);
    // Density in kg / m^3.
    this->config.PeekValue("Fixed_aerosol_density", "positive",
                           fixed_density_aer);

    // Options for radiative computation (photolysis rates="online").
    if (this->computed_photolysis == "online")
      {
        this->config.PeekValue("Wet_computation_option", option_wet_index);
        this->config.PeekValue("Well_mixed_computation_option",
                               option_well_mixed_index);
        this->config.PeekValue("Black_carbon_treatment",
                               option_black_carbon_treatment);
        this->config.PeekValue("Time_step_for_computing_photolysis_rates",
                               time_step_for_computing_photolysis_rates);
        iteration_radiatif = max((int)(time_step_for_computing_photolysis_rates /
                                       this->Delta_t), 1);
        this->config.PeekValue("Directory_efficiency_factor",
                               directory_efficiency_factor);

        this->config.PeekValue("Directory_OPAC", directory_OPAC);
        this->config.PeekValue("File_index_water", file_water_refractive_index);
        this->config.PeekValue("File_species_polyphemus_OPAC",
                               file_species_polyphemus_OPAC);

        this->config.PeekValue("FastJ_parameter_files",
                               fastJ_parameter_files);

        this->config.PeekValue("Tabulation_refractive_index_real",
                               tabulation_index_real);
        this->config.PeekValue("Tabulation_refractive_index_imaginary",
                               tabulation_index_imaginary);
        this->config.PeekValue("Ndiameter", index_diameter);

        this->config.PeekValue("N_OPAC_wavelength", N_OPAC_wavelength);
        this->config.PeekValue("N_water_wavelength", Nwater_wavelength);
        this->config.Find("FastJ_wavelength");
        split(this->config.GetLine(), FastJ_wavelength);
        Nwavelength = FastJ_wavelength.size();
      }

    /*** Disables management of useless fields ***/

    // Variables needed for computing dry deposition.
    this->option_manage["surface_temperature"] =
      this->option_manage["surface_temperature"]
      && this->option_process["compute_deposition_aer"];
    this->option_manage["surface_pressure"] =
      this->option_manage["surface_pressure"]
      && this->option_process["compute_deposition_aer"];
    this->option_manage["first_level_wind_module"] =
      this->option_manage["first_level_wind_module"]
      && this->option_process["compute_deposition_aer"];
    this->option_manage["snow_height"] = this->option_manage["snow_height"]
      && this->option_process["compute_deposition_aer"];
    this->option_manage["specific_humidity"] = true;

    // Other fields needed for computing dry deposition.
    this->option_manage["temperature"] =
      this->option_manage["temperature"]
      || this->option_process["compute_deposition_aer"];
    this->option_manage["pressure"] =
      this->option_manage["pressure"]
      || this->option_process["compute_deposition_aer"];
    this->option_manage["horizontal_wind"]
      = this->option_manage["horizontal_wind"]
      || this->option_process["compute_deposition_aer"];
    this->option_manage["zonal_wind"]
      = this->option_manage["zonal_wind"]
      || this->option_process["compute_deposition_aer"];

    // Rain and cloud height are only needed for scavenging.
    this->option_manage["rain"] = this->option_manage["rain"]
      && (this->scavenging_below_cloud_model != "none"
          || this->scavenging_in_cloud_model != "none"
          || this->option_process["with_scavenging_aer"])
      || this->Chemistry_.IsRequired("rain");
    this->option_manage["cloud_base_height"]
      = this->option_manage["cloud_base_height"]
      && (this->scavenging_below_cloud_model != "none"
          || this->scavenging_in_cloud_model != "none"
          || this->option_process["with_scavenging_aer"]);
    this->option_manage["cloud_top_height"] =
      this->option_manage["cloud_top_height"]
      && this->option_process["with_transport_in_cloud_scavenging"];

    // Other fields needed for computing scavenging.
    this->option_manage["temperature"] =
      this->option_manage["temperature"]
      || this->option_process["with_scavenging_aer"];

    this->option_manage["pressure"] =
      this->option_manage["pressure"]
      || this->option_process["with_scavenging_aer"];

    /*** Species and bins ***/

    // Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    // Section "[aerosol_species]" contains all aerosol species names.
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());
    this->list_aer = this->species_list_aer;
    // Add "Number" to the species list aerosol
    if (this->option_process["with_number_concentration"])
      this->list_aer.push_back("Number");

    // Reads bin bounds.
    this->config.SetSection("[domain]");
    this->config.Find("Bin_bounds");
    bin_list = split(this->config.GetLine());
    // this->Nbin_aer = int(bin_list.size()) - 1; // YK
    this->Nsize_section_aer = int(bin_list.size()) - 1;

    vector<string> species_bin;
    string species;
    int bin_index;

    if (this->option_process["with_external_composition"])
      {
        // Read aerosol groups.
        species_stream.SetSection("[aerosol_groups]");
        while (!species_stream.IsEmpty())
          this->groups_list_aer.push_back(species_stream.GetElement());
        this->Ngroup_aer = int(this->groups_list_aer.size());

        // Read relations between aerosol species and groups.
        species_stream.SetSection("[aerosol_species_group_relations]");
        map<string, string> parameter2;
        map<string, string>::iterator iter2;
        while (!species_stream.IsEmpty())
          {
            species = species_stream.GetElement();
            parameter2[species] = species_stream.GetElement();
          }
        this->aerosol_species_group_relation.resize(this->Ns_aer);
        for (int i = 0; i < this->Ns_aer; i++)
          {
            iter2 = parameter2.find(this->species_list_aer[i]);
            if(iter2 != parameter2.end())
              {
                int group_index = this->GetGroupIndex_aer(iter2->second);
                this->aerosol_species_group_relation(i) = group_index;
              }
            else
              {
                if (this->species_list_aer[i] == "PH2O")
                  this->aerosol_species_group_relation(i) = -1;
                else
                  throw string("The group relation is not available ")
                    + "for the species: "
                    + this->species_list_aer[i] + ".\n"
                    + "See [aerosol_species_group_relations] in the species file";
              }
          }
        parameter2.clear();
      }
    else//in case of internal mixing
      {
        this->groups_list_aer.push_back("None");
        this->Ngroup_aer=int(this->groups_list_aer.size());
        this->aerosol_species_group_relation.resize(this->Ns_aer);
        for (int i = 0; i < this->Ns_aer; i++)
          this->aerosol_species_group_relation(i)=0;
      }
      
    /*** Input files ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    if (this->option_process["with_external_composition"])
      {
        //compute composition discretisation based on number of groups and number of fraction sections
      
        // Reads fraction bounds.
      this->config.SetSection("[domain]");
      this->config.Find("Fraction_bounds");
      fraction_list = split(this->config.GetLine());
      Nfraction_aer = int(fraction_list.size()) - 1;      
      Fractionbound_aer.resize(Nfraction_aer + 1);
      for (int i = 0; i < Nfraction_aer + 1; i++)
        Fractionbound_aer(i) = convert<T>(fraction_list[i]);
      // this->config.Find("Secondary_Fraction_bounds");
      // fraction_list_2 = split(this->config.GetLine());
      // Nfraction_aer_2 = int(fraction_list_2.size()) - 1;
      // Fractionbound_aer_2.resize(Nfraction_aer_2 + 1);
      // for (int i = 0; i < Nfraction_aer_2 + 1; i++)
      //   Fractionbound_aer_2(i) = convert<T>(fraction_list_2[i]);      

        // this->config.SetSection("[options]");
        // this->config.PeekValue("Coagulation_coefficient_file", coagulation_coefficient_file);
        // if(!exists(coagulation_coefficient_file))
        //   throw string("Error, Coagulation Coefficient data is not found, please computed it first!");

        // NcFile fnc(coagulation_coefficient_file.c_str(), NcFile::ReadOnly);

      vector<T> compositions;
      T sumfrac = 0;      
      vector<int> counter(this->Ngroup_aer - 1);

      // Calculate the maximum fraction combinations
      for (int i = 0; i < Nfraction_aer; i++)
        {
          for(int g = 0; g < this->Ngroup_aer-1; g++)
            counter[g] = 0;//initial the counter


          // YK
	  //when the index counter of second group reaches its top,
	  //move to the next fraction bin of first group
	  if (this->Ngroup_aer > 2)
            {
              while (counter[1] <= Nfraction_aer - 1 )
                {
                  // Take the base fraction bounds of current section of 
                  // first group.
                  sumfrac = Fractionbound_aer(i);
                  for(int g = 1; g < this->Ngroup_aer-1; g++)
                    {
                      int j = counter[g];//the fraction list index for group g
                      sumfrac=sumfrac+Fractionbound_aer(j);//calculate one possible combination
                    }
                  if(sumfrac<1.0)
                    {//find a new possible combination
                      this->Ncomposition_aer++;
                      //for first group
                      compositions.push_back(Fractionbound_aer(i));
                      compositions.push_back(Fractionbound_aer(i+1));
                      for(int g=1; g<this->Ngroup_aer-1;g++)
                        {//save possible combinations
                          int j=counter[g];//the fraction list index for group g
                          compositions.push_back(Fractionbound_aer(j));
                          compositions.push_back(Fractionbound_aer(j+1));
                        }
                      //last group with default fraction section [0,1]
                      compositions.push_back(0.0);
                      compositions.push_back(1.0);
                    }

                  //when the second last group hasn't reaches i top,
                  if (counter[this->Ngroup_aer-2]<=Nfraction_aer)
                    counter[this->Ngroup_aer-2]++;//move the index of second last group
         
                  for(int g=3;g<this->Ngroup_aer;g++)
                    {//check every neighbor counter,[2,Ngroup_aer-2] form back to forward
                      int j=this->Ngroup_aer+1-g;
                      sumfrac=Fractionbound_aer(counter[j-1])+Fractionbound_aer(counter[j]);
                      //the bottom sum of two neighbor group fraction is already too big
                      if(sumfrac>=1.0)
                        {
                          //reset all j following counter
                          for(int s=j;s<this->Ngroup_aer-1;s++)
                            counter[s]=0;
                          //increase the
                          counter[j-1]++;
                        }
                    }
                }
            }
	  else
            {//n case of only two
              if(this->Ngroup_aer==2)
                {
                  this->Ncomposition_aer=Nfraction_aer;
                  //for first group
                  compositions.push_back(Fractionbound_aer(i));
                  compositions.push_back(Fractionbound_aer(i+1));
                  //for second group
                  compositions.push_back(0.0);
                  compositions.push_back(1.0);
                }
            }
        }

      this->composition_bounds.Resize(this->Ncomposition_aer,this->Ngroup_aer,2);
      int iter=0;
      for(int c=0; c<this->Ncomposition_aer; c++)
	for(int g=0; g<this->Ngroup_aer;g++)
	  for(int k=0; k<2; k++)
	  {
	    this->composition_bounds(c,g,k)=compositions[iter];
	    iter++;
	  }
      this->Nbin_aer=this->Nsize_section_aer*this->Ncomposition_aer;
    }
    else//in case of internal mixing
    {
      this->Ncomposition_aer=1;
      this->Nbin_aer=this->Nsize_section_aer;
      this->composition_bounds.Resize(1,1,2);
      this->composition_bounds(0,0,0)=0.0;
      this->composition_bounds(0,0,1)=1.0;
    }

    // Reads mass density aerosol
    map<string, T> parameter;
    typename map<string, T>::iterator iter;

    ConfigStream config_species(this->GetSpeciesFile());
    config_species.SetSection("[mass_density_aer]");
    Mass_Density_aer.resize(this->Ns_aer);
    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	config_species.GetNumber(parameter[species]);
      }
    
    for (int i = 0; i < this->Ns_aer; i++)
      {
	iter = parameter.find(this->species_list_aer[i]);
	if(iter != parameter.end())
	  Mass_Density_aer[i] = iter->second;
	else
	  throw string("Module: no mass density for ")
	    + string("aerosol species \"") + this->species_list_aer[i]
	    + string("\". Please provide one.");
      }
    parameter.clear();


    // Initial conditions.
    if (this->option_process["with_initial_condition_aer"])
    {
      data_description_stream.SetSection("[initial_condition_aerosol]");
      if (this->option_process["with_external_composition"])
        data_description_stream.PeekValue("Format", ic_format);
      else
        ic_format = "Internal";
      this->input_files["initial_condition_aer"]
	.ReadFiles(data_description_file, "initial_condition_aerosol");
    }
    else
      this->input_files["initial_condition_aer"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["initial_condition_aer"].Begin();
         i != this->input_files["initial_condition_aer"].End(); i++)
      {
	species_bin = split(i->first, "_");
	if (species_bin.size() != 2)
	  throw string("Species \"") + i->first + "\" is badly formatted.";
	species = species_bin[0];
	bin_index = convert<int>(species_bin[1]);
	if (bin_index < this->Nsize_section_aer)
	  // Built bin list with initial conditions => ic_bin_list_aer
	  this->InitialConditionBinList_aer(bin_index);
	else
	  throw string("initial_condition_aer: Bin index in data file are out of range");
	unsigned int j = 0;
	while (j < species_list_ic_aer.size()
	       && species_list_ic_aer[j].first != species)
	  j++;
	if (j == species_list_ic_aer.size())
	  species_list_ic_aer
	    .push_back(pair<string, vector<int> >(species, vector<int>()));
	species_list_ic_aer[j].second.push_back(bin_index);
      }
    Ns_ic_aer = int(species_list_ic_aer.size());
    Nb_ic_aer = int(ic_bin_list_aer.size());

    //add for number
    if (this->option_process["with_number_concentration"])//
      if (this->option_process["with_initial_condition_number_aer"])
	this->input_files["initial_condition_aer"]
	  .ReadFiles(data_description_file, "initial_number_aerosol");
    
    // Boundary conditions.
    if (this->option_process["with_boundary_condition_aer"])
    {
      data_description_stream.SetSection("[boundary_condition_aerosol]");
      if (this->option_process["with_external_composition"])
        data_description_stream.PeekValue("Format", bc_format);
      else
        bc_format = "Internal";
      this->input_files["boundary_condition_aer"]
	.Read(data_description_file, "boundary_condition_aerosol");
    }
    else
      this->input_files["boundary_condition_aer"].Empty();
    Nbin_bc_aer = 0;
    for (map<string, string>::iterator i
           = this->input_files["boundary_condition_aer"].Begin();
         i != this->input_files["boundary_condition_aer"].End(); i++)
      {
        species_bin = split(i->first, "_");
        if (species_bin.size() != 2)
          throw string("Species \"") + i->first + "\" is badly formatted.";
        species = species_bin[0];
        bin_index = convert<int>(species_bin[1]);
	    if (bin_index < this->Nsize_section_aer)
	       	// Built bin list with boundary conditions => bc_bin_list_aer
	  		this->BoundaryConditionBinList_aer(bin_index);
		else
	  		throw string("boundary_condition_aer: Bin index in data file are out of range");
        unsigned int j = 0;
        while (j < species_list_bc_aer.size()
               && species_list_bc_aer[j].first != species)
          j++;
        if (j == species_list_bc_aer.size())
          species_list_bc_aer
            .push_back(pair<string, vector<int> >(species, vector<int>()));
        species_list_bc_aer[j].second.push_back(bin_index);
        Nbin_bc_aer
          = max(Nbin_bc_aer, int(species_list_bc_aer[j].second.size()));
      }
    Ns_bc_aer = int(species_list_bc_aer.size());
    Nb_bc_aer = int(bc_bin_list_aer.size());
    //add for number
    if (this->option_process["with_number_concentration"])
      if (this->option_process["with_boundary_condition_number_aer"])
	this->input_files["boundary_condition_aer"]
	  .ReadFiles(data_description_file, "boundary_number_aerosol");
	
    // Deposition velocities.
    if (this->option_process["with_deposition_aer"])
      if (this->option_process["compute_deposition_aer"])
        this->input_files["deposition_velocity_aer"]
          .ReadFields(data_description_file, "deposition_velocity_aerosol");
      else
        this->input_files["deposition_velocity_aer"]
          .Read(data_description_file, "deposition_velocity_aerosol");
    else
      this->input_files["deposition_velocity_aer"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["deposition_velocity_aer"].Begin();
         i != this->input_files["deposition_velocity_aer"].End(); i++)
      {
        if (!is_integer(i->first))
          throw string("In deposition fields, \"") + i->first
            + "\" is not a bin index.";
        bin_index = convert<int>(i->first);
        bin_list_dep_aer.push_back(bin_index);
      }
    Nbin_dep_aer = int(bin_list_dep_aer.size());

    // Scavenging for gaseous species.
    if (this->scavenging_below_cloud_model == "microphysical-ph")
      {
        // Henry constants in mol / L / atm.
        species_stream.SetSection("[henry]");
        while (!species_stream.IsEmpty())
          {
            species = species_stream.GetElement();
            species_stream.GetNumber(this->henry_constant[species]);
          }

        // Gas phase diffusivity constants in cm^2 / s.
        species_stream.SetSection("[diffusivity]");
        while (!species_stream.IsEmpty())
          {
            species = species_stream.GetElement();
            species_stream.GetNumber(this->gas_phase_diffusivity[species]);
          }

        // Dissolution heat (kcal / mol at 298K).
        species_stream.SetSection("[dissolution_heat]");
        while (!species_stream.IsEmpty())
          {
            species = species_stream.GetElement();
            species_stream.GetNumber(dissolution_heat[species]);
          }
      }

    // Scavenging for aerosols.

    if (this->option_process["with_scavenging_aer"])
      {
        this->option_manage["rain"] = true;
        this->option_manage["cloud_height"] = true;
        this->input_files["scavenging_aer"]
          .ReadFields(data_description_file, "scavenging_aerosol");
      }
    else
      this->input_files["scavenging_aer"].Empty();

    for (map<string, string>::iterator i
           = this->input_files["scavenging_aer"].Begin();
         i != this->input_files["scavenging_aer"].End(); i++)
      {
        if (!is_integer(i->first))
          throw string("In scavenging fields, \"") + i->first
            + "\" is not a bin index.";
        bin_index = convert<int>(i->first);
        bin_list_scav_aer.push_back(bin_index);
      }
    Nbin_scav_aer = int(bin_list_scav_aer.size());

    // Point emissions.
    if (this->option_process["with_point_emission_aer"])
      {
        string point_emission_file;
        data_description_stream.SetSection("[point_emission_aerosol]");
        data_description_stream.PeekValue("file", point_emission_file);
        if (this->option_process["with_external_composition"])
          data_description_stream.PeekValue("Format", point_emis_format);
        else
          point_emis_format = "Internal";

        ConfigStream point_emission_stream(point_emission_file);

        // Parses the file that describes all point emissions.
        string type, str, line;
        vector<string> str_split;
        T x, y, z;
        int index;
        int Nsource = 0;
        while (!point_emission_stream.IsEmpty())
          {
            line = point_emission_stream.GetLine();
            // A new source should be added.
            if (split(line)[0] == "[source]")
              {
                point_emission_list_aer.push_back(map<string, string>());
                // Emitted species.
                point_emission_stream.PeekValue("Species", str);
                str_split = split(str, "_");
                if (str_split.size() != 2)
                  throw string("In point emissions, badly formatted ")
                    + string("species: \"") + str + "\".";
                point_emission_list_aer[Nsource]["species"] = str_split[0];
                point_emission_list_aer[Nsource]["bin"] = str_split[1];
                // Coordinates.
                point_emission_stream.PeekValue("Abscissa", x);
                point_emission_stream.PeekValue("Ordinate", y);
                point_emission_stream.PeekValue("Altitude", "positive", z);
                // Maps its coordinates to model grid indices.
                index = int((x - this->x_min + this->Delta_x / 2.)
                            / this->Delta_x);
                point_emission_list_aer[Nsource]["i"] = to_str(index);
                index = int((y - this->y_min + this->Delta_y / 2.)
                            / this->Delta_y);
                point_emission_list_aer[Nsource]["j"] = to_str(index);
                index = 0;
                while (index < this->Nz && z > this->LayerInterface(index + 1))
                  index++;
                point_emission_list_aer[Nsource]["k"] = to_str(index);
                // Source type (instantaneous or continuous).
                point_emission_stream.PeekValue("Type", "puff | continuous",
                                                type);
                point_emission_list_aer[Nsource]["type"] = type;
                if (type == "puff")
                  // Instantaneous emission.
                  {
                    point_emission_stream.PeekValue("Date", str);
                    point_emission_list_aer[Nsource]["date"] = str;
                    // Quantity emitted.
                    point_emission_stream.PeekValue("Quantity", "positive",
                                                    str);
                    point_emission_list_aer[Nsource]["quantity"] = str;
                  }
                else if (type == "continuous")
                  // Continuous emission.
                  {
                    point_emission_stream.PeekValue("Date_beg", str);
                    point_emission_list_aer[Nsource]["date_beg"] = str;
                    point_emission_stream.PeekValue("Date_end", str);
                    point_emission_list_aer[Nsource]["date_end"] = str;
                    // Emission rate.
                    point_emission_stream.PeekValue("Rate", "positive", str);
                    point_emission_list_aer[Nsource]["rate"] = str;
                  }
                Nsource++;
              }
          }
      }

    // Surface emissions.
    if (this->option_process["with_surface_emission_aer"])
    {
      data_description_stream.SetSection("[surface_emission_aerosol]");
      if (this->option_process["with_external_composition"])
        data_description_stream.PeekValue("Format", surface_emis_format);
      else
        surface_emis_format = "Internal";
      this->input_files["surface_emission_aer"]
	.Read(data_description_file, "surface_emission_aerosol");
    }
    else
      this->input_files["surface_emission_aer"].Empty();
    for (map<string, string>::iterator i
           = this->input_files["surface_emission_aer"].Begin();
         i != this->input_files["surface_emission_aer"].End(); i++)
      {
        species_bin = split(i->first, "_");
        if (species_bin.size() != 2)
          throw string("Species \"") + i->first + "\" is badly formatted.";
        species = species_bin[0];
        bin_index = convert<int>(species_bin[1]);
        if(bin_index < this->Nsize_section_aer)
          // Build bin list with surface emissions
          this->SurfaceEmissionBinList_aer(bin_index);
        else
          throw string("surface_emission_aer: Bin index in data file are out of range") + species;
        unsigned int j = 0;
        while (j < species_list_surf_emis_aer.size()
               && species_list_surf_emis_aer[j].first != species)
          j++;
        if (j == species_list_surf_emis_aer.size())
          species_list_surf_emis_aer
            .push_back(pair<string, vector<int> >(species, vector<int>()));
        species_list_surf_emis_aer[j].second.push_back(bin_index);
      }
    Ns_surf_emis_aer = int(species_list_surf_emis_aer.size());
    Nb_surf_emis_aer = int(surface_emission_bin_list_aer.size());

    //add for number
    if (this->option_process["with_surface_number_emission_aer"])//
      this->input_files["surface_emission_aer"]
	.ReadFiles(data_description_file, "surface_emission_number_aerosol");
	  
    // Volume emissions.
    if (this->option_process["with_volume_emission_aer"])
      this->input_files["volume_emission_aer"]
        .Read(data_description_file, "volume_emission_aerosol");
    else
      this->input_files["volume_emission_aer"].Empty();
    Nb_vol_emis_aer = 0;
    for (map<string, string>::iterator i
           = this->input_files["volume_emission_aer"].Begin();
         i != this->input_files["volume_emission_aer"].End(); i++)
      {
	species_bin = split(i->first, "_");
	if (species_bin.size() != 2)
	  throw string("Species \"") + i->first + "\" is badly formatted.";
	species = species_bin[0];
	bin_index = convert<int>(species_bin[1]);
	if(bin_index < this->Nsize_section_aer)
	  // Build bin list with volume emission
	  this->VolumeEmissionBinList_aer(bin_index);
	else
	  throw string("volume_emission_aer: Bin index in data file are out of range");
	unsigned int j = 0;
	while (j < species_list_vol_emis_aer.size()
	       && species_list_vol_emis_aer[j].first != species)
	  j++;
	if (j == species_list_vol_emis_aer.size())
	  species_list_vol_emis_aer
	    .push_back(pair<string, vector<int> >(species, vector<int>()));
	species_list_vol_emis_aer[j].second.push_back(bin_index);
	Nb_vol_emis_aer =
	  max(Nb_vol_emis_aer,
	      int(species_list_vol_emis_aer[j].second.size()));
      }
    Ns_vol_emis_aer = int(species_list_vol_emis_aer.size());
    Nb_vol_emis_aer = int(volume_emission_bin_list_aer.size());

    //add for number
    if (this->option_process["with_volume_number_emission_aer"])//
      this->input_files["volume_emission_aer"]
	.ReadFiles(data_description_file, "volume_emission_number_aerosol");
	
    // Additional variable: number of levels.
    if (this->option_process["with_volume_emission_aer"])
      {
	data_description_stream.SetSection("[volume_emission_aerosol]");
	data_description_stream.PeekValue("Nz", "> 0", Nz_vol_emis_aer);
        if (this->option_process["with_external_composition"])
          data_description_stream.PeekValue("Format", volume_emis_format);
        else
          volume_emis_format = "Internal";
      }
    else
      Nz_vol_emis_aer = 0;
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::CheckConfiguration()
  {
    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::CheckConfiguration();

    if (this->option_manage["snow_height"]
        && this->input_files["meteo"]("SnowHeight").empty())
      throw "Snow height is needed but no input data file was provided.";

    if (this->scavenging_below_cloud_model == "microphysical-ph"
        && !this->option_process["with_pH"])
      throw string("Scavenging model \"microphysical-ph\" cannot be used ")
        + "without pH (option \"With_pH\").";

    if (this->scavenging_below_cloud_model == "microphysical"
        || this->scavenging_below_cloud_model == "microphysical-ph")
      {
        // Checks if Henry constant is well initialized for each scavenged
        // species.
        for (int i = 0; i < this->Ns_scav; i++)
          if (this->henry_constant.find(this->ScavengingName(i))
              == this->henry_constant.end())
            throw string("ERROR! Henry constant for species \"")
              + this->ScavengingName(i)
              + string("\" not found in section [henry].");

        // Checks if gas phase diffusivity constant is well initialized for
        // all scavenged species.
        for (int i = 0; i < this->Ns_scav; i++)
          if (this->gas_phase_diffusivity.find(this->ScavengingName(i))
              == this->gas_phase_diffusivity.end())
            throw string("ERROR! Gas phase diffusivity constant")
              + string(" for species \"") + this->ScavengingName(i)
              + "\" not found in section [diffusivity].";
      }

    if (this->scavenging_below_cloud_model == "microphysical-ph")
      // Checks if dissolution heats are well initialized for all scavenged
      // species.
      for (int i = 0; i < this->Ns_scav; i++)
        if (dissolution_heat.find(this->ScavengingName(i))
            == dissolution_heat.end())
          throw string("ERROR! Dissolution heat for species \"")
            + this->ScavengingName(i)
            + string("\" not found in section [dissolution_heat].");

    if (!this->option_process["with_deposition_aer"]
        && this->option_process["collect_dry_flux_aer"])
      throw string("Aerosol dry deposition fluxes cannot be collected") +
        " without deposition.";

    if (this->option_process["collect_wet_flux_aer"]
        && (!this->option_process["with_scavenging_aer"]))
      throw string("Aerosol wet deposition fluxes cannot be collected") +
        " without scavenging.";
  }

  //! Checks whether the model deals with number concentration
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasNumberConcentration_aer()
  {
    if(this->option_process["with_number_concentration"])
      return true;
    else
      return false;
  }

  //! Checks whether an aerosol has initial conditions.
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has initial conditions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasInitialCondition_aer(int s, int b) const
  {
    return HasInitialCondition_aer(this->GetSpeciesName_aer(s), b);
  }


  //! Checks whether an aerosol has initial conditions.
  /*!
    \param name species name.
    \param b bin number.
    \return True if the aerosol has initial conditions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasInitialCondition_aer(string name, int x) const
  {
    unsigned int i;
    int b = Bin_to_size_index_aer(x);
    for (i = 0; i < species_list_ic_aer.size(); i++)
      if (species_list_ic_aer[i].first == name)
        return find(species_list_ic_aer[i].second.begin(),
                    species_list_ic_aer[i].second.end(), b)
          != species_list_ic_aer[i].second.end();
    return false;
  }

  //! Checks whether an aerosol has number initial conditions.
  /*!
    \param b bin number.
    \return True if the aerosol has number initial conditions, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasNumberInitialCondition_aer(int x) const
  {
    int b = Bin_to_size_index_aer(x);
    return find(ic_bin_list_aer.begin(),
		ic_bin_list_aer.end(), b)
      != ic_bin_list_aer.end();
  }


  //! Returns the index in number initial conditions of a given aerosol.
  /*!
    \param b bin number.
    \return The aerosol index in initial conditions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::NumberInitialConditionIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if(HasNumberInitialCondition_aer(x))
    {
      int binout= find(ic_bin_list_aer.begin(),
		  ic_bin_list_aer.end(), b)
	- ic_bin_list_aer.begin();
      return Bin_index_translate_aer(binout,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in Initial Conditions.";
  }

  
  //! Returns the bin index in initial conditions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in initial conditions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitialConditionBinList_aer(int b)
  {
    
    vector<int>::iterator pos = find(ic_bin_list_aer.begin(),
                                     ic_bin_list_aer.end(), b);
    if (pos == ic_bin_list_aer.end())
      ic_bin_list_aer.push_back(b);

    return ic_bin_list_aer;
  }
  

  //! Checks whether an aerosol has boundary conditions.
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has boundary conditions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasBoundaryCondition_aer(int s, int b) const
  {
    return HasBoundaryCondition_aer(this->GetSpeciesName_aer(s), b);
  }


  //! Checks whether an aerosol has boundary conditions.
  /*!
    \param name species name.
    \param b bin number.
    \return True if the aerosol has boundary conditions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasBoundaryCondition_aer(string name, int x) const
  {
    unsigned int i;
    int b = Bin_to_size_index_aer(x);
    for (i = 0; i < species_list_bc_aer.size(); i++)
      if (species_list_bc_aer[i].first == name)
        return find(species_list_bc_aer[i].second.begin(),
                    species_list_bc_aer[i].second.end(), b)
          != species_list_bc_aer[i].second.end();
    return false;
  }


  //! Checks whether an aerosol has number boundary conditions.
  /*!
    \param b bin number.
    \return True if the aerosol has number boundary conditions, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasNumberBoundaryCondition_aer(int x) const
  {
    int b = Bin_to_size_index_aer(x);
    return find(bc_bin_list_aer.begin(),
		bc_bin_list_aer.end(), b)
      != bc_bin_list_aer.end();
  }


  //! Returns the index in boundary conditions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in boundary conditions.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::BoundaryConditionIndex_aer(int s, int x) const
  {
    vector<int> output(2);
    unsigned int i;
    int b = Bin_to_size_index_aer(x);
    string name = this->GetSpeciesName_aer(s);
    for (i = 0; i < species_list_bc_aer.size(); i++)
      if (species_list_bc_aer[i].first == name)
	{
	  output[0] = int(i);
	  int binout= find(species_list_bc_aer[i].second.begin(),
			   species_list_bc_aer[i].second.end(), b)
	    - species_list_bc_aer[i].second.begin();
	  output[1]= Bin_index_translate_aer(binout,x);
	  return output;
	}
    throw string("Species \"") + name + string("_") + to_str(b)
      + "\" not found in boundary conditions.";
  }


  //! Returns the bin index in number boundary conditions.
  /*!
    \param b bin number.
    \return The aerosol index in number boundary conditions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::NumberBoundaryConditionIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if (HasNumberBoundaryCondition_aer(x))
    {
      int binout=find(bc_bin_list_aer.begin(),
		  bc_bin_list_aer.end(), b)
	- bc_bin_list_aer.begin();
	return Bin_index_translate_aer(binout,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in boundary conditions.";
  }


  //! Returns the bin index in boundary conditions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in boundary conditions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::BoundaryConditionBinList_aer(int b)
  {
    vector<int>::iterator pos = find(bc_bin_list_aer.begin(),
                                     bc_bin_list_aer.end(), b);
    if (pos == bc_bin_list_aer.end())
      bc_bin_list_aer.push_back(b);

    return bc_bin_list_aer;
  }


  //! Checks whether aerosols in a given bin have deposition velocities.
  /*!
    \param b bin number.
    \return True if the aerosols in bin \a b has deposition velocities, false
    otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasDepositionVelocity_aer(int x) const
  {
    int b = Bin_to_size_index_aer(x);
    return find(bin_list_dep_aer.begin(), bin_list_dep_aer.end(), b)
      != bin_list_dep_aer.end();
  }


//! Returns the index in deposition velocities of a given aerosol bin.
  /*!
    \param b bin number.
    \return The aerosol bin index in deposition velocities.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::DepositionVelocityIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    int binout=find(bin_list_dep_aer.begin(), bin_list_dep_aer.end(), b)
      - bin_list_dep_aer.begin();
    return Bin_index_translate_aer(binout,x);
  }


  //! Checks whether aerosols in a given bin have scavenging.
  /*!
    \param b bin number.
    \return True if the aerosols in bin \a b has deposition velocities, false
    otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasScavenging_aer(int x) const
  {
    int b = Bin_to_size_index_aer(x);
    return find(bin_list_scav_aer.begin(), bin_list_scav_aer.end(), b)
      != bin_list_scav_aer.end();
  }


  //! Returns the index in scavenging of a given aerosol bin.
  /*!
    \param b bin number.
    \return The aerosol bin index in scavenging.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ScavengingIndex_aer(int b) const
  {
    return find(bin_list_scav_aer.begin(), bin_list_scav_aer.end(), b)
      - bin_list_scav_aer.begin();
  }


  //! Checks whether an aerosol has surface emissions.
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has surface emissions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasSurfaceEmission_aer(int s, int b) const
  {
    return HasSurfaceEmission_aer(this->GetSpeciesName_aer(s), b);
  }


  //! Checks whether an aerosol has surface emissions.
  /*!
    \param name species name.
    \param b bin number.
    \return True if the aerosol has surface emissions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasSurfaceEmission_aer(string name, int x) const
  {
    unsigned int i;
    int b = Bin_to_size_index_aer(x);
    for (i = 0; i < species_list_surf_emis_aer.size(); i++)
      if (species_list_surf_emis_aer[i].first == name)
        return find(species_list_surf_emis_aer[i].second.begin(),
                    species_list_surf_emis_aer[i].second.end(), b)
          != species_list_surf_emis_aer[i].second.end();
    return false;
  }

  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 3>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetWetDepositionFluxNumber_aer()
  {
    return this->WetDepositionFluxNumber_aer;
  }
  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 3>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetDryDepositionFluxNumber_aer()
  {
    return this->DryDepositionFluxNumber_aer;
  }
  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 3>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetInCloudWetDepositionFluxNumber_aer()
  {
    return this->InCloudWetDepositionFluxNumber_aer;
  }

  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 4>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetWetDepositionFlux_aer()
  {
    return this->WetDepositionFlux_aer;
  }
  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 4>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetDryDepositionFlux_aer()
  {
    return this->DryDepositionFlux_aer;
  }
  //! Returns the WetDepositionFlux Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  Data<T, 4>& Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetInCloudWetDepositionFlux_aer()
  {
    return this->InCloudWetDepositionFlux_aer;
  }

  
  //! Checks whether an aerosol has number surface emissions.
  /*!
    \param b bin number.
    \return True if the aerosol has surface emissions, false otherwise.
  */
  
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasNumberSurfaceEmission_aer(int x) const
  {
    int b = Bin_to_size_index_aer(x);
    return find(surface_emission_bin_list_aer.begin(),
		surface_emission_bin_list_aer.end(), b)
      != surface_emission_bin_list_aer.end();
  }
  

  //! Returns the index in surface emissions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in surface emissions.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SurfaceEmissionIndex_aer(int s, int x) const
  {
    vector<int> output(2);
    unsigned int i;
    int b=Bin_to_size_index_aer(x);
    string name = this->GetSpeciesName_aer(s);
    for (i = 0; i < species_list_surf_emis_aer.size(); i++)
      if (species_list_surf_emis_aer[i].first == name)
	{
	  output[0] = int(i);
	  int outbin = find(species_list_surf_emis_aer[i].second.begin(),
			   species_list_surf_emis_aer[i].second.end(), b)
	    - species_list_surf_emis_aer[i].second.begin();
	  output[1] = Bin_index_translate_aer(outbin,x);
	  return output;
	}
    throw string("Species \"") + name + string("_") + to_str(b)
      + "\" not found in surface emissions.";
  }


  //! Returns the index in number surface emissions of a given aerosol.
  /*!
    \param b bin number.
    \return The aerosol index in surface emissions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::NumberSurfaceEmissionIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if(HasNumberSurfaceEmission_aer(x))
    {
      int outbin=find(surface_emission_bin_list_aer.begin(),
		  surface_emission_bin_list_aer.end(), b)
	- surface_emission_bin_list_aer.begin();
      return Bin_index_translate_aer(outbin,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in Surface Emissions.";
	  
  }
  
  //! Returns the bin index in surface emissions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in surface emissions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SurfaceEmissionBinList_aer(int b)
  {
    vector<int>::iterator pos = find(surface_emission_bin_list_aer.begin(),
                                     surface_emission_bin_list_aer.end(), b);
    if (pos == surface_emission_bin_list_aer.end())
      surface_emission_bin_list_aer.push_back(b);

    return surface_emission_bin_list_aer;
  }
  
  
  //! Checks whether an aerosol has volume emissions.
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasVolumeEmission_aer(int s, int b) const
  {
    return HasVolumeEmission_aer(this->GetSpeciesName_aer(s), b);
  }


  //! Checks whether an aerosol has volume emissions.
  /*!
    \param name species name.
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasVolumeEmission_aer(string name, int x) const
  {
    unsigned int i;
    int b = Bin_to_size_index_aer(x);
    for (i = 0; i < species_list_vol_emis_aer.size(); i++)
      if (species_list_vol_emis_aer[i].first == name)
        return find(species_list_vol_emis_aer[i].second.begin(),
                    species_list_vol_emis_aer[i].second.end(), b)
          != species_list_vol_emis_aer[i].second.end();
    return false;
  }

  //! Checks whether an aerosol has number volume emissions.
  /*!
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasNumberVolumeEmission_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    return find(this->volume_emission_bin_list_aer.begin(),
		this->volume_emission_bin_list_aer.end(), b)
      != this->volume_emission_bin_list_aer.end();
  }


  //! Returns the index in volume emissions of a given aerosol.
  /*!
    \param s aerosol global index.
    \return The aerosol index in volume emissions.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::VolumeEmissionIndex_aer(int s) const
  {
    string name = this->GetSpeciesName_aer(s);
    unsigned int i;
    for (i = 0; i < species_list_vol_emis_aer.size(); i++)
      if (species_list_vol_emis_aer[i].first == name)
        return int(i);
    throw string("Unable to find species \"") + name
      + "\" in aerosol volume emissions.";
  }


  //! Returns the name of an aerosol with volume emissions.
  /*!
    \param s species index in volume emissions (for aerosols).
    \return The species name.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  string Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::VolumeEmissionName_aer(int s) const
  {
    return species_list_vol_emis_aer.at(s).first;
  }


  //! Returns the bin index in volume emissions of a given aerosol.
  /*!
    \param s species global index.
    \param b bin number.
    \return The aerosol index in volume emissions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  vector<int>
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::VolumeEmissionBinList_aer(int b)
  {
    vector<int>::iterator pos = find(volume_emission_bin_list_aer.begin(),
                                     volume_emission_bin_list_aer.end(), b);
    if (pos == volume_emission_bin_list_aer.end())
      volume_emission_bin_list_aer.push_back(b);

    return volume_emission_bin_list_aer;
  }

  
  //! Returns the index in number volume emissions of a given aerosol.
  /*!
    \param b bin number.
    \return The aerosol index in volume emissions.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::NumberVolumeEmissionIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if(HasNumberVolumeEmission_aer(x))
    {
      int binout=find(volume_emission_bin_list_aer.begin(),
		  volume_emission_bin_list_aer.end(), b)
	- volume_emission_bin_list_aer.begin();
      return Bin_index_translate_aer(binout,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in Volume Emissions.";
  }

  //! Returns the index of size section based on the index of bin .
  /*!
    \param b bin number.
    \return index of size section.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Bin_to_size_index_aer(int b) const
  {
    if(b<this->Nbin_aer && b >=0)
      return (b/this->Ncomposition_aer);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }

  //! Returns the index of compsition section based on the index of bin .
  /*!
    \param b bin number.
    \return index of compsition section.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Bin_to_composition_index(int b) const
  {
    if(b<this->Nbin_aer && b >=0)
      return (b-Bin_to_size_index_aer(b)*this->Ncomposition_aer);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }

  //! translate bin b from current size section to new size section s.
  /*!
    \param b bin number. s the new size section number
    \return index of compsition section.
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Bin_index_translate_aer(int s, int b) const
  {
    int composition_id=Bin_to_composition_index(b);
    if(b<this->Nbin_aer && b >=0)
      return (s*this->Ncomposition_aer+composition_id);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }
  
  //! Returns the global index of a species (aerosol) with volume emissions.
  /*!
    \param s species index in volume emissions (for aerosols).
    \return The species global index (aerosol).
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::VolumeEmissionGlobalIndex_aer(int s) const
  {
    return this->GetSpeciesIndex_aer(VolumeEmissionName_aer(s));
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Allocate()
  {
    Polair3DChemistry < T, ClassAdvection, ClassDiffusion,
                        ClassChemistry >::Allocate();

    /*** Additional grids ***/

    GridZ5D = RegularGrid<T>(this->Nz);
    GridY5D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridX5D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);

    GridY5D_interf_bc = RegularGrid<T>(this->y_min - this->Delta_y,
                                       T(this->Ny + 1) * this->Delta_y, 2);
    GridX5D_interf_bc = RegularGrid<T>(this->x_min - this->Delta_x,
				       T(this->Nx + 1) * this->Delta_x, 2);

    GridY3D_interf_bc = RegularGrid<T>(this->y_min - this->Delta_y,
				       T(this->Ny + 1) * this->Delta_y, 2);
    GridX3D_interf_bc = RegularGrid<T>(this->x_min - this->Delta_x,
				       T(this->Nx + 1) * this->Delta_x, 2);				       

    BinBound_aer.resize(this->Nsize_section_aer + 1);
    // Reads bin bounds in micrometers and converts it to meters.
    for (int i = 0; i < this->Nsize_section_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]);
    
    GridG3D_aer = RegularGrid<T>(this->Ngroup_aer);
    GridG4D_aer = RegularGrid<T>(this->Ngroup_aer);
    GridG5D_aer = RegularGrid<T>(this->Ngroup_aer);

    GridC4D_aer = RegularGrid<T>(this->Ncomposition_aer);
    GridC5D_aer = RegularGrid<T>(this->Ncomposition_aer);
    
    GridB4D_aer = RegularGrid<T>(this->Nbin_aer);//Nsize_section_aer
    GridB5D_aer = RegularGrid<T>(this->Nbin_aer);
    GridB4D_aer_i = RegularGrid<T>(this->Nsize_section_aer);
    GridB5D_aer_i = RegularGrid<T>(this->Nsize_section_aer);
    GridS4D_aer = RegularGrid<T>(this->Ns_aer);
    GridS5D_aer = RegularGrid<T>(this->Ns_aer);
    GridS5D_aer_noH2O = RegularGrid<T>(this->Ns_aer - 1);

    GridS_bc_aer = RegularGrid<T>(Ns_bc_aer);
    GridB_bc_aer_i = RegularGrid<T>(Nb_bc_aer);
    GridB_bc_aer = RegularGrid<T>(Nb_bc_aer*this->Ncomposition_aer);

    GridB_dep_aer = RegularGrid<T>(Nbin_dep_aer*this->Ncomposition_aer);
    GridB_scav_aer = RegularGrid<T>(Nbin_scav_aer*this->Ncomposition_aer);

    GridS_surf_emis_aer = RegularGrid<T>(Ns_surf_emis_aer);
    GridB_surf_emis_aer = RegularGrid<T>(Nb_surf_emis_aer*this->Ncomposition_aer);

    GridS_vol_emis_aer = RegularGrid<T>(Ns_vol_emis_aer);
    GridB_vol_emis_aer = RegularGrid<T>(Nb_vol_emis_aer*this->Ncomposition_aer);   
    GridZ_vol_emis_aer = RegularGrid<T>(Nz_vol_emis_aer);

    GridZ5D.SetVariable(2);
    GridZ5D.SetDuplicate(false);

    GridY5D.SetVariable(3);
    GridY5D.SetDuplicate(false);
    GridY5D_interf_bc.SetVariable(3);
    GridY5D_interf_bc.SetDuplicate(false);

    GridX5D.SetVariable(4);
    GridX5D.SetDuplicate(false);
    GridX5D_interf_bc.SetVariable(4);
    GridX5D_interf_bc.SetDuplicate(false);

    /*** Meteorological fields ***/

    SurfaceTemperature_i.Resize(this->GridY2D, this->GridX2D);
    SurfaceTemperature_f.Resize(this->GridY2D, this->GridX2D);
    FileSurfaceTemperature_i.Resize(this->GridY2D, this->GridX2D);
    FileSurfaceTemperature_f.Resize(this->GridY2D, this->GridX2D);

    SurfacePressure_i.Resize(this->GridY2D, this->GridX2D);
    SurfacePressure_f.Resize(this->GridY2D, this->GridX2D);
    FileSurfacePressure_i.Resize(this->GridY2D, this->GridX2D);
    FileSurfacePressure_f.Resize(this->GridY2D, this->GridX2D);

    LiquidWaterContent_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileLiquidWaterContent_i.Resize(this->GridZ3D,
                                    this->GridY3D, this->GridX3D);
    FileLiquidWaterContent_f.Resize(this->GridZ3D,
                                    this->GridY3D, this->GridX3D);

    if (this->computed_photolysis == "online")
      {
        CloudOpticalDepth_i.Resize
          (this->GridZ3D, this->GridY3D, this->GridX3D);
        FileCloudOpticalDepth_i.Resize(this->GridZ3D,
                                       this->GridY3D, this->GridX3D);
        FileCloudOpticalDepth_f.Resize(this->GridZ3D,
                                       this->GridY3D, this->GridX3D);

        IceOpticalDepth_i.Resize
          (this->GridZ3D, this->GridY3D, this->GridX3D);
        FileIceOpticalDepth_i.Resize(this->GridZ3D,
                                     this->GridY3D, this->GridX3D);
        FileIceOpticalDepth_f.Resize(this->GridZ3D,
                                     this->GridY3D, this->GridX3D);
      }

    if (this->option_process["compute_deposition_aer"])
      {
        FirstLevelWindModule_i.Resize(this->GridY2D, this->GridX2D);
        FirstLevelWindModule_f.Resize(this->GridY2D, this->GridX2D);

        SnowHeight_i.Resize(this->GridY2D, this->GridX2D);
        SnowHeight_f.Resize(this->GridY2D, this->GridX2D);
        FileSnowHeight_i.Resize(this->GridY2D, this->GridX2D);
        FileSnowHeight_f.Resize(this->GridY2D, this->GridX2D);
      }

    RelativeHumidity_i.Resize
      (this->GridZ3D, this->GridY3D, this->GridX3D);
    RelativeHumidity_f.Resize
      (this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** Boundary conditions ***/

    // Along Z.
    BoundaryCondition_z_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
				     this->GridY4D, this->GridX4D);
    BoundaryCondition_z_aer_i.SetZero();
    FileBoundaryCondition_z_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
                                         this->GridY4D, this->GridX4D);
    FileBoundaryCondition_z_aer_f.Resize(GridS_bc_aer, GridB_bc_aer,
					 this->GridY4D, this->GridX4D);
    if (this->option_process["with_number_concentration"])
      {
		
	this->NumberBoundaryCondition_z_aer_i.Resize(GridB_bc_aer,
						     this->GridY3D,
						     this->GridX3D);
	this->NumberBoundaryCondition_z_aer_i.SetZero();
	FileNumberBoundaryCondition_z_aer_i.Resize(GridB_bc_aer,
						   this->GridY3D,
						   this->GridX3D);
	FileNumberBoundaryCondition_z_aer_f.Resize(GridB_bc_aer, 
						   this->GridY3D,
						   this->GridX3D);
      }

    // Along Y.
    BoundaryCondition_y_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
				     this->GridZ5D, this->GridY5D_interf_bc,
				     this->GridX5D);
    BoundaryCondition_y_aer_i.SetZero();
    FileBoundaryCondition_y_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
                                         this->GridZ5D,
                                         this->GridY5D_interf_bc,
                                         this->GridX5D);
    FileBoundaryCondition_y_aer_f.Resize(GridS_bc_aer, GridB_bc_aer,
					 this->GridZ5D,
					 this->GridY5D_interf_bc,
					 this->GridX5D);
    if (this->option_process["with_number_concentration"])
      {
	this->NumberBoundaryCondition_y_aer_i.Resize(GridB_bc_aer,
						     this->GridZ4D,
						     this->GridY4D_interf_bc,
						     this->GridX4D);
	this->NumberBoundaryCondition_y_aer_i.SetZero();
	FileNumberBoundaryCondition_y_aer_i.Resize(GridB_bc_aer,
						   this->GridZ4D,
						   this->GridY4D_interf_bc,
						   this->GridX4D);
	FileNumberBoundaryCondition_y_aer_f.Resize(GridB_bc_aer,
						   this->GridZ4D,
						   this->GridY4D_interf_bc,
						   this->GridX4D);

      }

    // Along X.
    BoundaryCondition_x_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
				     this->GridZ5D, this->GridY5D,
				     this->GridX5D_interf_bc);
    BoundaryCondition_x_aer_i.SetZero();
    FileBoundaryCondition_x_aer_i.Resize(GridS_bc_aer, GridB_bc_aer,
					 this->GridZ5D,
					 this->GridY5D,
					 this->GridX5D_interf_bc);
    FileBoundaryCondition_x_aer_i.SetZero();
    FileBoundaryCondition_x_aer_f.Resize(GridS_bc_aer, GridB_bc_aer,
					 this->GridZ5D,
					 this->GridY5D,
					 this->GridX5D_interf_bc);
    FileBoundaryCondition_x_aer_f.SetZero();
    if (this->option_process["with_number_concentration"])
      {
	this->NumberBoundaryCondition_x_aer_i.Resize(GridB_bc_aer, 
						     this->GridZ4D,
						     this->GridY4D,
						     this->GridX4D_interf_bc);
	this->NumberBoundaryCondition_x_aer_i.SetZero();
	FileNumberBoundaryCondition_x_aer_i.Resize(GridB_bc_aer,
						   this->GridZ4D,
						   this->GridY4D,
						   this->GridX4D_interf_bc);
	FileNumberBoundaryCondition_x_aer_f.Resize(GridB_bc_aer,
						   this->GridZ4D,
						   this->GridY4D,
						   this->GridX4D_interf_bc);
      }

    if(bc_format == "Internal" && 
       this->option_process["with_external_composition"])
      {
	BoundaryCondition_z_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridY4D, this->GridX4D);
	FileBoundaryCondition_z_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridY4D, this->GridX4D);
	FileBoundaryCondition_z_aer_f_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridY4D, this->GridX4D);
	BoundaryCondition_y_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D_interf_bc,
				     this->GridX5D);
	FileBoundaryCondition_y_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D_interf_bc,
				     this->GridX5D);
	FileBoundaryCondition_y_aer_f_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D_interf_bc,
				     this->GridX5D);
	BoundaryCondition_x_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D,
				     this->GridX5D_interf_bc);
	FileBoundaryCondition_x_aer_i_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D,
				     this->GridX5D_interf_bc);
	FileBoundaryCondition_x_aer_f_tmp.Resize(GridS_bc_aer, GridB_bc_aer_i,
				     this->GridZ5D, this->GridY5D,
				     this->GridX5D_interf_bc);

        if (this->option_process["with_number_concentration"])
	  {
	    NumberBoundaryCondition_z_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridY3D,this->GridX3D);
	    FileNumberBoundaryCondition_z_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridY3D,this->GridX3D);
	    FileNumberBoundaryCondition_z_aer_f_tmp.Resize(GridB_bc_aer_i,
						     this->GridY3D,this->GridX3D);
	    NumberBoundaryCondition_y_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D_interf_bc,this->GridX4D);
	    FileNumberBoundaryCondition_y_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D_interf_bc,this->GridX4D);
	    FileNumberBoundaryCondition_y_aer_f_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D_interf_bc,this->GridX4D);
	    NumberBoundaryCondition_x_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D,this->GridX4D_interf_bc);
	    FileNumberBoundaryCondition_x_aer_i_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D,this->GridX4D_interf_bc);
	    FileNumberBoundaryCondition_x_aer_f_tmp.Resize(GridB_bc_aer_i,
						     this->GridZ4D,this->GridY4D,this->GridX4D_interf_bc);
	  }
      }

    if (ic_format == "Internal" &&
        this->option_process["with_external_composition"])
      {
	Concentration_aer_i.Resize(GridS5D_aer, GridB5D_aer_i,
				   GridZ5D, GridY5D, GridX5D);
	Concentration_aer_i.SetZero();
	NumberConcentration_aer_i.Resize(GridB4D_aer_i,
				   this->GridZ4D, this->GridY4D, this->GridX4D);
	NumberConcentration_aer_i.SetZero();
      }

    /*** Loss terms ***/

    DepositionVelocity_aer_i.Resize(GridB_dep_aer,
				    this->GridY3D, this->GridX3D);
    DepositionVelocity_aer_i.SetZero();
    DepositionVelocity_aer_f.Resize(GridB_dep_aer,
				    this->GridY3D, this->GridX3D);
    DepositionVelocity_aer_f.SetZero();
    FileDepositionVelocity_aer_i.Resize(GridB_dep_aer,
                                        this->GridY3D, this->GridX3D);
    FileDepositionVelocity_aer_f.Resize(GridB_dep_aer,
                                        this->GridY3D, this->GridX3D);

    ScavengingCoefficient_aer_i.Resize(GridB_scav_aer, this->GridZ4D,
				       this->GridY4D, this->GridX4D);
    ScavengingCoefficient_aer_i.SetZero();
    ScavengingCoefficient_aer_f.Resize(GridB_scav_aer, this->GridZ4D,
				       this->GridY4D, this->GridX4D);
    ScavengingCoefficient_aer_f.SetZero();

    /*** Source terms ***/

    // Surface emissions.
    SurfaceEmission_aer_i.Resize(GridS_surf_emis_aer, GridB_surf_emis_aer,
				 this->GridY4D, this->GridX4D);
    SurfaceEmission_aer_i.SetZero();
    SurfaceEmission_aer_f.Resize(GridS_surf_emis_aer, GridB_surf_emis_aer,
				 this->GridY4D, this->GridX4D);
    SurfaceEmission_aer_f.SetZero();
    FileSurfaceEmission_aer_i.Resize(GridS_surf_emis_aer, GridB_surf_emis_aer,
				     this->GridY4D, this->GridX4D);
    FileSurfaceEmission_aer_i.SetZero();
    FileSurfaceEmission_aer_f.Resize(GridS_surf_emis_aer, GridB_surf_emis_aer,
				     this->GridY4D, this->GridX4D);
    FileSurfaceEmission_aer_f.SetZero();
    if (this->option_process["with_number_concentration"])
      {
		
	NumberSurfaceEmission_aer_i.Resize(GridB_surf_emis_aer,
					   this->GridY3D, this->GridX3D);
	NumberSurfaceEmission_aer_i.SetZero();
	NumberSurfaceEmission_aer_f.Resize(GridB_surf_emis_aer,
					   this->GridY3D, this->GridX3D);
	NumberSurfaceEmission_aer_f.SetZero();
	FileNumberSurfaceEmission_aer_i.Resize(GridB_surf_emis_aer,
					       this->GridY3D, this->GridX3D);
	FileNumberSurfaceEmission_aer_i.SetZero();
	FileNumberSurfaceEmission_aer_f.Resize(GridB_surf_emis_aer,
					       this->GridY3D, this->GridX3D);
	FileNumberSurfaceEmission_aer_f.SetZero();
      }
	
    // Volume emissions.
    VolumeEmission_aer_i.Resize(GridS_vol_emis_aer, GridB_vol_emis_aer,
				this->GridZ_vol_emis_aer, this->GridY5D,
				this->GridX5D);
    VolumeEmission_aer_i.SetZero();
    VolumeEmission_aer_f.Resize(GridS_vol_emis_aer, GridB_vol_emis_aer,
				this->GridZ_vol_emis_aer, this->GridY5D,
				this->GridX5D);
    VolumeEmission_aer_f.SetZero();

    FileVolumeEmission_aer_i.Resize(GridS_vol_emis_aer, GridB_vol_emis_aer,
				    this->GridZ_vol_emis_aer, this->GridY5D,
				    this->GridX5D);
    FileVolumeEmission_aer_i.SetZero();
    FileVolumeEmission_aer_f.Resize(GridS_vol_emis_aer, GridB_vol_emis_aer,
				    this->GridZ_vol_emis_aer, this->GridY5D,
				    this->GridX5D);
    FileVolumeEmission_aer_i.SetZero();

    if(this->option_process["with_number_concentration"])
      {
		
	NumberVolumeEmission_aer_i.Resize(GridB_vol_emis_aer, this->GridZ_vol_emis_aer,
					  this->GridY4D, this->GridX4D);
	NumberVolumeEmission_aer_i.SetZero();
	NumberVolumeEmission_aer_f.Resize(GridB_vol_emis_aer, this->GridZ_vol_emis_aer,
					  this->GridY4D,this->GridX4D);
	NumberVolumeEmission_aer_f.SetZero();
	FileNumberVolumeEmission_aer_i.Resize(GridB_vol_emis_aer, this->GridZ_vol_emis_aer,
					      this->GridY4D, this->GridX4D);
	FileNumberVolumeEmission_aer_i.SetZero();
	FileNumberVolumeEmission_aer_f.Resize(GridB_vol_emis_aer, this->GridZ_vol_emis_aer,
					      this->GridY4D, this->GridX4D);
	FileNumberVolumeEmission_aer_f.SetZero();
      }
	
    /*** Sources (source splitting) ***/

    if (this->source_splitting && this->option_process["with_chemistry"])
      {
        Source_aer_i.Resize(GridS5D_aer, GridB5D_aer,
                            this->GridZ5D, this->GridY5D, this->GridX5D);
        Source_aer_f.Resize(GridS5D_aer, GridB5D_aer,
                            this->GridZ5D, this->GridY5D, this->GridX5D);
      }

    /*** State ***/

    this->Concentration_aer.Resize(GridS5D_aer, GridB5D_aer,
				   GridZ5D, GridY5D, GridX5D);
    if (this->option_process["with_number_concentration"])
      {
	this->NumberConcentration_aer.Resize(GridB4D_aer, this->GridZ4D,
					     this->GridY4D, this->GridX4D);
      }

    /*** Deposition fluxes ***/

    if (this->option_process["collect_dry_flux_aer"])
      {
	DryDepositionFlux_aer.Resize(GridS4D_aer, GridB_dep_aer, this->GridY4D,
				     this->GridX4D);
	if (this->option_process["with_number_concentration"])
	  DryDepositionFluxNumber_aer.Resize(GridB_dep_aer, this->GridY3D,
					     this->GridX3D);
      }
	
    if (this->option_process["collect_wet_flux_aer"])
      {
	WetDepositionFlux_aer.Resize(GridS4D_aer, GridB_scav_aer,
				     this->GridY4D, this->GridX4D);
	InCloudWetDepositionFlux_aer.Resize(GridS4D_aer, GridB_scav_aer,
					    this->GridY4D, this->GridX4D);
	if (this->option_process["with_number_concentration"])
	  {
	    WetDepositionFluxNumber_aer.Resize(GridB_scav_aer,
					       this->GridY3D, this->GridX3D);
	    InCloudWetDepositionFluxNumber_aer.Resize(GridB_scav_aer,
						      this->GridY3D, this->GridX3D);
	  }
		
      }

    /*** Aerosol physical properties ***/

    // Density per bin in kg / m^3.
    Density_aer.resize(this->Nbin_aer);

    WetDiameter_aer.Resize(GridB4D_aer, this->GridZ4D, this->GridY4D,
                           this->GridX4D);

    pH.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    // Default cloud pH.
    pH.Fill(4.5);

    /*** Radiative computation (photolysis rates) ***/

    if (this->computed_photolysis == "online")
      {
        NLegendre = 8;

        GridWavelength = RegularGrid<T>(Nwavelength);
        for (int i = 0; i < Nwavelength; i++)
          GridWavelength(i) = FastJ_wavelength[i];
        GridIndexReal =  RegularGrid<T>(tabulation_index_real);
        GridIndexImaginary =  RegularGrid<T>(tabulation_index_imaginary);
        GridIndexDiameter =  RegularGrid<T>(index_diameter);
        GridLegendre =  RegularGrid<T>(NLegendre);
        PureSpeciesIndexReal.Resize(GridS5D_aer_noH2O, GridWavelength);
        PureSpeciesIndexImaginary.Resize(GridS5D_aer_noH2O, GridWavelength);
        WaterIndexReal.Resize(GridWavelength);
        WaterIndexImaginary.Resize(GridWavelength);
        // For radiative data calculation, H2O is treated apart.
        OPACNames.Resize(GridS5D_aer_noH2O);
        SpeciesNames_polyphemus_OPAC.Resize(GridS5D_aer_noH2O);
        OPACNames_tmp.Resize(GridS5D_aer_noH2O);
        SpeciesNames_polyphemus_OPAC_tmp.Resize(GridS5D_aer_noH2O);
        AbsorptionEfficiencyFactorTable.Resize(GridWavelength,
                                               GridIndexReal,
                                               GridIndexImaginary,
                                               GridIndexDiameter);
        ExtinctionEfficiencyFactorTable.Resize(GridWavelength,
                                               GridIndexReal,
                                               GridIndexImaginary,
                                               GridIndexDiameter);
        PhaseFunctionTable.Resize(GridWavelength, GridIndexReal,
                                  GridIndexImaginary, GridIndexDiameter,
                                  GridLegendre);
        OpticalDepthAerosol.Resize(GridWavelength, this->GridZ4D,
                                   this->GridY4D, this->GridX4D);
        SingleScatteringAlbedo.Resize(GridWavelength, this->GridZ4D,
                                      this->GridY4D, this->GridX4D);
        MeanExtinctionEfficiencyFactor.Resize(GridWavelength,
                                              this->GridZ4D,
                                              this->GridY4D,
                                              this->GridX4D);
        MeanAbsorbtionEfficiencyFactor.Resize(GridWavelength,
                                              this->GridZ4D,
                                              this->GridY4D,
                                              this->GridX4D);
        PhaseFunction.Resize(GridWavelength, this->GridZ4D,
                             this->GridY4D, this->GridX4D,
                             GridLegendre);
        this->Concentration_aer.Resize(GridS5D_aer, GridB5D_aer,
                                       GridZ5D, GridY5D, GridX5D);
        photolysis_reaction_index.resize(this->Nr_photolysis);
      }

  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Init()
  {
    Polair3DTransport<T, ClassAdvection, ClassDiffusion>::Init();

    /*** Aerosol density ***/

    // Density per bin in kg / m^3.
    Density_aer = fixed_density_aer;

    /*** Land data ***/

    if (this->option_process["compute_deposition_aer"])
      {
        string land_data_file, LUC_file;
        this->config.SetSection("[data]");
        this->config.PeekValue("Land_data", land_data_file);
        ConfigStream land_data_stream(land_data_file);
        // Land use cover.
        land_data_stream.SetSection("[LUC]");
        land_data_stream.PeekValue("Nland", "> 0", Nland);
        LUC.Resize(Nland, this->Ny, this->Nx);
        land_data_stream.PeekValue("LUC", LUC_file);
        FormatBinary<float>().Read(LUC_file, LUC);
        // Land data.
        // Reads roughness heights for all land categories and all five
        // seasons.
        land_data_stream.SetSection("[roughness]");
        LandRoughnessHeight.resize(Nland, 5);
        for (int c = 0; c < Nland; c++)
          for (int s = 0; s < 5; s++)
            land_data_stream.GetNumber(LandRoughnessHeight(c, s));
        // Reads characteristic radius of small receptors for all land
        // categories and all five seasons.
        land_data_stream.SetSection("[small_radius]");
        SmallRadius.resize(Nland, 5);
        for (int c = 0; c < Nland; c++)
          for (int s = 0; s < 5; s++)
            land_data_stream.GetNumber(SmallRadius(c, s));
        // Reads characteristic radius of large receptors for all land
        // categories and all five seasons.
        land_data_stream.SetSection("[large_radius]");
        LargeRadius.resize(Nland, 5);
        for (int c = 0; c < Nland; c++)
          for (int s = 0; s < 5; s++)
            land_data_stream.GetNumber(LargeRadius(c, s));
        // Reads Zhang coefficients alpha and gamma.
        land_data_stream.SetSection("[zhang_alpha]");
        ZhangAlpha.resize(Nland);
        for (int c = 0; c < Nland; c++)
          land_data_stream.GetNumber(ZhangAlpha(c));
        land_data_stream.SetSection("[zhang_gamma]");
        ZhangGamma.resize(Nland);
        for (int c = 0; c < Nland; c++)
          land_data_stream.GetNumber(ZhangGamma(c));
      }

    /*** Input fields ***/

    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::InitAllData();

    Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::InitAllData();

    /*** Initial conditions for aerosols ***/

    if (this->option_manage["initial_condition_aer"])
      {
        string species, species_bin;
        vector<int> isize_section;

	this->Concentration_aer.SetZero();
	if (ic_format=="Internal" && 
            this->option_process["with_external_composition"])
          {//necessary to translate into external composition
            NumberConcentration_aer_i.SetZero();
          }
	  
	for (int i = 0; i < Ns_ic_aer; i++)
	  {
	    species = species_list_ic_aer[i].first;
	    isize_section = species_list_ic_aer[i].second;
	    for (int j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["initial_condition_aer"](species_bin);
		int index = this->GetSpeciesIndex_aer(species);
		if(ic_format=="Internal" && 
                   this->option_process["with_external_composition"])
                  {
                    Data<T, 3> Concentration_tmp(&Concentration_aer_i(index, isize_section[j],0, 0, 0),
                                                 shape(this->Nz, this->Ny, this->Nx));
                    if (is_num(filename))
                      Concentration_tmp.Fill(to_num<T>(filename));
                    else
                      FormatBinary<float>().Read(filename, Concentration_tmp);
                  }
		else
                  {
                    Data<T, 4> Concentration_tmp(&this->Concentration_aer(index, isize_section[j]*this->Ncomposition_aer,0, 0, 0),
                                                 shape(this->Ncomposition_aer,this->Nz, this->Ny, this->Nx));
                    if (is_num(filename))
                      Concentration_tmp.Fill(to_num<T>(filename));
                    else
                      FormatBinary<float>().Read(filename, Concentration_tmp);
                  }
	      }
	  }

	if (this->option_process["with_initial_condition_number_aer"])
	  {
	    this->NumberConcentration_aer.SetZero();
	    for (int i = 0; i < Nb_ic_aer; i++)
	      {
		species_bin = string("Number_") + to_str(ic_bin_list_aer[i]);
		string filename= this->input_files["initial_condition_aer"](species_bin);//problem of None

		if (exists(filename))
		  {
		    if(ic_format=="Internal"&&this->option_process["with_external_composition"])
		    {
		      Data<T, 3> NumberConcentration_tmp(&NumberConcentration_aer_i(ic_bin_list_aer[i],0, 0, 0),
                                                         shape(this->Nz, this->Ny, this->Nx));
                      if (is_num(filename))
                        NumberConcentration_tmp.Fill(to_num<T>(filename));
                      else
                        FormatBinary<float>().Read(filename, NumberConcentration_tmp);
		    }
		    else
                      {
                        Data<T, 4> NumberConcentration_tmp(&this->NumberConcentration_aer(ic_bin_list_aer[i]*this->Ncomposition_aer,0, 0, 0),
                                                           shape(this->Ncomposition_aer,this->Nz, this->Ny, this->Nx));
                        if (is_num(filename))
                          NumberConcentration_tmp.Fill(to_num<T>(filename));
                        else
                          FormatBinary<float>().Read(filename, NumberConcentration_tmp);
                      }
		  }
	      }
	  }

        if (ic_format=="Internal"&&this->option_process["with_external_composition"])
	  {//for the computation of composition within each bins, loop must stated by bins
	    for (int j = 0; j < this->Nsize_section_aer; j++)
              {
                Data<T, 3> Conc_total(this->GridZ3D, this->GridY3D,this->GridX3D);
                Conc_total.SetZero();
                Data<T, 4> Conc_group(GridG4D_aer,this->GridZ4D, this->GridY4D,this->GridX4D);
                Conc_group.SetZero();
                //compute the total bin/group mass
                for (int i = 0; i < Ns_ic_aer; i++)
                  {
                    species = species_list_ic_aer[i].first;
                    if (species != "PH2O") // LWC is not used to determine aerosol composition index. (YK)
                      {
                        int index_species = this->GetSpeciesIndex_aer(species);
                        int index_group = this->aerosol_species_group_relation(index_species);
                        for(int z = 0; z < this->Nz; z++)
                          for(int y = 0; y < this->Ny; y++)
                            for(int x = 0; x < this->Nx; x++)
                              {
                                Conc_total(z,y,x)+=Concentration_aer_i(index_species,j,z,y,x);
                                Conc_group(index_group,z, y,x)+=Concentration_aer_i(index_species,j,z,y,x);
                              }
                      }
                  }
                //redistribute based on new composition id
                for(int z = 0; z < this->Nz; z++)
                  for(int y = 0; y < this->Ny; y++)
                    for(int x = 0; x < this->Nx; x++)
                      {//for each size section only one possible composition id
                        int composition_id=FindCompositionID(Conc_total,Conc_group,composition_bounds,
                                                             this->Ncomposition_aer,this->Ngroup_aer,z,y,x);
                        int id_bin=j*this->Ncomposition_aer+composition_id;
                        for (int i = 0; i < Ns_ic_aer; i++)
                          {
                            species = species_list_ic_aer[i].first;
                            int index_species = this->GetSpeciesIndex_aer(species);
                            this->Concentration_aer(index_species,id_bin,z,y,x)=Concentration_aer_i(index_species,j,z,y,x);
                          }
                        if (this->option_process["with_initial_condition_number_aer"])
                          {
                            this->NumberConcentration_aer(id_bin,z,y,x)=NumberConcentration_aer_i(j,z,y,x);
                          }
                      }
	      }
          }
	if (this->option_process["with_number_concentration"]&&
            !this->option_process["with_initial_condition_number_aer"])
	  for (int i = 0; i < Nb_ic_aer; i++)
            {
              this->ComputeNumberConcentration_forIC_aer(int(ic_bin_list_aer[i]));
            }
      }
    
    // if (this->option_process["with_number_concentration"]) 
    //   {
    //     this->NumberConcentration_aer.PrintInfo();
    //   }


    /*** Chemical mechanism ***/
    if (this->option_process["with_chemistry"])
      this->Chemistry_.Init(*this);

    /*** Optical properties and tabulation ***/

    // Reads correspondance between Polyphemus and OPAC species.
    if (this->computed_photolysis == "online")
      {
        black_carbon_index = -999;

        read_refractiveindex_tabulation(this->Ns_aer,
                                        file_species_polyphemus_OPAC,
                                        file_water_refractive_index,
                                        directory_OPAC,
                                        black_carbon_index,
                                        Nwater_wavelength,
                                        N_OPAC_wavelength,
                                        GridWavelength,
                                        PureSpeciesIndexReal,
                                        PureSpeciesIndexImaginary,
                                        WaterIndexReal,
                                        WaterIndexImaginary);

        read_Mie_tabulation(directory_efficiency_factor,
                            GridIndexReal, GridIndexImaginary,
                            GridIndexDiameter, GridWavelength,
                            AbsorptionEfficiencyFactorTable,
                            ExtinctionEfficiencyFactorTable,
                            PhaseFunctionTable);

        // Initializes photolysis rates.
        Radiatif(this->FilePhotolysis_i, this->current_date);

        previous_date_radiatif = this->current_date;
        next_date_radiatif = this->current_date;
        next_date_radiatif.AddSeconds
          (time_step_for_computing_photolysis_rates);
        Radiatif(this->FilePhotolysis_f, next_date_radiatif);

        this->PhotolysisRate_i = this->FilePhotolysis_i;

        this->Interpolate(previous_date_radiatif,
                          time_step_for_computing_photolysis_rates,
                          0, this->FilePhotolysis_i, this->FilePhotolysis_f,
                          this->next_date, this->PhotolysisRate_f);
      }

    /*** Deposition ***/

    int Nc=this->Ncomposition_aer;
    if (this->option_process["collect_dry_flux_aer"])
      for (int b = 0; b < Nbin_dep_aer; b++)
	for(int id = 0; id< Nc; id++)
	for (int j = 0; j < this->Ny; j++)
	  for (int i = 0; i < this->Nx; i++)
	    {
	      for (int s = 0; s < this->Ns_aer; s++)
		DryDepositionFlux_aer(s, b*Nc+id, j, i) =
		  DepositionVelocity_aer_f(b*Nc+id, j, i)
		  * this->Concentration_aer(s, bin_list_dep_aer[b]*Nc+id, 0, j, i);
	      if (this->option_process["with_number_concentration"])
		DryDepositionFluxNumber_aer(b*Nc+id, j, i) =
		  DepositionVelocity_aer_f(b*Nc+id, j, i)
		  * this->NumberConcentration_aer(bin_list_dep_aer[b]*Nc+id, 0, j, i);
	    }

    if (this->option_process["collect_wet_flux_aer"])
      {
	this->WetDepositionFlux_aer.SetZero();
	this->InCloudWetDepositionFlux_aer.SetZero();
	if (this->option_process["with_number_concentration"])
	  {
	    this->WetDepositionFluxNumber_aer.SetZero();
	    this->InCloudWetDepositionFlux_aer.SetZero();
	  }
	int Nc=this->Ncomposition_aer;
	for (int b = 0; b < Nbin_scav_aer; b++)
	  for(int id=0; id < Nc; id++)
	  for (int j = 0; j < this->Ny; j++)
	    for (int i = 0; i < this->Nx; i++)
	      if (this->Rain_f(j, i) > 0.)
		for (int k = 0; k < this->Nz; k++)
		  {
		    for (int s = 0; s < this->Ns_aer; s++)
		      if (this->CloudBaseHeight_f(j, i) < this->GridZ4D(k))
			WetDepositionFlux_aer(s, b *Nc +id, j, i) =
			  this->Concentration_aer(s, bin_list_scav_aer[b] * Nc + id, k, j, i)
			  / this->Delta_t
			  * (1.- exp(- this->Delta_t
				     * ScavengingCoefficient_aer_f(b*Nc+id, k, j, i)))
			  * (this->GridZ4D_interf(k+1) - this->GridZ4D_interf(k));

		    if (this->option_process["with_number_concentration"])
		      if (this->CloudBaseHeight_f(j, i) < this->GridZ4D(k))
			WetDepositionFluxNumber_aer(b * Nc + id, j, i) =
			  this->NumberConcentration_aer(bin_list_scav_aer[b] * Nc +id, k, j, i)
			  / this->Delta_t
			  * (1.- exp(- this->Delta_t
				     * ScavengingCoefficient_aer_f(b*Nc+id, k, j, i)))
			  * (this->GridZ4D_interf(k+1) - this->GridZ4D_interf(k));
		  }
      }
  }


  //! Model initialization for each step.
  /*! It reads on file the data needed for the current step.
   */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitStep()
  {
    if (this->data_gone_through_initstep
        && this->data_date == this->current_date)
      return;

    int i, j;
    string species, species_bin;
    vector<int> isize_section;
    string filename;
    Date date;
    T Delta_t;

    Polair3DChemistry < T, ClassAdvection, ClassDiffusion,
                        ClassChemistry >::InitStep();

    /*** Surface temperature ***/

    if (this->option_manage["surface_temperature"])
      this->UpdateData("meteo", "SurfaceTemperature",
                       FileSurfaceTemperature_i,
                       FileSurfaceTemperature_f, SurfaceTemperature_i,
                       SurfaceTemperature_f);

    /*** Surface pressure ***/

    if (this->option_manage["surface_pressure"])
      this->UpdateData("meteo", "SurfacePressure", FileSurfacePressure_i,
                       FileSurfacePressure_f, SurfacePressure_i,
                       SurfacePressure_f);

    /*** Snow height ***/

    if (this->option_manage["snow_height"])
      this->UpdateData("meteo", "SnowHeight", FileSnowHeight_i,
                       FileSnowHeight_f, SnowHeight_i, SnowHeight_f);

    /*** Liquid water content ***/

    if (this->option_manage["liquid_water_content"])
      this->UpdateData("meteo", "LiquidWaterContent",
                       FileLiquidWaterContent_i, FileLiquidWaterContent_f,
                       LiquidWaterContent_i);

    /*** CLoud and Ice optical depth ***/

    if (this->computed_photolysis == "online")
      {
        this->UpdateData("meteo", "CloudOpticalDepth",
                         FileCloudOpticalDepth_i, FileCloudOpticalDepth_f,
                         CloudOpticalDepth_i);
        this->UpdateData("meteo", "IceOpticalDepth",
                         FileIceOpticalDepth_i, FileIceOpticalDepth_f,
                         IceOpticalDepth_i);
      }


    /*** Boundary conditions ***/

    if (this->option_manage["boundary_condition_aer"])
		{		

	if(bc_format=="Internal"&&this->option_process["with_external_composition"])
	{
	  BoundaryCondition_z_aer_i.SetZero();
	  BoundaryCondition_y_aer_i.SetZero();
	  BoundaryCondition_x_aer_i.SetZero();
	  NumberBoundaryCondition_z_aer_i.SetZero();
	  NumberBoundaryCondition_z_aer_i_tmp.SetZero();
	  NumberBoundaryCondition_y_aer_i.SetZero();
	  NumberBoundaryCondition_y_aer_i_tmp.SetZero();
	  NumberBoundaryCondition_x_aer_i.SetZero();
	  NumberBoundaryCondition_x_aer_i_tmp.SetZero();
	}	

      for (i = 0; i < Ns_bc_aer; i++)
        {
          species = species_list_bc_aer[i].first;
          isize_section = species_list_bc_aer[i].second;
          for (j = 0; j < int(isize_section.size()); j++)
            {
              species_bin = species + string("_") + to_str(isize_section[j]);
              string filename
                = this->input_files["boundary_condition_aer"](species_bin);
              date = this->input_files["boundary_condition_aer"].GetDateMin();
              Delta_t
                = this->input_files["boundary_condition_aer"].GetDelta_t();
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->UpdateData(find_replace(filename, "&c", "z"), date,
				Delta_t, FileBoundaryCondition_z_aer_i_tmp,
				FileBoundaryCondition_z_aer_f_tmp, i, j,
				BoundaryCondition_z_aer_i_tmp);
		  this->UpdateData(find_replace(filename, "&c", "y"), date,
				Delta_t, FileBoundaryCondition_y_aer_i_tmp,
				FileBoundaryCondition_y_aer_f_tmp, i,  j,
				BoundaryCondition_y_aer_i_tmp);
		  this->UpdateData(find_replace(filename, "&c", "x"), date,
				Delta_t, FileBoundaryCondition_x_aer_i_tmp,
				FileBoundaryCondition_x_aer_f_tmp, i,  j,
				BoundaryCondition_x_aer_i_tmp);
		}
		else
		{
		    this->UpdateData(find_replace(filename, "&c", "z"), date,
				  Delta_t, FileBoundaryCondition_z_aer_i,
				  FileBoundaryCondition_z_aer_f, i, j,
				  BoundaryCondition_z_aer_i,this->Ncomposition_aer);
		    this->UpdateData(find_replace(filename, "&c", "y"), date,
				  Delta_t, FileBoundaryCondition_y_aer_i,
				  FileBoundaryCondition_y_aer_f, i,  j,
				  BoundaryCondition_y_aer_i, this->Ncomposition_aer);
		    this->UpdateData(find_replace(filename, "&c", "x"), date,
				  Delta_t, FileBoundaryCondition_x_aer_i,
				  FileBoundaryCondition_x_aer_f, i,  j,
				  BoundaryCondition_x_aer_i,this->Ncomposition_aer);
		}
	      }
	  }
	  
	if (this->option_process["with_number_concentration"])
	  for (i = 0; i < Nb_bc_aer; i++)
	    {
	      species_bin = "Number_" + to_str(bc_bin_list_aer[i]);
	      string filename
		= this->input_files["boundary_condition_aer"](species_bin);
	      date = this->input_files["boundary_condition_aer"].GetDateMin();
	      Delta_t
		= this->input_files["boundary_condition_aer"].GetDelta_t();
	      if (exists(find_replace(filename, "&c", "z")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->UpdateData(find_replace(filename, "&c", "z"), date,
				Delta_t, FileNumberBoundaryCondition_z_aer_i_tmp,
				FileNumberBoundaryCondition_z_aer_f_tmp, i,
				NumberBoundaryCondition_z_aer_i_tmp);
		}
		else
		{
		    this->UpdateData(find_replace(filename, "&c", "z"), date,
				  Delta_t, FileNumberBoundaryCondition_z_aer_i,
				  FileNumberBoundaryCondition_z_aer_f, i,
				  NumberBoundaryCondition_z_aer_i,this->Ncomposition_aer);
		}
	      }
			  
	      if (exists(find_replace(filename, "&c", "y")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->UpdateData(find_replace(filename, "&c", "y"), date,
				 Delta_t, FileNumberBoundaryCondition_y_aer_i_tmp,
				 FileNumberBoundaryCondition_y_aer_f_tmp, i,
				 NumberBoundaryCondition_y_aer_i_tmp);		  
		}
		else
		{
		    this->UpdateData(find_replace(filename, "&c", "y"), date,
				  Delta_t, FileNumberBoundaryCondition_y_aer_i,
				  FileNumberBoundaryCondition_y_aer_f, i,
				  NumberBoundaryCondition_y_aer_i,this->Ncomposition_aer);
		}
	      }

	      if (exists(find_replace(filename, "&c", "x")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->UpdateData(find_replace(filename, "&c", "x"), date,
				 Delta_t, FileNumberBoundaryCondition_x_aer_i_tmp,
				 FileNumberBoundaryCondition_x_aer_f_tmp, i,
				 NumberBoundaryCondition_x_aer_i_tmp);		  
		}
		else
		{
		    this->UpdateData(find_replace(filename, "&c", "x"), date,
				  Delta_t, FileNumberBoundaryCondition_x_aer_i,
				  FileNumberBoundaryCondition_x_aer_f, i,
				  NumberBoundaryCondition_x_aer_i,this->Ncomposition_aer);
		}
	      }
	    }
	    
	if(bc_format=="Internal"&&this->option_process["with_external_composition"])
	{//for the computation of composition within each bins, loop must stated by bins
	  BoundaryConditionTransformation();
	}
	if (this->option_process["with_number_concentration"])
	  for (i = 0; i < Nb_bc_aer; i++)
	    {
	      species_bin = "Number_" + to_str(bc_bin_list_aer[i]);
	      string filename
		= this->input_files["boundary_condition_aer"](species_bin);
	      if (!exists(find_replace(filename, "&c", "z")))
		this->ComputeNumberBoundaryCondition_z_aer(int(bc_bin_list_aer[i]));

	      if (!exists(find_replace(filename, "&c", "y")))
		this->ComputeNumberBoundaryCondition_y_aer(int(bc_bin_list_aer[i]));

	      if (!exists(find_replace(filename, "&c", "x")))
		this->ComputeNumberBoundaryCondition_x_aer(int(bc_bin_list_aer[i]));
	    }
      }

    /*** Surface emissions ***/

    if (this->option_manage["surface_emission_aer"])
      {
	date = this->input_files["surface_emission_aer"].GetDateMin();
	Delta_t
	  = this->input_files["surface_emission_aer"].GetDelta_t();	
	for (i = 0; i < Ns_surf_emis_aer; i++)
	  {
	    species = species_list_surf_emis_aer[i].first;
	    isize_section= species_list_surf_emis_aer[i].second;
	    int CompositionID=FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["surface_emission_aer"](species_bin);
		if(surface_emis_format=="Internal" &&
                   this->option_process["with_external_composition"])
		{		  
		  this->UpdateData(filename, date, Delta_t,
				    FileSurfaceEmission_aer_i,
				    FileSurfaceEmission_aer_f,
				    i, j,CompositionID,
				    SurfaceEmission_aer_i,
				    SurfaceEmission_aer_f);
		}
		else
		{
		  this->UpdateData(filename, date, Delta_t,
				  FileSurfaceEmission_aer_i,
				  FileSurfaceEmission_aer_f, i, j,
				  SurfaceEmission_aer_i,
				  SurfaceEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	  }

	if (this->option_process["with_number_concentration"])
	  for( j = 0; j < Nb_surf_emis_aer ; j++)
	    {
	      species_bin = "Number_"  + to_str(surface_emission_bin_list_aer[j]);
	      string filename = this->input_files["surface_emission_aer"](species_bin);
			
	      if (exists(filename))
	      {
		if(surface_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  TinyVector<int, 2> new_shape;
		  for (int i = 0; i < 2; i++)
		    new_shape(i) = SurfaceEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 2> FileData_extract_i(new_shape);
		  Data<T, 2> FileData_extract_f(new_shape);
		  Data<T, 2> CurrentData_extract_i(new_shape);
		  Data<T, 2> CurrentData_extract_f(new_shape);
		  this->UpdateData(filename, date, Delta_t,
				FileData_extract_i,
				FileData_extract_f,
				CurrentData_extract_i,
				CurrentData_extract_f);
		  for(int y=0; y<this->Ny; y++)
		    for(int x=0; x<this->Nx; x++)
		    {
		      T TotalMass=0;
		      Data<T, 1> CompositionMass(this->Ncomposition_aer);
		      CompositionMass.SetZero();
		      for(int s=0; s<Ns_surf_emis_aer; s++)
		      {
			species = species_list_surf_emis_aer[s].first;
			int CompositionID=FindExternalCompositionID(species);
			int RealID=this->Ncomposition_aer*j+CompositionID;
			TotalMass+=SurfaceEmission_aer_f(s,RealID,y,x);
			CompositionMass(CompositionID)+=SurfaceEmission_aer_f(s,RealID,y,x);
		      }
		      for(int id=0; id<this->Ncomposition_aer; id++)
		      {
			if(TotalMass>0)
			{
			  int RealID=this->Ncomposition_aer*j+id;
			  NumberSurfaceEmission_aer_i(RealID,y,x)=CurrentData_extract_i(y,x)*CompositionMass(id)/TotalMass;
			  NumberSurfaceEmission_aer_f(RealID,y,x)=CurrentData_extract_f(y,x)*CompositionMass(id)/TotalMass;
			  FileNumberSurfaceEmission_aer_i(RealID,y,x)=FileData_extract_i(y,x)*CompositionMass(id)/TotalMass;
			  FileNumberSurfaceEmission_aer_f(RealID,y,x)=FileData_extract_f(y,x)*CompositionMass(id)/TotalMass;
			}
		      }
		    }
		}
		else
		{
		  this->UpdateData(filename, date, Delta_t,
				  FileNumberSurfaceEmission_aer_i,
				  FileNumberSurfaceEmission_aer_f, j,
				  NumberSurfaceEmission_aer_i,
				  NumberSurfaceEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	      else
	      {
		 for(int id =0 ; id < this->Ncomposition_aer; id++ )
		  for (int l = 0; l < this->Ny; l++)
		    for (int m = 0; m < this->Nx; m++)
		      NumberSurfaceEmission_aer_i(j*this->Ncomposition_aer+id, l, m) =
			NumberSurfaceEmission_aer_f(j*this->Ncomposition_aer+id, l, m);
			
		 this->ComputeNumberSurfaceEmission_aer
		  (int(surface_emission_bin_list_aer[j]));
	      }
	    }


      }


    /*** Volume emissions ***/

    if (this->option_manage["volume_emission_aer"])
      {
	date = this->input_files["volume_emission_aer"].GetDateMin();
	Delta_t
	  = this->input_files["volume_emission_aer"].GetDelta_t();
	for (i = 0; i < Ns_vol_emis_aer; i++)
	  {
	    species = species_list_vol_emis_aer[i].first;
	    isize_section = species_list_vol_emis_aer[i].second;
	    int CompositionID=FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["volume_emission_aer"](species_bin);
		if(volume_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{	
		  this->UpdateData(filename, date, Delta_t,
				  FileVolumeEmission_aer_i,
				  FileVolumeEmission_aer_f,
				  i, j,CompositionID,
				  VolumeEmission_aer_i,
				  VolumeEmission_aer_f);
		}
		else
		{
		  this->UpdateData(filename, date, Delta_t,
				  FileVolumeEmission_aer_i,
				  FileVolumeEmission_aer_f, i, j,
				  VolumeEmission_aer_i,
				  VolumeEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	  }
	if (this->option_process["with_number_concentration"])
	  for( j = 0; j < Nb_vol_emis_aer; j++)
	    {
	      species_bin = "Number_"  + to_str(volume_emission_bin_list_aer[j]);
	      string filename = this->input_files["volume_emission_aer"](species_bin);

	      if (exists(filename))
	      {
		if(volume_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  TinyVector<int, 3> new_shape;
		  for (int i = 0; i < 3; i++)
		    new_shape(i) = VolumeEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 3> FileData_extract_i(new_shape);
		  Data<T, 3> FileData_extract_f(new_shape);
		  Data<T, 3> CurrentData_extract_i(new_shape);
		  Data<T, 3> CurrentData_extract_f(new_shape);
		  this->UpdateData(filename, date, Delta_t,
				FileData_extract_i,
				FileData_extract_f,
				CurrentData_extract_i,
				CurrentData_extract_f);//Nz_vol_emis_aer
		  for(int z=0; z<Nz_vol_emis_aer; z++)
		    for(int y=0; y<this->Ny; y++)
		      for(int x=0; x<this->Nx; x++)
		      {
			T TotalMass=0;
			Data<T, 1> CompositionMass(this->Ncomposition_aer);
			CompositionMass.SetZero();
			for(int s=0; s<Ns_surf_emis_aer; s++)
			{
			  species = species_list_surf_emis_aer[s].first;
			  int CompositionID=FindExternalCompositionID(species);
			  int RealID=this->Ncomposition_aer*j+CompositionID;
			  TotalMass+=VolumeEmission_aer_f(s,RealID,z,y,x);
			  CompositionMass(CompositionID)+=VolumeEmission_aer_f(s,RealID,z,y,x);
			}
			for(int id=0; id<this->Ncomposition_aer; id++)
			{
			  if(TotalMass>0)
			  {
			    int RealID=this->Ncomposition_aer*j+id;
			    NumberVolumeEmission_aer_i(RealID,z,y,x)=CurrentData_extract_i(z,y,x)*CompositionMass(id)/TotalMass;
			    NumberVolumeEmission_aer_f(RealID,z,y,x)=CurrentData_extract_f(z,y,x)*CompositionMass(id)/TotalMass;
			    FileNumberVolumeEmission_aer_i(RealID,z,y,x)=FileData_extract_i(z,y,x)*CompositionMass(id)/TotalMass;
			    FileNumberVolumeEmission_aer_f(RealID,z,y,x)=FileData_extract_f(z,y,x)*CompositionMass(id)/TotalMass;
			  }
			}
		      }
		}
		else
		{
		this->UpdateData(filename, date, Delta_t,
				 FileNumberVolumeEmission_aer_i,
				 FileNumberVolumeEmission_aer_f, j,
				 NumberVolumeEmission_aer_i,
				 NumberVolumeEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	      else
	      {
		for(int id =0 ; id < this->Ncomposition_aer; id++ )
		  for (int k = 0; k < this->Nz_vol_emis_aer; k++)
		    for (int l = 0; l < this->Ny; l++)
		      for (int m = 0; m < this->Nx; m++)
			NumberVolumeEmission_aer_i(j*this->Ncomposition_aer+id, k, l, m) =
			  NumberVolumeEmission_aer_f(j*this->Ncomposition_aer+id, k, l, m);
			
		this->ComputeNumberVolumeEmission_aer
		  (int(volume_emission_bin_list_aer[j]));
	      }
	    }
	    
      }

    if (this->computed_photolysis == "online")
      {
        if (this->step % iteration_radiatif == 0 && this->step != 0)
          // Update of FilePhotolysis: Photolysis rate calculated using the
          // radiative model Fast-J.
          {
            this->FilePhotolysis_i = this->FilePhotolysis_f;
            previous_date_radiatif = this->current_date;
            next_date_radiatif = this->current_date;
            next_date_radiatif.AddSeconds
              (time_step_for_computing_photolysis_rates);
            Radiatif(this->FilePhotolysis_f, next_date_radiatif);
          }

        // Interpolation of photolysis rates at the
        // current date from FilePhotolysis.
        this->PhotolysisRate_i.GetArray() = this->PhotolysisRate_f.GetArray();
        this->Interpolate(previous_date_radiatif,
                          time_step_for_computing_photolysis_rates,
                          0, this->FilePhotolysis_i, this->FilePhotolysis_f,
                          this->next_date, this->PhotolysisRate_f);
      }
  }


  //! Compute density of aerosol composition 
  //! from mass concentration, density of each species
  /*!
    \param b bin number
    \return rho density (ug.um-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  T Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeDensity(Data <T, 1> Conc_aer_tmp,
		   vector<T> Rho_species, T TotalMass, int Ns)
  {
    int s;
    T rho, subrho;

    rho = 0.0;
    subrho = 0.0;
	
    for (s = 0; s < Ns; s++)
      subrho += Conc_aer_tmp(s) / Rho_species[s];
	
    if (TotalMass == 0. or subrho == 0.)
      rho = 1.;
    else
      rho = 1.e-6 * TotalMass/subrho;
	
    return rho;
	
  }
  
  //! Compute density from mass concentration, and
  //! each species density
  /*!
    \ param b bin number
    \ return density 
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberConcentration_forIC_aer(int b)
  {
    int k, j, i;
    T TotalMass, Rho_aer;
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species;
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");

    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_ic_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s =0; s < Ns_ic_aer; s++)
      {
	rho_stream_aer.PeekValue(this->species_list_ic_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }

    Data<T,1> Mean_domain;
    Mean_domain.Resize(this->Ns_ic_aer);
    Mean_domain.SetZero();
    T Mean_domain_number = 0.0;

    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
	for (i = 0; i < this->Nx; i++)
	  {
	    TotalMass = 0.0;
	    Rho_aer = 0.0;
	    for (int s = 0; s < Ns_ic_aer; s++)
		{
		  TotalMass += this->Concentration_aer
		    (this->GetSpeciesIndex_aer(species_list_ic_aer[s].first),b,k,j,i);

		  Mean_domain(s) += TotalMass / (this->Nz*this->Ny*this->Nx);
		  Conc_aer_tmp(s) = this->Concentration_aer
		    (this->GetSpeciesIndex_aer(species_list_ic_aer[s].first),b,k,j,i);

		}

	    Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_ic_aer);
	    Mean_domain_number +=
	      TotalMass/Rho_aer/PI*6.
	      /(MeanDiameter*MeanDiameter*MeanDiameter)/(this->Nz*this->Ny*this->Nx);
	  }
   T Mean_domain_number_out = 0.0;

    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
	for (i = 0; i < this->Nx; i++)
	  for(int id=0;id<this->Ncomposition_aer;id++)
	    {//KS: add a loop on Ncomposition_aer	     
	      TotalMass = 0.0;
	      Rho_aer = 0.0;
			 
	      for (int s = 0; s < Ns_ic_aer; s++)
		{
		  TotalMass += this->Concentration_aer
		    (this->GetSpeciesIndex_aer(species_list_ic_aer[s].first),b*this->Ncomposition_aer+ id,k,j,i);
		  Conc_aer_tmp(s) = this->Concentration_aer
		    (this->GetSpeciesIndex_aer(species_list_ic_aer[s].first),b*this->Ncomposition_aer+ id,k,j,i);
		}

	      Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_ic_aer);
			
	      this->NumberConcentration_aer(b*this->Ncomposition_aer+ id,k,j,i) =
		TotalMass/Rho_aer/PI*6.
		/(MeanDiameter*MeanDiameter*MeanDiameter);
	      Mean_domain_number_out += this->NumberConcentration_aer(b*this->Ncomposition_aer+ id,k,j,i)/(this->Nz*this->Ny*this->Nx);
	  }
  }
  
  
  //! Compute number boundary condition from mass concentration, diameter and
  //! density along z
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberBoundaryCondition_z_aer(int b)
  {
    
    int id_b = this->NumberBoundaryConditionIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);
    int i, j;
    T TotalMass, Rho_aer;
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species;

    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_bc_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_bc_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_bc_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
	
    for (j = 0; j < this->Ny; j++)
      for (i = 0; i < this->Nx; i++)
	{
	  TotalMass = 0.0;
	  Rho_aer = 0.0;
	  for (int ic=0; ic<Nc;ic++) //ZS
	  {
	    for (int s = 0; s < Ns_bc_aer; s++)
	      {
		it_begin = species_list_bc_aer[s].second.begin();
		it_end = species_list_bc_aer[s].second.end();

		pos = find(it_begin,it_end,b);
		dist = distance(it_begin,pos);

		if (dist < int(species_list_bc_aer[s].second.size()))
		  {
		      TotalMass = TotalMass + this->BoundaryCondition_z_aer_i(s,dist*Nc+ic,j,i);
		      Conc_aer_tmp(s) = this->BoundaryCondition_z_aer_i(s,dist*Nc+ic,j,i);
		  }
	      }
	    Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_bc_aer);
	      NumberBoundaryCondition_z_aer_i(index_b*Nc+ic,j,i) =
		TotalMass/Rho_aer/PI*6.
		/(MeanDiameter*MeanDiameter*MeanDiameter);	      
	  }
	}
  }


  //! Compute number boundary condition from mass concentration, diameter and
  //! density along y
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberBoundaryCondition_y_aer(int b)
  {
    int id_b = this->NumberBoundaryConditionIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);
    int i, j, k;
    T TotalMass, Rho_aer;
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species;
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_bc_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_bc_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_bc_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < 2; j++)
	for (i = 0; i < this->Nx; i++)
	  {
	    TotalMass = 0.0;
	    Rho_aer = 0.0;
	    for (int ic=0; ic<Nc;ic++) //ZS
	    {			 
	      for (int s = 0; s < Ns_bc_aer; s++)
		{
		  it_begin = species_list_bc_aer[s].second.begin();
		  it_end = species_list_bc_aer[s].second.end();
		  pos = find(it_begin,it_end,b);
		  dist = distance(it_begin,pos);

		  if (dist < int(species_list_bc_aer[s].second.size()))
		    {
			TotalMass = TotalMass + BoundaryCondition_y_aer_i(s,dist*Nc+ic,k,j,i);
			Conc_aer_tmp(s) = BoundaryCondition_y_aer_i(s,dist*Nc+ic,k,j,i);
		    }
		}
	      Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_bc_aer);
		NumberBoundaryCondition_y_aer_i(index_b*Nc+ic,k,j,i) =
		  TotalMass/Rho_aer/PI*6.
		  /(MeanDiameter*MeanDiameter*MeanDiameter);
	    }
	  }
  }
  

  //! Compute number boundary condition from mass concentration, diameter and
  //! density along x
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberBoundaryCondition_x_aer(int b)
  {
    int id_b = this->NumberBoundaryConditionIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);
    int i, j, k;
    T TotalMass, Rho_aer;
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species;
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;

    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_bc_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_bc_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_bc_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
	
    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
	for (i = 0; i < 2; i++)
	  {
	    TotalMass = 0.0;
	    Rho_aer = 0.0;
	    for (int ic=0; ic<Nc;ic++) //ZS
	    {
	      for (int s = 0; s < Ns_bc_aer; s++)
		{
		  it_begin = species_list_bc_aer[s].second.begin();
		  it_end = species_list_bc_aer[s].second.end();
		  pos = find(it_begin,it_end,b);
		  dist = distance(it_begin,pos);

		  if (dist < int(species_list_bc_aer[s].second.size()))
		    {
			TotalMass = TotalMass + BoundaryCondition_x_aer_i(s,dist*Nc+ic,k,j,i);
			Conc_aer_tmp(s) = BoundaryCondition_x_aer_i(s,dist*Nc+ic,k,j,i);
		    }
		}
	      Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_bc_aer);

		NumberBoundaryCondition_x_aer_i(index_b*Nc+ic,k,j,i) =
		  TotalMass/Rho_aer/PI*6.
		  /(MeanDiameter*MeanDiameter*MeanDiameter);		
	    }
	  }
  }
  

  //! Compute number surface emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberSurfaceEmission_aer(int b)
  {
    int id_b = this->NumberSurfaceEmissionIndex_aer(b*this->Ncomposition_aer);

    int index_b=Bin_to_size_index_aer(id_b);
    int i, j;
    T TotalMass; // ug.m-3
    T Rho_aer; 
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
	
    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_surf_emis_aer);
    Conc_aer_tmp.SetZero();

    for (int s = 0; s < Ns_surf_emis_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_surf_emis_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
    for (j = 0; j < this->Ny; j++)
      for (i = 0; i < this->Nx; i++)
	{
	  for (int ic=0; ic<Nc;ic++) //ZS
	  {
	    TotalMass = 0.0;
	    Rho_aer = 0.0;
	    for (int s = 0; s < Ns_surf_emis_aer; s++)
	      {
		it_begin = species_list_surf_emis_aer[s].second.begin();
		it_end = species_list_surf_emis_aer[s].second.end();
		pos = find(it_begin,it_end,b);
		dist = distance(it_begin,pos);

		if (dist < int(species_list_surf_emis_aer[s].second.size()))
		  {
		    TotalMass = TotalMass + this->SurfaceEmission_aer_f(s,dist*Nc+ic,j,i);
		    Conc_aer_tmp(s) = this->SurfaceEmission_aer_f(s,dist*Nc+ic,j,i);
		  }

	      }
	    Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass,  Ns_surf_emis_aer);

	    this->NumberSurfaceEmission_aer_f(index_b*Nc+ic,j,i) =
	      TotalMass/Rho_aer/PI*6.
	      /(MeanDiameter*MeanDiameter*MeanDiameter);
	      T tmp_n=this->NumberSurfaceEmission_aer_f(index_b*Nc+ic,j,i);
	      if(TotalMass*tmp_n==0&&TotalMass!=tmp_n)
		cout<<j<<" , "<<i<<" m:"<<TotalMass<<" n:"<<tmp_n
		<<" Rho:"<<Rho_aer<<" MeanDiameter:"<<MeanDiameter<<endl;
	  }
	}
  }

  //! Compute number volume emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::ComputeNumberVolumeEmission_aer(int b)
  {
    //here b is the index of size section
    //id_b is the index of bin
    int id_b = this->NumberVolumeEmissionIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);	
    int i, j, k;
    T TotalMass; // ug.m-3
    T Rho_aer; 
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); 
    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this-> Ns_vol_emis_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_vol_emis_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_vol_emis_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
	
    for (k = 0; k < this->Nz_vol_emis_aer; k++)
      for (j = 0; j < this->Ny; j++)
	for (i = 0; i < this->Nx; i++)
	  {
	    for (int ic=0; ic<Nc;ic++) //ZS
	    {
	      TotalMass = 0.0;
	      Rho_aer = 0.0;
	      for (int s = 0; s < Ns_vol_emis_aer; s++)
		{
		  it_begin = species_list_vol_emis_aer[s].second.begin();
		  it_end = species_list_vol_emis_aer[s].second.end();
		  pos = find(it_begin,it_end,b);
		  dist = distance(it_begin, pos);

		  if (dist < int(species_list_vol_emis_aer[s].second.size()))
		    {
		      TotalMass = TotalMass + this->VolumeEmission_aer_f(s,dist*Nc+ic,k,j,i);
		      Conc_aer_tmp(s) = this->VolumeEmission_aer_f(s,dist*Nc+ic,k,j,i);
		    }

		}

	      Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_vol_emis_aer);
	      this->NumberVolumeEmission_aer_f(index_b*Nc+ic,k,j,i)=
		TotalMass/Rho_aer/PI*6.
		/(MeanDiameter*MeanDiameter*MeanDiameter);
	      T tmp_n=this->NumberVolumeEmission_aer_f(index_b*Nc+ic,k,j,i);
	      if(TotalMass*tmp_n==0&&TotalMass!=tmp_n)
		cout<<k<<" , "<<j<<" , "<<i<<" m:"<<TotalMass<<" n:"<<tmp_n
		<<" Rho:"<<Rho_aer<<" MeanDiameter:"<<MeanDiameter<<endl;
		
	    }
	  }
  }
  
 //! Initializes wet aerosol diameters.
  /*! Computes wet aerosol diameters for aerosols using Gerber formula.
    \param RelativeHumidity_ relative humidity.
    \param Temperature_ temperature (K).
    \param WetDiameter_aer_ (output) wet aerosol diameter (m).
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitWetDiameter_aer(Data<T, 3>& RelativeHumidity_,
			Data<T, 3>& Temperature_,
			Data<T, 4>& WetDiameter_aer_)
  {	
    int i, j, k, b, id;

    Array<T, 1> MeanDiameter(this->Nsize_section_aer);
    for (b = 0; b < this->Nsize_section_aer; b++)
      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    for (b = 0; b < this->Nsize_section_aer; b++)
      for (id = 0; id < this->Ncomposition_aer; id ++)
	for (k = 0; k < this->Nz; k++)
	  for (j = 0; j < this->Ny; j++)
	    for (i = 0; i < this->Nx; i++)
	    {
	      int bin_id=b*this->Ncomposition_aer+id;
	      if(RelativeHumidity_(k, j, i)>0.0)
		_gerber_wet_diameter(&RelativeHumidity_(k, j, i),
				    &Temperature_(k, j, i), &MeanDiameter(b),
				    &WetDiameter_aer_(bin_id, k, j, i));
	      else
		WetDiameter_aer_(bin_id, k, j, i)=MeanDiameter(b);
	    }

    // Back to meters.
    WetDiameter_aer_() *= 1.e-6;
  }


  
  //! Initializes wet aerosol diameters.
  /*! Computes wet aerosol diameters for aerosols using Gerber formula.
    \param RelativeHumidity_ relative humidity.
    \param Temperature_ temperature (K).
    \param WetDiameter_aer_ (output) wet aerosol diameter (m).
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitWetDiameter_aer(Data<T, 3>& RelativeHumidity_,
			Data<T, 3>& Temperature_,
			Data<T, 5>& Concentration_aer_,
		        Data<T, 4>& NumberConcentration_aer_,
			Data<T, 4>& WetDiameter_aer_,
		        bool with_computation_drydiameter)
  {	
    int i, j, k, b, id;
    T TotalMass, Rho_aer;
    // KS Number and mass concentrations are input to the subroutine
    // in order to compute the MeanDryDiameter

    Array<T, 1> MeanDiameter(this->Nsize_section_aer);

    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_aer);
    Conc_aer_tmp.SetZero();

    for (b = 0; b < this->Nsize_section_aer; b++)
      for (id = 0; id < this->Ncomposition_aer; id ++)
        for (k = 0; k < this->Nz; k++)
          for (j = 0; j < this->Ny; j++)
            for (i = 0; i < this->Nx; i++)
              {
                int bin_id = b * this->Ncomposition_aer + id;
                if (with_computation_drydiameter 
                    and this->option_process["with_number_concentration"])
                  {

                    TotalMass = 0.0;
                    Rho_aer = 0.0;
		  
                    for (int s = 0; s < this->Ns_aer-1; s++)
                      {
                        TotalMass += this->Concentration_aer 
                          (this->GetSpeciesIndex_aer(this->species_list_aer[s]),bin_id,k,j,i);
                        Conc_aer_tmp(s) = this->Concentration_aer
                          (this->GetSpeciesIndex_aer(this->species_list_aer[s]),bin_id,k,j,i);
                      }
                    if (this->option_process["with_fixed_density"]) 
                      Rho_aer = 1.e-9*fixed_density_aer;  // convert to [kg/m3] -> [ug/um3]
                    else
                      Rho_aer = ComputeDensity(Conc_aer_tmp, Mass_Density_aer, 
                                               TotalMass, this->Ns_aer-1);
                    if (this->NumberConcentration_aer(bin_id,k,j,i) > 1.0 and TotalMass > 0.)
                      MeanDiameter(b) = pow((TotalMass/Rho_aer/PI*6.
                                             /this->NumberConcentration_aer(bin_id,k,j,i)),(1./3.)) ;
                    else
                      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));
                  }
                else
                  MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

                _gerber_wet_diameter(&RelativeHumidity_(k, j, i),
                                     &Temperature_(k, j, i), &MeanDiameter(b),
                                     &WetDiameter_aer_(bin_id, k, j, i));
              }
    // Back to meters.
    WetDiameter_aer_() *= 1.e-6;
  }


  //! Initializes deposition velocities for aerosols.
  /*! Computes deposition velocities for aerosols on the basis of Zhang
    parameterization.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitDepositionVelocity(Data<T, 3>& Temperature_,
                           Data<T, 3>& Pressure_,
                           Data<T, 2>& SurfaceTemperature_,
                           Data<T, 2>& SurfacePressure_,
                           Data<T, 2>& FirstLevelWindModule_,
                           Data<T, 4>& WetDiameter_aer_,
                           Data<T, 2>& SnowHeight_,
                           Data<T, 3>& DepositionVelocity_aer_)
  {
    int i, j, c, b;
    T seconds(this->current_date.GetNumberOfSeconds());
    int Nbin_dep_=this->Nbin_dep_aer*this->Ncomposition_aer;

    // 1D variables.
    Array<T, 1> LUC_1D(Nland);
    Array<T, 1> DepositionVelocity_1D(Nbin_dep_);
    Array<T, 1> Diameter_aer(Nbin_dep_);

    DepositionVelocity_aer_.SetZero();
    for (j = 0; j < this->Ny; j++)
      for (i = 0; i < this->Nx; i++)
	{
	  for (c = 0; c < Nland; c++)
	    LUC_1D(c) = LUC(c, j, i);

	  for (b = 0; b < this->Nbin_dep_aer; b++)
	    for (int id=0; id < this->Ncomposition_aer; id++)
	    {
	      int RealID=b*this->Ncomposition_aer+id;
	      Diameter_aer(RealID) = WetDiameter_aer_(
			bin_list_dep_aer[b]*this->Ncomposition_aer+id, 0, j, i);
	    }

	  _compute_dry_deposition(&seconds, &Nbin_dep_, &Nland,
				  Diameter_aer.data(), Density_aer.data(),
				  LUC_1D.data(),
				  &(this->GridZ3D(0)),
				  &FirstLevelWindModule_(j, i),
				  &Temperature_(0, j, i),
				  &SurfaceTemperature_(j, i),
				  &Pressure_(0, j, i),
				  &SurfacePressure_(j, i), &SnowHeight_(j, i),
				  LandRoughnessHeight.data(),
				  ZhangGamma.data(), ZhangAlpha.data(),
				  SmallRadius.data(), LargeRadius.data(),
				  DepositionVelocity_1D.data());
				  
	  for (b = 0; b < Nbin_dep_; b++)
	  {
	    DepositionVelocity_aer_(b, j, i) =DepositionVelocity_1D(b);
	  }
	}
  }

  //! Initializes scavenging coefficients.
  /*! Computes below-cloud scavenging coefficients of a reversibly soluble
    gas.
    \param Temperature_ temperature (K).
    \param Pressure_ pressure (Pa).
    \param LiquidWaterContent_ liquid water content (kg/kg).
    \param pH_ pH of cloud droplet.
    \param CloudBaseHeight_ cloud basis height (m).
    \param Rain_ rain intensity (mm/h).
    \param ScavengingCoefficient_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitScavengingCoefficient(Data<T, 3>& Temperature_,
                              Data<T, 3>& Pressure_,
                              Data<T, 3>& LiquidWaterContent_,
                              Data<T, 3>& pH_,
                              Data<T, 2>& CloudBaseHeight_,
                              Data<T, 2>& Rain_,
                              Data<T, 4>& ScavengingCoefficient_)
  {
    int i, j, k, s;

    int Ns_scav = ScavengingCoefficient_.GetLength(0);
    int Nz = ScavengingCoefficient_.GetLength(1);
    int Ny = ScavengingCoefficient_.GetLength(2);
    int Nx = ScavengingCoefficient_.GetLength(3);

    T concH, pH_loc(0.), lwc_loc, factor;
    T H1, H2, H3, H4, H5, H6, Kw;
    T Pr = 8.2e-2;

    bool has_scavenging_SO2 = this->HasScavenging("SO2");
    bool has_scavenging_NH3 = this->HasScavenging("NH3");
    bool has_scavenging_HNO3 = this->HasScavenging("HNO3");
    bool has_scavenging_HONO = this->HasScavenging("HONO");
    bool has_scavenging_HNO2 = this->HasScavenging("HNO2");
    bool has_scavenging_HCL = this->HasScavenging("HCL");

    int index_SO2 = this->ScavengingIndex("SO2");
    int index_NH3 = this->ScavengingIndex("NH3");
    int index_HNO3 = this->ScavengingIndex("HNO3");
    int index_HONO = this->ScavengingIndex("HONO");
    int index_HNO2 = this->ScavengingIndex("HNO2");
    int index_HCL = this->ScavengingIndex("HCL");

    // Henry constants of scavenged species.
    Array<T, 4> henry_constant_(Ns_scav, Nz, Ny, Nx);

    // Heat of dissolution.
    Array<T, 1> dissolution_heat_(Ns_scav);
    // Gas phase diffusivities of scavenged species.
    Array<T, 1> gas_phase_diffusivity_(Ns_scav);

    // Extracts species data for scavenged species only.
    for (s = 0; s < this->Ns_scav; s++)
      {
        dissolution_heat_(s) = dissolution_heat[this->ScavengingName(s)];
        gas_phase_diffusivity_(s) =
          this->gas_phase_diffusivity[this->ScavengingName(s)];
      }

    // Extracts species data for scavenged species only.
    for (s = 0; s < Ns_scav; s++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            henry_constant_(s, k, j, i) =
              this->henry_constant[this->ScavengingName(s)];

    // Modifies Henry constants according to pH.
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        if (Rain_(j, i) > 0.0)
          {
            for (k = 1; k < Nz; k++)
              {
                // Converts from kg / kg to g / m^3.
                lwc_loc = LiquidWaterContent_(k, j, i) * 1000.0
                  * Pressure_(k, j, i) / 101325.0 * 28.97 / Pr
                  / Temperature_(k, j, i);

                // Default cloud pH.
                pH_loc = 4.5;

                if (CloudBaseHeight_(j, i) >= this->GridZ4D(k - 1) &&
                    CloudBaseHeight_(j, i) < this->GridZ4D(k) &&
                    lwc_loc >= lwc_cloud_threshold
                    && pH_(k, j, i) > 0.)
                  // Assumed bounds on H+ concentration implies that pH is
                  // between 4 and 6.
                  pH_loc = min(max(pH_(k, j, i), 4.), 6.);
              }

            // H+ concentration.
            concH = pow(10., -pH_loc);

            for (k = 0; k < Nz; k++)
              if (CloudBaseHeight_(j, i) > this->GridZ4D(k))
                {
                  factor = 1. / 298. - 1. / Temperature_(k, j, i);

                  for (s = 0; s < this->Ns_scav; s++)
                    henry_constant_(s, k, j, i) *=
                      exp(dissolution_heat_(s) * factor / 0.001985);

                  H1 = 0.0130 * exp(-2000. * factor);
                  H2 = 6.3e-8 * exp(-1495. * factor);
                  H3 = 1.7e-5 * exp(4325. * factor);
                  H4 = 15.1 * exp(-8700. * factor);
                  H5 = 5.9e-4 * exp(1760. * factor);
                  H6 = 1.7e6 * exp(-6900. * factor);
                  Kw = 1.e-14 * exp(6700. * factor);

                  if (has_scavenging_SO2)
                    henry_constant_(index_SO2, k, j, i) *=
                      1. + H1 / concH + H1 * H2 / (concH * concH);
                  if (has_scavenging_NH3)
                    henry_constant_(index_NH3, k, j, i) *= 1.
                      + H3 * concH / Kw;
                  if (has_scavenging_HNO3)
                    henry_constant_(index_HNO3, k, j, i) *= 1. + H4 / concH;
                  if (has_scavenging_HONO)
                    henry_constant_(index_HONO, k, j, i) *= 1. + H5 / concH;
                  if (has_scavenging_HNO2)
                    henry_constant_(index_HNO2, k, j, i) *= 1. + H5 / concH;
                  if (has_scavenging_HCL)
                    henry_constant_(index_HCL, k, j, i) *= 1. + H6 / concH;
                }
          }

    _compute_scavenging_coefficient(&Nx, &Ny, &Nz, &Ns_scav,
                                    &Nx, &Ny, &Nz,
                                    this->GridZ4D.GetArray().data(),
                                    gas_phase_diffusivity_.data(),
                                    henry_constant_.data(),
                                    Temperature_.GetData(),
                                    Pressure_.GetData(),
                                    Rain_.GetData(),
                                    CloudBaseHeight_.GetData(),
                                    ScavengingCoefficient_.GetData());
  }


  //! Initializes scavenging coefficients for aerosols.
  /*! Computes below-cloud scavenging coefficients for aerosols.
    \param Temperature_ temperature (K).
    \param Pressure_ pressure (Pa).
    \param WetDiameter_aer_ (output) wet aerosol diameter (m).
    \param CloudBaseHeight_ cloud basis height (m).
    \param Rain_ rain intensity (mm/h).
    \param ScavengingCoefficient_aer_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitScavengingCoefficient_aer(Data<T, 3>& Temperature_,
                                  Data<T, 3>& Pressure_,
                                  Data<T, 4>& WetDiameter_aer_,
                                  Data<T, 2>& CloudBaseHeight_,
                                  Data<T, 2>& Rain_,
                                  Data<T, 4>& ScavengingCoefficient_aer_)
  {
    int i, j, k, b;

    // 1D variables.
    int Nbin_scav_=this->Nbin_scav_aer*this->Ncomposition_aer;
    Array<T, 1> Coefficient(Nbin_scav_);
    Array<T, 1> wet_diameter_loc(Nbin_scav_);

    ScavengingCoefficient_aer_.SetZero();
    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
	for (i = 0; i < this->Nx; i++)
	  {
	    for (b = 0; b < this->Nbin_scav_aer; b++)
	      for(int id=0; id< this->Ncomposition_aer; id++)
		wet_diameter_loc(b*this->Ncomposition_aer+id) =
		  WetDiameter_aer_(bin_list_scav_aer[b]*this->Ncomposition_aer+id, k, j, i);

            _compute_scavenging_coefficient_aer(&Nbin_scav_,
                                                &Temperature_(k, j, i),
                                                &Pressure_(k, j, i),
                                                &Rain_(j, i),
                                                wet_diameter_loc.data(),
                                                Density_aer.data(),
                                                &(this->GridZ3D(k)),
                                                &CloudBaseHeight_(j, i),
                                                Coefficient.data());
            for (b = 0; b < Nbin_scav_; b++)
              ScavengingCoefficient_aer_(b, k, j, i) =
                Coefficient(b);
          }
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SetDate(Date date)
  {
    BaseModel<T>::SetDate(date);

    if (date != this->data_date)
      {
        Polair3DTransport<T, ClassAdvection, ClassDiffusion>::InitAllData();
        Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
          ::InitAllData();
        Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
          ::InitAllData();
      }
  }

  //////////////////////////////////
  // External composition methods //
  //////////////////////////////////

  //! Return external composition id based on species name
  //! (where this specie related group is pure)
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::FindExternalCompositionID(string species)
  {
    int index_species = this->GetSpeciesIndex_aer(species);
    int index_group = this->aerosol_species_group_relation(index_species);
    if(index_group==(this->Ngroup_aer-1))
      return 0;
    else
      for(int i = 1; i < this->Ncomposition_aer; i++)
      {
	if(composition_bounds(i,index_group,1)==1.0)
	  return i;
      }
    cout<<"Error!:can not find related composition ID,";
    cout<<"please check species-group relations"<<endl;
    return 0;
  }

  //! Perform the transformation from internal bc data to external bc data 
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::BoundaryConditionTransformation()
  {
    for (int j = 0; j < this->Nsize_section_aer; j++)
    {
      Data<T, 4> Conc_group_x(GridG4D_aer, this->GridZ4D, this->GridY4D,this->GridX4D_interf_bc);
      Data<T, 4> Conc_group_y(GridG4D_aer, this->GridZ4D, this->GridY4D_interf_bc,this->GridX4D);
      Data<T, 3> Conc_group_z(GridG3D_aer, this->GridY3D, this->GridX3D);
      Conc_group_x.SetZero();
      Conc_group_y.SetZero();
      Conc_group_z.SetZero();
      Data<T, 3> Conc_total_x(this->GridZ3D, this->GridY3D,GridX3D_interf_bc);
      Data<T, 3> Conc_total_y(this->GridZ3D, GridY3D_interf_bc,this->GridX3D);
      Data<T, 2> Conc_total_z(this->GridY2D, this->GridX2D);
      Conc_total_x.SetZero();
      Conc_total_y.SetZero();
      Conc_total_z.SetZero();

      //compute the total bin/group mass
      for (int i = 0; i < Ns_bc_aer; i++)
      {
	string species = species_list_bc_aer[i].first;
	int index_species = this->GetSpeciesIndex_aer(species);
	int index_group = this->aerosol_species_group_relation(index_species);
	for(int y = 0; y < this->Ny; y++)
	  for(int x = 0; x < this->Nx; x++)
	  {
	    Conc_total_z(y,x)+=BoundaryCondition_z_aer_i_tmp(i,j,y,x);
	    Conc_group_z(index_group,y,x)+=BoundaryCondition_z_aer_i_tmp(i,j,y,x);
	  }

	for(int z = 0; z < this->Nz; z++)
	  for(int x = 0; x < this->Nx; x++)
	  {
	    Conc_total_y(z,0,x)+=BoundaryCondition_y_aer_i_tmp(i,j,z,0,x);
	    Conc_group_y(index_group,z,0,x)+=BoundaryCondition_y_aer_i_tmp(i,j,z,0,x);
	    Conc_total_y(z,1,x)+=BoundaryCondition_y_aer_i_tmp(i,j,z,1,x);
	    Conc_group_y(index_group,z,1,x)+=BoundaryCondition_y_aer_i_tmp(i,j,z,1,x);
	  }

	for(int z = 0; z < this->Nz; z++)
	  for(int y = 0; y < this->Ny; y++)
	  {
	    Conc_total_x(z,y,0)+=BoundaryCondition_x_aer_i_tmp(i,j,z,y,0);
	    Conc_group_x(index_group,z,y,0)+=BoundaryCondition_x_aer_i_tmp(i,j,z,y,0);
	    Conc_total_x(z,y,1)+=BoundaryCondition_x_aer_i_tmp(i,j,z,y,1);
	    Conc_group_x(index_group,z,y,1)+=BoundaryCondition_x_aer_i_tmp(i,j,z,y,1);
	  }
      }
      //redistribute based on new composition id
      for(int y = 0; y < this->Ny; y++)
	for(int x = 0; x < this->Nx; x++)
	{
	  int composition_id=FindCompositionID(Conc_total_z,Conc_group_z,composition_bounds,
			this->Ncomposition_aer,this->Ngroup_aer,y,x);
	  int id_bin=j*this->Ncomposition_aer+composition_id;
	  for (int i = 0; i < Ns_bc_aer; i++)
	  {
	    string species = species_list_bc_aer[i].first;
	    this->BoundaryCondition_z_aer_i(i,id_bin,y,x)=BoundaryCondition_z_aer_i_tmp(i,j,y,x);
	  }
	  if (this->option_process["with_number_concentration"])
	  {
	    this->NumberBoundaryCondition_z_aer_i(id_bin,y,x)=NumberBoundaryCondition_z_aer_i_tmp(j,y,x);
	  }
	}
      for(int z = 0; z < this->Nz; z++)
	for(int x = 0; x < this->Nx; x++)
	{
	  int composition_id_0=FindCompositionID(Conc_total_y,Conc_group_y,composition_bounds,
			this->Ncomposition_aer,this->Ngroup_aer,z,0,x);
	  int id_bin_0=j*this->Ncomposition_aer+composition_id_0;
	  int composition_id_1=FindCompositionID(Conc_total_y,Conc_group_y,composition_bounds,
			this->Ncomposition_aer,this->Ngroup_aer,z,1,x);
	  int id_bin_1=j*this->Ncomposition_aer+composition_id_1;
	  for (int i = 0; i < Ns_bc_aer; i++)
	  {
	    string species = species_list_bc_aer[i].first;
	    this->BoundaryCondition_y_aer_i(i,id_bin_0,z,0,x)=BoundaryCondition_y_aer_i_tmp(i,j,z,0,x);
	    this->BoundaryCondition_y_aer_i(i,id_bin_1,z,1,x)=BoundaryCondition_y_aer_i_tmp(i,j,z,1,x);
	  }
	  if (this->option_process["with_number_concentration"])
	  {
	    this->NumberBoundaryCondition_y_aer_i(id_bin_0,z,0,x)=NumberBoundaryCondition_y_aer_i_tmp(j,z,0,x);
	    this->NumberBoundaryCondition_y_aer_i(id_bin_1,z,1,x)=NumberBoundaryCondition_y_aer_i_tmp(j,z,1,x);
	  }
	}
      for(int z = 0; z < this->Nz; z++)
	for(int y = 0; y < this->Ny; y++)
	{
	  int composition_id_0=FindCompositionID(Conc_total_x,Conc_group_x,composition_bounds,
			this->Ncomposition_aer,this->Ngroup_aer,z,y,0);
	  int id_bin_0=j*this->Ncomposition_aer+composition_id_0;
	  int composition_id_1=FindCompositionID(Conc_total_x,Conc_group_x,composition_bounds,
			this->Ncomposition_aer,this->Ngroup_aer,z,y,1);
	  int id_bin_1=j*this->Ncomposition_aer+composition_id_1;
	  for (int i = 0; i < Ns_bc_aer; i++)
	  {
	    string species = species_list_bc_aer[i].first;
	    this->BoundaryCondition_x_aer_i(i,id_bin_0,z,y,0)=BoundaryCondition_x_aer_i_tmp(i,j,z,y,0);
	    this->BoundaryCondition_x_aer_i(i,id_bin_1,z,y,1)=BoundaryCondition_x_aer_i_tmp(i,j,z,y,1);
	  }
	  if (this->option_process["with_number_concentration"])
	  {
	    this->NumberBoundaryCondition_x_aer_i(id_bin_0,z,y,0)=NumberBoundaryCondition_x_aer_i_tmp(j,z,y,0);
	    this->NumberBoundaryCondition_x_aer_i(id_bin_1,z,y,1)=NumberBoundaryCondition_x_aer_i_tmp(j,z,y,1);
	  }
	}
    }
    return 1;
  }

  //! Return composition id based on group mass fractions 
  //! for 2 D data
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::FindCompositionID(Data<T,2>& total_mass_,
			    Data<T,3>& group_mass_,
			    Data<T,3>& composition_bounds_,
			    int Nc ,int Ng, int y, int x)
  {
    T frac_group=0;
    Data<T,1> fgroup(Ng);
    for(int i = 0; i < Nc; i++)
    {
      int counter=0;
      if(Ng>1)
      {
	for(int g = 0; g < Ng; g++)
	{
	  if(total_mass_(y,x)>0.0)
	  {
	    frac_group=group_mass_(g,y,x)/total_mass_(y,x);
	    fgroup(g)=frac_group;
	  }
	  else
	  frac_group=0.0;
	  
	  if(composition_bounds_(i,g,0)>0)
	  {
	    if(frac_group>composition_bounds_(i,g,0)&&frac_group<=composition_bounds_(i,g,1))
	      counter++;
	  }
	  else
	  {
	    if(frac_group>=composition_bounds_(i,g,0)&&frac_group<=composition_bounds_(i,g,1))
	      counter++;
	  }
	}
	if(counter==Ng)
	  return i;
      }
      else
      return 0;
    }
    cout<<"Error missing composition_bounds!"<<endl;
    cout<<"("<<y<<","<<x<<")"<<endl;
    for(int g = 0; g < Ng; g++)
      cout<<"g("<<g<<")="<<fgroup(g)<<endl;
    abort();
    return Nc;    
  }

  //! Return composition id based on group mass fractions
  //! for 3 D data
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassAdvection,
	   class ClassDiffusion, class ClassChemistry>
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::FindCompositionID(Data<T,3>& total_mass_,
			    Data<T,4>& group_mass_,
			    Data<T,3>& composition_bounds_,
			    int Nc ,int Ng, int z, int y, int x)
  {
    T frac_group=0;
    Data<T,1> fgroup(Ng);
    for(int i = 0; i < Nc; i++)
    {
      int counter=0;
      if(Ng>1)
      {
	for(int g = 0; g < Ng; g++)
	{
	  if(total_mass_(z,y,x)>0.0)
	  {
	    frac_group=group_mass_(g, z, y,x)/total_mass_(z,y,x);
	    fgroup(g)=frac_group;
	  }
	  else
	  frac_group=0.0;

	  if(composition_bounds_(i,g,0)>0)
	  {
	    if(frac_group>composition_bounds_(i,g,0)&&frac_group<=composition_bounds_(i,g,1))
	      counter++;
	  }
	  else
	  {
	    if(frac_group>=composition_bounds_(i,g,0)&&frac_group<=composition_bounds_(i,g,1))
	      counter++;
	  }
	}
	if(counter==Ng)
	  return i;
      }
      else
      return 0;
    }
    cout<<"Error missing composition_bounds!"<<endl;
    cout<<" "<<z<<" "<<y<<" "<<x<<endl;
    for(int g = 0; g < Ng; g++)
      cout<<"g("<<g<<")="<<fgroup(g)<<endl;
    abort();
    return Nc;        
  }

  /////////////////
  // INTEGRATION //
  /////////////////


  //! Performs one advection step.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Advection()
  {
    this->Advection_.Forward(*this);
    this->Advection_.Forward_aer(*this);
  }


  //! Performs one diffusion step.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Diffusion()
  {
    this->Diffusion_.Forward(*this);
    this->Diffusion_.Forward_aer(*this);
  }


  //! Performs the time integration of aerosol point emissions.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::PointEmission_aer()
  {
    int e, i, j, k, s, b;
    int Nsource(point_emission_list_aer.size());
    T quantity;
    // Loop over all point sources.
    for (e = 0; e < Nsource; e++)
      {
	i = convert<int>(point_emission_list_aer[e]["i"]);
	j = convert<int>(point_emission_list_aer[e]["j"]);
	k = convert<int>(point_emission_list_aer[e]["k"]);
	string species=point_emission_list_aer[e]["species"];
	s = this->GetSpeciesIndex_aer(species);
	b = convert<int>(point_emission_list_aer[e]["bin"]);
	if (point_emission_list_aer[e]["type"] == "puff")
	  {
	    if(point_emis_format=="Internal"&&this->option_process["with_external_composition"])
	    {
	      int CompositionID=FindExternalCompositionID(species);
	      int RealID=this->Ncomposition_aer*b+CompositionID;
	      Date date(point_emission_list_aer[e]["date"]);
	      if (this->current_date <= date && date < this->next_date)
		this->Concentration_aer(s, RealID, k, j, i) +=
		  convert<T>(point_emission_list_aer[e]["quantity"])
		  / (this->CellWidth_x(i) * this->CellWidth_y(j)
		    * this->CellWidth_z(k));
	    }
	    else
	    {
	      Date date(point_emission_list_aer[e]["date"]);
	      if (this->current_date <= date && date < this->next_date)
		this->Concentration_aer(s, b, k, j, i) +=
		  convert<T>(point_emission_list_aer[e]["quantity"])
		  / (this->CellWidth_x(i) * this->CellWidth_y(j)
		    * this->CellWidth_z(k));	      
	    }
	  }
	else if (point_emission_list_aer[e]["type"] == "continuous")
	  {
	    Date date_beg(point_emission_list_aer[e]["date_beg"]);
	    Date date_end(point_emission_list_aer[e]["date_end"]);
	    if (this->next_date > date_beg && this->current_date < date_end)
	      {
		quantity = this->Delta_t
		  - max(0., date_beg.GetSecondsFrom(this->current_date))
		  - max(0., this->next_date.GetSecondsFrom(date_end));
		quantity *= convert<T>(point_emission_list_aer[e]["rate"]);
		if(point_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  int CompositionID=FindExternalCompositionID(species);
		  int RealID=this->Ncomposition_aer*b+CompositionID;
		  this->Concentration_aer(s, RealID, k, j, i) += quantity
		    / (this->CellWidth_x(i) * this->CellWidth_y(j)
		      * this->CellWidth_z(k));
		}
		else
		{
		  this->Concentration_aer(s, b, k, j, i) += quantity
		    / (this->CellWidth_x(i) * this->CellWidth_y(j)
		      * this->CellWidth_z(k));
		}
	      }
	  }
      }
  }


  //! Performs one chemistry step.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Chemistry()
  {
    this->Chemistry_.Forward(*this);

    if (this->option_process["with_volume_emission_aer"])
	{
	int b, k, j, i, emis_s, emis_b, emis_b_number;
	for (int s = 0; s < this->Ns_aer; s++)
	  {
	    emis_b=0;
	    for (b = 0; b < this->Nbin_aer; b++)
	      if (this->HasVolumeEmission_aer(s,b))
		{
		  emis_s = VolumeEmissionIndex_aer(s);
		  for (k = 0; k < Nz_vol_emis_aer; k++)
		    for (j = 0; j < this->Ny; j++)
		      for (i = 0; i < this->Nx; i++)
			this->Concentration_aer(s, b, k, j, i) +=
			  this->Delta_t *
			  VolumeEmission_aer_i(emis_s, emis_b, k, j, i);
		  emis_b++;
		}
	  }
		
	if (this->option_process["with_number_concentration"])
	  {
	    for (b = 0; b < this->Nbin_aer; b++)
	      if (this->HasNumberVolumeEmission_aer(b))
		{
				 
		  emis_b_number = this->NumberVolumeEmissionIndex_aer(b);
		  for (k = 0; k < Nz_vol_emis_aer; k++)
		    for (j = 0; j < this->Ny; j++)
		      for (i = 0; i < this->Nx; i++)
			this->NumberConcentration_aer(b, k, j, i) +=
			  this->Delta_t *
			  NumberVolumeEmission_aer_i(emis_b_number, k, j, i);
		}
	  }
	}

    this->Chemistry_.Forward_aer(*this);
  }


  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Radiatif(Data<T, 4>& PhotolysisRate, Date date)
  {
    Data<T, 4> ConcentrationWater(this->GridB4D_aer, this->GridZ4D,
                                  this->GridY4D, this->GridX4D);

    Data<T, 5> Concentration_aer_WithoutWater(GridS5D_aer_noH2O, GridB5D_aer,
                                              GridZ5D, GridY5D, GridX5D);

    for (int b = 0; b < this->Nbin_aer; b++)
      for (int k = 0; k < this->Nz; k++)
        for (int j = 0; j < this->Ny; j++)
          for (int i = 0; i < this->Nx; i++)
            {
              ConcentrationWater(b, k, j, i) = this->Concentration_aer
                (this->GetSpeciesIndex_aer("PH2O"), b, k, j, i);
              for (int w = 0; w < (this->Ns_aer - 1); w++)
                // Water is the last specie stored in Concentration_aer
                Concentration_aer_WithoutWater(w, b, k, j, i)
                  = this->Concentration_aer(w, b, k, j, i);
            }


    Data<T, 1> MeanDiameter(this->Nbin_aer);
    for (int b = 0; b < this->Nbin_aer; b++)
      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));

    WetDiameter_aer.Mlt(1e6);
    fixed_density_aer = fixed_density_aer * 1.e3;

    // Calcul of Aerosol Optical Properties
    cout << "Calcul of Aerosol Optical Properties ..." ;
    cout.flush();

    compute_optical_properties(OpticalDepthAerosol,
                               SingleScatteringAlbedo,
                               MeanExtinctionEfficiencyFactor,
                               MeanAbsorbtionEfficiencyFactor,
                               PhaseFunction,
                               this->GridZ3D_interf,
                               GridIndexReal,
                               GridIndexImaginary,
                               GridIndexDiameter,
                               AbsorptionEfficiencyFactorTable,
                               ExtinctionEfficiencyFactorTable,
                               PhaseFunctionTable,
                               Concentration_aer_WithoutWater,
                               ConcentrationWater,
                               PureSpeciesIndexReal,
                               PureSpeciesIndexImaginary,
                               WaterIndexReal,
                               WaterIndexImaginary,
                               MeanDiameter,
                               WetDiameter_aer,
                               fixed_density_aer,
                               black_carbon_index,
                               option_black_carbon_treatment,
                               option_wet_index,
                               option_well_mixed_index);

    WetDiameter_aer.Mlt(1e-6);
    fixed_density_aer = fixed_density_aer / 1.e3;

    int Year = date.GetYear();
    int Day = date.GetOrdinalDay();

    float FloatHour = (date.GetHour() * 3600 + date.GetMinutes() * 60) / 3600.;

    cout << "done" << endl;

    // Pour le moment on utilise SurfacePressure_i et IceOpticalDepth_i,
    // CloudOpticalDepth_i. A terme, voir si il est judicieux de les lire
    // au bon pas de temps.

    float x_min_f = this->x_min;
    float Delta_x_f = this->Delta_x;
    float y_min_f = this->y_min;
    float Delta_y_f = this->Delta_y;
    RegularGrid<float> GridZ3D_interf_fl(this->Nz + 1);
    GridZ3D_interf_fl.GetArray() = cast<float>(this->GridZ3D_interf.GetArray());


    Data<float, 3, T> CloudOpticalDepth_i_fl
      (this->GridZ3D, this->GridY3D, this->GridX3D);
    Data<float, 3, T> IceOpticalDepth_i_fl
      (this->GridZ3D, this->GridY3D, this->GridX3D);
    Data<float, 2, T> SurfacePressure_i_fl
      (this->GridY2D, this->GridX2D);
    Data<float, 2, T> SurfacePressure_f_fl
      (this->GridY2D, this->GridX2D);
    Data<float, 3, T> OpticalDepthAerosol600_fl
      (this->GridZ3D, this->GridY3D, this->GridX3D);
    Data<float, 4, T> SingleScatteringAlbedo_fl
      (GridWavelength, this->GridZ4D, this->GridY4D, this->GridX4D);
    Data<float, 4, T> MeanExtinctionEfficiencyFactor_fl
      (GridWavelength, this->GridZ4D, this->GridY4D, this->GridX4D);
    Data<float, 5, T> PhaseFunction_fl
      (GridWavelength, this->GridZ4D, this->GridY4D,
       this->GridX4D, GridLegendre);
    Data<float, 4, T> PhotolysisRate_fl
      (this->GridR_photolysis, this->GridZ4D,
       this->GridY4D, this->GridX4D);

    // Blitz (at least 0.9) does not support casting views,
    // so a manual conversion is done in this case.
    for (int z = 0; z < this->Nz; z++)
      for (int j = 0; j < this->Ny; j++)
        for (int i = 0; i < this->Nx; i++)
          OpticalDepthAerosol600_fl(z, j, i) = float(OpticalDepthAerosol(2, z, j, i));

    CloudOpticalDepth_i_fl.GetArray() = cast<float>(CloudOpticalDepth_i.GetArray());
    IceOpticalDepth_i_fl.GetArray() = cast<float>(IceOpticalDepth_i.GetArray());
    SingleScatteringAlbedo_fl.GetArray() = cast<float>(SingleScatteringAlbedo.GetArray());
    MeanExtinctionEfficiencyFactor_fl.GetArray() = cast<float>(MeanExtinctionEfficiencyFactor.GetArray());
    PhaseFunction_fl.GetArray() = cast<float>(PhaseFunction.GetArray());
    SurfacePressure_i_fl.GetArray() = cast<float>(SurfacePressure_i.GetArray());
    SurfacePressure_f_fl.GetArray() = cast<float>(SurfacePressure_f.GetArray());

    const char* DirectoryParameter = fastJ_parameter_files.c_str();
    int DirectoryParameter_len = strlen(DirectoryParameter);

    Array<char, 2> photolysis_specie_name(this->Nr_photolysis, 10);
    photolysis_specie_name = ' ';
    for (int i = 0; i < this->Nr_photolysis; ++i)
      {
        const string& specie_name = this->GetPhotolysisReactionList()[i];
        for (int j = 0; j < (int) specie_name.size(); ++j)
          photolysis_specie_name(i, j) = specie_name[j];
      }
    cout << "Photolysis rates computation..." ;
    cout.flush();

    _fastjx(&this->Nx, &this->Ny, &this->Nz,
            &x_min_f, &Delta_x_f,
            &y_min_f, &Delta_y_f,
            GridZ3D_interf_fl.GetArray().data(),
            CloudOpticalDepth_i_fl.GetData(),
            IceOpticalDepth_i_fl.GetData(),
            SurfacePressure_f_fl.GetData(),
            OpticalDepthAerosol600_fl.GetData(),
            SingleScatteringAlbedo_fl.GetData(),
            MeanExtinctionEfficiencyFactor_fl.GetData(),
            PhaseFunction_fl.GetData(),
            PhotolysisRate_fl.GetData(),
            photolysis_specie_name.dataFirst(),
            &this->Nr_photolysis, &Year, &Day, &FloatHour,
            DirectoryParameter, &DirectoryParameter_len);

    cout << "done" << endl;

    PhotolysisRate.GetArray() = cast<T>(PhotolysisRate_fl.GetArray());
  }


  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::Forward()
  {
    int i, j, k;

    /*** Air density ***/

    if (!this->data_gone_through_forward
        && this->option_process["with_air_density"])
      {
        this->InterpolateInterface_z(this->CellCenterDistance_z,
                                     this->CellWidth_z,
                                     this->AirDensity_f,
                                     this->AirDensity_interf_z_f);
        this->InterpolateInterface_y(this->CellCenterDistance_y,
                                     this->CellWidth_y,
                                     this->AirDensity_f,
                                     this->AirDensity_interf_y_f);
        this->InterpolateInterface_x(this->CellCenterDistance_x,
                                     this->CellWidth_x,
                                     this->AirDensity_f,
                                     this->AirDensity_interf_x_f);
      }

    /*** Wind ***/

    if (!this->data_gone_through_forward
        && this->option_manage["horizontal_wind"])
      if (!this->option_cartesian)
        {
          this->TransformZonalWind(this->ZonalWind_i);
          this->TransformMeridionalWind(this->MeridionalWind_i);
        }

    if (!this->data_gone_through_forward
        && this->option_manage["vertical_wind"])
      if (this->option_process["with_air_density"])
        this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
                                  this->CellWidth_z,
                                  this->AirDensity_interf_x_i,
                                  this->ZonalWind_i,
                                  this->AirDensity_interf_y_i,
                                  this->MeridionalWind_i,
                                  this->AirDensity_interf_z_i,
                                  this->VerticalWind_i);
      else
        this->ComputeVerticalWind(this->CellWidth_x, this->CellWidth_y,
                                  this->CellWidth_z, this->ZonalWind_i,
                                  this->MeridionalWind_i,
                                  this->VerticalWind_i);

    /*** Diffusion coefficients ***/

    if (!this->data_gone_through_forward
        && this->option_manage["vertical_diffusion"])
      {
        // Computes rho * Kz.
        if (this->option_process["with_air_density"])
          this->VerticalDiffusionCoefficient_f.GetArray() =
            this->AirDensity_interf_z_f.GetArray()
            * this->VerticalDiffusionCoefficient_f.GetArray();
      }

    if (!this->data_gone_through_forward
        && this->option_manage["horizontal_diffusion"]
        && this->option_isotropic_diffusion)
      {
        LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->ZonalDiffusionCoefficient_i);
        this->ZonalDiffusionCoefficient_i.ThresholdMin(0.);
        LinearInterpolationRegular(this->VerticalDiffusionCoefficient_i,
                                   this->MeridionalDiffusionCoefficient_i);
        this->MeridionalDiffusionCoefficient_i.ThresholdMin(0.);

        if (!this->option_cartesian)
          {
            this->TransformZonalDiffusion(this->GridY3D_interf.GetArray(),
                                          this->ZonalDiffusionCoefficient_i);
            this->TransformMeridionalDiffusion
              (this->MeridionalDiffusionCoefficient_i);
          }
      }

    /*** First level wind module ***/

    if (!this->data_gone_through_forward)
      FirstLevelWindModule_i.GetArray() = FirstLevelWindModule_f.GetArray();
    T zonal_wind, meridional_wind;
    if (!this->data_gone_through_forward
        && this->option_manage["first_level_wind_module"])
      for (j = 0; j < this->Ny; j++)
        for (i = 0; i < this->Nx; i++)
          {
            zonal_wind = 0.5 * (this->ZonalWind_i(0, j, i)
                                + this->ZonalWind_i(0, j, i + 1));
            meridional_wind = 0.5 * (this->MeridionalWind_i(0, j, i)
                                     + this->MeridionalWind_i(0, j + 1, i));
            FirstLevelWindModule_f(j, i) =
              sqrt(zonal_wind * zonal_wind
                   + meridional_wind * meridional_wind);
          }

    /*** Relative humidity ***/

    if (!this->data_gone_through_forward)
      RelativeHumidity_i.GetArray() = RelativeHumidity_f.GetArray();
    if (!this->data_gone_through_forward)
      for (k = 0; k < this->Nz; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              _compute_relative_humidity(&(this->SpecificHumidity_f(k, j, i)),
                                         &(this->Temperature_f(k, j, i)),
                                         &(this->Pressure_f(k, j, i)),
                                         &RelativeHumidity_f(k, j, i));
              RelativeHumidity_f(k, j, i) =
                min(max(RelativeHumidity_f(k, j, i), 0.05), 0.95);
            }

    /*** Wet aerosol diameter ***/

    if (!this->data_gone_through_forward)
      this->InitWetDiameter_aer(RelativeHumidity_i, this->Temperature_i,
                                this->Concentration_aer,this->NumberConcentration_aer,
                                WetDiameter_aer, true);

    /*** Deposition velocities ***/

    if (!this->data_gone_through_forward
        && this->option_manage["deposition_velocity_aer"])
      if (!this->option_process["compute_deposition_aer"])
        {
          for (i = 0; i < Nbin_dep_aer; i++)
            {
              Date date=this->input_files["deposition_velocity_aer"].GetDateMin();
              T Delta_t = this->input_files["deposition_velocity_aer"].GetDelta_t();
              string filename = this->input_files["deposition_velocity_aer"](to_str(bin_list_dep_aer[i]));
              TinyVector<int, 2> new_shape;
              for (int j = 0; j < 2; j++)
                new_shape(j) = DepositionVelocity_aer_f.GetArray().shape()(j + 1);
	  
              Data<T, 2> FileData_extract_i(new_shape);
              Data<T, 2> FileData_extract_f(new_shape);
              Data<T, 2> CurrentData_extract_i(new_shape);
              Data<T, 2> CurrentData_extract_f(new_shape);

              this->UpdateData(filename, date, Delta_t,
                               FileData_extract_i,
                               FileData_extract_f,
                               CurrentData_extract_i,
                               CurrentData_extract_f);
              
              for(int y=0; y<this->Ny; y++)
                for(int x=0; x<this->Nx; x++)
                  for(int id=0; id< this->Ncomposition_aer; id++)
                    {
                      int RealID=this->Ncomposition_aer*i+id;
                      FileDepositionVelocity_aer_i(RealID,y,x)=FileData_extract_i(y,x);
                      FileDepositionVelocity_aer_f(RealID,y,x)=FileData_extract_f(y,x);
                      DepositionVelocity_aer_i(RealID,y,x)=CurrentData_extract_i(y,x);
                      DepositionVelocity_aer_f(RealID,y,x)=CurrentData_extract_f(y,x);
                    }
            }
        }
      else
        {
          DepositionVelocity_aer_i.GetArray()
            = DepositionVelocity_aer_f.GetArray();
          this->InitDepositionVelocity(this->Temperature_f, this->Pressure_f,
                                       this->SurfaceTemperature_f,
                                       this->SurfacePressure_f,
                                       FirstLevelWindModule_f,
                                       WetDiameter_aer, SnowHeight_f,
                                       DepositionVelocity_aer_f);
        }

    /*** Scavenging coefficients ***/

    if (!this->data_gone_through_forward
        && this->option_manage["scavenging_below_cloud_coefficient"]
        && this->scavenging_below_cloud_model == "microphysical-ph")
      InitScavengingCoefficient(this->Temperature_f, this->Pressure_f,
                                LiquidWaterContent_i, pH,
                                this->CloudBaseHeight_f, this->Rain_f,
                                this->ScavengingBelowCloudCoefficient_f);

    if (!this->data_gone_through_forward
        && this->option_manage["scavenging_coefficient_aer"]
        && this->option_process["with_scavenging_aer"])
      {
        ScavengingCoefficient_aer_i.GetArray()
          = ScavengingCoefficient_aer_f.GetArray();
        InitScavengingCoefficient_aer(this->Temperature_f,
                                      this->Pressure_f,
                                      WetDiameter_aer,
                                      this->CloudBaseHeight_f,
                                      this->Rain_f,
                                      ScavengingCoefficient_aer_f);
      }

    this->data_gone_through_forward = true;

    /*** Time integration ***/

    if (this->source_splitting && this->option_process["with_chemistry"])
      {
        this->Source_f.GetArray() = this->Concentration.GetArray();
        this->Source_aer_f.GetArray() = this->Concentration_aer.GetArray();
      }

    if (this->option_process["with_advection"])
      this->Advection();
    //! Dry fluxes are collected before dry deposition is performed
    //! in diffusion.
    if (this->option_process["collect_dry_flux"])
      for (int s = 0; s < this->Ns_dep; s++)
        for (int j = 0; j < this->Ny; j++)
          for (int i = 0; i < this->Nx; i++)
            this->DryDepositionFlux(s, j, i) = 0.5
              * (this->DepositionVelocity_i(s, j, i)
                 + this->DepositionVelocity_f(s, j, i))
              * this->Concentration(this->DepositionVelocityGlobalIndex(s),
                                    0, j, i);
    if (this->option_process["collect_dry_flux_aer"])
      for (int b = 0; b < Nbin_dep_aer; b++)
	for (int id=0; id < this->Ncomposition_aer; id++)
	for (int j = 0; j < this->Ny; j++)
	  for (int i = 0; i < this->Nx; i++)
	    {
	      int Nc= this->Ncomposition_aer;
	      for (int s = 0; s < this->Ns_aer; s++)
		DryDepositionFlux_aer(s, b*Nc+id, j, i) = 0.5
		  * (DepositionVelocity_aer_i(b*Nc+id, j, i)
		     + DepositionVelocity_aer_f(b*Nc+id, j, i))
		  * this->Concentration_aer(s, bin_list_dep_aer[b]*Nc+id, 0, j, i);
	      if(this->option_process["with_number_concentration"])
		DryDepositionFluxNumber_aer(b*Nc+id, j, i) = 0.5
		  * (DepositionVelocity_aer_i(b*Nc+id, j, i)
		     + DepositionVelocity_aer_f(b*Nc+id, j, i))
		  * this->NumberConcentration_aer(bin_list_dep_aer[b]*Nc+id, 0, j, i);
	      
	    }
		
    if (this->option_process["with_diffusion"]) 
      this->Diffusion();

    if (this->option_process["with_point_emission"])
      this->PointEmission();
    if (this->option_process["with_point_emission_aer"])
      this->PointEmission_aer();
    if (this->option_process["with_volume_emission"]
        && !this->option_process["with_chemistry"])
      for (int s = 0; s < this->Ns; s++)
        {
          int k, j, i, emis_s;
          if (this->HasVolumeEmission(s))
            {
              emis_s = this->VolumeEmissionIndex(s);
              for (k = 0; k < this->Nz_vol_emis; k++)
                for (j = 0; j < this->Ny; j++)
                  for (i = 0; i < this->Nx; i++)
                    this->Concentration(s, k, j, i) +=
                      this->Delta_t * this->VolumeEmission_i(emis_s, k, j, i);
            }
        }
    if (this->option_process["with_volume_emission_aer"]
        && !this->option_process["with_chemistry"])
	{
          int b, k, j, i, emis_s, emis_b, emis_b_number;
      for (int s = 0; s < this->Ns_aer; s++)
        {
          emis_b = 0;
          for (b = 0; b < this->Nbin_aer; b++)
            if (this->HasVolumeEmission_aer(s, b))
              {
                emis_s = VolumeEmissionIndex_aer(s);
                for (k = 0; k < Nz_vol_emis_aer; k++)
                  for (j = 0; j < this->Ny; j++)
                    for (i = 0; i < this->Nx; i++)
                      this->Concentration_aer(s, b, k, j, i) +=
                        this->Delta_t *
                        VolumeEmission_aer_i(emis_s, emis_b, k, j, i);
                emis_b++;
              }
        }

		if (this->option_process["with_number_concentration"])
		  {
	    emis_b_number=0;
	    for (b = 0; b < this->Nbin_aer; b++)
	      if (this->HasNumberVolumeEmission_aer(b))
			{
		  for (k = 0; k < Nz_vol_emis_aer; k++)
		    for (j = 0; j < this->Ny; j++)
		      for (i = 0; i < this->Nx; i++)
			this->NumberConcentration_aer(b, k, j, i) +=
			  this->Delta_t *
			  NumberVolumeEmission_aer_i(emis_b_number, k, j, i); 
		  emis_b_number++;
			}
		  }
      }	    


    if (this->option_process["collect_wet_flux"])
      {
        if (this->scavenging_below_cloud_model != "none")
          this->WetDepositionFlux.SetZero();
        if (this->scavenging_in_cloud_model != "none")
          this->InCloudWetDepositionFlux.SetZero();
      }
    if (this->scavenging_below_cloud_model != "none")
      for (int s = 0; s < this->Ns_scav; s++)
        {
          T scavenging_ratio;
          T cloud_height_mean;
          // Depth of the layer that is below the cloud.
          T depth_below;
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  cloud_height_mean = 0.5
                    * (this->CloudBaseHeight_i(j, i)
                       + this->CloudBaseHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (cloud_height_mean > this->GridZ4D_interf(k))
                      {
                        depth_below = min(cloud_height_mean,
                                          this->GridZ4D_interf(k + 1))
                          - this->GridZ4D_interf(k);
                        scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (this->ScavengingBelowCloudCoefficient_i
                                 (s, k, j, i)
                                 + this->ScavengingBelowCloudCoefficient_f
                                 (s, k, j, i)));
                        if (this->option_process["collect_wet_flux"])
                          this->WetDepositionFlux(s, j, i) +=
                            this->Concentration
                            (this->ScavengingGlobalIndex(s), k, j, i)
                            / this->Delta_t * (1. - scavenging_ratio)
                            * depth_below;
                        this->Concentration
                          (this->ScavengingGlobalIndex(s), k, j, i)
                          *= 1. - (1. - scavenging_ratio)
                          * depth_below / (this->GridZ4D_interf(k + 1)
                                           - this->GridZ4D_interf(k));
                      }
                }
        }

    if (this->option_process["with_transport_in_cloud_scavenging"])
      for (int s = 0; s < this->Ns_scav; s++)
        {
          T scavenging_ratio;
          T cloud_base_height_mean, cloud_top_height_mean;
          // Depth of the layer that is inside the cloud.
          T depth_inside;
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              if (this->Rain_i(j, i) > 0.)
                {
                  cloud_base_height_mean = 0.5
                    * (this->CloudBaseHeight_i(j, i)
                       + this->CloudBaseHeight_f(j, i));
                  cloud_top_height_mean = 0.5
                    * (this->CloudTopHeight_i(j, i)
                       + this->CloudTopHeight_f(j, i));
                  for (int k = 0; k < this->Nz; k++)
                    if (cloud_base_height_mean < this->GridZ4D_interf(k + 1)
                        && cloud_top_height_mean > this->GridZ4D_interf(k))
                      {
                        depth_inside = min(cloud_top_height_mean,
                                           this->GridZ4D_interf(k + 1))
                          - max(cloud_base_height_mean,
                                this->GridZ4D_interf(k));
                        scavenging_ratio =
                          exp(- 0.5 * this->Delta_t
                              * (this->ScavengingInCloudCoefficient_i
                                 (s, k, j, i)
                                 + this->ScavengingInCloudCoefficient_f
                                 (s, k, j, i)));
                        if (this->option_process["collect_wet_flux"])
                          this->InCloudWetDepositionFlux(s, j, i) +=
                            this->Concentration
                            (this->ScavengingGlobalIndex(s), k, j, i)
                            / this->Delta_t * (1. - scavenging_ratio)
                            * depth_inside;
                        this->Concentration(this->ScavengingGlobalIndex(s),
                                            k, j, i)
                          *= 1. - (1. - scavenging_ratio)
                          * depth_inside / (this->GridZ4D_interf(k + 1)
                                            - this->GridZ4D_interf(k));
                      }
                }
        }

    if (this->option_process["collect_wet_flux_aer"])
	{
 	   this->WetDepositionFlux_aer.SetZero();
		if (this->option_process["with_number_concentration"])
		  this->WetDepositionFluxNumber_aer.SetZero();
      }


    if (this->option_process["with_scavenging_aer"])
      {
	int Nc=this->Ncomposition_aer;
	for (int b = 0; b < Nbin_scav_aer; b++)
	  {
	    T scavenging_ratio;
	    for(int id=0; id < Nc; id++)
	    for (int j = 0; j < this->Ny; j++)
	      for (int i = 0; i < this->Nx; i++)
		if (this->Rain_i(j, i) > 0.) 
		  {
		    for (int s = 0; s < this->Ns_aer; s++)
		      {
			for (int k = 0; k < this->Nz; k++)
			  {
			    scavenging_ratio =
			      exp(- 0.5 * this->Delta_t
				  * (ScavengingCoefficient_aer_i(b * Nc + id, k, j, i)
				     + ScavengingCoefficient_aer_f(b * Nc + id, k, j, i)));
			    if (this->option_process["collect_wet_flux_aer"])
			      WetDepositionFlux_aer(s * Nc + id, b, j, i) +=
				this->Concentration_aer(s, bin_list_scav_aer[b] * Nc + id,
							k, j, i)
				/ this->Delta_t * (1.- scavenging_ratio)
				* (this->GridZ4D_interf(k+1)
				   - this->GridZ4D_interf(k));
			    this->Concentration_aer(s, bin_list_scav_aer[b] * Nc + id,
						    k, j, i) *= scavenging_ratio;
			  }
			if (this->option_process["collect_wet_flux_aer"])
			  WetDepositionFlux_aer(s, b * Nc +id, j, i) +=
			    InCloudWetDepositionFlux_aer(s, b * Nc + id, j, i);
		      }
		  
		    if  (this->option_process["with_number_concentration"])
		      {
			for (int k = 0; k < this->Nz; k++)
			  {
			    scavenging_ratio =
			      exp(- 0.5 * this->Delta_t
				  * (ScavengingCoefficient_aer_i(b*Nc+id, k, j, i)
				     + ScavengingCoefficient_aer_f(b*Nc+id, k, j, i)));
			    if (this->option_process["collect_wet_flux_aer"])
			      WetDepositionFluxNumber_aer(b*Nc+id, j, i) +=
				this->NumberConcentration_aer(bin_list_scav_aer[b*Nc+id],
							      k, j, i)
				/ this->Delta_t * (1.- scavenging_ratio)
				* (this->GridZ4D_interf(k+1) - this->GridZ4D_interf(k));
							
						  
			    this->NumberConcentration_aer(bin_list_scav_aer[b]*Nc+id,
							  k, j, i) *= scavenging_ratio;
			  }
			if (this->option_process["collect_wet_flux_aer"])
			  WetDepositionFluxNumber_aer(b*Nc+id, j, i) +=
			    InCloudWetDepositionFluxNumber_aer(b*Nc+id, j, i);
		      }
					
		  }
	  }
      }

    if (this->source_splitting && this->option_process["with_chemistry"])
      {
	this->Source_i.GetArray() =
	  (this->Concentration.GetArray() - this->Source_f.GetArray())
	  / this->Delta_t;
	this->Concentration.GetArray() = this->Source_f.GetArray();
	this->Source_f.GetArray() = this->Source_i.GetArray();
		
	if (this->option_process["with_volume_emission"])
	  for (int s = 0; s < this->Ns_vol_emis; s++)
	    {
	      int gs = this->VolumeEmissionGlobalIndex(s);
	      for (k = 0; k < this->Nz_vol_emis; k++)
		for (j = 0; j < this->Ny; j++)
		  for (i = 0; i < this->Nx; i++)
		    {
		      this->Source_i(gs, k, j, i)
			+= this->VolumeEmission_i(s, k, j, i);
		      this->Source_f(gs, k, j, i)
			+= this->VolumeEmission_f(s, k, j, i);
		    }
	    }
		
	this->Source_aer_i.GetArray() =
	  (this->Concentration_aer.GetArray() - this->Source_aer_f.GetArray())
	  / this->Delta_t;
	this->Concentration_aer.GetArray() = this->Source_aer_f.GetArray();
	this->Source_aer_f.GetArray() = this->Source_aer_i.GetArray();
	
	if (this->option_process["with_volume_emission_aer"])
	  for (int s = 0; s < Ns_vol_emis_aer; s++)
	    {
	      int gs = this->VolumeEmissionGlobalIndex_aer(s);
	      for (int b = 0; b < Nb_vol_emis_aer; b++)
		{
		  int Nc=this->Ncomposition_aer;
		  int gb = species_list_vol_emis_aer[s].second[b];
		  for(int id=0; id < Nc; id++)
		  for (k = 0; k < Nz_vol_emis_aer; k++)
		    for (j = 0; j < this->Ny; j++)
		      for (i = 0; i < this->Nx; i++)
			{
			  this->Source_aer_i(gs, gb*Nc+id, k, j, i)
			    += this->VolumeEmission_aer_i(s, b*Nc+id, k, j, i);
			  this->Source_aer_f(gs, gb*Nc+id, k, j, i)
			    += this->VolumeEmission_aer_f(s, b*Nc+id, k, j, i);
			}
		}
	    }
      }

    if (this->option_process["with_chemistry"])
      Chemistry();

    this->AddTime(this->Delta_t);
    this->step++;
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////

  //! Returns the concentration array of aerosol species at a given point.
  /*!
    \return The concentration at the point.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  T Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetConcentration_aer(int s, int b, T z, T y, T x)
  {
    T concentration;
    Data<T, 3>
      Concentration_tmp(&this->Concentration_aer(s, b, 0, 0, 0),
                        shape(this->Nz, this->Ny, this->Nx));
    Concentration_tmp.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    Array<T, 1> Coord(3);
    Coord(0) = z;
    Coord(1) = y;
    Coord(2) = x;
    LinearInterpolationPoint(Concentration_tmp, Coord, concentration);
    return concentration;
  }



  //! Returns the number of aerosol species with volume sources.
  /*! If source splitting is used, the number of species with volume sources
    is the total number of species; otherwise, it is the number of volume
    emissions.
    \return The number of species with volume sources.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNs_source_aer()
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return this->Ns_aer;
    else
      return Ns_vol_emis_aer;
  }


  //! Returns the number of levels of aerosol volume sources.
  /*! If source splitting is used, the number of levels of volume sources
    is the total number of vertical layers; otherwise, it is the number of
    levels in volume emissions.
    \return The number of levels of volume sources.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNz_source_aer()
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return this->Nz;
    else
      return Nz_vol_emis_aer;
  }


  //! Returns the number of bins in aerosol volume sources.
  /*! \return The maximum number of bins of aerosol volume sources.
   */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetNbinMax_source_aer()
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return this->Nbin_aer;
    else
      return Nb_vol_emis_aer;
  }


  //! Checks whether an aerosol has sources.
  /*! If source splitting is used, any aerosol has sources; otherwise, it
    checks whether the aerosol has volume emissions.
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has sources, false otherwise.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  bool Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::HasSource_aer(int s, int b)
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return true;
    else
      return HasVolumeEmission_aer(s, b);
  }


  //! Returns the global index of an aerosol with volume sources.
  /*!
    \param s species index in volume sources.
    \return The species global index.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SourceGlobalIndex_aer(int s)
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return s;
    else
      return VolumeEmissionGlobalIndex_aer(s);
  }


  /*! \brief Returns the global index of an aerosol bin with volume sources
    for a given species. */
  /*!
    \param s species index in volume sources.
    \param b bin index in volume sources.
    \return The bin global index for the given species.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::SourceGlobalBinIndex_aer(int s, int b)
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return b;
    else
      return this->species_list_vol_emis_aer[s].second[b];
  }


  //! Returns the aerosol sources at current date.
  /*!
    \return The aerosol sources at current date.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Data<T, 5>&
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetSource_aer_i()
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return Source_aer_i;
    else
      return VolumeEmission_aer_i;
  }


  //! Returns the aerosol sources at next date.
  /*!
    \return The aerosol sources at next date.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Data<T, 5>&
  Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::GetSource_aer_f()
  {
    if (this->source_splitting && this->option_process["with_chemistry"])
      return Source_aer_f;
    else
      return VolumeEmission_aer_f;
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Moves model input-data to the current date.
  /*! This method prepares the model for a time integration from the current
    date. It reads input data to related to chemistry be read before InitStep
    and Forward.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DAerosol<T, ClassAdvection, ClassDiffusion, ClassChemistry>
  ::InitAllData()
  {
    int k, i, j;
    string species, species_bin;
    vector<int> isize_section;
    string filename;
    Date date;
    T Delta_t;

    /*** Surface temperature ***/

    if (this->option_manage["surface_temperature"])
      this->InitData("meteo", "SurfaceTemperature",
                     FileSurfaceTemperature_i, FileSurfaceTemperature_f,
                     this->current_date, SurfaceTemperature_f);

    /*** Surface pressure ***/

    if (this->option_manage["surface_pressure"])
      this->InitData("meteo", "SurfacePressure",
                     FileSurfacePressure_i, FileSurfacePressure_f,
                     this->current_date, SurfacePressure_f);

    /*** First level wind module ***/

    T zonal_wind, meridional_wind;
    if (this->option_manage["first_level_wind_module"])
      for (j = 0; j < this->Ny; j++)
        for (i = 0; i < this->Nx; i++)
          {
            zonal_wind = 0.5 * (this->ZonalWind_i(0, j, i)
                                + this->ZonalWind_i(0, j, i + 1));
            meridional_wind = 0.5 * (this->MeridionalWind_i(0, j, i)
                                     + this->MeridionalWind_i(0, j + 1, i));
            FirstLevelWindModule_f(j, i) =
              sqrt(zonal_wind * zonal_wind
                   + meridional_wind * meridional_wind);
          }

    /*** Liquid water content ***/

    if (this->option_manage["liquid_water_content"])
      this->InitData("meteo", "LiquidWaterContent",
                     FileLiquidWaterContent_i, FileLiquidWaterContent_f,
                     this->current_date, LiquidWaterContent_i);

    if (this->computed_photolysis == "online")
      {
        this->InitData("meteo", "CloudOpticalDepth",
                       FileCloudOpticalDepth_i, FileCloudOpticalDepth_f,
                       this->current_date, CloudOpticalDepth_i);
        this->InitData("meteo", "IceOpticalDepth",
                       FileIceOpticalDepth_i, FileIceOpticalDepth_f,
                       this->current_date, IceOpticalDepth_i);
      }

    /*** Relative humidity ***/

    for (k = 0; k < this->Nz; k++)
      for (j = 0; j < this->Ny; j++)
        for (i = 0; i < this->Nx; i++)
          {
            _compute_relative_humidity(&(this->SpecificHumidity_f(k, j, i)),
                                       &(this->Temperature_f(k, j, i)),
                                       &(this->Pressure_f(k, j, i)),
                                       &RelativeHumidity_f(k, j, i));
            RelativeHumidity_f(k, j, i) =
              min(max(RelativeHumidity_f(k, j, i), 0.05), 0.95);
          }

    /*** Aerosol wet diameter ***/

    this->InitWetDiameter_aer(RelativeHumidity_f, this->Temperature_f,
			      this->Concentration_aer,this->NumberConcentration_aer,
			      WetDiameter_aer,false);

    /*** Snow height ***/

    if (this->option_manage["snow_height"])
      this->InitData("meteo", "SnowHeight", FileSnowHeight_i,
                     FileSnowHeight_f, this->current_date, SnowHeight_f);

    /*** Boundary conditions ***/

    if (this->option_manage["boundary_condition_aer"])
      {
	date = this->input_files["boundary_condition_aer"].GetDateMin();
	Delta_t
	  = this->input_files["boundary_condition_aer"].GetDelta_t();
	if(bc_format=="Internal"&&this->option_process["with_external_composition"])
	{
	  BoundaryCondition_z_aer_i.SetZero();
	  BoundaryCondition_y_aer_i.SetZero();
	  BoundaryCondition_x_aer_i.SetZero();
	  NumberBoundaryCondition_z_aer_i.SetZero();
	  NumberBoundaryCondition_z_aer_i_tmp.SetZero();
	  NumberBoundaryCondition_y_aer_i.SetZero();
	  NumberBoundaryCondition_y_aer_i_tmp.SetZero();
	  NumberBoundaryCondition_x_aer_i.SetZero();
	  NumberBoundaryCondition_x_aer_i_tmp.SetZero();
	}

	for (i = 0; i < Ns_bc_aer; i++)
	  {
	    species = species_list_bc_aer[i].first;
	    isize_section = species_list_bc_aer[i].second;
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["boundary_condition_aer"](species_bin);
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->InitData(find_replace(filename, "&c", "z"), date,
				Delta_t, FileBoundaryCondition_z_aer_i_tmp,
				FileBoundaryCondition_z_aer_f_tmp,
				this->current_date, i, j,
				BoundaryCondition_z_aer_i_tmp);
		  this->InitData(find_replace(filename, "&c", "y"), date,
				Delta_t, FileBoundaryCondition_y_aer_i_tmp,
				FileBoundaryCondition_y_aer_f_tmp,
				this->current_date, i,  j,
				BoundaryCondition_y_aer_i_tmp);
		  this->InitData(find_replace(filename, "&c", "x"), date,
				Delta_t, FileBoundaryCondition_x_aer_i_tmp,
				FileBoundaryCondition_x_aer_f_tmp,
				this->current_date, i,  j,
				BoundaryCondition_x_aer_i_tmp);
		}
		else
		{
		    this->InitData(find_replace(filename, "&c", "z"), date,
				  Delta_t, FileBoundaryCondition_z_aer_i,
				  FileBoundaryCondition_z_aer_f,
				  this->current_date, i, j,
				  BoundaryCondition_z_aer_i,this->Ncomposition_aer);
		    this->InitData(find_replace(filename, "&c", "y"), date,
				  Delta_t, FileBoundaryCondition_y_aer_i,
				  FileBoundaryCondition_y_aer_f,
				  this->current_date, i,  j,
				  BoundaryCondition_y_aer_i,this->Ncomposition_aer);
		    this->InitData(find_replace(filename, "&c", "x"), date,
				  Delta_t, FileBoundaryCondition_x_aer_i,
				  FileBoundaryCondition_x_aer_f,
				  this->current_date, i,  j,
				  BoundaryCondition_x_aer_i,this->Ncomposition_aer);
		}
	      }
	  }
	if (this->option_process["with_number_concentration"])
	  for (i = 0; i < Nb_bc_aer; i++)
	    {
	      species_bin = "Number_" + to_str(bc_bin_list_aer[i]);
	      string filename
		= this->input_files["boundary_condition_aer"](species_bin);
	      if (exists(find_replace(filename, "&c", "z")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->InitData(find_replace(filename, "&c", "z"), date,
				Delta_t, FileNumberBoundaryCondition_z_aer_i_tmp,
				FileNumberBoundaryCondition_z_aer_f_tmp,
				this->current_date, i,
				NumberBoundaryCondition_z_aer_i_tmp);
		}
		else
		{
		    this->InitData(find_replace(filename, "&c", "z"), date,
				  Delta_t, FileNumberBoundaryCondition_z_aer_i,
				  FileNumberBoundaryCondition_z_aer_f,
				  this->current_date, i,
				  NumberBoundaryCondition_z_aer_i,this->Ncomposition_aer);
		}
	      }


	      if (exists(find_replace(filename, "&c", "y")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->InitData(find_replace(filename, "&c", "y"), date,
				 Delta_t, FileNumberBoundaryCondition_y_aer_i_tmp,
				 FileNumberBoundaryCondition_y_aer_f_tmp,
				 this->current_date, i,
				 NumberBoundaryCondition_y_aer_i_tmp);
		}
		else
		{
		    this->InitData(find_replace(filename, "&c", "y"), date,
				  Delta_t, FileNumberBoundaryCondition_y_aer_i,
				  FileNumberBoundaryCondition_y_aer_f,
				  this->current_date, i,
				  NumberBoundaryCondition_y_aer_i,this->Ncomposition_aer);
		}
	      }

	      if (exists(find_replace(filename, "&c", "x")))
	      {
		if(bc_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->InitData(find_replace(filename, "&c", "x"), date,
				 Delta_t, FileNumberBoundaryCondition_x_aer_i_tmp,
				 FileNumberBoundaryCondition_x_aer_f_tmp,
				 this->current_date, i,
				 NumberBoundaryCondition_x_aer_i_tmp);
		}
		else
		{
		    this->InitData(find_replace(filename, "&c", "x"), date,
				  Delta_t, FileNumberBoundaryCondition_x_aer_i,
				  FileNumberBoundaryCondition_x_aer_f,
				  this->current_date, i,
				  NumberBoundaryCondition_x_aer_i,this->Ncomposition_aer);
		}
	      }
	    }
	if(bc_format=="Internal"&&this->option_process["with_external_composition"])
	{//for the computation of composition within each bins, loop must stated by bins
	  BoundaryConditionTransformation();
	}
	if (this->option_process["with_number_concentration"])
	  for (i = 0; i < Nb_bc_aer; i++)
	    {
	      species_bin = "Number_" + to_str(bc_bin_list_aer[i]);
	      string filename
		= this->input_files["boundary_condition_aer"](species_bin);
	      if (!exists(find_replace(filename, "&c", "z")))
		this->ComputeNumberBoundaryCondition_z_aer(int(bc_bin_list_aer[i]));

	      if (!exists(find_replace(filename, "&c", "y")))
		this->ComputeNumberBoundaryCondition_y_aer(int(bc_bin_list_aer[i]));

	      if (!exists(find_replace(filename, "&c", "x")))
		this->ComputeNumberBoundaryCondition_x_aer(int(bc_bin_list_aer[i]));
	    }	
      }

    /*** Deposition ***/

    if (this->option_manage["deposition_velocity_aer"])
      if (!this->option_process["compute_deposition_aer"])
	for (i = 0; i < Nbin_dep_aer; i++)
	{
	  date=this->input_files["deposition_velocity_aer"].GetDateMin();
	  Delta_t = this->input_files["deposition_velocity_aer"].GetDelta_t();
	  string filename = this->input_files["deposition_velocity_aer"](to_str(bin_list_dep_aer[i]));
	  TinyVector<int, 2> new_shape;
	  for (int j = 0; j < 2; j++)
	    new_shape(j) = DepositionVelocity_aer_f.GetArray().shape()(j + 1);

	  Data<T, 2> FileData_extract_i(new_shape);
	  Data<T, 2> FileData_extract_f(new_shape);
	  Data<T, 2> CurrentData_extract(new_shape);
	  this->InitData(filename, date, Delta_t,
			 FileData_extract_i,
			 FileData_extract_f,
			 this->current_date,
			 CurrentData_extract);
			 
	  for(int y=0; y<this->Ny; y++)
	    for(int x=0; x<this->Nx; x++)
	      for(int id=0; id< this->Ncomposition_aer; id++)
	      {
		int RealID=this->Ncomposition_aer*i+id;
		FileDepositionVelocity_aer_i(RealID,y,x)=FileData_extract_i(y,x);
		FileDepositionVelocity_aer_f(RealID,y,x)=FileData_extract_f(y,x);
		DepositionVelocity_aer_f(RealID,y,x)=CurrentData_extract(y,x);
	      }	
	}
      else 
	this->InitDepositionVelocity(this->Temperature_f, this->Pressure_f,
				     this->SurfaceTemperature_f,
				     this->SurfacePressure_f,
				     FirstLevelWindModule_f,
				     WetDiameter_aer, SnowHeight_f,
				     DepositionVelocity_aer_f);
	  
    /*** Scavenging coefficients ***/

    if (this->option_manage["scavenging_below_cloud_coefficient"]
        && this->scavenging_below_cloud_model == "microphysical-ph")
      InitScavengingCoefficient(this->Temperature_f, this->Pressure_f,
                                LiquidWaterContent_i, pH,
                                this->CloudBaseHeight_f, this->Rain_f,
                                this->ScavengingBelowCloudCoefficient_f);

    if (this->option_manage["scavenging_coefficient_aer"]
        && this->option_process["with_scavenging_aer"])
      this->InitScavengingCoefficient_aer(this->Temperature_f,
                                          this->Pressure_f,
                                          WetDiameter_aer,
                                          this->CloudBaseHeight_f,
                                          this->Rain_f,
                                          ScavengingCoefficient_aer_f);

    /*** Surface emissions ***/

    if (this->option_manage["surface_emission_aer"])
      {
	date = this->input_files["surface_emission_aer"].GetDateMin();
	Delta_t
	  = this->input_files["surface_emission_aer"].GetDelta_t();
	if (surface_emis_format == "Internal" &&
            this->option_process["with_external_composition"])
          {
            SurfaceEmission_aer_i.SetZero();
            SurfaceEmission_aer_f.SetZero();
            FileSurfaceEmission_aer_i.SetZero();
            FileSurfaceEmission_aer_f.SetZero();
            NumberSurfaceEmission_aer_i.SetZero();
            NumberSurfaceEmission_aer_f.SetZero();
            FileNumberSurfaceEmission_aer_i.SetZero();
            FileNumberSurfaceEmission_aer_f.SetZero();
          }
	for (i = 0; i < Ns_surf_emis_aer; i++)
	  {
	    species = species_list_surf_emis_aer[i].first;
	    isize_section = species_list_surf_emis_aer[i].second;
	    int CompositionID = FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["surface_emission_aer"](species_bin);
		if (surface_emis_format == "Internal" &&
                    this->option_process["with_external_composition"])
                  { 
                    this->InitData(filename, date, Delta_t,
                                   FileSurfaceEmission_aer_i,
                                   FileSurfaceEmission_aer_f,
                                   this->current_date, i, j,CompositionID,
                                   SurfaceEmission_aer_f);
                  }
		else
                  {
                    this->InitData(filename, date, Delta_t,
                                   FileSurfaceEmission_aer_i,
                                   FileSurfaceEmission_aer_f,
                                   this->current_date, i, j,
                                   SurfaceEmission_aer_f, this->Ncomposition_aer);
                  }
	      }
	  }
	if (this->option_process["with_number_concentration"])
	  for( j = 0; j < Nb_surf_emis_aer; j++)
	    {
	      species_bin = "Number_"  + to_str(surface_emission_bin_list_aer[j]);
	      string filename = this->input_files["surface_emission_aer"](species_bin);
	      if (exists(filename))
	      {
		if(surface_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  TinyVector<int, 2> new_shape;
		  for (int i = 0; i < 2; i++)
		    new_shape(i) = SurfaceEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 2> FileData_extract_i(new_shape);
		  Data<T, 2> FileData_extract_f(new_shape);
		  Data<T, 2> CurrentData_extract(new_shape);
		  this->InitData(filename, date, Delta_t,
				FileData_extract_i,
				FileData_extract_f,
				this->current_date,
				CurrentData_extract);
		  for(int y=0; y<this->Ny; y++)
		    for(int x=0; x<this->Nx; x++)
		    {
		      T TotalMass=0;
		      Data<T, 1> CompositionMass(this->Ncomposition_aer);
		      CompositionMass.SetZero();
		      for(int s=0; s<Ns_surf_emis_aer; s++)
		      {
			species = species_list_surf_emis_aer[s].first;
			int CompositionID=FindExternalCompositionID(species);
			int RealID=this->Ncomposition_aer*j+CompositionID;
			TotalMass+=SurfaceEmission_aer_f(s,RealID,y,x);
			CompositionMass(CompositionID)+=SurfaceEmission_aer_f(s,RealID,y,x);
		      }
		      for(int id=0; id<this->Ncomposition_aer; id++)
		      {
			if(TotalMass>0)
			{
			  int RealID=this->Ncomposition_aer*j+id;
			  NumberSurfaceEmission_aer_f(RealID,y,x)=CurrentData_extract(y,x)*CompositionMass(id)/TotalMass;
			  FileNumberSurfaceEmission_aer_i(RealID,y,x)=FileData_extract_i(y,x)*CompositionMass(id)/TotalMass;
			  FileNumberSurfaceEmission_aer_f(RealID,y,x)=FileData_extract_f(y,x)*CompositionMass(id)/TotalMass;
			}
		      }
		    }
		}
		else
		{
		  this->InitData(filename, date, Delta_t,
				  FileNumberSurfaceEmission_aer_i,
				  FileNumberSurfaceEmission_aer_f,
				  this->current_date, j,
				  NumberSurfaceEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	      else //need to be calculate at last
	      {
		this->ComputeNumberSurfaceEmission_aer
		  (int(surface_emission_bin_list_aer[j]));
	      }
	    }
      }

    /*** Volume emissions ***/

    if (this->option_manage["volume_emission_aer"])
      {
	date = this->input_files["volume_emission_aer"].GetDateMin();
	Delta_t
	  = this->input_files["volume_emission_aer"].GetDelta_t();
	for (i = 0; i < Ns_vol_emis_aer; i++)
	  {
	    species = species_list_vol_emis_aer[i].first;
	    isize_section= species_list_vol_emis_aer[i].second;
	    int CompositionID=FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["volume_emission_aer"](species_bin);
		if(volume_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  this->InitData(filename, date, Delta_t,
				FileVolumeEmission_aer_i,
				FileVolumeEmission_aer_f,
				this->current_date,
				i, j, CompositionID,
				VolumeEmission_aer_f);
		}
		else
		{
		  this->InitData(filename, date, Delta_t,
				FileVolumeEmission_aer_i,
				FileVolumeEmission_aer_f,
				this->current_date, i, j,
				VolumeEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	  }
	if (this->option_process["with_number_concentration"])
	  for( j = 0; j < Nb_vol_emis_aer; j++)
	    {
	      species_bin = "Number_"  + to_str(volume_emission_bin_list_aer[j]);
	      string filename = this->input_files["volume_emission_aer"](species_bin);

	      if (exists(filename))
	      {
		if(volume_emis_format=="Internal"&&this->option_process["with_external_composition"])
		{
		  TinyVector<int, 3> new_shape;
		  for (int i = 0; i < 3; i++)
		    new_shape(i) = VolumeEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 3> FileData_extract_i(new_shape);
		  Data<T, 3> FileData_extract_f(new_shape);
		  Data<T, 3> CurrentData_extract(new_shape);
		  this->InitData(filename, date, Delta_t,
				FileData_extract_i,
				FileData_extract_f,
				this->current_date,
				CurrentData_extract);//Nz_vol_emis_aer
		  for(int z=0; z<Nz_vol_emis_aer; z++)
		    for(int y=0; y<this->Ny; y++)
		      for(int x=0; x<this->Nx; x++)
		      {
			T TotalMass=0;
			Data<T, 1> CompositionMass(this->Ncomposition_aer);
			CompositionMass.SetZero();
			for(int s=0; s<Ns_surf_emis_aer; s++)
			{
			  species = species_list_surf_emis_aer[s].first;
			  int CompositionID=FindExternalCompositionID(species);
			  int RealID=this->Ncomposition_aer*j+CompositionID;
			  TotalMass+=VolumeEmission_aer_f(s,RealID,z,y,x);
			  CompositionMass(CompositionID)+=VolumeEmission_aer_f(s,RealID,z,y,x);
			}
			for(int id=0; id<this->Ncomposition_aer; id++)
			{
			  if(TotalMass>0)
			  {
			    int RealID=this->Ncomposition_aer*j+id;
			    NumberVolumeEmission_aer_f(RealID,z,y,x)=CurrentData_extract(z,y,x)*CompositionMass(id)/TotalMass;
			    FileNumberVolumeEmission_aer_i(RealID,z,y,x)=FileData_extract_i(z,y,x)*CompositionMass(id)/TotalMass;
			    FileNumberVolumeEmission_aer_f(RealID,z,y,x)=FileData_extract_f(z,y,x)*CompositionMass(id)/TotalMass;
			  }
			}
		      }
		}
		else
		{
		this->InitData(filename, date, Delta_t,
			       FileNumberVolumeEmission_aer_i,
			       FileNumberVolumeEmission_aer_f,
			       this->current_date, j,
			       NumberVolumeEmission_aer_f, this->Ncomposition_aer);
		}
	      }
	      else
		this->ComputeNumberVolumeEmission_aer
		  (int(volume_emission_bin_list_aer[j]));
	    }
      }

    this->data_date = this->current_date;
    this->data_gone_through_initstep = false;
    this->data_gone_through_forward = false;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DAEROSOL_CXX
#endif
