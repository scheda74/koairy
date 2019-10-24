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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMEINGRIDAEROSOL_CXX


#include "PlumeInGridAerosol.hxx"
#include <time.h>

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*! Builds the model and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::PlumeInGridAerosol(string config_file):
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>(config_file)
  {
    this->D3_map["DryDepositionFluxNumber_aer"] = &DryDepositionFluxNumber_aer;
    this->D3_map["WetDepositionFluxNumber_aer"] = &WetDepositionFluxNumber_aer;
    this->D3_map["InCloudWetDepositionFluxNumber_aer"]
      = &InCloudWetDepositionFluxNumber_aer;
    this->D4_map["DryDepositionFlux_aer"] = &DryDepositionFlux_aer;
    this->D4_map["WetDepositionFlux_aer"] = &WetDepositionFlux_aer;
    this->D4_map["InCloudWetDepositionFlux_aer"] = &InCloudWetDepositionFlux_aer;
    this->field_bins["DepositionVelocity_aer"] = &bin_list_dep_aer;
    this->field_bins["ScavengingCoefficient_aer"] = &bin_list_scav_aer;

  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::~PlumeInGridAerosol()
  {
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    the species list and display options.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
      ::ReadConfiguration();

    /*** Species and bins ***/

    // Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    // Section "[aerosol_species]" contains all aerosol species names.
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());

    // Reads bin bounds.
    this->config.SetSection("[domain]");
    this->config.Find("Bin_bounds");
    bin_list = split(this->config.GetLine());
    this->Nbin_aer = int(bin_list.size()) - 1;

    this->config.SetSection("[options]");
    //For number computation.
    this->config.PeekValue("With_number_concentration",
			   this->option_process["with_number_concentration"]);
    this->config.PeekValue("With_deposition_aerosol",
                           this->option_process["with_deposition_aer"]);
   this->config.PeekValue("Collect_dry_flux_aerosol",
			   this->option_process["collect_dry_flux_aer"]);
    this->config.PeekValue("Collect_wet_flux_aerosol",
			   this->option_process["collect_wet_flux_aer"]);

    string data_file;
    int bin_index;
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", data_file);

    this->config.SetSection("[options]");
    if (this->option_process["with_deposition_aer"])
      this->config.PeekValue("Compute_deposition_aerosol",
                             this->option_process["compute_deposition_aer"]);
    else
      this->option_process["compute_deposition_aer"] = false;

    // Deposition velocities.
    if (this->option_process["with_deposition_aer"])
      if (this->option_process["compute_deposition_aer"])
	this->input_files["deposition_velocity_aer"]
	  .ReadFields(data_file, "deposition_velocity_aerosol");
      else
	this->input_files["deposition_velocity_aer"]
	  .Read(data_file, "deposition_velocity_aerosol");
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

    this->config.PeekValue("With_scavenging_aerosol",
			   this->option_process["with_scavenging_aer"]);
    if (this->option_process["with_scavenging_aer"])
      {
	this->input_files["scavenging_aer"]
	  .ReadFields(data_file, "scavenging_aerosol");
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
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::CheckConfiguration()
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>::CheckConfiguration();
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.

  }

  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::Allocate()
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>::Allocate();

    /*** Additional grids ***/

    GridZ5D = RegularGrid<T>(this->Nz);
    GridY5D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridX5D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);

    GridZ5D.SetVariable(2);
    GridZ5D.SetDuplicate(false);

    GridY5D.SetVariable(3);
    GridY5D.SetDuplicate(false);

    GridX5D.SetVariable(4);
    GridX5D.SetDuplicate(false);

    // test //
    GridZ4D = RegularGrid<T>(this->Nz);
    GridY4D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridX4D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);

    GridZ4D.SetVariable(1);
    GridZ4D.SetDuplicate(false);

    GridY4D.SetVariable(2);
    GridY4D.SetDuplicate(false);

    GridX4D.SetVariable(3);
    GridX4D.SetDuplicate(false);
    GridS4D = RegularGrid<T>(this->Ns);


    // test //

    GridB5D_aer = RegularGrid<T>(this->Nbin_aer);
    GridS5D_aer = RegularGrid<T>(this->Ns_aer);
    GridS4D_number = RegularGrid<T>(this->Nbin_aer);

    /*** State ***/

    this->Concentration_aer.Resize(GridS5D_aer, GridB5D_aer,
				   GridZ5D, GridY5D, GridX5D);
    this->delta_number_puff.Resize(GridS4D_number,
				   GridZ4D, GridY4D, GridX4D);

    // this->Concentration_pre_aer.Resize(GridS5D_aer, GridB5D_aer,
    // 				   GridZ5D, GridY5D, GridX5D);

    // this->Concentration_pre.Resize(GridS4D,
    // 				   GridZ4D, GridY4D, GridX4D);
    this->GasCompensation.Resize(GridS4D,
				   GridZ4D, GridY4D, GridX4D); 
    if (this->option_process["with_number_concentration"])
      {
	this->NumberConcentration_aer.Resize(GridS4D_number, 
					     GridZ4D, GridY4D, GridX4D);
      }
  }
  //! Initialization.
  /*! Initializes the Eulerian model and the meteorological conditions. Reads
    the sources and creates the corresponding Gaussian model instances.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
      ::Init();

   /*! Concentration_aer is initialized during the GaussianPuff model 
     initialization. This is needed for Saver initialization.
   */ 
    int l,s,b,k,i;

    this->Concentration_aer.Copy(this->Model.GetConcentration_aer());
    
    if (this->option_process["collect_dry_flux"])
      {
	this->DryDepositionFlux_aer.Copy(this->Model.GetDryDepositionFlux_aer());
	this->DryDepositionFluxNumber_aer.Copy(this->Model.GetDryDepositionFluxNumber_aer());
      }
    if (this->option_process["collect_wet_flux"])
      {
	this->WetDepositionFlux_aer.Copy(this->Model.GetWetDepositionFlux_aer());
	this->InCloudWetDepositionFlux_aer.Copy(this->Model.GetInCloudWetDepositionFlux_aer());
      
	this->WetDepositionFluxNumber_aer.Copy(this->Model.GetWetDepositionFluxNumber_aer());
	this->InCloudWetDepositionFluxNumber_aer.Copy(this->Model.GetInCloudWetDepositionFluxNumber_aer());
      }

    if (this->option_process["with_number_concentration"])
      {
	this->NumberConcentration_aer.Copy(this->Model.GetNumberConcentration_aer());
      }
  
    // Meteorological data. 
    LiquidWaterContent_i.Copy(this->Model.D3("LiquidWaterContent_i"));
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::InitStep()
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
      ::InitStep();

    LiquidWaterContent_i.Copy(this->Model.D3("LiquidWaterContent_i"));
  }


//   //! Performs one step forward.
//   /*! Calls the Eulerian model. For each major point emission, the
//     corresponding Gaussian Model instance is called.
//   */
//   template<class T, class ClassEulerianModel, class ClassLocalModel>
//   void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
//   ::Forward()
//   {
//     int rank, rank_size;
// #ifdef POLYPHEMUS_PARALLEL_WITH_MPI
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
// #else
//     rank = 0;
// #endif

//     cout << "PinG model forward " << rank_size << " " << rank<<  endl;
    
//     int l, k, j, i, s, b;

//     this->Concentration.Copy(this->Model.GetConcentration());
//     this->Concentration_aer.Copy(this->Model.GetConcentration_aer());
   
//     if (this->option_process["collect_dry_flux"])      
//       {
// 	this->DryDepositionFlux_aer.Copy
//           (this->Model.GetDryDepositionFlux_aer());
// 	this->DryDepositionFluxNumber_aer.Copy
//           (this->Model.GetDryDepositionFluxNumber_aer());
//       }
    
//     if (this->option_process["collect_wet_flux"])
//       {
// 	this->WetDepositionFlux_aer.Copy
//           (this->Model.GetWetDepositionFlux_aer());
// 	this->InCloudWetDepositionFlux_aer.Copy
//           (this->Model.GetInCloudWetDepositionFlux_aer());

// 	this->WetDepositionFluxNumber_aer.Copy
//           (this->Model.GetWetDepositionFluxNumber_aer());
// 	this->InCloudWetDepositionFluxNumber_aer.Copy
//           (this->Model.GetInCloudWetDepositionFluxNumber_aer());
//       }
    
//     if (this->option_process["with_number_concentration"])
//       {
// 	this->NumberConcentration_aer.Copy
//           (this->Model.GetNumberConcentration_aer());
//       }

//     int Npuff, puff_index;
//     bool isday;
//     if (rank == 0)
//       {

// 	int index_x, index_y, index_z;
// 	bool puff_transfer;

// 	Npuff = this->GaussianPuffModel->GetPuffNumber();
// 	/*** Tests if puff has to be transfered to Eulerian model. ***/
// 	for (puff_index = 0; puff_index < Npuff; puff_index++)
// 	  {
// 	    // Puff center coordinates.
// 	    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
// 	    this->GaussianPuffModel->
// 	      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
// 			      puff_time);

// 	    // Conversion to longitude/latitude (degrees).
// 	    this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

// 	    // Night or day.
// 	    isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->
// 			  GetCurrentDate());

// 	    // Gets corresponding cell in eulerian grid.
// 	    this->GetCellIndices(lon_c, lat_c, z_c,
// 				 index_z, index_y, index_x);

// 	    // Computing cell width along y.
// 	    T cell_width_z, cell_width_y, cell_width_x, cell_volume;
// 	    this->ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
// 				   cell_width_y, cell_width_x, cell_volume);

// 	    // Getting puff standard deviations.
// 	    T sigma_x, sigma_y, sigma_z;
// 	    this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
// 						  sigma_y, sigma_z);
// 	    puff_transfer = 0;

	
// 	    T puff_release_time = this->GaussianPuffModel->
// 	      GetPuffReleaseTime(puff_index); // YK
	    
// 	    // Tests if puff has reached the end of the domain.
// 	    if (index_x == 0 || index_x == this->Nx - 1
// 		|| index_y == 0 || index_y == this->Ny - 1
// 		|| index_z == this->Nz)
// 	      this->GaussianPuffModel->ErasePuff(puff_index);

// 	    // Tests if puff has reached the injection time.
// 	    else if (this->option_time)
// 	      {
// 		if (puff_time >= this->reinjection_time)
// 		  puff_transfer = 1;
// 		// Tests if puff size has reached the cell width.
// 		else if (this->coefficient_y * sigma_y >= cell_width_y)	
// 		  puff_transfer = 1;
// 	      }

// 	    // Tests if puff size has reached the cell width.
// 	    else if (this->coefficient_y * sigma_y >= cell_width_y)
// 	      puff_transfer = 1;
	    
// 	    // Puff transfer.
// 	    if (puff_transfer)
// 	      {
// 		if (this->injection_integrated)
// 		  {
// 		    this->delta_number_puff.SetZero();
// 		    this->PuffIntegratedTransfer
// 		      (puff_index, this->Model.GetConcentration());
		    
// 		    this->PuffIntegratedTransfer_aer
// 		      (puff_index, this->Model.GetConcentration_aer(), 
// 		       this->Model.GetConcentration());
// 		    if (this->option_process["with_number_concentration"])
// 		      this->PuffIntegratedTransfer_number
// 			(puff_index, this->Model.GetConcentration_aer(), 
// 			 this->Model.GetNumberConcentration_aer());
// 		  }
// 		else
// 		  {
// 		    for (int s = 0; s < this->Ns; s++)
// 		      this->PuffTransfer(this->GaussianPuffModel->
// 					 GetPuffQuantity(puff_index, s),
// 					 sigma_z, s, z_c, lat_c, lon_c,
// 					 isday, this->Model.GetConcentration());
// 		    for (int s = 0; s < this->Ns_aer; s++)
// 		      for (int b = 0; b < this->Nbin_aer; b++)
// 			this->PuffTransfer_aer(this->GaussianPuffModel->
// 					       GetPuffQuantity_aer(puff_index, s, b), 
// 					       sigma_z, 
// 					       s,  z_c, lat_c, lon_c,
// 					       isday, this->Model.GetConcentration_aer(), b);
// 		    if (this->option_process["with_number_concentration"])
// 		      for (int b = 0; b < this->Nbin_aer; b++)
// 			this->PuffTransfer_number(this->GaussianPuffModel->
// 						  GetPuffQuantity_number(puff_index, b), 
// 						  sigma_z, 
// 						  z_c, lat_c, lon_c,
// 						  isday, this->Model.GetNumberConcentration_aer(), b);
// 		  }
// 		this->GaussianPuffModel->ErasePuff(puff_index);
// 	      }
// 	  }
	
// 	/*** Inner time-loop for Gaussian models. ***/
	
// 	Array<T, 4> PuffConcentration(this->Ns, this->Nz, this->Ny, this->Nx);
// 	Array<T, 5> PuffConcentration_aer(this->Ns_aer, this->Nbin_aer, 
// 					  this->Nz, this->Ny, this->Nx);
// 	PuffConcentration = 0.;
// 	PuffConcentration_aer = 0.;
	
// 	if (this->option_process["with_number_concentration"])
// 	  {
// 	    Array<T, 4> PuffConcentration_number(this->Nbin_aer, this->Nz, 
// 						 this->Ny, this->Nx);
// 	    PuffConcentration_number = 0.;
// 	  }
	
// 	for (j = 0; j < this->Nt_local; j++)
// 	  {
// 	    this->GaussianPuffModel->InitStep();
// 	    Npuff = this->GaussianPuffModel->GetPuffNumber();
// 	    Array<int, 2> PuffCellList(Npuff, 3);
	    
// 	    // Loop on puffs to get the meteorological data.
// 	    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time, volume_puff;
// 	    T CellWidth_y_tmp, CellWidth_x_tmp, CellWidth_z_tmp, CellVolume_tmp;
// 	    int jj;
// 	    for (puff_index = 0; puff_index < Npuff; puff_index++)
// 	      {
// 		// Puff center coordinates.
// 		this->GaussianPuffModel->
// 		  GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
// 				  puff_time);
		
// 		// Conversion to longitude/latitude (degrees).
// 		this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
		
// 		// Gets corresponding cell in Eulerian grid.
// 		this->GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);
// 		PuffCellList(puff_index, 0) = index_z;
// 		PuffCellList(puff_index, 1) = index_y;
// 		PuffCellList(puff_index, 2) = index_x;

// 		// Updating the puff meteorological data.
// 		if (j == 0 || puff_time == 0.)
// 		  UpdateMeteo(puff_index);
// 	      }
	    
// 	    // Puff effective height.
// 	    if (!this->evolutive_plume_rise)
// 	      this->GaussianPuffModel->ComputePlumeRise();
	    
// 	    // Advection and diffusion for puffs.
// 	    this->GaussianPuffModel->Advection();
	    
// 	    if (this->evolutive_plume_rise)
// 	      this->GaussianPuffModel->ComputeEvolutivePlumeRise();
// 	    this->GaussianPuffModel->Diffusion();
// 	    this->GaussianPuffModel->ComputeLossFactor();
	    
// 	    if (this->merge_puff)
// 	      {
// 		//Combine overlaping puff 
// 		int puff_index_alpha=0;
// 		int puff_index_beta;
// 		T puff_alpha_volume;
// 		T puff_beta_volume;
// 		T puff_alpha_beta_volume;
// 		T overlap_tmp;
// 		T tmp2;
		
		
// 		Array<vector<int>, 1> PuffInteractionList(Npuff);
// 		vector<int> PuffList_tmp;
// 		vector<int> Puff_erased;
// 		int combine_puff;
// 		int Npuff_combine;
// 		int Npuff_erase;
// 		string puff_alpha_id;
// 		string puff_beta_id;
		
// 		// Loop to define which puff will be merged
// 		for (puff_index_alpha = 0; puff_index_alpha < Npuff; puff_index_alpha++)
// 		  {
// 		    puff_index_beta = 0;
// 		    PuffList_tmp.clear();
// 		    puff_alpha_id = this->GaussianPuffModel
// 		      ->GetPuffSourceId(puff_index_alpha);
// 		    for (puff_index_beta = 0; puff_index_beta < Npuff; puff_index_beta++)
// 		      {
// 			puff_beta_id = this->GaussianPuffModel
// 			  ->GetPuffSourceId(puff_index_beta);
			
// 			puff_alpha_volume = this->GaussianPuffModel
// 			  ->ComputePuffOverlap(puff_index_alpha, puff_index_alpha);
// 			puff_beta_volume = this->GaussianPuffModel
// 			  ->ComputePuffOverlap(puff_index_beta, puff_index_beta);
// 			puff_alpha_beta_volume = this->GaussianPuffModel
// 			  ->ComputePuffOverlap(puff_index_alpha, puff_index_beta);
			
// 			if (puff_index_alpha == puff_index_beta)
// 			  puff_alpha_beta_volume = 0.;
			
// 			if (puff_alpha_beta_volume == 0.)
// 			  overlap_tmp = 0.;
// 			else
// 			  overlap_tmp = puff_alpha_beta_volume / (puff_alpha_volume * puff_beta_volume);
// 			// If overlap between puff A and B is greater than 80% 
// 			if (overlap_tmp > 0.8 *  (1. / puff_alpha_volume))
// 			  if ((1. / puff_alpha_volume) > (1. / puff_beta_volume))
// 			    {		  
// 			      if (this->merge_source_id)
// 				{
// 				  if (puff_alpha_id == puff_beta_id)
// 				    PuffList_tmp.push_back(puff_index_beta);
// 				}
// 			      else
// 				PuffList_tmp.push_back(puff_index_beta);
// 			    }
// 		      }
// 		    //List contains puff to be merged with puff A
// 		    PuffInteractionList(puff_index_alpha) = PuffList_tmp;
// 		  }
		
// 		int j_tmp, i_tmp;
// 		//Combine overlaping puff
// 		for (puff_index_alpha = 0; puff_index_alpha < Npuff; puff_index_alpha++)
// 		  {
// 		    PuffList_tmp = PuffInteractionList(puff_index_alpha);
// 		    Npuff_combine = PuffList_tmp.size();
// 		    for (j_tmp = 0; j_tmp < Npuff_combine; j_tmp++)
// 		      {
// 			combine_puff = 0;
// 			puff_index_beta = PuffList_tmp[j_tmp];
// 			Npuff_erase = Puff_erased.size();
// 			//Check if puff has already been merged
// 			for (i_tmp=0; i_tmp<Npuff_erase; i_tmp++)
// 			  if (puff_index_beta == Puff_erased[i_tmp])
// 			    combine_puff += 1;
// 			if (combine_puff == 0)
// 			  {
// 			    puff_alpha_volume = 1. / this->GaussianPuffModel
// 			      ->ComputePuffOverlap(puff_index_alpha, puff_index_alpha);
// 			    puff_beta_volume = 1. / this->GaussianPuffModel
// 			      ->ComputePuffOverlap(puff_index_beta, puff_index_beta);
// 			    this->GaussianPuffModel
// 			      ->CombineOverlappingPuff(puff_index_alpha, 
// 						       puff_index_beta, puff_alpha_volume,
// 						       puff_beta_volume);
// 			    Puff_erased.push_back(puff_index_beta);
// 			  }
// 		      }
// 		  }
		
// 		// Erase merged puffs
// 		Npuff_erase = Puff_erased.size();
// 		sort(Puff_erased.begin(), Puff_erased.end());
// 		for (i_tmp = 0; i_tmp < Npuff_erase; i_tmp ++)
// 		  this->GaussianPuffModel->ErasePuff(Puff_erased[Npuff_erase - 1 - i_tmp]);
// 		Puff_erased.clear();
// 	      }
// 	    Npuff = this->GaussianPuffModel->GetPuffNumber();
	    
// 	    for (puff_index = 0; puff_index < Npuff; puff_index++)
// 	      {
// 		T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
// 		this->GaussianPuffModel->
// 		  GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
// 				  puff_time);
// 		// Conversion to longitude/latitude (degrees).
// 		this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
		
// 		// Gets corresponding cell in eulerian grid.
// 		this->GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);
// 		UpdateMeteo(puff_index);
// 		PuffCellList(puff_index,0) = index_z;
// 		PuffCellList(puff_index,1) = index_y;
// 		PuffCellList(puff_index,2) = index_x;
// 	      }
	    
// 	    // Chemistry for puffs.
// 	    this->GaussianPuffModel->Chemistry();

// 	    // Save the quantities of the species in the puff 
// 	    // if (rank == 0)
// 	    this->GaussianPuffModel->PuffOutputSaver();

// 	    this->GaussianPuffModel->AddTime(this->delta_t_local);
// 	  } // GaussianPuffModel iteration 
//       }

//     // Eulerian model.
//     this->Model.Forward();

//     this->Concentration.Copy(this->Model.GetConcentration());
//     this->Concentration_aer.Copy(this->Model.GetConcentration_aer());

//     if (this->option_process["collect_dry_flux"])      
//       {
// 	this->DryDepositionFlux_aer.Copy(this->Model.GetDryDepositionFlux_aer());
// 	this->DryDepositionFluxNumber_aer.Copy(this->Model.GetDryDepositionFluxNumber_aer());
//       }

//     if (this->option_process["collect_wet_flux"])
//       {
// 	this->WetDepositionFlux_aer.Copy(this->Model.GetWetDepositionFlux_aer());
// 	this->InCloudWetDepositionFlux_aer.Copy(this->Model.GetInCloudWetDepositionFlux_aer());
	
// 	this->WetDepositionFluxNumber_aer.Copy(this->Model.GetWetDepositionFluxNumber_aer());
// 	this->InCloudWetDepositionFluxNumber_aer.Copy(this->Model.GetInCloudWetDepositionFluxNumber_aer());
//       }

//     if (this->option_process["with_number_concentration"])
//       this->NumberConcentration_aer.Copy(this->Model.GetNumberConcentration_aer());


//     if (rank == 0)
//       {
// 	/* Adding all puffs to concentration Data in case concentrations are saved
// 	   on the domain.*/
// 	this->GaussianPuffModel->SubtractTime(this->delta_t_local);
// 	Npuff = this->GaussianPuffModel->GetPuffNumber();
// 	for (puff_index = 0; puff_index < Npuff; puff_index++)
// 	  {
// 	    this->GasCompensation.SetZero();
// 	    T puff_release_time = this->GaussianPuffModel->
// 	      GetPuffReleaseTime(puff_index);
// 	    if (this->GaussianPuffModel->GetCurrentTime() >= puff_release_time)
// 	      {
// 		// Puff center coordinates.
// 		T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
// 		this->GaussianPuffModel->
// 		  GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
// 				  puff_time);
		
// 		// Conversion to longitude/latitude (degrees).
// 		this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
		
// 		//Get indices of the Eulerian cell
// 		int ind_x, ind_y, ind_z;
// 		this->GetCellIndices(lon_c,lat_c,z_c, ind_z, ind_y, ind_x);
		
// 		T sigma_x, sigma_y, sigma_z;
// 		this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
// 						      sigma_y, sigma_z);
		
// 		isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->
// 			      GetCurrentDate());
// 		// Adding puff concentration to the Concentration Data.
// 		if (this->injection_integrated)
// 		  {
// 		    this->delta_number_puff.SetZero();
// 		    this->PuffIntegratedTransfer(puff_index, this->Concentration);
// 		    this->PuffIntegratedTransfer_aer(puff_index, this->Concentration_aer, this->Concentration);
// 		    if (this->option_process["with_number_concentration"])
// 		      this->PuffIntegratedTransfer_number(puff_index, this->Concentration_aer, 
// 							  this->NumberConcentration_aer);
// 		  }
// 		else
// 		  {
// 		    for (int s = 0; s < this->Ns; s++)
// 		      this->PuffTransfer(this->GaussianPuffModel->
// 					 GetPuffQuantity(puff_index, s),
// 					 sigma_z, s, z_c, lat_c, lon_c,
// 					 isday, this->Concentration);
// 		    for (int s = 0; s < this->Ns_aer; s++)
// 		      for (int b = 0; b < this->Nbin_aer; b++)
// 			this->PuffTransfer_aer(this->GaussianPuffModel->
// 					       GetPuffQuantity_aer(puff_index, s, b),
// 					       sigma_z,
// 					       s,  z_c, lat_c, lon_c,
// 					       isday, this->Concentration_aer, b);
// 		    if (this->option_process["with_number_concentration"])
// 		      for (int b = 0; b < this->Nbin_aer; b++)
// 			this->PuffTransfer_number(this->GaussianPuffModel->
// 						  GetPuffQuantity_number(puff_index, b), 
// 						  sigma_z, 
// 						  z_c, lat_c, lon_c,
// 						  isday, this->NumberConcentration_aer, b);

// 		  }
// 	      }
// 	  }

// 	this->GaussianPuffModel->AddTime(this->delta_t_local);
//       }
//     this->AddTime(this->delta_t_eulerian);
//     this->step++;
//   }

  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    int rank;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif    
    //Gaussian model is not parrallelized;
    if (rank == 0)
      this->Forward_Gaussian();
    
    this->Forward_Euler();   

    if (rank == 0 && this->save_gaussian_domain)
      this->AddGaussianConcentrationsDomain();
  }

  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::Forward_Gaussian()
  {
    int l, k, j, i, s, b;
    this->Concentration.Copy(this->Model.GetConcentration());
    this->Concentration_aer.Copy(this->Model.GetConcentration_aer());
   
    if (this->option_process["collect_dry_flux"])      
      {
	this->DryDepositionFlux_aer.Copy
          (this->Model.GetDryDepositionFlux_aer());
	this->DryDepositionFluxNumber_aer.Copy
          (this->Model.GetDryDepositionFluxNumber_aer());
      }
    
    if (this->option_process["collect_wet_flux"])
      {
	this->WetDepositionFlux_aer.Copy
          (this->Model.GetWetDepositionFlux_aer());
	this->InCloudWetDepositionFlux_aer.Copy
          (this->Model.GetInCloudWetDepositionFlux_aer());

	this->WetDepositionFluxNumber_aer.Copy
          (this->Model.GetWetDepositionFluxNumber_aer());
	this->InCloudWetDepositionFluxNumber_aer.Copy
          (this->Model.GetInCloudWetDepositionFluxNumber_aer());
      }
    
    if (this->option_process["with_number_concentration"])
      {
	this->NumberConcentration_aer.Copy
          (this->Model.GetNumberConcentration_aer());
      }
    
    int Npuff;
    int index_x, index_y, index_z, puff_index;
    bool puff_transfer;
    bool isday;

    Npuff = this->GaussianPuffModel->GetPuffNumber();
    /*** Tests if puff has to be transfered to Eulerian model. ***/
    for (puff_index = 0; puff_index < Npuff; puff_index++)
      {
        // Puff center coordinates.
        T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
        this->GaussianPuffModel->
          GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                          puff_time);

        // Conversion to longitude/latitude (degrees).
        this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

        // Night or day.
        isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->
                      GetCurrentDate());

        // Gets corresponding cell in eulerian grid.
        this->GetCellIndices(lon_c, lat_c, z_c,
                             index_z, index_y, index_x);

        // Computing cell width along y.
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        this->ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
                               cell_width_y, cell_width_x, cell_volume);

        // Getting puff standard deviations.
        T sigma_x, sigma_y, sigma_z;
        this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                              sigma_y, sigma_z);
        puff_transfer = 0;

	
	T puff_release_time = this->GaussianPuffModel->
	  GetPuffReleaseTime(puff_index); // YK
	    
        // Tests if puff has reached the end of the domain.
        if (index_x == 0 || index_x == this->Nx - 1
            || index_y == 0 || index_y == this->Ny - 1
            || index_z == this->Nz)
          this->GaussianPuffModel->ErasePuff(puff_index);

        // Tests if puff has reached the injection time.
        else if (this->option_time)
          {
            if (puff_time >= this->reinjection_time)
              puff_transfer = 1;
	    // Tests if puff size has reached the cell width.
	    else if (this->coefficient_y * sigma_y >= cell_width_y)	
	      puff_transfer = 1;
          }

        // Tests if puff size has reached the cell width.
        else if (this->coefficient_y * sigma_y >= cell_width_y)
          puff_transfer = 1;

        // Puff transfer.
        if (puff_transfer)
          {
            if (this->injection_integrated)
              {
		this->delta_number_puff.SetZero();
		this->PuffIntegratedTransfer
                  (puff_index, this->Model.GetConcentration());

		this->PuffIntegratedTransfer_aer
                  (puff_index, this->Model.GetConcentration_aer(), 
                   this->Model.GetConcentration());
		if (this->option_process["with_number_concentration"])
		  this->PuffIntegratedTransfer_number
                    (puff_index, this->Model.GetConcentration_aer(), 
                     this->Model.GetNumberConcentration_aer());
              }
            else
              {
		for (int s = 0; s < this->Ns; s++)
		  this->PuffTransfer(this->GaussianPuffModel->
				     GetPuffQuantity(puff_index, s),
				     sigma_z, s, z_c, lat_c, lon_c,
				     isday, this->Model.GetConcentration());
		for (int s = 0; s < this->Ns_aer; s++)
		  for (int b = 0; b < this->Nbin_aer; b++)
		    this->PuffTransfer_aer(this->GaussianPuffModel->
					   GetPuffQuantity_aer(puff_index, s, b), 
					   sigma_z, 
					   s,  z_c, lat_c, lon_c,
					   isday, this->Model.GetConcentration_aer(), b);
		if (this->option_process["with_number_concentration"])
		  for (int b = 0; b < this->Nbin_aer; b++)
		    this->PuffTransfer_number(this->GaussianPuffModel->
					      GetPuffQuantity_number(puff_index, b), 
					      sigma_z, 
					      z_c, lat_c, lon_c,
					      isday, this->Model.GetNumberConcentration_aer(), b);
              }
            this->GaussianPuffModel->ErasePuff(puff_index);
          }
      }

    /*** Inner time-loop for Gaussian models. ***/

    Array<T, 4> PuffConcentration(this->Ns, this->Nz, this->Ny, this->Nx);
    Array<T, 5> PuffConcentration_aer(this->Ns_aer, this->Nbin_aer, 
                                      this->Nz, this->Ny, this->Nx);
    PuffConcentration = 0.;
    PuffConcentration_aer = 0.;
   
    if (this->option_process["with_number_concentration"])
      {
	Array<T, 4> PuffConcentration_number(this->Nbin_aer, this->Nz, 
                                             this->Ny, this->Nx);
	PuffConcentration_number = 0.;
      }

    for (j = 0; j < this->Nt_local; j++)
      {
        this->GaussianPuffModel->InitStep();
        Npuff = this->GaussianPuffModel->GetPuffNumber();
        Array<int, 2> PuffCellList(Npuff, 3);

        // Loop on puffs to get the meteorological data.
	T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time, volume_puff;
	T CellWidth_y_tmp, CellWidth_x_tmp, CellWidth_z_tmp, CellVolume_tmp;
	int jj;
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            // Puff center coordinates.
            this->GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

            // Gets corresponding cell in Eulerian grid.
            this->GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);
            PuffCellList(puff_index, 0) = index_z;
            PuffCellList(puff_index, 1) = index_y;
            PuffCellList(puff_index, 2) = index_x;

            // Updating the puff meteorological data.
            if (j == 0 || puff_time == 0.)
              UpdateMeteo(puff_index);
          }

        // Puff effective height.
	if (!this->evolutive_plume_rise)
	  this->GaussianPuffModel->ComputePlumeRise();

        // Advection and diffusion for puffs.
        this->GaussianPuffModel->Advection();

	if (this->evolutive_plume_rise)
	  this->GaussianPuffModel->ComputeEvolutivePlumeRise();
        this->GaussianPuffModel->Diffusion();
        this->GaussianPuffModel->ComputeLossFactor();

	if (this->merge_puff)
	  {
	    //Combine overlaping puff 
	    int puff_index_alpha=0;
	    int puff_index_beta;
	    T puff_alpha_volume;
	    T puff_beta_volume;
	    T puff_alpha_beta_volume;
	    T overlap_tmp;
	    T tmp2;


	    Array<vector<int>, 1> PuffInteractionList(Npuff);
	    vector<int> PuffList_tmp;
	    vector<int> Puff_erased;
	    int combine_puff;
	    int Npuff_combine;
	    int Npuff_erase;
	    string puff_alpha_id;
	    string puff_beta_id;
	   
	    // Loop to define which puff will be merged
	    for (puff_index_alpha = 0; puff_index_alpha < Npuff; puff_index_alpha++)
	      {
		puff_index_beta = 0;
		PuffList_tmp.clear();
		puff_alpha_id = this->GaussianPuffModel
		  ->GetPuffSourceId(puff_index_alpha);
		for (puff_index_beta = 0; puff_index_beta < Npuff; puff_index_beta++)
		  {
		    puff_beta_id = this->GaussianPuffModel
		      ->GetPuffSourceId(puff_index_beta);
		    
		    puff_alpha_volume = this->GaussianPuffModel
		      ->ComputePuffOverlap(puff_index_alpha, puff_index_alpha);
		    puff_beta_volume = this->GaussianPuffModel
		      ->ComputePuffOverlap(puff_index_beta, puff_index_beta);
		    puff_alpha_beta_volume = this->GaussianPuffModel
		      ->ComputePuffOverlap(puff_index_alpha, puff_index_beta);
		    
		    if (puff_index_alpha == puff_index_beta)
		      puff_alpha_beta_volume = 0.;
		    
		    if (puff_alpha_beta_volume == 0.)
		      overlap_tmp = 0.;
		    else
		      overlap_tmp = puff_alpha_beta_volume / (puff_alpha_volume * puff_beta_volume);
		    // If overlap between puff A and B is greater than 80% 
		    if (overlap_tmp > 0.8 *  (1. / puff_alpha_volume))
		      if ((1. / puff_alpha_volume) > (1. / puff_beta_volume))
			{		  
			  if (this->merge_source_id)
			    {
			      if (puff_alpha_id == puff_beta_id)
				PuffList_tmp.push_back(puff_index_beta);
			    }
			  else
			    PuffList_tmp.push_back(puff_index_beta);
			}
		  }
		//List contains puff to be merged with puff A
		PuffInteractionList(puff_index_alpha) = PuffList_tmp;
	      }
	    
	    int j_tmp, i_tmp;
	    //Combine overlaping puff
	    for (puff_index_alpha = 0; puff_index_alpha < Npuff; puff_index_alpha++)
	      {
		PuffList_tmp = PuffInteractionList(puff_index_alpha);
		Npuff_combine = PuffList_tmp.size();
		for (j_tmp = 0; j_tmp < Npuff_combine; j_tmp++)
		  {
		    combine_puff = 0;
		    puff_index_beta = PuffList_tmp[j_tmp];
		    Npuff_erase = Puff_erased.size();
		    //Check if puff has already been merged
		    for (i_tmp=0; i_tmp<Npuff_erase; i_tmp++)
		      if (puff_index_beta == Puff_erased[i_tmp])
			combine_puff += 1;
		    if (combine_puff == 0)
		      {
			puff_alpha_volume = 1. / this->GaussianPuffModel
			  ->ComputePuffOverlap(puff_index_alpha, puff_index_alpha);
			puff_beta_volume = 1. / this->GaussianPuffModel
			  ->ComputePuffOverlap(puff_index_beta, puff_index_beta);
			this->GaussianPuffModel
			  ->CombineOverlappingPuff(puff_index_alpha, 
						   puff_index_beta, puff_alpha_volume,
						   puff_beta_volume);
			Puff_erased.push_back(puff_index_beta);
		      }
		  }
	      }
	    
	    // Erase merged puffs
	    Npuff_erase = Puff_erased.size();
	    sort(Puff_erased.begin(), Puff_erased.end());
	    for (i_tmp = 0; i_tmp < Npuff_erase; i_tmp ++)
	      this->GaussianPuffModel->ErasePuff(Puff_erased[Npuff_erase - 1 - i_tmp]);
	    Puff_erased.clear();
	  }
	Npuff = this->GaussianPuffModel->GetPuffNumber();
	
	for (puff_index = 0; puff_index < Npuff; puff_index++)
	  {
	    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
	    this->GaussianPuffModel->
	      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
			      puff_time);
	    // Conversion to longitude/latitude (degrees).
	    this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
	    
	    // Gets corresponding cell in eulerian grid.
	    this->GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);
	    UpdateMeteo(puff_index);
	    PuffCellList(puff_index,0) = index_z;
	    PuffCellList(puff_index,1) = index_y;
	    PuffCellList(puff_index,2) = index_x;
	  }

	// Chemistry for puffs.
	this->GaussianPuffModel->Chemistry();

	// Save the quantities of the species in the puff 
	// if (rank == 0)
	this->GaussianPuffModel->PuffOutputSaver();
	
	this->GaussianPuffModel->AddTime(this->delta_t_local);
      } // GaussianPuffModel iteration 
  }


  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::Forward_Euler()
  {
    int l, k, j, i, s, b;
    // Eulerian model.
    this->Model.Forward();

    this->Concentration.Copy(this->Model.GetConcentration());
    this->Concentration_aer.Copy(this->Model.GetConcentration_aer());

    if (this->option_process["collect_dry_flux"])      
      {
	this->DryDepositionFlux_aer.Copy(this->Model.GetDryDepositionFlux_aer());
	this->DryDepositionFluxNumber_aer.Copy(this->Model.GetDryDepositionFluxNumber_aer());
      }

    if (this->option_process["collect_wet_flux"])
      {
	this->WetDepositionFlux_aer.Copy(this->Model.GetWetDepositionFlux_aer());
	this->InCloudWetDepositionFlux_aer.Copy(this->Model.GetInCloudWetDepositionFlux_aer());
	
	this->WetDepositionFluxNumber_aer.Copy(this->Model.GetWetDepositionFluxNumber_aer());
	this->InCloudWetDepositionFluxNumber_aer.Copy(this->Model.GetInCloudWetDepositionFluxNumber_aer());
      }

    if (this->option_process["with_number_concentration"])
      this->NumberConcentration_aer.Copy(this->Model.GetNumberConcentration_aer());

    this->AddTime(this->delta_t_eulerian);
    this->step++;
  }

  //! Add gaussian model concentration.
  /*! Add concentrations in the gaussian model
    in case concentrations are saved on the domain
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::AddGaussianConcentrationsDomain()
  {
    int Npuff, puff_index;
    bool isday;

    this->GaussianPuffModel->SubtractTime(this->delta_t_local);
    Npuff = this->GaussianPuffModel->GetPuffNumber();
    for (puff_index = 0; puff_index < Npuff; puff_index++)
      {
	this->GasCompensation.SetZero();
        T puff_release_time = this->GaussianPuffModel->
          GetPuffReleaseTime(puff_index);
        if (this->GaussianPuffModel->GetCurrentTime() >= puff_release_time)
          {
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            this->GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
	    
	    //Get indices of the Eulerian cell
	    int ind_x, ind_y, ind_z;
	    this->GetCellIndices(lon_c,lat_c,z_c, ind_z, ind_y, ind_x);
	    
	    T sigma_x, sigma_y, sigma_z;
	    this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
						  sigma_y, sigma_z);
	    
	    isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->
			  GetCurrentDate());
	    // Adding puff concentration to the Concentration Data.
	    if (this->injection_integrated)
	      {
		this->delta_number_puff.SetZero();
		this->PuffIntegratedTransfer(puff_index, this->Concentration);
		this->PuffIntegratedTransfer_aer(puff_index, this->Concentration_aer, this->Concentration);
		if (this->option_process["with_number_concentration"])
		  this->PuffIntegratedTransfer_number(puff_index, this->Concentration_aer, 
						      this->NumberConcentration_aer);
	      }
	    else
	      {
		for (int s = 0; s < this->Ns; s++)
		  this->PuffTransfer(this->GaussianPuffModel->
				     GetPuffQuantity(puff_index, s),
				     sigma_z, s, z_c, lat_c, lon_c,
				     isday, this->Concentration);
		for (int s = 0; s < this->Ns_aer; s++)
		  for (int b = 0; b < this->Nbin_aer; b++)
		    this->PuffTransfer_aer(this->GaussianPuffModel->
					   GetPuffQuantity_aer(puff_index, s, b),
					   sigma_z,
					   s,  z_c, lat_c, lon_c,
					   isday, this->Concentration_aer, b);
		if (this->option_process["with_number_concentration"])
		  for (int b = 0; b < this->Nbin_aer; b++)
		    this->PuffTransfer_number(this->GaussianPuffModel->
					      GetPuffQuantity_number(puff_index, b), 
					      sigma_z, 
					      z_c, lat_c, lon_c,
					      isday, this->NumberConcentration_aer, b);
		
	      }
	  }
      }

    this->GaussianPuffModel->AddTime(this->delta_t_local);
  }


  //! Checks whether the model deals with number concentration
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  bool PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::HasNumberConcentration_aer()
  {
    if(this->option_process["with_number_concentration"])
      return true;
    else
      return false;
  }


  //! Updates the puff meteorological data.
  /*! Updates the puff meteorological data taking them from the center cell.
    \param puff_index index of the puff.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateMeteo(int puff_index)
  {

    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    this->GaussianPuffModel->
      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance, puff_time);

    // Conversion to longitude/latitude (degrees).
    this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Night or day.
    bool isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->GetCurrentDate());

    // Extracting meteorological data.
    bool rural;
    string stability;
    map<string, T> meteorological_data;
    bool option_similarity = this->GaussianPuffModel->WithSimilarity();

    ExtractMeteo(z_c, lat_c, lon_c, isday, option_similarity,
                 meteorological_data, stability, rural);

    if (!option_similarity && !this->evolutive_plume_rise)
      this->GaussianPuffModel->
        SetPuffMeteo(puff_index,
                     meteorological_data["temperature"],
                     meteorological_data["wind_angle"],
                     meteorological_data["wind"],
                     stability, lon_c, lat_c, isday, rural,
                     meteorological_data["boundaryheight"]);
    else
      this->GaussianPuffModel->
        SetPuffMeteo(puff_index,
                     meteorological_data["temperature"],
                     meteorological_data["wind_angle"],
                     meteorological_data["wind"],
                     meteorological_data["frictionmodule"],
                     meteorological_data
                     ["convectivevelocity"],
                     meteorological_data["boundaryheight"],
                     meteorological_data["lmo"],
                     meteorological_data["coriolis"],
                     lon_c, lat_c, isday, rural);

    if (this->option_chemistry)
      {
        this->GaussianPuffModel->SetPuffAdditionalMeteo
          (puff_index,
           meteorological_data["attenuation"],
           meteorological_data["pressure"],
           meteorological_data["specific_humidity"]);

        this->GaussianPuffModel->SetPuffAdditionalMeteo_aer
          (puff_index,
           meteorological_data["liquid_water_content"]);
      }

    // Data depending on the puff species.
    Array<T, 1> deposition_velocity;
    deposition_velocity.resize(this->Ns);
    deposition_velocity = 0.;

    Array<T, 1> scavenging_coefficient;
    scavenging_coefficient.resize(this->Ns);
    scavenging_coefficient = 0.;

    Array<T, 1> photolysis_rate;
    photolysis_rate.resize(this->Nr_photolysis);

    Array<T, 1> background_concentration;
    background_concentration.resize(this->Ns);

    Array<T, 1> previous_background_concentration;
    previous_background_concentration.resize(this->Ns);

    Array<T, 2> background_concentration_aer;
    background_concentration_aer.resize(this->Ns_aer, this->Nbin_aer);

    Array<T, 2> previous_background_concentration_aer;
    previous_background_concentration_aer.resize(this->Ns_aer, this->Nbin_aer);

    Array<T, 1> background_concentration_number;
    background_concentration_number.resize(this->Nbin_aer);

    Array<T, 1> previous_background_concentration_number;
    previous_background_concentration_number.resize(this->Nbin_aer);

    Array<T, 1> deposition_velocity_aer;
    deposition_velocity_aer.resize(this->Ns_aer);
    deposition_velocity_aer = 0.;

    Array<T, 1> scavenging_coefficient_aer;
    scavenging_coefficient_aer.resize(this->Ns_aer);
    scavenging_coefficient_aer = 0.;

    ExtractSpeciesData(z_c, lat_c, lon_c,
                       deposition_velocity,
                       deposition_velocity_aer,
                       scavenging_coefficient,
                       scavenging_coefficient_aer,
                       photolysis_rate,
                       background_concentration,
                       background_concentration_number,
                       background_concentration_aer);

    if (this->option_deposition)
      this->GaussianPuffModel
        ->SetPuffDepositionVelocity(puff_index,
                                    deposition_velocity);

    if (this->option_scavenging)
      this->GaussianPuffModel
        ->SetPuffScavengingCoefficient(puff_index,
                                       scavenging_coefficient);

    if (this->option_chemistry)
      {
        this->GaussianPuffModel
          ->SetPuffPhotolysisRate(puff_index,
                                  photolysis_rate);
        this->GaussianPuffModel->SetPuffBackgroundConcentration
          (puff_index, background_concentration);

        this->GaussianPuffModel->SetPuffBackgroundConcentration_aer
          (puff_index, background_concentration_aer);

		this->GaussianPuffModel->SetPuffBackgroundConcentration_number
		  (puff_index, background_concentration_number);

      }

  }


  //! Extracts meteorological data.
  /*! Extracts meteorological data from 3D fields at a given point.
    \param D3_map pointer to map containing 3D meteorological data.
    \param index_z cell index along z.
    \param index_y cell index along y.
    \param index_x cell index along x.
    \param met_data map containing meteorological data at the given point.
    \param stability stability class at the given point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractMeteo(T height, T lat, T lon, bool isday,
                 bool option_similarity, map<string, T> &met_data,
                 string &stability, bool &rural)
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
      ::ExtractMeteo(height, lat, lon, isday,
                     option_similarity, met_data,
                     stability, rural);

    int index_x, index_y, index_z;
    this->GetCellIndices(lon, lat, height, index_z, index_y, index_x);

    Array<T, 1> Coord3D(3);
    Coord3D(0) = height;
    Coord3D(1) = lat;
    Coord3D(2) = lon;

    if (this->option_chemistry)
      {
        T liquid_water_content;
        if (!this->interpolated)
          liquid_water_content = LiquidWaterContent_i(index_z, index_y, index_x);
        else
          LinearInterpolationPoint(LiquidWaterContent_i, Coord3D, liquid_water_content);
        met_data["liquid_water_content"] = liquid_water_content;
      }
  }


  //! Extracts data that depend on meteorology and species (photolysis rates).
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractSpeciesData(T height, T lat, T lon,
                       Array<T, 1>& deposition_velocity,
                       Array<T, 1>& deposition_velocity_aer,
                       Array<T, 1>& scavenging_coefficient,
                       Array<T, 1>& scavenging_coefficient_aer,
                       Array<T, 1>& photolysis_rate,
                       Array<T, 1>& background_concentration,
                       Array<T, 1>& background_concentration_number,
                       Array<T, 2>& background_concentration_aer)
  {
    PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
      ::ExtractSpeciesData(height, lat, lon,
                           deposition_velocity,
                           scavenging_coefficient,
                           photolysis_rate,
                           background_concentration);


    int index_x, index_y, index_z;
    this->GetCellIndices(lon, lat, height, index_z, index_y, index_x);

    if (this->option_chemistry){
      if (this->option_process["with_number_concentration"])
	{
	  Array<T, 1> Coord4D(4);
	  Coord4D(0) = 0;
	  Coord4D(1) = height;
	  Coord4D(2) = lat;
	  Coord4D(3) = lon;
	  // Aerosol background number concentrations.
	  for (int b = 0; b < this->Nbin_aer; b++)
	      {
		Coord4D(0) = b;
		if (!this->interpolated)
		  background_concentration_number(b) =
		    this->Model.GetNumberConcentration_aer()(b, index_z, index_y, index_x);
		else
		  LinearInterpolationPoint(this->Model.GetNumberConcentration_aer(),
					   Coord4D, background_concentration_number(b));
		// To correct extrapolation that could give concentrations < 0.
		background_concentration_number(b) =
		  max(background_concentration_number(b), 0.); 
	      }
	}
      else
	{
	  for (int b = 0; b < this->Nbin_aer; b++)
	    {
	      background_concentration_number(b) = 0.0;
	    }
	}
    }
    if (this->option_chemistry)
      {
        Array<T, 1> Coord5D(5);
        Coord5D(0) = 0;
        Coord5D(1) = 0;
        Coord5D(2) = height;
        Coord5D(3) = lat;
        Coord5D(4) = lon;
        // Aerosol background concentrations.
        for (int s = 0; s < this->Ns_aer; s++)
          for (int b = 0; b < this->Nbin_aer; b++)
            {
              Coord5D(0) = s;
              Coord5D(1) = b;
              if (!this->interpolated)
                background_concentration_aer(s, b) =
                  this->Model.GetConcentration_aer()(s, b, index_z, index_y, index_x);
              else
                LinearInterpolationPoint(this->Model.GetConcentration_aer(),
                                         Coord5D, background_concentration_aer(s, b));
              // To correct extrapolation that could give concentrations < 0.
              background_concentration_aer(s, b) =
                max(background_concentration_aer(s, b), 0.);
            }
      }
  }

  /*! Computes the Gaussian chemistry and the subsequent feedback to
    background cells concentrations (in case there is feedback with
    chemistry).*/
  /*!
    \param PuffCellList list of coordinates of the cells containing puffs.
    \param PuffConcentration Concentration perturbation to be added
    to each cell after puff chemistry.
  */ 
  // template<class T, class ClassEulerianModel, class ClassLocalModel>
  // void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  // ::ComputeChemistry(Array<int, 2> PuffCellList,
  // 		     Array<T, 4>& PuffConcentration,
  //                    Array<T, 5>& PuffConcentration_aer)
  // {
  //   int l, k, i, s, puff_index, b;
  //   int index_x = 0;
  //   int index_y = 0;
  //   int index_z = 0;
  //   int Npuff = this->GaussianPuffModel->GetPuffNumber();
  //   // Loop on cells to take them into account for background chemistry.
  //   bool is_in_list;
  //   list<int> pufflist_tmp;
  //   vector<list<int> > PuffList;
  //   vector<T> PuffCellVolume;
  //   Array<int, 1> Coord(3);
  //   vector<Array<int, 1> > PuffCellCoordinates;
  //   Array<vector<T>, 1 > PuffCellConcentration;
  //   Array<vector<T>, 2 > PuffCellConcentration_aer;
  //   PuffCellConcentration.resize(this->Ns);
  //   PuffCellConcentration_aer.resize(this->Ns_aer, this->Nbin_aer);
  //   int Ncell = 0;
  //   for (l = 0; l < this->Nz; l++)
  //     for (k = 0; k < this->Ny; k++)
  // 	for (i = 0; i < this->Nx; i++)
  // 	  {
  // 	    is_in_list = false;
  // 	    pufflist_tmp.clear();
  // 	    // Is there at least one puff in the cell?
  // 	    for (puff_index = 0; puff_index < Npuff; puff_index++)
  // 	      if (PuffCellList(puff_index,0) == l
  // 		  && PuffCellList(puff_index,1) == k
  // 		  && PuffCellList(puff_index,2) == i)
  // 		{
  // 		  index_z = l;
  // 		  index_y = k;
  // 		  index_x = i;
  // 		  is_in_list = true;
  // 		  pufflist_tmp.push_back(puff_index);
  // 		}
  // 	    if (is_in_list)
  // 	      {
  // 		Ncell++;
  // 		// Computing cell volume.
  // 		T cell_width_z, cell_width_y, cell_width_x, cell_volume;
  // 		this->ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
  // 				 cell_width_y, cell_width_x, cell_volume);
  // 		PuffCellVolume.push_back(cell_volume);
  // 		Coord(0) = index_z;
  // 		Coord(1) = index_y;
  // 		Coord(2) = index_x;
  // 		PuffCellCoordinates.push_back(Coord);
  // 		for (s = 0; s < this->Ns; s++)
  // 		  PuffCellConcentration(s).push_back
  // 		    (this->Model.GetConcentration()(s, index_z, index_y, index_x));
  //               for (s = 0; s < this->Ns_aer; s++)
  //                 for (b = 0; b < this->Nbin_aer; b++)
  //                   PuffCellConcentration_aer(s, b).push_back
  //                     (this->Model.GetConcentration_aer()(s, b, index_z, index_y, index_x));

  // 		PuffList.push_back(pufflist_tmp);
  // 	      }
  // 	  }
  //   this->GaussianPuffModel->Chemistry(PuffList, PuffCellVolume,
  //                                      PuffCellConcentration,
  //                                      PuffCellConcentration_aer);
  //   // Adding the remaining concentrations in a buffer array.
  //   for (i = 0; i < Ncell; i++)
  //     {
  // 	index_z = PuffCellCoordinates[i](0);
  // 	index_y = PuffCellCoordinates[i](1);
  // 	index_x = PuffCellCoordinates[i](2);
	
  // 	for (s = 0; s < this->Ns; s++)
  // 	  {
  // 	    PuffConcentration(s, index_z, index_y, index_x)
  // 	      += PuffCellConcentration(s)[i];
  // 	  }
  //       for (s = 0; s < this->Ns_aer; s++)
  //         for (b = 0; b < this->Nbin_aer; b++)
  // 	    PuffConcentration_aer(s, b, index_z, index_y, index_x)
  // 	      += PuffCellConcentration_aer(s, b)[i];
  //     }
  // }

  //! Transfers a given puff to the Eulerian model.
  /*! Transfers a given puff to the Eulerian model on one vertical column.
    \param quantity mass of the puff.
    \param sigma_z vertical standard deviation of the puff.
    \param index_s index of the puff species.
    \param z_c puff center coordinate along z.
    \param lat_c puff center coordinate along y.
    \param lon_c puff center coordinate along x.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::PuffTransfer_aer(T quantity, T sigmaz, int s, T z, T lat, T lon, bool isday,
                     Data<T, 5>& Concentration_out_aer, int b)
  {
    int ibz, iby, ibx, itz, ity, itx;
    // Computes puff vertical extent.
    T puff_bottom = max(z - (this->coefficient_z * sigmaz) / 2, 0.);
    this->GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    T puff_top = z + (this->coefficient_z * sigmaz) / 2;
    T boundary_height = this->BoundaryHeight_i(iby, ibx);

    // Takes into account the inversion height.
    if (isday && z < boundary_height)
      puff_top = min(puff_top, boundary_height);
    if (isday && z > boundary_height)
      puff_bottom = max(puff_bottom, boundary_height);
    this->GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    this->GetCellIndices(lon, lat, puff_top, itz, ity, itx);

    // Number of cells vertically covered by the puff.
    int Ncell = (itz - ibz) + 1;

    // Computing vertical extent of cells where puff will be transfered.
    T vertical_extent = 0.;
    for (int k = 0; k < Ncell; k++)
      vertical_extent += (this->GridZ3D_interf(ibz + k + 1)
                          - this->GridZ3D_interf(ibz + k));

    // Concentration to be added to each cell.
    T cell_width_z, cell_width_y, cell_width_x, cell_volume;
    this->ComputeCellWidth(ibz, iby, ibx, cell_width_z, cell_width_y,
                           cell_width_x, cell_volume);
    T concentration = quantity /
      (vertical_extent * cell_width_y * cell_width_x);

    //    Adding puff aerosol concentration to the Eulerian concentration.
    for (int k = 0; k < Ncell; k++)
      Concentration_out_aer(s, b, ibz + k, iby, ibx)
        = max(Concentration_out_aer(s, b, ibz + k, iby, ibx) + concentration, 0.);

  }

  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::PuffTransfer_number(T quantity, T sigmaz, T z, T lat, T lon, bool isday,
                        Data<T, 4>& Concentration_out_number, int b)
  {
    int ibz, iby, ibx, itz, ity, itx;
    // Computes puff vertical extent.
    T puff_bottom = max(z - (this->coefficient_z * sigmaz) / 2, 0.);
    this->GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    T puff_top = z + (this->coefficient_z * sigmaz) / 2;
    T boundary_height = this->BoundaryHeight_i(iby, ibx);

    // Takes into account the inversion height.
    if (isday && z < boundary_height)
      puff_top = min(puff_top, boundary_height);
    if (isday && z > boundary_height)
      puff_bottom = max(puff_bottom, boundary_height);
    this->GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    this->GetCellIndices(lon, lat, puff_top, itz, ity, itx);

    // Number of cells vertically covered by the puff.
    int Ncell = (itz - ibz) + 1;

    // Computing vertical extent of cells where puff will be transfered.
    T vertical_extent = 0.;
    for (int k = 0; k < Ncell; k++)
      vertical_extent += (this->GridZ3D_interf(ibz + k + 1)
			  - this->GridZ3D_interf(ibz + k));
    
    // Concentration to be added to each cell.
    T cell_width_z, cell_width_y, cell_width_x, cell_volume;
    this->ComputeCellWidth(ibz, iby, ibx, cell_width_z, cell_width_y,
                           cell_width_x, cell_volume);
    T concentration = quantity /
      (vertical_extent * cell_width_y * cell_width_x);

    //    Adding puff aerosol concentration to the Eulerian concentration.
    for (int k = 0; k < Ncell; k++)
      Concentration_out_number(b, ibz + k, iby, ibx)
	= max(Concentration_out_number(b, ibz + k, iby, ibx) + concentration, 0.);
    
  }



  //! Transfers a given puff to the Eulerian model.
  /*! Transfers a puff to the Eulerian model with the integrated method.
    \param puff_index index of the puff.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::PuffIntegratedTransfer_aer(int puff_index, 
                               Data<T, 5>& Concentration_out_aer, 
                               Data<T, 4>& Concentration_out)
  {
    T sigma_x, sigma_y, sigma_z;
    this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                          sigma_y, sigma_z);

    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    this->GaussianPuffModel->
      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                      puff_time);
    // Conversion to longitude/latitude (degrees).
    this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Computes puff vertical extent.
    int ibz, iby, ibx, itz, ity, itx;
    T puff_bottom = max(z_c - (this->coefficient_z * sigma_z) / 2, 0.);
    this->GetCellIndices(lon_c, lat_c, puff_bottom, ibz, iby, ibx);
    T puff_top = z_c + (this->coefficient_z * sigma_z) / 2;
    this->GetCellIndices(lon_c, lat_c, puff_top, itz, ity, itx);

    // Computes puff horizontal extent.
    int ilyz, ily, ilyx, iryz, iry, iryx;
    T puff_left = y_c - (this->coefficient_y * sigma_y) / 2;
    T lon_left, lat_left;
    this->CartesianToLatLon(x_c, puff_left, lon_left, lat_left);
    this->GetCellIndices(lon_c, lat_left, z_c, ilyz, ily, ilyx);
    T puff_right = y_c + (this->coefficient_y * sigma_y) / 2;
    T lon_right, lat_right;
    this->CartesianToLatLon(x_c, puff_right, lon_right, lat_right);
    this->GetCellIndices(lon_c, lat_right, z_c, iryz, iry, iryx);

    // Computes puff horizontal extent.
    int ilxz, ilxy, ilx, irxz, irxy, irx;
    puff_left = x_c - (this->coefficient_y * sigma_x) / 2;
    this->CartesianToLatLon(puff_left, y_c, lon_left, lat_left);
    this->GetCellIndices(lon_left, lat_c, z_c, ilxz, ilxy, ilx);
    puff_right = x_c + (this->coefficient_y * sigma_x) / 2;
    this->CartesianToLatLon(puff_right, y_c, lon_right, lat_right);
    this->GetCellIndices(lon_right, lat_c, z_c, irxz, irxy, irx);

    // Number of cells covered by the puff.
    int Nz = (itz - ibz) + 1;
    int Ny = (iry - ily) + 1;
    int Nx = (irx - ilx) + 1;
    int Ns_aer = this->Ns_aer;
    int Nbin_aer = this->Nbin_aer;

    T concentration;
    T gas_compensation;
    int b,s, b_bis;

    //Total gas/particle mass conservation ...
    ConfigStream config_species(this->GetSpeciesFile());
    map<string, string> gas_aer_interaction;
    map<string, string>::iterator iter_interact_tmp;
    string species_interact;
    Array<int, 1> species_index_aerosol_interact_tmp;
    int species_comp;

    species_index_aerosol_interact_tmp.resize(this->Ns_aer);
    config_species.SetSection("[gas_species_aerosol_interact]");
    while (!config_species.IsEmpty())
      {
    	species_interact = config_species.GetElement();
    	gas_aer_interaction[species_interact] = config_species.GetElement();
      }
    
    for (int i = 0; i < this->Ns_aer; i++)
      {
    	iter_interact_tmp = gas_aer_interaction.find(this->species_list_aer[i]);
    	if (iter_interact_tmp != gas_aer_interaction.end())
    	  {
    	    int gas_index = this->GetSpeciesIndex(iter_interact_tmp->second);
    	    species_index_aerosol_interact_tmp(i) = gas_index;
    	  }
    	else
    	  species_index_aerosol_interact_tmp(i) = -1;
      }
    gas_aer_interaction.clear();
    
    // If the puff is in only one cell.
    if (Nz * Ny * Nx == 1)
      {
        // Puff center cell.
        int icz, icy, icx;
        this->GetCellIndices(lon_c, lat_c, z_c, icz, icy, icx);
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        this->ComputeCellWidth(icz, icy, icx,
                               cell_width_z, cell_width_y,
                               cell_width_x, cell_volume);

        // Concentration for all species.
        for (int s = 0; s < Ns_aer; s++)
          for (int b = 0; b < Nbin_aer; b++)
            {
              concentration = this->GaussianPuffModel
                ->GetPuffQuantity_aer(puff_index, s, b)/ cell_volume;

              Concentration_out_aer(s, b, icz, icy, icx) =
                max(Concentration_out_aer(s, b, icz, icy, icx) + concentration,
                    0.);
            }	
      }
    else
      {
        Array<T, 2> total_mass(Ns_aer, Nbin_aer);
        total_mass = 0.;
        Array<T, 5> PuffConcentration_aer(Ns_aer, Nbin_aer, Nz, Ny, Nx);
        PuffConcentration_aer = 0.;
        int s, k, j, i, index_z, index_y, index_x, b;
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        T concentration_erf;
        Array<T, 2> quantity(Ns_aer, Nbin_aer);
        T x_c, y_c, z_c, lat_c, lon_c;
        
        // Computing concentration to be added to each cell.
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                index_z = ibz + k;
                index_y = ily + j;
                index_x = ilx + i;
                this->ComputeCellWidth(index_z, index_y, index_x,
                                       cell_width_z, cell_width_y,
                                       cell_width_x, cell_volume);

                // Cell center.
                lon_c = this->x_min + index_x * this->Delta_x;
                lat_c = this->y_min + index_y * this->Delta_y;
                z_c = this->GridZ3D(index_z);
                this->LatLonToCartesian(lon_c, lat_c, x_c, y_c);

                // Computing concentrations.
                for (s = 0; s < Ns_aer; s++)
                  for (b = 0; b < Nbin_aer; b++)
                    {
                      quantity(s, b) = this->GaussianPuffModel
                        ->GetPuffQuantity_aer(puff_index, s, b);
                      if (quantity(s, b) != 0.)
                        {
                          concentration_erf =
                            this->GaussianPuffModel->
                            ComputePuffIntegral_aer(puff_index, s, b, 
                                                    x_c, y_c, z_c,
                                                    cell_width_x, cell_width_y,
                                                    cell_width_z);
                          PuffConcentration_aer(s, b, k, j, i) += concentration_erf;
                          
                          total_mass(s, b) += concentration_erf * cell_volume;
                        }
                    }
              }
        
        // Adding concentration to each cell (with mass conservation).
        for (s = 0; s < Ns_aer; s++)
          for (k = 0; k < Nz; k++)
            for (j = 0; j < Ny; j++)
              for (i = 0; i < Nx; i++)
                {
                  index_z = ibz + k;
                  index_y = ily + j;
                  index_x = ilx + i;

		  for (b = 0; b < Nbin_aer; b++)
		    {
		      if (quantity(s, b) != 0. && total_mass(s, b) != 0.)
			{
			  concentration = (quantity(s, b) / total_mass(s, b))
			    * PuffConcentration_aer(s, b, k, j, i);
			}
		      else
			concentration = 0.;

                      Concentration_out_aer(s, b, index_z, index_y, index_x)
                        = max(Concentration_out_aer(s, b, 
                                                    index_z, index_y, index_x)
                              + concentration, 0.);
                    }
                }
        
      }
  }

  //! Transfers a given puff to the Eulerian model.
  /*! Transfers a puff to the Eulerian model with the integrated method.
    \param puff_index index of the puff.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridAerosol<T, ClassEulerianModel, ClassLocalModel>
  ::PuffIntegratedTransfer_number(int puff_index, 
                                  Data<T, 5> Concentration_out_aer, 
				  Data<T, 4>& Concentration_out_number)
  {
    T duree;
    time_t debut, fin;
    time_t inter;
    debut = time(NULL);
    T sigma_x, sigma_y, sigma_z;
    this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
				    sigma_y, sigma_z);

    Array<T, 2> concentration_aer_tmp(this->Ns_aer, this->Nbin_aer);
    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    this->GaussianPuffModel->
      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
		      puff_time);
    // Conversion to longitude/latitude (degrees).
    this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);
    
    // Computes puff vertical extent.
    int ibz, iby, ibx, itz, ity, itx;
    T puff_bottom = max(z_c - (this->coefficient_z * sigma_z) / 2, 0.);
    this->GetCellIndices(lon_c, lat_c, puff_bottom, ibz, iby, ibx);
    T puff_top = z_c + (this->coefficient_z * sigma_z) / 2;
    this->GetCellIndices(lon_c, lat_c, puff_top, itz, ity, itx);

    // Computes puff horizontal extent.
    int ilyz, ily, ilyx, iryz, iry, iryx;
    T puff_left = y_c - (this->coefficient_y * sigma_y) / 2;
    T lon_left, lat_left;
    this->CartesianToLatLon(x_c, puff_left, lon_left, lat_left);
    this->GetCellIndices(lon_c, lat_left, z_c, ilyz, ily, ilyx);
    T puff_right = y_c + (this->coefficient_y * sigma_y) / 2;
    T lon_right, lat_right;
    this->CartesianToLatLon(x_c, puff_right, lon_right, lat_right);
    this->GetCellIndices(lon_c, lat_right, z_c, iryz, iry, iryx);
    
    // Computes puff horizontal extent.
    int ilxz, ilxy, ilx, irxz, irxy, irx;
    puff_left = x_c - (this->coefficient_y * sigma_x) / 2;
    this->CartesianToLatLon(puff_left, y_c, lon_left, lat_left);
    this->GetCellIndices(lon_left, lat_c, z_c, ilxz, ilxy, ilx);
    puff_right = x_c + (this->coefficient_y * sigma_x) / 2;
    this->CartesianToLatLon(puff_right, y_c, lon_right, lat_right);
    this->GetCellIndices(lon_right, lat_c, z_c, irxz, irxy, irx);

    // Number of cells covered by the puff.
    int Nz = (itz - ibz) + 1;
    int Ny = (iry - ily) + 1;
    int Nx = (irx - ilx) + 1;
    int Ns_aer = this->Ns_aer;
    int Nbin_aer = this->Nbin_aer;

    int b,s;
    T puff_number_new;
    int b_number;
    T total_number;
    T total_number_positive;
    T delta_number;
    T mass_ratio;
    T concentration_tmp;
    T delta_concentration_number;

    Array<T, 1> MassDensity_aer(Nbin_aer);
    MassDensity_aer = 0.;
    this->GaussianPuffModel->GetMassDensity_aer(MassDensity_aer);
    
    // Array<T, 1> delta_concentration_list(Nbin_aer);
    // delta_concentration_list = 0.;
    
    // If the puff is in only one cell.
    if (Nz * Ny * Nx == 1)
      {
	// Puff center cell.
	int icz, icy, icx;
	this->GetCellIndices(lon_c, lat_c, z_c, icz, icy, icx);
	T cell_width_z, cell_width_y, cell_width_x, cell_volume;
	this->ComputeCellWidth(icz, icy, icx,
			 cell_width_z, cell_width_y,
			 cell_width_x, cell_volume);
	
	total_number =0.;
	total_number_positive = 0.;
        
	for (b = 0; b < Nbin_aer; b++)
	  {
	    delta_concentration_number = this->GaussianPuffModel
	      ->GetPuffQuantity_number(puff_index, b) / cell_volume;
            
	    if (Concentration_out_number(b, icz, icy, icx)
		+ delta_concentration_number < 0.)
	      {
		for (s=0; s<Ns_aer-1; s++)
		  concentration_aer_tmp(s, b) = 
                    Concentration_out_aer(s, b, icz, icy, icx);  
		puff_number_new = 
		  this->GaussianPuffModel->
                  ComputeNumberPuff(b, concentration_aer_tmp);
		
		Concentration_out_number(b, icz, icy, icx) = puff_number_new;
	      }
	    else
	      {
		Concentration_out_number(b, icz, icy, icx) =
		  max(Concentration_out_number(b, icz, icy, icx) + 
                      delta_concentration_number, 0.);
	      }
	  } 
      }
    else
      {
	Array<T, 1> total_mass(Nbin_aer);
	total_mass = 0.;
	Array<T, 4> PuffConcentration_number(Nbin_aer, Nz, Ny, Nx);
	PuffConcentration_number = 0.;
	int s, k, j, i, index_z, index_y, index_x, b;
	T cell_width_z, cell_width_y, cell_width_x, cell_volume;

	Array<T, 1> quantity(Nbin_aer);
	T concentration_erf;
	T concentration;
	T x_c, y_c, z_c, lat_c, lon_c;

	// Computing concentration to be added to each cell.
	for (k = 0; k < Nz; k++)
	  for (j = 0; j < Ny; j++)
	    for (i = 0; i < Nx; i++)
	      {
		index_z = ibz + k;
		index_y = ily + j;
		index_x = ilx + i;
		this->ComputeCellWidth(index_z, index_y, index_x,
                                       cell_width_z, cell_width_y,
                                       cell_width_x, cell_volume);

		// Cell center.
		lon_c = this->x_min + index_x * this->Delta_x;
		lat_c = this->y_min + index_y * this->Delta_y;
		z_c = this->GridZ3D(index_z);
		this->LatLonToCartesian(lon_c, lat_c, x_c, y_c);

		// Computing concentrations.
		for (b = 0; b < Nbin_aer; b++)
		  {
		    quantity(b) = this->GaussianPuffModel
		      ->GetPuffQuantity_number(puff_index, b);
		    if (quantity(b) != 0.)
		      {
			concentration_erf =
			  this->GaussianPuffModel->
			  ComputePuffIntegral_number(puff_index, b, 
                                                     x_c, y_c, z_c,
                                                     cell_width_x, cell_width_y,
                                                     cell_width_z);
			PuffConcentration_number(b, k, j, i) += 
                          concentration_erf;

			total_mass(b) += concentration_erf * cell_volume;
		      }
		  }
              }

	for (k = 0; k < Nz; k++)
	  for (j = 0; j < Ny; j++)
	    for (i = 0; i < Nx; i++)
	      {     
		index_z = ibz + k;
		index_y = ily + j;
		index_x = ilx + i;
		
		total_number = 0.;
		total_number_positive = 0.;
		
		for (b = 0; b < Nbin_aer; b++)
		  if (quantity(b) != 0. && total_mass(b) != 0.)
		    {
		      delta_concentration_number = (quantity(b) / total_mass(b))
			* PuffConcentration_number(b, k, j, i);
		      if (Concentration_out_number(b, index_z, index_y, index_x)
                          + delta_concentration_number < 0.)
		      	{
		      	  for (s=0; s<Ns_aer; s++)
		      	    concentration_aer_tmp(s, b) = 
		      	      Concentration_out_aer(s, b,
                                                    index_z, index_y, index_x);
		      	  puff_number_new = 
		      	    this->GaussianPuffModel->
                            ComputeNumberPuff(b, concentration_aer_tmp);
		      	  Concentration_out_number(b, index_z, index_y, index_x) = puff_number_new;
		      	}
		      else
		      	{
			  Concentration_out_number(b, index_z, index_y, index_x)
                            = max(Concentration_out_number(b, index_z, index_y,
                                                           index_x) + 
                                  delta_concentration_number, 0.);
			}
		    }
	      }
      }
  } 
  
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDAEROSOL_CXX
#endif
