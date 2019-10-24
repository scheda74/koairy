// Copyright (C) 2016, CEREA - ENPC - EDF R&D
// Author(s): Youngseob Kim
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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


#ifndef POLYPHEMUS_FILE_MODELS_STREETINGRIDTRANSPORT_CXX


#include "StreetInGridTransport.hxx"


namespace Polyphemus
{

  template<class T, class ClassEulerianModel, class ClassLocalModel>
  const T StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>::pi = acos(-1);

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*! Builds the model and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::StreetInGridTransport(string config_file): BaseModel<T>(config_file),
                                               Model(config_file)
  {
    /*** Pointers to 3D and 2D data ***/

    this->D3_map["MeridionalWind_i"] = &MeridionalWind_i;
    this->D3_map["ZonalWind_i"] = &ZonalWind_i;
    this->D3_map["Temperature_i"] = &Temperature_i;

    this->D2_map["FrictionModule_i"] = &FrictionModule_i;
    this->D2_map["BoundaryHeight_i"] = &BoundaryHeight_i;
    this->D2_map["LMO_i"] = &LMO_i;

    this->D2_map["StreetConcentration"] = &StreetConcentration;
  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::~StreetInGridTransport()
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
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    this->config.SetSection("[domain]");
    this->config.PeekValue("Cartesian", option_cartesian);

    this->config.SetSection("[options]");
    string scavenging_model;
    this->config.PeekValue("With_deposition", option_deposition);
    this->config.PeekValue("Scavenging_model", scavenging_model);
    option_scavenging = (scavenging_model != "none");

    /*** Street-network configuration ***/

    this->config.SetSection("[street]");
    this->config.PeekValue("Street_configuration", street_config);
    this->config.PeekValue("With_interpolation", interpolated);

    /*** LUC ***/
    this->config.SetSection("[data]");
    string data_file;
    this->config.PeekValue("Data_description", data_file);
    ConfigStream data(data_file);
    data.SetSection("[ground]");
    data.PeekValue("LUC_file", LUC_file);
    if (!exists(LUC_file))
      throw "Unable to open land use cover file \"" + LUC_file + "\".";
    Nc = int(file_size(LUC_file)) / sizeof(float) / (this->Ny * this->Nx);
    data.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                   urban_index);

  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::CheckConfiguration()
  {
    BaseModel<T>::CheckConfiguration();
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    ConfigStream data_description_stream(data_description_file);
    this->input_files["meteo"].Read(data_description_file, "meteo");
  }

  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Allocate()
  {
    BaseModel<T>::Allocate();

    /*** Winds ***/

    ZonalWind_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D_interf);
    FileZonalWind_i.Resize(this->GridZ3D, this->GridY3D,
			   this->GridX3D_interf);
    FileZonalWind_f.Resize(this->GridZ3D, this->GridY3D,
			   this->GridX3D_interf);

    MeridionalWind_i.Resize(this->GridZ3D, this->GridY3D_interf,
			    this->GridX3D);
    FileMeridionalWind_i.Resize(this->GridZ3D, this->GridY3D_interf,
				this->GridX3D);
    FileMeridionalWind_f.Resize(this->GridZ3D, this->GridY3D_interf,
				this->GridX3D);

    /*** Temperature ***/

    Temperature_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_i.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileTemperature_f.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);

    /*** FrictionModule ***/

    FrictionModule_i.Resize( this->GridY2D, this->GridX2D);
    FileFrictionModule_i.Resize( this->GridY2D, this->GridX2D);
    FileFrictionModule_f.Resize( this->GridY2D, this->GridX2D);

    /*** BoundaryHeight ***/

    BoundaryHeight_i.Resize( this->GridY2D, this->GridX2D);
    FileBoundaryHeight_i.Resize( this->GridY2D, this->GridX2D);
    FileBoundaryHeight_f.Resize( this->GridY2D, this->GridX2D);

    /*** LMO ***/

    LMO_i.Resize( this->GridY2D, this->GridX2D);
    FileLMO_i.Resize( this->GridY2D, this->GridX2D);
    FileLMO_f.Resize( this->GridY2D, this->GridX2D);

    /*** LUC ***/

    RegularGrid<T> GridC(Nc);
    LUC.Resize(GridC, this->GridY3D, this->GridX3D);

    /* */
    ConcentrationOverCanopy.Resize(this->GridS3D, this->GridY3D, this->GridX3D);
    ConcentrationOverCanopy.SetZero();
  }


  //! Initialization.
  /*! Initializes the Eulerian model and the meteorological conditions. Reads
    the sources and creates the corresponding Street-network model instances.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    BaseModel<T>::Init();

    /*** Initializations ***/
    Model.Init();
    this->Concentration.Copy(Model.GetConcentration());

    // Extracts simulation information from Eulerian model.
    delta_t_eulerian = Model.GetDelta_t();
    if (this->Delta_t != delta_t_eulerian
	|| this->Date_min != Model.GetDate_min())
      throw string("Time steps and date for Eulerian model") +
	" and street-in-grid model must be equal.";

    if (option_deposition)
      {
        if (Model.HasField("DepositionVelocity"))
          DepositionVelocity_i.Copy(Model.D3("DepositionVelocity_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"DepositionVelocity\".";
      }

    if (option_scavenging)
      {
        if (Model.HasField("ScavengingCoefficient"))
          ScavengingCoefficient_i.Copy(Model.D4("ScavengingCoefficient_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"ScavengingCoefficient\".";
      }

    // Is there a field "Altitude" to manage vertical levels?
    altitude_field = Model.HasField("Altitude");

    // Read LUC data.
    FormatBinary<float> InputLUC;
    InputLUC.Read(LUC_file, LUC);

    /***  Meteorological data ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");

    ConfigStream data_description_stream(data_description_file);
    this->input_files["meteo"].Read(data_description_file, "meteo");

    /*** Winds ***/
    this->InitData("meteo", "ZonalWind", FileZonalWind_i, FileZonalWind_f,
		   this->current_date, ZonalWind_i);

    this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
		   FileMeridionalWind_f, this->current_date,
		   MeridionalWind_i);

    /*** Temperature ***/
    this->InitData("meteo", "Temperature", FileTemperature_i,
		   FileTemperature_f, this->current_date, Temperature_i);

    /*** FrictionModule ***/
    this->InitData("meteo", "FrictionModule",
        	   FileFrictionModule_i, FileFrictionModule_f,
        	   this->current_date, FrictionModule_i);

    /*** BoundaryHeight ***/
    this->InitData("meteo", "BoundaryHeight",
        	   FileBoundaryHeight_i, FileBoundaryHeight_f,
        	   this->current_date, BoundaryHeight_i);

    /*** LMO ***/
    this->InitData("meteo", "LMO", FileLMO_i,
        	   FileLMO_f, this->current_date, LMO_i);

     if (option_deposition)
       DepositionVelocity_i.Copy(Model.D3("DepositionVelocity_i"));

     if (option_scavenging)
       ScavengingCoefficient_i.Copy(Model.D4("ScavengingCoefficient_i"));

    /*** Street-Network model ***/

    // Street-Network Model creation.
    StreetNetworkModel = new ClassLocalModel(street_config);
  
    // Read the configuration for Street-Network model..
    StreetNetworkModel->ReadConfiguration();

    delta_t_local = StreetNetworkModel->GetDelta_t();
    if (delta_t_local != delta_t_eulerian)
      {
        cout << "=== Warning: the time step in the street model should be equal to that in the Eulerian model." << endl; 
        cout << "=== Warning: the time step in the street model is set to " << delta_t_eulerian << endl;

        delta_t_local = delta_t_eulerian;
        StreetNetworkModel->SetDelta_t(delta_t_eulerian);
      }
    
    if (StreetNetworkModel->GetDateMin() != this->Date_min)
      {
        cout << "=== Warning: the starting date in the street model should be same to that in the Polair3D model." << endl;       
        StreetNetworkModel->SetDateMin(this->Date_min);
      }

    // Street-Network model initialization.
    StreetNetworkModel->Init();

    StreetConcentration.Copy(StreetNetworkModel->GetStreetConcentration()); // YK

    // Initialize the input data 
    StreetNetworkModel->InitData();

    // Compute the urban fraction in the grid cells.
    ComputeUrbanFraction();

    for (int s = 0; s < this->Ns; s++)
      for (int j = 0; j < this->Ny; j++)
        for (int i = 0; i < this->Nx; i++)
          ConcentrationOverCanopy(s, j, i) = 
            this->Concentration(s, 0, j, i); 

  }

  //! Compute the mass flux between the street network and the atmosphere.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeMassFlux()
  {
    mass_flux_grid.resize(this->Ns, this->Ny, this->Nx);
    mass_flux_grid = 0.0;

    Nstreet = StreetNetworkModel->GetNStreet();
    for (int street_index = 0; street_index < Nstreet; street_index++)
      {
        // Gets the coordinate of the street.
        T longitude, latitude;
        StreetNetworkModel->GetStreetCoordinate(street_index, longitude, latitude); 

        // Gets corresponding cell in Eulerian grid.
        int index_x, index_y, index_z;
        GetCellIndices(longitude, latitude, 0.0, index_z, index_y, index_x);

        for (int s = 0; s < this->Ns; ++s)
          {
            // Gets the mass flux. 
            T mass_flux = StreetNetworkModel->
              GetMassTransferBackground(street_index, s);
            // Adds the mass flux to the correspoding grid cell.
            mass_flux_grid(s, index_y, index_x) += mass_flux;
          }
      }
  }


  //! Background concentrations are corrected accounting for the urban fraction.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::AddMassFluxToBackgroundConcentration()
  {
    for (int j = 0; j < this->Ny; ++j)
      for (int i = 0; i < this->Nx; ++i)
        {        
          if (is_urban(j, i))
            {
              // Computes ratio of the urban volume and the volume of the grid cell.
              T cell_width_z, cell_width_y, cell_width_x, cell_volume;
              ComputeCellWidth(0, j, i, cell_width_z, cell_width_y,
                               cell_width_x, cell_volume);

              // Applies the ratio to compute the new concentrations for each species.
              for (int s = 0; s < this->Ns; ++s)
                {
                  T old_concentration = Model.GetConcentration()(s, 0, j, i);
                  T added_mass = mass_flux_grid(s, j, i) * this->Delta_t;
                  T new_concentration = (old_concentration * cell_volume + 
                                         added_mass) / cell_volume;
                  if (new_concentration < 0.0)
                    new_concentration = 0.0;
                  
                  Model.GetConcentration()(s, 0, j, i) = new_concentration;
                }
            }
        }
  }


  //! Compute urban fractions in the grid cells.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeUrbanFraction()
  {
     is_urban.resize(this->Ny, this->Nx);
     is_urban = false;

     urban_fraction.resize(this->Ny, this->Nx);
     urban_fraction = 0.0;

     weighted_height.resize(this->Ny, this->Nx);
     weighted_height = 0.0;

     sum_height_length.resize(this->Ny, this->Nx);
     sum_height_length = 0.0;

     sum_length.resize(this->Ny, this->Nx);
     sum_length = 0.0;

     street_number.resize(this->Ny, this->Nx);
     street_number = 0.0;

     sum_street_volume.resize(this->Ny, this->Nx);
     sum_street_volume = 0.0;

     // Checks if street canyons exist in the grid cell and
     // Gets the canyon heights.
     Nstreet = StreetNetworkModel->GetNStreet();
     for (int street_index = 0; street_index < Nstreet; street_index++)
       {
         T street_height = StreetNetworkModel->GetStreetHeight(street_index);
         T street_length = StreetNetworkModel->GetStreetLength(street_index);
         T street_volume = StreetNetworkModel->GetStreetVolume(street_index);

         int s_id =  StreetNetworkModel->GetStreetID(street_index);
         T longitude, latitude;
         StreetNetworkModel->GetStreetCoordinate(street_index, longitude, latitude); 

         // Gets corresponding cell in Eulerian grid.
         int index_x, index_y, index_z;
         GetCellIndices(longitude, latitude, 0.0, index_z, index_y, index_x);

         is_urban(index_y, index_x) = true; 
         sum_height_length(index_y, index_x) += (street_height * street_length);
         sum_length(index_y, index_x) += street_length; 
         street_number(index_y, index_x) += 1; 

         sum_street_volume(index_y, index_x) += street_volume;
       }

     // Gets the urban fraction from the LUC data for the identified grid cells.
     for (int j = 0; j < this->Ny; ++j)
       for (int i = 0; i < this->Nx; ++i)
         if (is_urban(j, i))
           {
             // Gets the urban fraction from the given binary LUC file.
             urban_fraction(j, i) = LUC(urban_index, j, i);

             // Computes the mean height of the urban canopy in the grid cell.
             weighted_height(j, i) = sum_height_length(j, i) / sum_length(j, i);
           }
  }

  //! Background concentrations are corrected accounting for the urban fraction.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::CorrectBackgroundConcentration()
  {

    for (int i = 0; i < this->Nx; ++i)
      for (int j = 0; j < this->Ny; ++j)
        {        
          if (is_urban(j, i))
            {
              // Computes ratio of the urban volume and the volume of the grid cell.
              T cell_width_z, cell_width_y, cell_width_x, cell_volume;
              ComputeCellWidth(0, j, i, cell_width_z, cell_width_y,
                               cell_width_x, cell_volume);

              // Computes ratio of the non-urban volume over the total volume.
              T factor = cell_volume / (cell_volume - sum_street_volume(j, i));

              // Applies the ratio to compute the new concentrations for each species.
              for (int s = 0; s < this->Ns; ++s)
                {
                  T old_concentration = Model.GetConcentration()(s, 0, j, i);
                  T new_concentration = old_concentration * factor;
                  Model.GetConcentration()(s, 0, j, i) = new_concentration;
                }
            }
        }
  }

  //! Background concentrations are corrected accounting for the urban fraction.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::RestoreBackgroundConcentration()
  {
    for (int i = 0; i < this->Nx; ++i)
      for (int j = 0; j < this->Ny; ++j)
        {        
          if (is_urban(j, i))
            {
              // Computes ratio of the urban volume and the volume of the grid cell.
              T cell_width_z, cell_width_y, cell_width_x, cell_volume;
              ComputeCellWidth(0, j, i, cell_width_z, cell_width_y,
                               cell_width_x, cell_volume);

              // Computes ratio of the non-urban volume over the total volume.
              T factor = cell_volume / (cell_volume - sum_street_volume(j, i));

              // Applies the ratio to compute the new concentrations for each species.
              for (int s = 0; s < this->Ns; ++s)
                {
                  T old_concentration = Model.GetConcentration()(s, 0, j, i);
                  T new_concentration = old_concentration / factor;
                  if (new_concentration < 0.0)
                      new_concentration = 0.0;
                  Model.GetConcentration()(s, 0, j, i) = new_concentration;
                }
            }
        }
  }

  //! Background concentrations are corrected accounting for the urban fraction.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeConcentrationOverCanopy()
  {
    for (int j = 0; j < this->Ny; ++j)
      for (int i = 0; i < this->Nx; ++i)
        {        
          if (is_urban(j, i))
            {
              // Computes ratio of the urban volume and the volume of the grid cell.
              T cell_width_z, cell_width_y, cell_width_x, cell_volume;
              ComputeCellWidth(0, j, i, cell_width_z, cell_width_y,
                               cell_width_x, cell_volume);

              // Computes ratio of the non-urban volume over the total volume.
              T H = weighted_height(j, i);
              T z = cell_width_z;
              T alpha = urban_fraction(j, i);
              T factor = 1.0;

              // Applies the ratio to compute the new concentrations for each species.
              for (int s = 0; s < this->Ns; ++s)
                {
                  T old_concentration = Model.GetConcentration()(s, 0, j, i);
                  T added_mass = mass_flux_grid(s, j, i) * this->Delta_t;
                  T new_concentration = (old_concentration * factor * cell_volume + 
                                         added_mass) / (factor * cell_volume);
                  if (new_concentration < 0.0)
                      new_concentration = 0.0;
                  ConcentrationOverCanopy(s, j, i) = new_concentration;
                }
            }
          else //! Non-urban
            {
              for (int s = 0; s < this->Ns; ++s)
                ConcentrationOverCanopy(s, j, i) = Model.GetConcentration()(s, 0, j, i);
            }
        }
  }


  //! Returns the concentrations Data.
  /*!
    \return The concentrations Data.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  Data<T, 3>& StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetConcentrationOverCanopy()
  {
    return ConcentrationOverCanopy;
  }



  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */  
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::InitStep()
  {
    //! InitStep of Polair3DTransport
    Model.InitStep();

    UpdateMeteo();


     /*** Meteo data on the streets***/

    Nstreet = StreetNetworkModel->GetNStreet();
    // Updating the meteorological data for the streets.
    for (int street_index = 0; street_index < Nstreet; street_index++)
        UpdateStreetData(street_index);

    Nintersection = StreetNetworkModel->GetNumberIntersection();
    // Updating the meteorological data for the intersections.
    for (int intersection_index = 0; intersection_index < Nintersection;
         intersection_index++)
      UpdateIntersectionData(intersection_index);

    StreetNetworkModel->InitStep();

    bool emission_save = true;
    if (emission_save)
      EmissionRateSave();

  }

  //! 
  /*! 
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateMeteo()
  {

    /*** Winds ***/
    
    this->UpdateData("meteo", "ZonalWind", FileZonalWind_i,
                     FileZonalWind_f, ZonalWind_i);
    
    this->UpdateData("meteo", "MeridionalWind", FileMeridionalWind_i,
                     FileMeridionalWind_f, MeridionalWind_i);

    /*** Temperature ***/
    
    this->UpdateData("meteo", "Temperature", FileTemperature_i,
                     FileTemperature_f, Temperature_i);
    
    /*** FrictionModule ***/
    
    this->UpdateData("meteo", "FrictionModule",
                     FileFrictionModule_i, FileFrictionModule_f,
                     FrictionModule_i);
    
    /*** BoundaryHeight ***/
    
    this->UpdateData("meteo", "BoundaryHeight",
                     FileBoundaryHeight_i, FileBoundaryHeight_f,
                     BoundaryHeight_i);
    
    /*** LMO ***/
    
    this->UpdateData("meteo", "LMO", FileLMO_i,
                     FileLMO_f, LMO_i);


     if (option_deposition)
       DepositionVelocity_i.Copy(Model.D3("DepositionVelocity_i"));

     if (option_scavenging)
       ScavengingCoefficient_i.Copy(Model.D4("ScavengingCoefficient_i"));


  }



  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    int rank;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    if (rank == 0)
      {
        cell_quantity.resize(this->Ns, this->Ny, this->Nx);
        cell_quantity = 0.0;

        CorrectBackgroundConcentration();

        UpdateBackgroundConcentration();

        // Street-network model.
        StreetNetworkModel->Transport();

        StreetNetworkModel->SetStreetConcentration();

        // Compute the mass flux.
        ComputeMassFlux();

        // Correct the background concentration by the ratio of the urban volume.
        RestoreBackgroundConcentration();

        ComputeConcentrationOverCanopy();

      }

    // Eulerian model.
    Model.Forward();

    if (rank == 0)
      {
    /* Adding all streets to concentration Data in case concentrations are saved
       on the domain.*/

        StreetConcentration.Copy(StreetNetworkModel->GetStreetConcentration());

        this->Concentration.Copy(Model.GetConcentration());
        Nstreet = StreetNetworkModel->GetNStreet();
        for (int street_index = 0; street_index < Nstreet; street_index++)
          {
            // Street center coordinates.
            T z_c, lon_c, lat_c;
            z_c = 2.0;
            StreetNetworkModel->
              GetStreetCoordinate(street_index, lon_c, lat_c);
            int iz, iy, ix;
            GetCellIndices(lon_c, lat_c, z_c, iz, iy, ix);

            for (int s = 0; s < this->Ns; ++s)
              {
                T quantity_street = StreetNetworkModel->
                  GetStreetQuantity(street_index, s);
                cell_quantity(s, iy, ix) += quantity_street;
                
              }
          }
		
        // Adding street concentration to the Concentration Data.
        StreetTransfer(this->Concentration);

        /* Increase the time step */
        StreetNetworkModel->AddTime(delta_t_local);
        this->AddTime(delta_t_eulerian);
        this->step++;
      }
  }

  //! Updates the street meteorological data.
  /*! Updates the street meteorological data taking them from the center cell.
    \param street_index index of the street.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateBackgroundConcentration()
  {
    Nstreet = StreetNetworkModel->GetNStreet();
    // Updating the meteorological data for the streets.
    for (int street_index = 0; street_index < Nstreet; street_index++)
      {
        // Street center coordinates.
        T lon_c, lat_c;
        T z_c = 2.0; 
        StreetNetworkModel->
          GetStreetCoordinate(street_index, lon_c, lat_c);
    
        Array<T, 1> background_concentration;
        background_concentration.resize(this->Ns);
        ExtractSpeciesData(z_c, lat_c, lon_c,
                           background_concentration);

        StreetNetworkModel
          ->SetStreetBackgroundConcentration(street_index,
                                             background_concentration);
      }
  }



  //! Updates the street meteorological data.
  /*! Updates the street meteorological data taking them from the center cell.
    \param street_index index of the street.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateStreetData(int street_index)
  {
    // Street center coordinates.
    T lon_c, lat_c;
    T z_c = 2.0; 
    StreetNetworkModel->
      GetStreetCoordinate(street_index, lon_c, lat_c);
    
    // Extracting meteorological data.
    map<string, T> meteorological_data;
    ExtractMeteo(z_c, lat_c, lon_c, meteorological_data);

    T temp = meteorological_data["wind_angle"];
    T wind_angle;
    // temp 
    // wind_angle = 0 for the wind from south to north.
    // wind_angle = pi for the wind from north to south.
    if (temp <= (pi / 2.0))
      wind_angle = pi / 2.0 - temp;
    else
      wind_angle = pi / 2.0 - temp + pi * 2.0;
    meteorological_data["wind_angle"] = wind_angle;

    StreetNetworkModel->
      SetStreetMeteo(street_index,
		     meteorological_data["wind_angle"],
		     meteorological_data["wind"],
		     meteorological_data["boundaryheight"],
		     meteorological_data["frictionmodule"],
		     meteorological_data["lmo"]);

    Array<T, 1> background_concentration;
    background_concentration.resize(this->Ns);
    ExtractSpeciesData(z_c, lat_c, lon_c,
        	       background_concentration);

    StreetNetworkModel
      ->SetStreetBackgroundConcentration(street_index,
                                         background_concentration);

  }

  //! Updates the intersection meteorological data.
  /*! Updates the intersection meteorological data taking them from the center cell.
    \param intersection_index index of the intersection.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateIntersectionData(int intersection_index)
  {
    // Street center coordinates.
    T lon_c, lat_c;
    T z_c = 2.0;
    StreetNetworkModel->
      GetIntersectionCoordinate(intersection_index, lon_c, lat_c);
    
    // Extracting meteorological data.
    map<string, T> meteorological_data;

    bool check = false;
    int id = StreetNetworkModel->GetIntersectionID(intersection_index);
    ExtractMeteo(z_c, lat_c, lon_c, meteorological_data, check);

    T temp = meteorological_data["wind_angle"];
    T wind_angle;
    if (temp <= (pi / 2.0))
      wind_angle = pi / 2.0 - temp;
    else
      wind_angle = pi / 2.0 - temp + pi * 2.0;
    meteorological_data["wind_angle"] = wind_angle;

    StreetNetworkModel->
      SetIntersectionMeteo(intersection_index,
                           meteorological_data["wind_angle"],
                           meteorological_data["wind"],
                           meteorological_data["boundaryheight"],
                           meteorological_data["frictionmodule"],
                           meteorological_data["lmo"]);
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
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractMeteo(T height, T lat, T lon, map<string, T> &met_data, bool check)
  {
    const T pi = 3.14159265358979323846264;
    int index_x, index_y, index_z;
    GetCellIndices(lon, lat, height, index_z, index_y, index_x);
    Array<T, 1> Coord3D(3);
    Array<T, 1> Coord2D(2);
    Coord3D(0) = height;
    Coord3D(1) = lat;
    Coord3D(2) = lon;
    Coord2D(0) = lat;
    Coord2D(1) = lon;

    T temperature, wind_angle, wind, meridional, zonal, boundary_height;

    // Interpolation for winds on street center.
    LinearInterpolationPoint(MeridionalWind_i, Coord3D, meridional);
    LinearInterpolationPoint(ZonalWind_i, Coord3D, zonal);

    // Test if other fields have to be interpolated.
    if (!interpolated)
      {
	temperature = Temperature_i(index_z, index_y, index_x);
	boundary_height = BoundaryHeight_i(index_y, index_x);
      }
    else
      {
	LinearInterpolationPoint(Temperature_i, Coord3D, temperature);
	LinearInterpolationPoint(BoundaryHeight_i,
				 Coord2D, boundary_height);
      }
    wind_angle = atan2(meridional, zonal);
    wind = sqrt(zonal * zonal + meridional * meridional);

    met_data["temperature"] = temperature;
    met_data["wind_angle"] = wind_angle;
    met_data["wind"] = wind;
    met_data["boundaryheight"] = boundary_height;

    T friction_velocity, lmo;
    if (!interpolated)
      {
        friction_velocity = FrictionModule_i(index_y, index_x);
        lmo = LMO_i(index_y, index_x);
      }
    else
      {
        LinearInterpolationPoint(FrictionModule_i,
                                 Coord2D, friction_velocity);
        LinearInterpolationPoint(LMO_i, Coord2D, lmo);
      }
    met_data["frictionmodule"] = friction_velocity;
    met_data["lmo"] = lmo;
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
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractMeteo(T height, T lat, T lon, map<string, T> &met_data)
  {
    const T pi = 3.14159265358979323846264;
    int index_x, index_y, index_z;
    GetCellIndices(lon, lat, height, index_z, index_y, index_x);
    Array<T, 1> Coord3D(3);
    Array<T, 1> Coord2D(2);
    Coord3D(0) = height;
    Coord3D(1) = lat;
    Coord3D(2) = lon;
    Coord2D(0) = lat;
    Coord2D(1) = lon;

    T temperature, wind_angle, wind, meridional, zonal, boundary_height;

    // Interpolation for winds on street center.
    LinearInterpolationPoint(MeridionalWind_i, Coord3D, meridional);
    LinearInterpolationPoint(ZonalWind_i, Coord3D, zonal);

    // Test if other fields have to be interpolated.
    if (!interpolated)
      {
	temperature = Temperature_i(index_z, index_y, index_x);
	boundary_height = BoundaryHeight_i(index_y, index_x);
      }
    else
      {
	LinearInterpolationPoint(Temperature_i, Coord3D, temperature);
	LinearInterpolationPoint(BoundaryHeight_i,
				 Coord2D, boundary_height);
      }
    wind_angle = atan2(meridional, zonal);
    wind = sqrt(zonal * zonal + meridional * meridional);

    met_data["temperature"] = temperature;
    met_data["wind_angle"] = wind_angle;
    met_data["wind"] = wind;
    met_data["boundaryheight"] = boundary_height;

    T friction_velocity, lmo;
    if (!interpolated)
      {
        friction_velocity = FrictionModule_i(index_y, index_x);
        lmo = LMO_i(index_y, index_x);
      }
    else
      {
        LinearInterpolationPoint(FrictionModule_i,
                                 Coord2D, friction_velocity);
        LinearInterpolationPoint(LMO_i, Coord2D, lmo);
      }
    met_data["frictionmodule"] = friction_velocity;
    met_data["lmo"] = lmo;
  }

  //! Extracts data that depend on meteorology and species.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractSpeciesData(T height, T lat, T lon,
		       Array<T, 1>& background_concentration)
  {
    int r, s, index_x, index_y, index_z;
    GetCellIndices(lon, lat, height, index_z, index_y, index_x);

    Array<T, 1> Coord4D(4);
    Coord4D(0) = 0;
    Coord4D(1) = height;
    Coord4D(2) = lat;
    Coord4D(3) = lon;
    
    // Background concentrations.
    for (s = 0; s < this->Ns; s++)
      {
        Coord4D(0) = s;
        if (!interpolated)
          background_concentration(s) =
            Model.GetConcentration()(s, index_z, index_y, index_x);
        else
          LinearInterpolationPoint(Model.GetConcentration(),
                                   Coord4D, background_concentration(s));

        // To correct extrapolation that could give concentrations < 0.
        background_concentration(s) =
          max(background_concentration(s), 0.);
      }
  }

  /*! Transfers street concentrations to the Eulerian model.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::StreetTransfer(Data<T, 4>& Concentration_out)
  {
    // Concentration to be added to each cell.
    for (int i = 0; i < this->Nx; ++i)
      for (int j = 0; j < this->Ny; ++j)
         if (is_urban(j, i))
           {   
             T cell_width_z, cell_width_y, cell_width_x, cell_volume;
             ComputeCellWidth(0, j, i, cell_width_z, cell_width_y,
                              cell_width_x, cell_volume);

             // Computes ratio of the non-urban volume over the total volume.
             T factor = cell_volume / (cell_volume - sum_street_volume(j, i));

             for (int s = 0; s < this->Ns; ++s)
               {
                 T quantity_background = Concentration_out(s, 0, j, i) *
                   cell_volume; // in mass (ug)
                 T new_concentration = (quantity_background + 
                                        cell_quantity(s, j, i)) / cell_volume;

                 Concentration_out(s, 0, j, i)
                   = max(new_concentration, 0.);
               }
           }
  }



  //! Computes cell width (in meters) and cell volume at given indices.
  /*!
    \param k cell index along z.
    \param j cell index along y.
    \param i cell index along x.
    \param CellWidth_x (output) cell width along x in meters.
    \param CellWidth_y (output) cell width along y in meters.
    \param CellWidth_z (output) cell width along z in meters.
    \param CellVolume (output) cell volume in cubic meters.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeCellWidth(int k, int j, int i, T& CellWidth_z,
		     T& CellWidth_y, T& CellWidth_x, T& CellVolume)
  {

    const T pi = 3.14159265358979323846264;
    const T earth_radius = 6371229.;

    if (!option_cartesian)
      {
	T lat = this->y_min + j * this->Delta_y;
	CellWidth_y = earth_radius * this->Delta_y * pi / 180.;
	CellWidth_x = earth_radius * this->Delta_x * pi / 180.
	  * cos(lat * pi / 180.);
      }
    else
      {
	CellWidth_x = this->Delta_x;
	CellWidth_y = this->Delta_y;
      }
    if (altitude_field)
      CellWidth_z = Model.D3("Altitude")(k + 1, j, i)
	- Model.D3("Altitude")(k, j, i);
    else
      CellWidth_z = this->GridZ3D_interf(k + 1) - this->GridZ3D_interf(k);
    CellVolume = CellWidth_x * CellWidth_y * CellWidth_z;
  }


  //! Computes indices of the cell containing a given point in Eulerian grid.
  /*!
    \param lon longitude of the point (degrees).
    \param lat latitude of the point (degrees).
    \param height of the point (meters).
    \param k (output) cell index along z.
    \param j (output) cell index along y.
    \param i (output) cell index along x.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetCellIndices(T lon, T lat, T height,
		   int& index_z, int& index_y, int& index_x)
  {
    index_x = max(int((lon - this->x_min + this->Delta_x / 2.)
		      / this->Delta_x), 0);
    index_x = min(index_x, this->Nx - 1);
    index_y = max(int((lat - this->y_min + this->Delta_y / 2.)
		      / this->Delta_y), 0);
    index_y = min(index_y, this->Ny - 1);
    if (altitude_field)
      for (int k=0; k < this->Nz; k++)
	{
	  index_z = k;
	  if (Model.D3("Altitude")(k+1, index_y, index_x) > height)
	    break;
	}
    else
      for (int k=0; k < this->Nz; k++)
	{
	  index_z = k;
	  if (this->GridZ3D_interf(k+1) > height)
	    break;
	}
  }

  ///////////////////////////////////////////
  // Additional functions are needed       //
  // to save concentrations at the streets //
  ///////////////////////////////////////////
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  int StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetNStreet() const
  {
    return StreetNetworkModel->GetNStreet();
  }

  template<class T, class ClassEulerianModel, class ClassLocalModel>
  Data<T, 2>& StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetStreetConcentration()
  {
    return StreetNetworkModel->GetStreetConcentration();
  }

  //! Make binaray files for traffic emissions in the corresponding corrgrid cells. Make sure the path to the output files.
  /*!
    \param k cell index along z.
    \param j cell index along y.
    \param i cell index along x.
    \param CellWidth_x (output) cell width along x in meters.
    \param CellWidth_y (output) cell width along y in meters.
    \param CellWidth_z (output) cell width along z in meters.
    \param CellVolume (output) cell volume in cubic meters.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::EmissionRateSave()
  {
    emission_rate.resize(this->Ns, this->Ny, this->Nx);
    emission_rate = 0.0;

    Array<T, 2> rate_NO2(this->Ny, this->Nx);
    rate_NO2 = 0.0;
    Array<T, 2> rate_NO(this->Ny, this->Nx);
    rate_NO = 0.0;
    for (int street_index = 0; street_index < Nstreet; street_index++)
      {
        // Street center coordinates.
        T z_c, lon_c, lat_c;
        z_c = 2.0;
        StreetNetworkModel->
          GetStreetCoordinate(street_index, lon_c, lat_c);
        int iz, iy, ix;
        GetCellIndices(lon_c, lat_c, z_c, iz, iy, ix);

        // Computes ratio of the urban volume and the volume of the grid cell.
        T cell_width_z, cell_width_y, cell_width_x, cell_volume, surface;
        ComputeCellWidth(iz, iy, ix, cell_width_z, cell_width_y,
                         cell_width_x, cell_volume);
        surface = cell_width_x * cell_width_y;
        for (int s = 0; s < this->Ns; ++s)
          {
            emission_rate(s, iy, ix) += (StreetNetworkModel->GetStreetEmissionRate(street_index, s) / surface);
            if (s == 51)
              {
                rate_NO2(iy, ix) += (StreetNetworkModel->GetStreetEmissionRate(street_index, s) / surface);
              }
            if (s == 47)
              {
                rate_NO(iy, ix) += (StreetNetworkModel->GetStreetEmissionRate(street_index, s) / surface);
              }
          }
      }

    string filename_no2 = "/net/libre/yomi/kimy/StreetInGrid/Trafipollu/sing-voc/nox20-streetvolume/reference/results_test/emission_collect/NO2.bin";
    FormatBinary<float>().Append(rate_NO2, filename_no2);
    string filename_no = "/net/libre/yomi/kimy/StreetInGrid/Trafipollu/sing-voc/nox20-streetvolume/reference/results_test/emission_collect/NO.bin";
    FormatBinary<float>().Append(rate_NO, filename_no);
  }   
    

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETINGRIDTRANSPORT_CXX
#endif
