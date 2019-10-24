// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMEINGRIDTRANSPORT_CXX


#include "PlumeInGridTransport.hxx"


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
  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::PlumeInGridTransport(string config_file): BaseModel<T>(config_file),
    Model(config_file)
  {
    /*** Pointers to 3D and 2D data ***/

    this->D3_map["MeridionalWind_i"] = &MeridionalWind_i;
    this->D3_map["ZonalWind_i"] = &ZonalWind_i;
    this->D3_map["Temperature_i"] = &Temperature_i;
    this->D3_map["DryDepositionFlux"] = &DryDepositionFlux;
    this->D3_map["WetDepositionFlux"] = &WetDepositionFlux;
    this->D3_map["InCloudWetDepositionFlux"]
      = &InCloudWetDepositionFlux;

    this->D2_map["FirstLevelWind_i"] = &FirstLevelWind_i;
    this->D2_map["LowCloudiness_i"] = &LowCloudiness_i;
    this->D2_map["MediumCloudiness_i"] = &MediumCloudiness_i;
    this->D2_map["HighCloudiness_i"] = &HighCloudiness_i;
    this->D2_map["Insolation_i"] = &Insolation_i;
    this->D2_map["FrictionModule_i"] = &FrictionModule_i;
    this->D2_map["BoundaryHeight_i"] = &BoundaryHeight_i;
    this->D2_map["LMO_i"] = &LMO_i;
    this->field_species["ScavengingCoefficient"] = &species_list_scav;
    this->field_species["ScavengingBelowCloudCoefficient"] = &species_list_scav;
    this->field_species["ScavengingInCloudCoefficient"] = &species_list_scav;
    //    this->field_species["ScavengingCoefficient"] = &species_list_scav;
    this->field_species["DepositionVelocity"] = &species_list_dep;
  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::~PlumeInGridTransport()
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
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    this->config.SetSection("[domain]");
    this->config.PeekValue("Cartesian", option_cartesian);

    this->config.SetSection("[options]");
    this->config.PeekValue("With_deposition", option_deposition);

    string scavenging_model;
    this->config.PeekValue("Scavenging_below_cloud_model",
                           "none|constant|belot|microphysical|microphysical-ph",
                           scavenging_model);
    // in-cloud scavenging need to be taken into account in the future version.
    // (YK: 2018/04/06)
    // this->config.PeekValue("Scavenging_in_cloud_model",
    //                        "none|belot|pudykiewicz|aqueous",
    //                        this->scavenging_in_cloud_model);
    option_scavenging = (scavenging_model != "none");

    this->config.SetSection("[data]");

    string data_file, emission_file;
    this->config.PeekValue("Data_description", data_file);


    // Scavenging.
    this->config.SetSection("[options]");
    this->config.PeekValue("With_deposition",
			   this->option_process["with_deposition"]);
    if (this->option_process["with_deposition"])
      this->input_files["deposition_velocity"].Read(data_file,
						    "deposition");
    else
      this->input_files["deposition_velocity"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["deposition_velocity"].Begin();
	 i != this->input_files["deposition_velocity"].End(); i++)
      species_list_dep.push_back(i->first);
    this->config.PeekValue("Collect_dry_flux",
			   this->option_process["collect_dry_flux"]);
    this->config.PeekValue("Collect_wet_flux",
			   this->option_process["collect_wet_flux"]);

    this->config.SetSection("[data]");
    if (scavenging_model != "none")
      this->input_files["scavenging_coefficient"]
	.ReadFields(data_file, "scavenging");
    else
      this->input_files["scavenging_coefficient"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["scavenging_coefficient"].Begin();
	 i != this->input_files["scavenging_coefficient"].End(); i++)
      species_list_scav.push_back(i->first);
    int Ns_scav = int(species_list_scav.size());



    // Plume-in-grid general data.
    ConfigStream data(data_file);

    data.SetSection("[ground]");

    data.PeekValue("LUC_file", LUC_file);
    if (LUC_file == "rural" || LUC_file == "urban")
      {
        land_type = LUC_file;
        Nc = 1;
      }
    else
      {
        if (!exists(LUC_file))
          throw "Unable to open land use cover file \"" + LUC_file + "\".";
        land_type = "LUC";
        Nc = int(file_size(LUC_file)) / sizeof(float) / (this->Ny * this->Nx);
        data.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                       urban_index);
        data.PeekValue("Urban_proportion", ">=0 | <= 1", percent_urban);
      }

    data.SetSection("[plume-in-grid]");
    data.PeekValue("With_interpolation", interpolated);
    data.PeekValue("Horizontal_coefficient", "positive", coefficient_y);
    data.PeekValue("Vertical_coefficient", "positive", coefficient_z);
    data.PeekValue("Save_Gaussian_concentrations", save_gaussian_domain);

    // Injection data.
    string injection_method;
    string merging_criteria;

    data.PeekValue("With_reinjection_time", option_time);
    data.PeekValue("Reinjection_time", "positive", reinjection_time);
    data.PeekValue("Injection_method", "integrated | column",
                   injection_method);
    if (injection_method == "integrated")
      injection_integrated = 1;
    else if (injection_method == "column")
      injection_integrated = 0;

    data.PeekValue("Merge_overlaping_puff", merge_puff);
    data.PeekValue("Merging_criteria", merging_criteria);
    if (merging_criteria == "source_id")
      merge_source_id = 1;
    else
      merge_source_id = 0;

    /*** Gaussian configuration ***/

    this->config.SetSection("[gaussian]");
    this->config.PeekValue("file_gaussian", gaussian_file);
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::CheckConfiguration()
  {
    BaseModel<T>::CheckConfiguration();
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    ConfigStream data_description_stream(data_description_file);
    this->input_files["meteo"].Read(data_description_file, "meteo");
    this->input_files["gaussian_meteo"]
      .Read(data_description_file, "gaussian_meteo");

    if (this->input_files["gaussian_meteo"]("FirstLevelWindModule").empty())
      throw "FirstLevelWindModule is needed but no input data file was provided.";

  }

  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Allocate()
  {
    BaseModel<T>::Allocate();

    /*** Ground ***/

    if (land_type == "LUC")
      {
        RegularGrid<T> GridC(Nc);
        LUC.Resize(GridC, this->GridY3D, this->GridX3D);
      }

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

    /*** First level wind module ***/

    FirstLevelWind_i.Resize(this->GridY2D, this->GridX2D);
    FileFirstLevelWind_i.Resize(this->GridY2D, this->GridX2D);
    FileFirstLevelWind_f.Resize(this->GridY2D, this->GridX2D);

    /*** Cloudiness ***/

    LowCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileLowCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileLowCloudiness_f.Resize(this->GridY2D, this->GridX2D);

    MediumCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileMediumCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileMediumCloudiness_f.Resize(this->GridY2D, this->GridX2D);

    HighCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileHighCloudiness_i.Resize(this->GridY2D, this->GridX2D);
    FileHighCloudiness_f.Resize(this->GridY2D, this->GridX2D);

    /*** Insolation ***/

    Insolation_i.Resize(this->GridY2D, this->GridX2D);
    FileInsolation_i.Resize(this->GridY2D, this->GridX2D);
    FileInsolation_f.Resize(this->GridY2D, this->GridX2D);

    /*** FrictionModule ***/

    FrictionModule_i.Resize(this->GridY2D, this->GridX2D);
    FileFrictionModule_i.Resize(this->GridY2D, this->GridX2D);
    FileFrictionModule_f.Resize(this->GridY2D, this->GridX2D);

    /*** BoundaryHeight ***/

    BoundaryHeight_i.Resize(this->GridY2D, this->GridX2D);
    FileBoundaryHeight_i.Resize(this->GridY2D, this->GridX2D);
    FileBoundaryHeight_f.Resize(this->GridY2D, this->GridX2D);

    /*** LMO ***/

    LMO_i.Resize(this->GridY2D, this->GridX2D);
    FileLMO_i.Resize(this->GridY2D, this->GridX2D);
    FileLMO_f.Resize(this->GridY2D, this->GridX2D);

  }


  //! Initialization.
  /*! Initializes the Eulerian model and the meteorological conditions. Reads
    the sources and creates the corresponding Gaussian model instances.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    BaseModel<T>::Init();

    /*** Initializations ***/
    Model.Init();
    this->Concentration.Copy(Model.GetConcentration());
    if (this->option_process["collect_dry_flux"])      
      this->DryDepositionFlux.Copy(Model.GetDryDepositionFlux());
    if (this->option_process["collect_wet_flux"])
      {
	this->WetDepositionFlux.Copy(Model.GetWetDepositionFlux());
	this->InCloudWetDepositionFlux.Copy(Model.GetInCloudWetDepositionFlux());
      }
    // Extracts simulation information from Eulerian model.
    delta_t_eulerian = Model.GetDelta_t();
    if (this->Delta_t != delta_t_eulerian
        || this->Date_min != Model.GetDate_min())
      throw string("Time steps and date for Eulerian model") +
        " and plume-in-grid model must be equal.";

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
        if (Model.HasField("ScavengingBelowCloudCoefficient"))
          ScavengingCoefficient_i.Copy(Model.D4("ScavengingBelowCloudCoefficient_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"ScavengingCoefficient\".";
      }

    // Is there a field "Altitude" to manage vertical levels?
    altitude_field = Model.HasField("Altitude");

    /*** Ground data ***/

    if (land_type == "LUC")
      {
        FormatBinary<float> InputLUC;
        InputLUC.Read(LUC_file, LUC);
      }

    /***  Meteorological data ***/

    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    ConfigStream data_description_stream(data_description_file);
    this->input_files["meteo"].Read(data_description_file, "meteo");
    this->input_files["gaussian_meteo"]
      .Read(data_description_file, "gaussian_meteo");

    /*** Winds ***/
    this->InitData("meteo", "ZonalWind", FileZonalWind_i, FileZonalWind_f,
                   this->current_date, ZonalWind_i);
    this->InitData("meteo", "MeridionalWind", FileMeridionalWind_i,
                   FileMeridionalWind_f, this->current_date,
                   MeridionalWind_i);

    /*** Temperature ***/
    this->InitData("meteo", "Temperature", FileTemperature_i,
                   FileTemperature_f, this->current_date, Temperature_i);

    /*** FirstLevelWind ***/
    this->InitData("gaussian_meteo", "FirstLevelWindModule",
                   FileFirstLevelWind_i, FileFirstLevelWind_f,
                   this->current_date, FirstLevelWind_i);

    /*** LowCloudiness ***/
    this->InitData("gaussian_meteo", "LowCloudiness", FileLowCloudiness_i,
                   FileLowCloudiness_f, this->current_date, LowCloudiness_i);

    /*** MediumCloudiness ***/
    this->InitData("gaussian_meteo", "MediumCloudiness",
                   FileMediumCloudiness_i, FileMediumCloudiness_f,
                   this->current_date, MediumCloudiness_i);

    /*** HighCloudiness ***/
    this->InitData("gaussian_meteo", "HighCloudiness",
                   FileHighCloudiness_i, FileHighCloudiness_f,
                   this->current_date, HighCloudiness_i);

    /*** Insolation ***/
    this->InitData("gaussian_meteo", "SolarRadiation", FileInsolation_i,
                   FileInsolation_f, this->current_date, Insolation_i);

    /*** FrictionModule ***/
    this->InitData("gaussian_meteo", "FrictionModule",
                   FileFrictionModule_i, FileFrictionModule_f,
                   this->current_date, FrictionModule_i);

    /*** BoundaryHeight ***/
    this->InitData("gaussian_meteo", "BoundaryHeight",
                   FileBoundaryHeight_i, FileBoundaryHeight_f,
                   this->current_date, BoundaryHeight_i);

    /*** LMO ***/
    this->InitData("gaussian_meteo", "LMO", FileLMO_i,
                   FileLMO_f, this->current_date, LMO_i);

    /*** Emissions ***/

    // Gaussian Model creation.
    GaussianPuffModel = new ClassLocalModel(gaussian_file);

    // Gaussian model initialization.
    GaussianPuffModel->Init();
    delta_t_local = min(GaussianPuffModel->GetDelta_t(), delta_t_eulerian);
    Nt_local = int(delta_t_eulerian / delta_t_local);
    this->evolutive_plume_rise = this->GaussianPuffModel->ComputeVerticalTrajectory();

    GaussianPuffModel->SetDateMin(this->Date_min);
    GaussianPuffModel->SetTimeStep(delta_t_local);
    GaussianPuffModel->SetNt(this->Nt * Nt_local);

    // Is there puff deposition?
    if (!GaussianPuffModel->WithDeposition())
      option_deposition = false;
    else if (!option_deposition)
      throw string("Depostion cannot be activated in GaussianPuff ")
        + "model without deposition "
        + "in Eulerian model.";

    // Is there puff scavenging?
    if (!GaussianPuffModel->WithScavenging())
      option_scavenging = false;
    else if (!option_scavenging)
      throw string("Scavenging cannot be activated in GaussianPuff ")
        + "model without scavenging "
        + "in Eulerian model.";

    // Sources initialization in Gaussian Puff model.
    GaussianPuffModel->InitPuffSource();

    // Conversion of sources coordinates from lat/lon to cartesian.
    Array<T, 2> source_coordinates;
    GaussianPuffModel->GetSourceCoordinates(source_coordinates);
    int Nemis = source_coordinates.extent(0);

    // Array containing the sources cartesian coordinates.
    Array<T, 2> source_coordinates_tmp = source_coordinates;
    for (int index = 0; index < Nemis; index++)
      {
        T abscissa, ordinate;
        LatLonToCartesian(source_coordinates(index, 0),
                          source_coordinates(index, 1),
                          abscissa, ordinate);
        source_coordinates_tmp(index, 0) = abscissa;
        source_coordinates_tmp(index, 1) = ordinate;
      }
    GaussianPuffModel->SetSourceCoordinates(source_coordinates_tmp);

    // Clears the list of puffs and initializes several variables.
    GaussianPuffModel->InitPosition();
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */  template<class T, class ClassEulerianModel, class ClassLocalModel>
   void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
   ::InitStep()
   {
     Model.InitStep();

     /*** Winds ***/

     this->UpdateData("meteo", "ZonalWind", FileZonalWind_i,
                      FileZonalWind_f, ZonalWind_i);

     this->UpdateData("meteo", "MeridionalWind", FileMeridionalWind_i,
                      FileMeridionalWind_f, MeridionalWind_i);

     /*** Temperature ***/

     this->UpdateData("meteo", "Temperature", FileTemperature_i,
                      FileTemperature_f, Temperature_i);

     /*** FirstLevelWind ***/

     this->UpdateData("gaussian_meteo", "FirstLevelWindModule",
                      FileFirstLevelWind_i, FileFirstLevelWind_f,
                      FirstLevelWind_i);

     /*** LowCloudiness ***/

     this->UpdateData("gaussian_meteo", "LowCloudiness", FileLowCloudiness_i,
                      FileLowCloudiness_f,  LowCloudiness_i);

     /*** MediumCloudiness ***/

     this->UpdateData("gaussian_meteo", "MediumCloudiness",
                      FileMediumCloudiness_i, FileMediumCloudiness_f,
                      MediumCloudiness_i);

     /*** HighCloudiness ***/

     this->UpdateData("gaussian_meteo", "HighCloudiness",
                      FileHighCloudiness_i, FileHighCloudiness_f,
                      HighCloudiness_i);

     /*** Insolation ***/

     this->UpdateData("gaussian_meteo", "SolarRadiation", FileInsolation_i,
                      FileInsolation_f, Insolation_i);

     /*** FrictionModule ***/

     this->UpdateData("gaussian_meteo", "FrictionModule",
                      FileFrictionModule_i, FileFrictionModule_f,
                      FrictionModule_i);

     /*** BoundaryHeight ***/

     this->UpdateData("gaussian_meteo", "BoundaryHeight",
                      FileBoundaryHeight_i, FileBoundaryHeight_f,
                      BoundaryHeight_i);

     /*** LMO ***/

     this->UpdateData("gaussian_meteo", "LMO", FileLMO_i,
                      FileLMO_f, LMO_i);

     if (option_deposition)
       DepositionVelocity_i.Copy(Model.D3("DepositionVelocity_i"));

     if (option_scavenging)
       ScavengingCoefficient_i.Copy(Model.D4("ScavengingBelowCloudCoefficient_i"));

   }


  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    this->Concentration.Copy(Model.GetConcentration());
    if (this->option_process["collect_dry_flux"])      
      this->DryDepositionFlux.Copy(Model.GetDryDepositionFlux());
    if (this->option_process["collect_wet_flux"])
      {
	this->WetDepositionFlux.Copy(Model.GetWetDepositionFlux());
	this->InCloudWetDepositionFlux.Copy(Model.GetInCloudWetDepositionFlux());
      }

    int Npuff;
    int index_x, index_y, index_z, puff_index;
    bool puff_transfer;
    bool isday;

    /*** Tests if puff has to be transfered to Eulerian model. ***/

    Npuff = GaussianPuffModel->GetPuffNumber();
    for (puff_index = 0; puff_index < Npuff; puff_index++)
      {
        // Puff center coordinates.
        T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
        GaussianPuffModel->
          GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                          puff_time);

        // Conversion to longitude/latitude (degrees).
        CartesianToLatLon(x_c, y_c, lon_c, lat_c);

        // Night or day.
        isday = IsDay(lon_c, lat_c, GaussianPuffModel->
                      GetCurrentDate());

        // Gets corresponding cell in eulerian grid.
        GetCellIndices(lon_c, lat_c, z_c,
                       index_z, index_y, index_x);

        // Computing cell width along y.
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
                         cell_width_y, cell_width_x, cell_volume);

        // Getting puff standard deviations.
        T sigma_x, sigma_y, sigma_z;
        GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                        sigma_y, sigma_z);
        puff_transfer = 0;

        // Tests if puff has reached the end of the domain.
        if (index_x == 0 || index_x == this->Nx - 1
            || index_y == 0 || index_y == this->Ny - 1
            || index_z == this->Nz)
          GaussianPuffModel->ErasePuff(puff_index);

        // Tests if puff has reached the injection time.
        else if (option_time)
          {
            if (puff_time >= reinjection_time)
              puff_transfer = 1;
          }

        // Tests if puff size has reached the cell width.
        else if (coefficient_y * sigma_y >= cell_width_y)
          puff_transfer = 1;

        // Puff transfer.
        if (puff_transfer)
          {
            if (injection_integrated)
              PuffIntegratedTransfer(puff_index, Model.GetConcentration());
            else
              for (int s = 0; s < this->Ns; s++)
                PuffTransfer(GaussianPuffModel->
                             GetPuffQuantity(puff_index, s),
                             sigma_z, s, z_c, lat_c, lon_c,
                             isday, Model.GetConcentration());
            GaussianPuffModel->ErasePuff(puff_index);
          }
      }

    /*** Inner time-loop for Gaussian models. ***/

    Array<T, 4> PuffConcentration(this->Ns, this->Nz, this->Ny, this->Nx);
    PuffConcentration = 0.;
    int l, k, j, i, s;
    for (j = 0; j < Nt_local; j++)
      {
        GaussianPuffModel->InitStep();
        Npuff = GaussianPuffModel->GetPuffNumber();
        Array<int, 2> PuffCellList(Npuff, 3);

        // Loop on puffs to get the meteorological data.
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            CartesianToLatLon(x_c, y_c, lon_c, lat_c);

            // Gets corresponding cell in Eulerian grid.
            GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);
            PuffCellList(puff_index, 0) = index_z;
            PuffCellList(puff_index, 1) = index_y;
            PuffCellList(puff_index, 2) = index_x;

            // Updating the puff meteorological data.
            if (j == 0 || puff_time == 0.)
              UpdateMeteo(puff_index);
          }

        // Puff effective height.
        GaussianPuffModel->ComputePlumeRise();

        // Advection and diffusion for puffs.
        GaussianPuffModel->Advection();
        GaussianPuffModel->Diffusion();
        GaussianPuffModel->ComputeLossFactor();

        // Loop on puffs to update data if necessary.
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            bool has_changed = false;
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            CartesianToLatLon(x_c, y_c, lon_c, lat_c);

            // Gets corresponding cell in eulerian grid.
            GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);

            has_changed = (PuffCellList(puff_index, 0) != index_z
                           || PuffCellList(puff_index, 1) != index_y
                           || PuffCellList(puff_index, 2) != index_x);

            // Updating the meteorological data if necessary.
            if (has_changed)
              {
                UpdateMeteo(puff_index);
                PuffCellList(puff_index, 0) = index_z;
                PuffCellList(puff_index, 1) = index_y;
                PuffCellList(puff_index, 2) = index_x;
              }
          }

        GaussianPuffModel->AddTime(delta_t_local);
      }

    // Eulerian model.
    Model.Forward();

    this->Concentration.Copy(Model.GetConcentration());
    if (this->option_process["collect_dry_flux"])
      this->DryDepositionFlux.Copy(Model.GetDryDepositionFlux());
    if (this->option_process["collect_wet_flux"])
      {    
	this->WetDepositionFlux.Copy(Model.GetWetDepositionFlux());
	this->InCloudWetDepositionFlux.Copy(Model.GetInCloudWetDepositionFlux());
      }

    /* Adding all puffs to concentration Data in case concentrations are saved
       on the domain.*/

    GaussianPuffModel->SubtractTime(delta_t_local);
    Npuff = GaussianPuffModel->GetPuffNumber();
    for (puff_index = 0; puff_index < Npuff; puff_index++)
      {
        T puff_release_time = GaussianPuffModel->
          GetPuffReleaseTime(puff_index);
        if (GaussianPuffModel->GetCurrentTime() >= puff_release_time)
          {
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            CartesianToLatLon(x_c, y_c, lon_c, lat_c);
            T sigma_x, sigma_y, sigma_z;
            GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                            sigma_y, sigma_z);

            isday = IsDay(lon_c, lat_c, GaussianPuffModel->
                          GetCurrentDate());
            // Adding puff concentration to the Concentration Data.
            if (injection_integrated)
              PuffIntegratedTransfer(puff_index, this->Concentration);
            else
              for (int s = 0; s < this->Ns; s++)
                PuffTransfer(GaussianPuffModel->
                             GetPuffQuantity(puff_index, s),
                             sigma_z, s, z_c, lat_c, lon_c,
                             isday, this->Concentration);
          }
      }
    GaussianPuffModel->AddTime(delta_t_local);
    this->AddTime(delta_t_eulerian);
    this->step++;
  }



  //! Updates the puff meteorological data.
  /*! Updates the puff meteorological data taking them from the center cell.
    \param puff_index index of the puff.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateMeteo(int puff_index)
  {
    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    GaussianPuffModel->
      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance, puff_time);

    // Conversion to longitude/latitude (degrees).
    CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Night or day.
    bool isday = IsDay(lon_c, lat_c, GaussianPuffModel->GetCurrentDate());

    // Extracting meteorological data.
    bool rural;
    string stability;
    map<string, T> meteorological_data;
    bool option_similarity = GaussianPuffModel->WithSimilarity();

    ExtractMeteo(z_c, lat_c, lon_c, isday, option_similarity,
                 meteorological_data, stability, rural);

    if (!option_similarity && !this->evolutive_plume_rise)
      GaussianPuffModel->
        SetPuffMeteo(puff_index,
                     meteorological_data["temperature"],
                     meteorological_data["wind_angle"],
                     meteorological_data["wind"],
                     stability, lon_c, lat_c, isday, rural,
                     meteorological_data["boundaryheight"]);
    else
      GaussianPuffModel->
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

    // Data depending on the puff species.
    Array<T, 1> deposition_velocity;
    deposition_velocity.resize(this->Ns);
    deposition_velocity = 0.;

    Array<T, 1> scavenging_coefficient;
    scavenging_coefficient.resize(this->Ns);
    scavenging_coefficient = 0.;

    ExtractSpeciesData(z_c, lat_c, lon_c,
                       deposition_velocity,
                       scavenging_coefficient);

    if (option_deposition)
      GaussianPuffModel
        ->SetPuffDepositionVelocity(puff_index,
                                    deposition_velocity);

    if (option_scavenging)
      GaussianPuffModel
        ->SetPuffScavengingCoefficient(puff_index,
                                       scavenging_coefficient);

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
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractMeteo(T height, T lat, T lon, bool isday,
                 bool option_similarity, map<string, T> &met_data,
                 string &stability, bool& rural)
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

    /*** Base data for Gaussian models. ***/

    T temperature, wind_angle, wind, meridional, zonal, boundary_height;

    // Interpolation for winds on puff center.
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

    if (land_type == "LUC")
      rural = (LUC(urban_index, index_y, index_x) < percent_urban);
    else
      rural = (land_type == "rural");

    if (!option_similarity && !this->evolutive_plume_rise)
      {
        T surface_wind, solar_radiation;
        T cloudiness, low_cloudiness, medium_cloudiness, high_cloudiness;

        if (!interpolated)
          {
            surface_wind = FirstLevelWind_i(index_y, index_x);
            solar_radiation = Insolation_i(index_y, index_x);
            low_cloudiness = LowCloudiness_i(index_y, index_x);
            medium_cloudiness = MediumCloudiness_i(index_y, index_x);
            high_cloudiness = HighCloudiness_i(index_y, index_x);
          }
        else
          {
            LinearInterpolationPoint(FirstLevelWind_i, Coord2D, surface_wind);
            LinearInterpolationPoint(Insolation_i, Coord2D, solar_radiation);
            LinearInterpolationPoint(LowCloudiness_i,
                                     Coord2D, low_cloudiness);
            LinearInterpolationPoint(MediumCloudiness_i,
                                     Coord2D, medium_cloudiness);
            LinearInterpolationPoint(HighCloudiness_i,
                                     Coord2D, high_cloudiness);
          }
        // Computing total cloudiness.
        cloudiness = ComputeTotalCloudiness<T, T>(low_cloudiness,
                                                  medium_cloudiness,
                                                  high_cloudiness);
        // Computing Pasquill stability class.
        stability = ComputePasquillStabilityClass(surface_wind,
                                                  solar_radiation,
                                                  cloudiness, isday);
      }
    else
      {
        T friction_velocity, convective_velocity;
        T lmo, coriolis;

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

        // Computing convective velocity.
        T friction_velocity_3 = friction_velocity * friction_velocity
          * friction_velocity;
        T convective_velocity_3 = - friction_velocity_3 *
          boundary_height / (0.4 * lmo);
        if (convective_velocity_3 < 0.)
          convective_velocity = - pow(- convective_velocity_3, 1. / 3.);
        else
          convective_velocity = pow(convective_velocity_3, 1. / 3.);

        // Computing coriolis parameter.
        coriolis = max(4 * pi  * sin(lat * pi / 180.) / 86400., 5.e-5);

        met_data["frictionmodule"] = friction_velocity;
        met_data["lmo"] = lmo;
        met_data["convectivevelocity"] = convective_velocity;
        met_data["coriolis"] = coriolis;
        if (boundary_height / lmo < - 0.3)
          stability = "unstable";
        else if (boundary_height / lmo > 1.)
          stability = "stable";
        else
          stability = "neutral";
      }
  }


  //! Extracts data that depend on meteorology and species (photolysis rates).
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractSpeciesData(T height, T lat, T lon,
                       Array<T, 1>& deposition_velocity,
                       Array<T, 1>& scavenging_coefficient)
  {
    int r, s, index_x, index_y, index_z;
    GetCellIndices(lon, lat, height, index_z, index_y, index_x);

    if (option_deposition)
      {
        Array<T, 1> Coord3D(3);
        Coord3D(0) = 0;
        Coord3D(1) = lat;
        Coord3D(2) = lon;
        for (s = 0; s < this->Ns; s++)
          {
            if (Model.HasDepositionVelocity(s))
              {
                int dep_s = Model.DepositionVelocityIndex(s);
                Coord3D(0) = dep_s;
                if (!interpolated)
                  deposition_velocity(s) =
                    DepositionVelocity_i(dep_s, index_y, index_x);
                else
                  LinearInterpolationPoint(DepositionVelocity_i, Coord3D,
                                           deposition_velocity(s));
              }
          }
      }

    if (option_scavenging)
      {
        Array<T, 1> Coord4D(4);
        Coord4D(0) = 0;
        Coord4D(1) = height;
        Coord4D(2) = lat;
        Coord4D(3) = lon;
        for (s = 0; s < this->Ns; s++)
          {
            if (Model.HasScavenging(s))
              {
                int scav_s = Model.ScavengingIndex(s);
                Coord4D(0) = scav_s;
                if (!interpolated)
                  scavenging_coefficient(s) =
                    ScavengingCoefficient_i(scav_s, index_z,
                                            index_y, index_x);
                else
                  LinearInterpolationPoint(ScavengingCoefficient_i, Coord4D,
                                           scavenging_coefficient(s));
              }
          }
      }

  }


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
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::PuffTransfer(T quantity, T sigmaz, int s, T z, T lat, T lon, bool isday,
                 Data<T, 4>& Concentration_out)
  {
    int ibz, iby, ibx, itz, ity, itx;
    // Computes puff vertical extent.
    T puff_bottom = max(z - (coefficient_z * sigmaz) / 2, 0.);
    GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    T puff_top = z + (coefficient_z * sigmaz) / 2;
    T boundary_height = BoundaryHeight_i(iby, ibx);

    // Takes into account the inversion height.
    if (isday && z < boundary_height)
      puff_top = min(puff_top, boundary_height);
    if (isday && z > boundary_height)
      puff_bottom = max(puff_bottom, boundary_height);
    GetCellIndices(lon, lat, puff_bottom, ibz, iby, ibx);
    GetCellIndices(lon, lat, puff_top, itz, ity, itx);

    // Number of cells vertically covered by the puff.
    int Ncell = (itz - ibz) + 1;

    // Computing vertical extent of cells where puff will be transfered.
    T vertical_extent = 0.;
    for (int k = 0; k < Ncell; k++)
      vertical_extent += (this->GridZ3D_interf(ibz + k + 1)
                          - this->GridZ3D_interf(ibz + k));

    // Concentration to be added to each cell.
    T cell_width_z, cell_width_y, cell_width_x, cell_volume;
    ComputeCellWidth(ibz, iby, ibx, cell_width_z, cell_width_y,
                     cell_width_x, cell_volume);
    T concentration = quantity /
      (vertical_extent * cell_width_y * cell_width_x);

    //    Adding puff concentration to the Eulerian concentration.
    for (int k = 0; k < Ncell; k++)
      Concentration_out(s, ibz + k, iby, ibx)
        = max(Concentration_out(s, ibz + k, iby, ibx) + concentration, 0.);
  }


  //! Transfers a given puff to the Eulerian model.
  /*! Transfers a puff to the Eulerian model with the integrated method.
    \param puff_index index of the puff.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::PuffIntegratedTransfer(int puff_index, Data<T, 4>& Concentration_out)
  {
    T sigma_x, sigma_y, sigma_z;
    GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                    sigma_y, sigma_z);

    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    GaussianPuffModel->
      GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                      puff_time);
    // Conversion to longitude/latitude (degrees).
    CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Computes puff vertical extent.
    int ibz, iby, ibx, itz, ity, itx;
    T puff_bottom = max(z_c - (coefficient_z * sigma_z) / 2, 0.);
    GetCellIndices(lon_c, lat_c, puff_bottom, ibz, iby, ibx);
    T puff_top = z_c + (coefficient_z * sigma_z) / 2;
    GetCellIndices(lon_c, lat_c, puff_top, itz, ity, itx);

    // Computes puff horizontal extent along y.
    int ilyz, ily, ilyx, iryz, iry, iryx;
    T puff_left = y_c - (coefficient_y * sigma_y) / 2;
    T lon_left, lat_left;
    CartesianToLatLon(x_c, puff_left, lon_left, lat_left);
    GetCellIndices(lon_c, lat_left, z_c, ilyz, ily, ilyx);
    T puff_right = y_c + (coefficient_y * sigma_y) / 2;
    T lon_right, lat_right;
    CartesianToLatLon(x_c, puff_right, lon_right, lat_right);
    GetCellIndices(lon_c, lat_right, z_c, iryz, iry, iryx);

    // Computes puff horizontal extent along x.
    int ilxz, ilxy, ilx, irxz, irxy, irx;
    puff_left = x_c - (coefficient_y * sigma_x) / 2;
    CartesianToLatLon(puff_left, y_c, lon_left, lat_left);
    GetCellIndices(lon_left, lat_c, z_c, ilxz, ilxy, ilx);
    puff_right = x_c + (coefficient_y * sigma_x) / 2;
    CartesianToLatLon(puff_right, y_c, lon_right, lat_right);
    GetCellIndices(lon_right, lat_c, z_c, irxz, irxy, irx);

    T gas_compensation;

    // Number of cells covered by the puff.
    int Nz = (itz - ibz) + 1;
    int Ny = (iry - ily) + 1;
    int Nx = (irx - ilx) + 1;
    int Ns = this->Ns;

    // If the puff is in only one cell.
    if (Nz * Ny * Nx == 1)
      {
        // Puff center cell.
        int icz, icy, icx;
        GetCellIndices(lon_c, lat_c, z_c, icz, icy, icx);
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        ComputeCellWidth(icz, icy, icx,
                         cell_width_z, cell_width_y,
                         cell_width_x, cell_volume);

        // Concentration for all species.
        for (int s = 0; s < this->Ns; s++)
          Concentration_out(s, icz, icy, icx) =
            max(Concentration_out(s, icz, icy, icx) + GaussianPuffModel
                ->GetPuffQuantity(puff_index, s) / cell_volume, 0.);
      }
    else
      {
        Array<T, 1> total_mass(Ns);
        total_mass = 0.;
        Array<T, 4> PuffConcentration(Ns, Nz, Ny, Nx);
        PuffConcentration = 0.;
        int s, k, j, i, index_z, index_y, index_x;
        T cell_width_z, cell_width_y, cell_width_x, cell_volume;
        T concentration;
        Array<T, 1> quantity(Ns);
        T concentration_erf;
        T x_c, y_c, z_c, lat_c, lon_c;

        // Computing concentration to be added to each cell.
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                index_z = ibz + k;
                index_y = ily + j;
                index_x = ilx + i;

                ComputeCellWidth(index_z, index_y, index_x,
                                 cell_width_z, cell_width_y,
                                 cell_width_x, cell_volume);

                // Cell center.
                lon_c = this->x_min + index_x * this->Delta_x;
                lat_c = this->y_min + index_y * this->Delta_y;
                z_c = this->GridZ3D(index_z);
                LatLonToCartesian(lon_c, lat_c, x_c, y_c);

                // Computing concentrations.
                for (s = 0; s < Ns; s++)
                  {
                    quantity(s) = GaussianPuffModel
                      ->GetPuffQuantity(puff_index, s);
                    if (quantity(s) != 0.)
                      {
                        concentration_erf =
                          GaussianPuffModel->
                          ComputePuffIntegral(puff_index, s, x_c, y_c, z_c,
                                              cell_width_x, cell_width_y,
                                              cell_width_z);
                        PuffConcentration(s, k, j, i) += concentration_erf;

                        total_mass(s) += concentration_erf * cell_volume;
                      }
                  }
              }

        T total_background = 0.0; // YK
        T total_puff = 0.0;
        // Adding concentration to each cell (with mass conservation).
        for (s = 0; s < Ns; s++)
          if (quantity(s) != 0. && total_mass(s) != 0.)
            for (k = 0; k < Nz; k++)
              for (j = 0; j < Ny; j++)
                for (i = 0; i < Nx; i++)
                  {
                    index_z = ibz + k;
                    index_y = ily + j;
                    index_x = ilx + i;
                    ComputeCellWidth(index_z, index_y, index_x,
                                     cell_width_z, cell_width_y,
                                     cell_width_x, cell_volume);
                    concentration = (quantity(s) / total_mass(s))
                      * PuffConcentration(s, k, j, i);

                    Concentration_out(s, index_z, index_y, index_x)
                      = max(Concentration_out(s, index_z, index_y, index_x)
                            + concentration, 0.);
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
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
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
  void PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
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
      for (int k = 0; k < this->Nz; k++)
        {
          index_z = k;
          if (Model.D3("Altitude")(k + 1, index_y, index_x) > height)
            break;
        }
    else
      for (int k = 0; k < this->Nz; k++)
        {
          index_z = k;
          if (this->GridZ3D_interf(k + 1) > height)
            break;
        }
  }


  //! Converts longitude/latitude to cartesian coordinates.
  /*!
    \param lon longitude (degrees).
    \param lat latitude (degrees).
    \param x abscissa (meters).
    \param y ordinate (meters).
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::LatLonToCartesian(T lon, T lat, T& x, T& y)
  {
    const T pi = 3.14159265358979323846264;
    const T earth_radius = 6371229.;
    if (!option_cartesian)
      {
        x = earth_radius * cos(lat * pi / 180.) * (lon * pi / 180.);
        y = earth_radius * (lat * pi / 180.);
      }
    else
      {
        x = lon;
        y = lat;
      }
  }


  //! Converts cartesian coordinates to longitude/latitude.
  /*!
    \param x abscissa (meters).
    \param y ordinate (meters).
    \param lon longitude (degrees).
    \param lat latitude (degrees).
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::CartesianToLatLon(T x, T y, T& lon, T& lat)
  {
    const T pi = 3.14159265358979323846264;
    const T earth_radius = 6371229.;
    if (!option_cartesian)
      {
        lat = (y / earth_radius) * (180. / pi);
        lon = x / (earth_radius * cos(lat * pi / 180.)) * (180. / pi);
      }
    else
      {
        lon = x;
        lat = y;
      }
  }


  //! Returns the concentrations Data.
  /*!
    \return The concentrations Data given by the Eulerian model.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  Data<T, 4>& PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetConcentration()
  {
    return this->Concentration;
  }


  //! Computes concentration at a given point.
  /*!
    \return The concentration at the point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  T  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetConcentration(int species, T z, T y, T x)
  {
    T concentration;
    concentration = Model.GetConcentration(species, z, y, x);
    // Conversion to cartesian coordinates.
    T abscissa, ordinate;
    LatLonToCartesian(x, y, abscissa, ordinate);
    concentration +=
      GaussianPuffModel->GetConcentration(species, z, ordinate , abscissa);
    return concentration;
  }

  /*!
    \return The gas compensation at the point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  T  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetGasCompensation(int species, int z, int y, int x)
  {
    return GasCompensation(species, z, y, x);
  }

  /*!
    \Sets The gas compensation at the point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::SetGasCompensation(int species, int z, int y, int x, T gas_compensation)
  {
    
    GasCompensation(species, z, y, x) = gas_compensation;
  }


  //! Returns the name of the model.
  /*!
  */
   template<class T, class ClassEulerianModel, class ClassLocalModel>
  string  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetModelName() const
  {
    return "PlumeInGrid";
  }

  //! Computes the mean concentration over a given volume.
  /*!
    \return The concentration over the given volume.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  T  PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  ::GetIntegratedConcentration(int species, T z, T y, T x, T lz, T ly, T lx)
  {
    T concentration;
    concentration = Model.GetConcentration(species, z, y, x);
    // Conversion to cartesian coordinates.
    T abscissa, ordinate;
    LatLonToCartesian(x, y, abscissa, ordinate);
    concentration +=
      GaussianPuffModel->GetIntegratedConcentration(species, z, ordinate,
                                                    abscissa, lz, ly, lx);
    return concentration;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDTRANSPORT_CXX
#endif
