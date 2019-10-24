// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMEINGRID_CXX


#include "PlumeInGrid.hxx"


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
  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::PlumeInGrid(string config_file): BaseModel<T>(config_file),
    Model(config_file)
  {
    /*** Pointers to 3D and 2D data ***/

    this->D3_map["MeridionalWind_i"] = &MeridionalWind_i;
    this->D3_map["ZonalWind_i"] = &ZonalWind_i;
    this->D3_map["Temperature_i"] = &Temperature_i;

    this->D2_map["FirstLevelWind_i"] = &FirstLevelWind_i;
    this->D2_map["LowCloudiness_i"] = &LowCloudiness_i;
    this->D2_map["MediumCloudiness_i"] = &MediumCloudiness_i;
    this->D2_map["HighCloudiness_i"] = &HighCloudiness_i;
    this->D2_map["Insolation_i"] = &Insolation_i;
    this->D2_map["FrictionModule_i"] = &FrictionModule_i;
    this->D2_map["BoundaryHeight_i"] = &BoundaryHeight_i;
    this->D2_map["LMO_i"] = &LMO_i;
  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::~PlumeInGrid()
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    BaseModel<T>::ReadConfiguration();

    this->config.SetSection("[domain]");
    this->config.PeekValue("Cartesian", option_cartesian);

    this->config.SetSection("[options]");
    bool option_photolysis_tabulation;
    string scavenging_model;
    this->config.PeekValue("With_deposition", option_deposition);
    this->config.PeekValue("Scavenging_model", scavenging_model);
    this->config.PeekValue("With_chemistry", option_chemistry);
    this->config.PeekValue("With_photolysis", option_photolysis_tabulation);

    option_scavenging = (scavenging_model != "none");
    option_chemistry = (option_chemistry && option_photolysis_tabulation);

    this->config.SetSection("[data]");

    string data_file, emission_file;
    this->config.PeekValue("Data_description", data_file);

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

    /*** Gaussian configuration ***/

    this->config.SetSection("[gaussian]");
    this->config.PeekValue("gaussian_type", gaussian_type);
    this->config.PeekValue("file_gaussian", gaussian_file);
    this->config.PeekValue("With_temporal_profile", option_temporal_profile);
    this->config.PeekValue("Proportional_emission",
                           option_proportional_emission);
    if (option_proportional_emission)
      {
        emission_coefficient.resize(this->Ns);
        emission_coefficient = 0.;
        this->config.Find("Emitted_species");
        const vector<string>& token =  split(this->config.GetLine());
        for (int i = 0; i + 1 < int(token.size()); i += 2)
          {
            int species = BaseModel<T>::GetSpeciesIndex(token[i]);
            emission_coefficient(species) = to_num<T>(token[i + 1]);
          }
      }

    data.SetSection("[plume-in-grid]");
    data.PeekValue("With_interpolation", interpolated);
    if (gaussian_type == "puff")
      {
        if (option_chemistry)
          data.PeekValue("With_chemistry_feedback", option_chemical_feedback);
        data.PeekValue("Horizontal_coefficient", "positive", coefficient_y);
        data.PeekValue("Vertical_coefficient", "positive", coefficient_z);

        data.PeekValue("With_reinjection_time", option_time);
        data.PeekValue("Reinjection_time", "positive", reinjection_time);
      }

    // Injection data.
    string injection_method;

    data.PeekValue("Injection_method", "integrated | column",
                   injection_method);
    if (injection_method == "integrated")
      injection_integrated = 1;
    else if (injection_method == "column")
      injection_integrated = 0;
  }


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    BaseModel<T>::Init();

    /*** Initializations ***/

    Model.Init();
    this->Concentration.Copy(Model.GetConcentration());
    EulerianConcentration.Copy(this->Concentration);
    GaussianConcentration.Copy(this->Concentration);
    GaussianConcentration.SetZero();
    PlumeConcentrationGrid.Copy(this->Concentration);
    PlumeConcentrationGrid.SetZero();

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
          ScavengingBelowCloudCoefficient_i.Copy(Model.D4("ScavengingBelowCloudCoefficient_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"ScavengingCoefficient\".";
      }

    Nr_photolysis = 0;
    vector<string> photolysis_reaction_list;

    if (option_chemistry)
      {
        // Extracts information useful for chemistry from Eulerian model.
        Nr_photolysis = Model.GetNr_photolysis();
        photolysis_reaction_list = Model.GetPhotolysisReactionList();

        // Photolysis rates.
        if (Model.HasField("PhotolysisRate"))
          PhotolysisRate_i.Copy(Model.D4("PhotolysisRate_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"PhotolysisRate\".";

        // Meteorological data.
        Attenuation_i.Copy(Model.D3("Attenuation_i"));
        Pressure_i.Copy(Model.D3("Pressure_i"));
        SpecificHumidity_i.Copy(Model.D3("SpecificHumidity_i"));
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

    InitGaussianModel(photolysis_reaction_list);

    /*** Latitude and longitude ***/

    latitude.resize(this->Ny, this->Nx);
    longitude.resize(this->Ny, this->Nx);
    T delta_x = this->GetDelta_x();
    T delta_y = this->GetDelta_y();
    T x_min = this->GetX_min();
    T y_min = this->GetY_min();
    for (int j = 0; j < this->Ny; j++)
      for (int i = 0; i < this->Nx; i++)
        if (option_cartesian)
          CartesianToLatLon(x_min + (i + 0.5) * delta_x,
                            y_min + (j + 0.5) * delta_y,
                            longitude(j, i), latitude(j, i));
        else
          {
            longitude(j, i) = x_min + (i + 0.5) * delta_x;
            latitude(j, i) = y_min + (j + 0.5) * delta_y;
          }
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */  template<class T, class ClassEulerianModel, class ClassLocalModel>
   void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
       ScavengingBelowCloudCoefficient_i.Copy(Model.D4("ScavengingBelowCloudCoefficient_i"));

     if (option_chemistry)
       {
         Attenuation_i.Copy(Model.D3("Attenuation_i"));
         Pressure_i.Copy(Model.D3("Pressure_i"));
         SpecificHumidity_i.Copy(Model.D3("SpecificHumidity_i"));

         // Photolysis rates.
         PhotolysisRate_i.Copy(Model.D4("PhotolysisRate_i"));
       }
   }


  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    if (gaussian_type == "puff")
      {
        this->Concentration.Copy(Model.GetConcentration());


        int Npuff;
        int index_x, index_y, index_z, puff_index;
        bool puff_transfer;
        bool isday;

        /*** Tests if puff has to be transfered to Eulerian model. ***/

        Npuff = GaussianModel->GetGaussianSourceCount();
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            GaussianModel->
              GetGaussianPosition(puff_index, x_c, y_c, z_c, puff_distance,
                                  puff_time);

            // Conversion to longitude/latitude (degrees).
            CartesianToLatLon(x_c, y_c, lon_c, lat_c);

            // Night or day.
            isday = IsDay(lon_c, lat_c, GaussianModel->
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
            GaussianModel->GetGaussianSigma(puff_index, sigma_x,
                                            sigma_y, sigma_z);
            puff_transfer = 0;

            // Tests if puff has reached the end of the domain.
            if (index_x == 0 || index_x == this->Nx - 1
                || index_y == 0 || index_y == this->Ny - 1
                || index_z == this->Nz)
              GaussianModel->ErasePuff(puff_index);

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
                  PuffIntegratedTransfer(puff_index,
                                         Model.GetConcentration());
                else
                  for (int s = 0; s < this->Ns; s++)
                    GaussianTransfer(GaussianModel->
                                     GetPuffQuantity(puff_index, s),
                                     sigma_z, s, z_c, lat_c, lon_c,
                                     isday, Model.GetConcentration());
                GaussianModel->ErasePuff(puff_index);
              }
          }

        /*** Inner time-loop for Gaussian models. ***/

        Array<T, 4> PuffConcentration(this->Ns, this->Nz, this->Ny, this->Nx);
        PuffConcentration = 0.;
        int l, k, j, i, s;
        for (j = 0; j < Nt_local; j++)
          {
            GaussianModel->InitStep();
            Npuff = GaussianModel->GetGaussianSourceCount();
            Array<int, 2> PuffCellList(Npuff, 3);

            // Loop on puffs to get the meteorological data.
            for (puff_index = 0; puff_index < Npuff; puff_index++)
              {
                // Puff center coordinates.
                T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
                GaussianModel->
                  GetGaussianPosition(puff_index, x_c, y_c, z_c,
                                      puff_distance, puff_time);

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
            GaussianModel->ComputePlumeRise();

            // Advection and diffusion for puffs.
            GaussianModel->Advection();
            GaussianModel->Diffusion();
            GaussianModel->ComputeLossFactor();

            // Loop on puffs to update data if necessary.
            for (puff_index = 0; puff_index < Npuff; puff_index++)
              {
                bool has_changed = false;
                // Puff center coordinates.
                T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
                GaussianModel->
                  GetGaussianPosition(puff_index, x_c, y_c, z_c,
                                      puff_distance, puff_time);

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

            // Chemistry for puffs.
            if (option_chemistry)
              if (option_chemical_feedback)
                ComputePuffChemistry(PuffCellList, PuffConcentration);
              else
                GaussianModel->ComputeChemistry();
            GaussianModel->AddTime(delta_t_local);
          }

        // Eulerian model.
        Model.Forward();

        if (option_chemistry)
          if (option_chemical_feedback)
            Model.GetConcentration().GetArray() += PuffConcentration;

        this->Concentration.Copy(Model.GetConcentration());


        /* Adding all puffs to concentration Data in case concentrations are
           saved on the domain.*/

        GaussianModel->SubtractTime(delta_t_local);
        Npuff = GaussianModel->GetGaussianSourceCount();
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            T puff_release_time = GaussianModel->
              GetPuffReleaseTime(puff_index);
            if (GaussianModel->GetCurrentTime() >= puff_release_time)
              {
                // Puff center coordinates.
                T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
                GaussianModel->
                  GetGaussianPosition(puff_index, x_c, y_c, z_c,
                                      puff_distance, puff_time);

                // Conversion to longitude/latitude (degrees).
                CartesianToLatLon(x_c, y_c, lon_c, lat_c);
                T sigma_x, sigma_y, sigma_z;
                GaussianModel->GetGaussianSigma(puff_index, sigma_x,
                                                sigma_y, sigma_z);

                isday = IsDay(lon_c, lat_c, GaussianModel->
                              GetCurrentDate());
                // Adding puff concentration to the Concentration Data.
                if (injection_integrated)
                  PuffIntegratedTransfer(puff_index,
                                         this->Concentration);
                else
                  for (int s = 0; s < this->Ns; s++)
                    GaussianTransfer(GaussianModel->
                                     GetPuffQuantity(puff_index, s),
                                     sigma_z, s, z_c, lat_c, lon_c,
                                     isday, this->Concentration);
              }
          }
        GaussianModel->AddTime(delta_t_local);
      }
    else if (gaussian_type == "plume-line")
      {
        Model.GetConcentration().GetArray() = EulerianConcentration.GetArray()
          + GaussianConcentration.GetArray();

        // Resets to zero negative concentrations.
        // Negative concentration can occur because Gaussian chemistry scheme
        // can produce a deficit of O3.
        Array<T, 4>& concentration_array = Model.GetConcentration().GetArray();
        for (typename Array<T, 4>::iterator it = concentration_array.begin();
             it != concentration_array.end(); ++it)
          {
            T& concentration = *it;
            if (concentration < 0.)
              concentration = 0.;
          }

        if (this->step * int(delta_t_eulerian) % int(delta_t_local) == 0)
          {
            if (option_temporal_profile)
              {
                T profile_coeficient = temporal_profile_coefficients(this->step);
                GaussianModel->MultiplySourcesRate(profile_coeficient);
              }
            PlumeConcentrationGrid.SetZero();

            GaussianConcentration.SetZero();
            GaussianModel->SetZeroConcentrationList();

            Array<T, 2> GaussianCoordinateList = GaussianModel
              ->GetCoordinateList();

            // Set Background Concentration.
            Data<T, 4> BackgroundConcentration(this->Ns, this->Nz, this->Ny,
                                               this->Nx);
            Data<T, 2> BackgroundConcentration_list(GaussianCoordinateList
                                                    .extent(0), this->Ns);

            if (rank == 0)
              {
                BackgroundConcentration.Copy(Model.GetConcentration());

                // Set Background Concentration for the coordinate list.
                int index_z, index_y, index_x;
                T lat, lon;
                for (int i = 0; i < GaussianCoordinateList.extent(0); i++)
                  for (int j = 0; j < this->Ns; j++)
                    {
                      CartesianToLatLon(GaussianCoordinateList(i, 2),
                                        GaussianCoordinateList(i, 1),
                                        lon, lat);

                      GetCellIndices(lon, lat, GaussianCoordinateList(i, 0),
                                     index_z, index_y, index_x);

                      BackgroundConcentration_list(i, j)
                        = BackgroundConcentration(j, index_z, index_y, index_x);
                    }

                Array<T, 3> SpecificHumidity_grid(this->Nz, this->Ny, this->Nx);
                Array<T, 3> Temperature_grid(this->Nz, this->Ny, this->Nx);
                Array<T, 3> Pressure_grid(this->Nz, this->Ny, this->Nx);
                Array<int, 3> Stability_grid(this->Nz, this->Ny, this->Nx);
                Array<int, 3> Rural_grid(this->Nz, this->Ny, this->Nx);
                Array<T, 3> Wind_grid(this->Nz, this->Ny, this->Nx);
                Array<T, 3> Attenuation_grid(this->Nz, this->Ny, this->Nx);

                for (int k = 0; k < this->Nz; k++)
                  for (int j = 0; j < this->Ny; j++)
                    for (int i = 0; i < this->Nx; i++)
                      GetGaussianMeteo(this->GridZ4D(k), this->GridY4D(j),
                                       this->GridX4D(i),
                                       SpecificHumidity_grid(k, j, i),
                                       Temperature_grid(k, j, i),
                                       Pressure_grid(k, j, i),
                                       Stability_grid(k, j, i),
                                       Rural_grid(k, j, i),
                                       Wind_grid(k, j, i),
                                       Attenuation_grid(k, j, i));
                GaussianModel->
                  SetChemistryMeteorologicalParameters(SpecificHumidity_grid,
                                                       Temperature_grid,
                                                       Pressure_grid,
                                                       Stability_grid,
                                                       Rural_grid, Wind_grid,
                                                       Attenuation_grid);

                Array<T, 1>
                  SpecificHumidity_list(GaussianCoordinateList.extent(0));
                Array<T, 1>
                  Temperature_list(GaussianCoordinateList.extent(0));
                Array<T, 1> Pressure_list(GaussianCoordinateList.extent(0));
                Array<int, 1> Stability_list(GaussianCoordinateList.
                                             extent(0));
                Array<int, 1> Rural_list(GaussianCoordinateList.extent(0));
                Array<T, 1> Wind_list(GaussianCoordinateList.extent(0));
                Array<T, 1>
                  Attenuation_list(GaussianCoordinateList.extent(0));

                for (int i = 0; i < GaussianCoordinateList.extent(0); i++)
                  {
                    CartesianToLatLon(GaussianCoordinateList(i, 2)
                                      , GaussianCoordinateList(i, 1),
                                      lon, lat);
                    GetGaussianMeteo(GaussianCoordinateList(0, i), lat, lon,
                                     SpecificHumidity_list(i),
                                     Temperature_list(i),
                                     Pressure_list(i),
                                     Stability_list(i),
                                     Rural_list(i), Wind_list(i),
                                     Attenuation_list(i));
                  }

                GaussianModel->
                  SetChemistryMeteorologicalParameters(SpecificHumidity_list,
                                                       Temperature_list,
                                                       Pressure_list,
                                                       Stability_list,
                                                       Rural_list, Wind_list,
                                                       Attenuation_list);
              }

            /*** Inner time-loop for Gaussian models. ***/

            int plume_index, Nplume;
            for (plume_index = 0;
                 plume_index < GaussianModel->GetGaussianSourceCount();
                 plume_index++)
              {
                GaussianModel->RestrictComputationToSource(plume_index);

                // Updating the plume meteorological data.
                UpdateMeteo(plume_index);

                GaussianModel->InitStep();

                // Compute source contribution for the Eulerian grid.
                ComputeGaussianSourceContribution(plume_index,
                                                  GaussianConcentration);

                // Compute concentration for the list of ouput point.
                GaussianModel->Compute_list(plume_index);

                // Loops for all points sources added for the
                // line source / point source discretization.
                Nplume = GaussianModel->GetGaussianSourceCount();
                int inputSourceCount = GaussianModel->GetInputSourceCount();

                if (inputSourceCount < Nplume)
                  {
                    GaussianModel->SaveCurrentPlume();
                    for (int i = inputSourceCount; i < Nplume; i++)
                      {
                        GaussianModel->RestrictComputationToSource(i);
                        ComputeGaussianSourceContribution(i,
                                                          GaussianConcentration);
                      }
                    GaussianModel->RestoreCurrentPlume();

                    GaussianModel->FillPointSourceList(inputSourceCount,
                                                       Nplume);
                    GaussianModel->Compute_list(inputSourceCount);
                    GaussianModel->EmptyPointSourceList();
                  }

                // Compute Grid concentration that will be saved.
                int species = GaussianModel->GetSpecies(plume_index);
                T x, y, Concentration_temp;
                GaussianModel->FillPointSourceList(plume_index);
                if (inputSourceCount < Nplume)
                  GaussianModel->FillPointSourceList(inputSourceCount,
                                                     Nplume);

                for (int k = 0; k < this->Nz; k++)
                  for (int j = 0; j < this->Ny; j++)
                    for (int i = 0; i < this->Nx; i++)
                      {
                        LatLonToCartesian(this->GridX4D(i),
                                          this->GridY4D(j), x, y);

                        PlumeConcentrationGrid(species, k, j, i) += GaussianModel
                          ->ComputeGaussianConcentration(species,
                                                         this->GridZ4D(k)
                                                         , y, x);
                      }
                GaussianModel->EmptyPointSourceList();

                GaussianModel->SaveCurrentPlume();
                if (inputSourceCount < Nplume)
                  {
                    GaussianModel->EraseDiscretizedSource();
                    GaussianModel->RestoreRate();
                  }
              }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
            int code;
            MPI_Datatype Array_1;
            MPI_Type_contiguous(this->Ns * this->Nz * this->Ny * this->Nx,
                                MPI_DOUBLE, &Array_1);
            MPI_Type_commit(&Array_1);

            T tmp_1[this->Ns][this->Nz][this->Ny][this->Nx];
            T tmp_2[Nproc][this->Ns][this->Nz][this->Ny][this->Nx];
            MPI_Barrier(MPI_COMM_WORLD);

            // PlumeConcentrationGrid.
            if (rank > 0)
              for (int s = 0; s < this->Ns; s++)
                for (int k = 0; k < this->Nz; k++)
                  for (int j = 0; j < this->Ny; j++)
                    for (int i = 0; i < this->Nx; i++)
                      tmp_1[s][k][j][i] = PlumeConcentrationGrid()(s, k, j, i);

            code = MPI_Gather(tmp_1, 1, Array_1 , tmp_2, 1, Array_1 , 0,
                              MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (code != MPI_SUCCESS)
              throw "MPI communication error";

            if (rank == 0)
              for (int np = 1; np < Nproc; np++)
                for (int s = 0; s < this->Ns; s++)
                  for (int k = 0; k < this->Nz; k++)
                    for (int j = 0; j < this->Ny; j++)
                      for (int i = 0; i < this->Nx; i++)
                        PlumeConcentrationGrid(s, k, j, i) +=
                          tmp_2[np][s][k][j][i];

            // GaussianConcentration.
            if (rank > 0)
              for (int s = 0; s < this->Ns; s++)
                for (int k = 0; k < this->Nz; k++)
                  for (int j = 0; j < this->Ny; j++)
                    for (int i = 0; i < this->Nx; i++)
                      tmp_1[s][k][j][i] = GaussianConcentration()(s, k, j, i);

            code = MPI_Gather(tmp_1, 1, Array_1 , tmp_2, 1, Array_1 , 0,
                              MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (code != MPI_SUCCESS)
              throw "MPI communication error";

            if (rank == 0)
              for (int np = 1; np < Nproc; np++)
                for (int s = 0; s < this->Ns; s++)
                  for (int k = 0; k < this->Nz; k++)
                    for (int j = 0; j < this->Ny; j++)
                      for (int i = 0; i < this->Nx; i++)
                        GaussianConcentration(s, k, j, i) +=
                          tmp_2[np][s][k][j][i];

            // GaussianModel->Concentration_list_.
            Array<T, 2> GaussianConcentrationList(GaussianCoordinateList.
                                                  extent(0), this->Ns);
            GaussianConcentrationList = GaussianModel->ConcentrationList();
            int N = GaussianCoordinateList.extent(0);

            MPI_Datatype Array_2;
            MPI_Type_contiguous(N * this->Ns, MPI_DOUBLE, &Array_2);
            MPI_Type_commit(&Array_2);

            T tmp_3[N][this->Ns];
            T tmp_4[Nproc][N][this->Ns];

            if (rank > 0)
              for (int s = 0; s < this->Ns; s++)
                for (int n = 0; n < N; n++)
                  tmp_3[n][s] = GaussianConcentrationList(n, s);

            code = MPI_Gather(tmp_3, 1, Array_2 , tmp_4, 1, Array_2 , 0,
                              MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (code != MPI_SUCCESS)
              throw "MPI communication error";

            if (rank == 0)
              for (int nproc = 1; nproc < Nproc; nproc++)
                for (int s = 0; s < this->Ns; s++)
                  for (int n = 0; n < N; n++)
                    GaussianConcentrationList(n, s) += tmp_4[nproc][n][s];

            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0)
              GaussianModel->ConcentrationList() = GaussianConcentrationList;
#endif
            if (rank == 0)
              {
                // Computes emission for all species.
                if (option_proportional_emission)
                  {
                    int s_ref = GaussianModel->GetSpecies(0);
                    for (int s = 0; s < this->Ns; ++s)
                      {
                        T emission_coef = emission_coefficient(s);
                        if (emission_coef <= 0.)
                          continue;
                        Array<T, 4>& PlumeConcentrationArray =
                          PlumeConcentrationGrid.GetArray();
                        PlumeConcentrationArray(s) =
                          PlumeConcentrationArray(s_ref) * emission_coef;

                        Array<T, 4>& GaussianConcentrationArray =
                          GaussianConcentration.GetArray();
                        GaussianConcentrationArray(s) =
                          GaussianConcentrationArray(s_ref) * emission_coef;

                        Array<T, 2>& concentration_list =
                          GaussianModel->ConcentrationList();
                        concentration_list(Range::all(), s) =
                          concentration_list(Range::all(), s_ref)
                          * emission_coef;
                      }
                  }

                // Chemistry for plumes.
                if (option_chemistry)
                  {
                    Array<T, 3> temperature(this->Nz, this->Ny, this->Nx);
                    temperature = Temperature_i.GetArray();
                    GaussianModel->
                      SetModelParameter("temperature",
                                        mean(temperature(0, Range::all(),
                                                         Range::all())));
                    GaussianModel->Chemistry_Plume(this->GetCurrentDate(),
                                                   latitude, longitude,
                                                   BackgroundConcentration,
                                                   GaussianConcentration);
                    GaussianModel->Chemistry_Plume(this->GetCurrentDate(),
                                                   latitude, longitude,
                                                   BackgroundConcentration,
                                                   PlumeConcentrationGrid);
                    GaussianModel->
                      Chemistry_Plume_list(this->GetCurrentDate(),
                                           BackgroundConcentration_list);
                  }
              }

            if (option_temporal_profile)
              {
                T profile_coefficient = temporal_profile_coefficients(this->step);
                GaussianModel->MultiplySourcesRate(1. / profile_coefficient);
              }

            GaussianModel->AddTime(delta_t_local);
          }

        // Eulerian model.
        Model.Forward();

        if (rank == 0)
          {
            EulerianConcentration.Copy(Model.GetConcentration());

            Data<T, 4> concentration_data = Model.GetConcentration();
            Array<T, 4>& concentration_array = concentration_data.GetArray();
            concentration_array += PlumeConcentrationGrid.GetArray();

            // Check if there are some negative value. This can occur because
            // Gaussian chemistry scheme can produce a deficit of O3.
            for (typename Array<T, 4>::iterator it =
                   concentration_array.begin();
                 it != concentration_array.end(); ++it)
              {
                T& concentration = *it;
                if (concentration < 0.)
                  concentration = 0.;
              }

            // Transfer to the Eulerien model.
            this->Concentration.Copy(concentration_data);
          }
      }

    this->AddTime(delta_t_eulerian);
    this->step++;
  }


  //! Updates the gaussian meteorological data.
  /*! Updates the gaussian meteorological data taking them from the center
    cell.
    \param gaussian_index index of the gaussian source.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateMeteo(int gaussian_index)
  {
    // Source center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, gaussian_distance, gaussian_time;
    GaussianModel->
      GetGaussianPosition(gaussian_index, x_c, y_c, z_c, gaussian_distance,
                          gaussian_time);

    // Conversion to longitude/latitude (degrees).
    CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Night or day.
    bool isday = IsDay(lon_c, lat_c, GaussianModel->GetCurrentDate());

    // Extracting meteorological data.
    bool rural;
    string stability;
    map<string, T> meteorological_data;
    bool option_similarity = GaussianModel->WithSimilarity();

    ExtractMeteo(z_c, lat_c, lon_c, isday, option_similarity,
                 meteorological_data, stability, rural);

    if (!option_similarity)
      GaussianModel->
        SetGaussianMeteo(gaussian_index,
                         meteorological_data["temperature"],
                         meteorological_data["wind_angle"],
                         meteorological_data["wind"],
                         stability, lon_c, lat_c, isday, rural,
                         meteorological_data["boundaryheight"]);
    else
      GaussianModel->
        SetGaussianMeteo(gaussian_index,
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

    if (GaussianModel->IsOption("plume_rise_breakup")
        && !option_similarity)
      {
        GaussianModel->
          SetModelParameter("friction_velocity",
                            meteorological_data["frictionmodule"]);
        GaussianModel->
          SetModelParameter("convective_velocity",
                            meteorological_data["convectivevelocity"]);
      }

    // Data depending on the species.
    Array<T, 1> deposition_velocity;
    deposition_velocity.resize(this->Ns);
    deposition_velocity = 0.;

    Array<T, 1> scavenging_coefficient;
    scavenging_coefficient.resize(this->Ns);
    scavenging_coefficient = 0.;

    Array<T, 1> photolysis_rate;
    photolysis_rate.resize(Nr_photolysis);

    Array<T, 1> background_concentration;
    background_concentration.resize(this->Ns);

    ExtractSpeciesData(z_c, lat_c, lon_c,
                       deposition_velocity,
                       scavenging_coefficient,
                       photolysis_rate,
                       background_concentration);

    if (option_deposition)
      GaussianModel
        ->SetGaussianDepositionVelocity(gaussian_index,
                                        deposition_velocity);

    if (option_scavenging)
      GaussianModel
        ->SetGaussianScavengingCoefficient(gaussian_index,
                                           scavenging_coefficient);

    if (gaussian_type == "puff" && option_chemistry)
      {
        GaussianModel->SetPuffAdditionalMeteo
          (gaussian_index,
           meteorological_data["attenuation"],
           meteorological_data["pressure"],
           meteorological_data["specific_humidity"]);

        GaussianModel
          ->SetGaussianPhotolysisRate(gaussian_index,
                                      photolysis_rate);
        GaussianModel->SetGaussianBackgroundConcentration
          (gaussian_index, background_concentration);
      }
    else if (gaussian_type == "plume-line")
      GaussianModel->InitCorrectionCoefficients();
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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

    // Interpolation for winds on gaussian source center.
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

    if (!option_similarity)
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

    if (option_chemistry)
      {
        T attenuation, pressure, specific_humidity;
        if (!interpolated)
          {
            attenuation = Attenuation_i(index_z, index_y, index_x);
            pressure = Pressure_i(index_z, index_y, index_x);
            specific_humidity =
              SpecificHumidity_i(index_z, index_y, index_x);
          }
        else
          {
            LinearInterpolationPoint(Attenuation_i, Coord3D, attenuation);
            LinearInterpolationPoint(Pressure_i, Coord3D, pressure);
            LinearInterpolationPoint(SpecificHumidity_i, Coord3D,
                                     specific_humidity);
          }
        met_data["attenuation"] = attenuation;
        met_data["pressure"] = pressure;
        met_data["specific_humidity"] = specific_humidity;
      }
  }


  //! Extracts data that depend on meteorology and species (photolysis rates).
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractSpeciesData(T height, T lat, T lon,
                       Array<T, 1>& deposition_velocity,
                       Array<T, 1>& scavenging_coefficient,
                       Array<T, 1>& photolysis_rate,
                       Array<T, 1>& background_concentration)
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
                    ScavengingBelowCloudCoefficient_i(scav_s, index_z,
                                                      index_y, index_x);
                else
                  LinearInterpolationPoint(ScavengingBelowCloudCoefficient_i, Coord4D,
                                           scavenging_coefficient(s));
              }
          }
      }

    if (option_chemistry)
      {
        Array<T, 1> Coord4D(4);
        Coord4D(0) = 0;
        Coord4D(1) = height;
        Coord4D(2) = lat;
        Coord4D(3) = lon;
        // Photolysis rates.
        for (r = 0; r < Nr_photolysis; r++)
          {
            Coord4D(0) = r;
            if (!interpolated)
              photolysis_rate(r) =
                PhotolysisRate_i(r, index_z, index_y, index_x);
            else
              LinearInterpolationPoint(PhotolysisRate_i, Coord4D,
                                       photolysis_rate(r));
          }

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
  }


  //! Transfers a given quantity to the Eulerian model.
  /*! Transfers a given quantity to the Eulerian model on one vertical column.
    \param quantity mass.
    \param sigma_z vertical standard deviation of the gaussian source.
    \param index_s index of the source species.
    \param z_c source center coordinate along z.
    \param lat_c source center coordinate along y.
    \param lon_c source center coordinate along x.
    \param Concentration_out concentrations Data to be modified.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::GaussianTransfer(T quantity, T sigmaz, int s, T z, T lat, T lon,
                     bool isday, Data<T, 4>& Concentration_out)
  {
    int ibz, iby, ibx, itz, ity, itx;
    // Computes source vertical extent.
    T source_bottom = max(z - (coefficient_z * sigmaz) / 2, 0.);
    GetCellIndices(lon, lat, source_bottom, ibz, iby, ibx);
    T source_top = z + (coefficient_z * sigmaz) / 2;
    T boundary_height = BoundaryHeight_i(iby, ibx);

    // Takes into account the inversion height.
    if (isday && z < boundary_height)
      source_top = min(source_top, boundary_height);
    if (isday && z > boundary_height)
      source_bottom = max(source_bottom, boundary_height);
    GetCellIndices(lon, lat, source_bottom, ibz, iby, ibx);
    GetCellIndices(lon, lat, source_top, itz, ity, itx);

    // Number of cells vertically covered by the source.
    int Ncell = (itz - ibz) + 1;

    // Computing vertical extent of cells where source will be transfered.
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

    // Adding source concentration to the Eulerian concentration.
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::PuffIntegratedTransfer(int puff_index, Data<T, 4>& Concentration_out)
  {
    T sigma_x, sigma_y, sigma_z;
    GaussianModel->GetGaussianSigma(puff_index, sigma_x,
                                    sigma_y, sigma_z);

    // Puff center coordinates.
    T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
    GaussianModel->
      GetGaussianPosition(puff_index, x_c, y_c, z_c, puff_distance,
                          puff_time);
    // Conversion to longitude/latitude (degrees).
    CartesianToLatLon(x_c, y_c, lon_c, lat_c);

    // Computes puff vertical extent.
    int ibz, iby, ibx, itz, ity, itx;
    T puff_bottom = max(z_c - (coefficient_z * sigma_z) / 2, 0.);
    GetCellIndices(lon_c, lat_c, puff_bottom, ibz, iby, ibx);
    T puff_top = z_c + (coefficient_z * sigma_z) / 2;
    GetCellIndices(lon_c, lat_c, puff_top, itz, ity, itx);

    // Computes puff horizontal extent.
    int ilyz, ily, ilyx, iryz, iry, iryx;
    T puff_left = y_c - (coefficient_y * sigma_y) / 2;
    T lon_left, lat_left;
    CartesianToLatLon(x_c, puff_left, lon_left, lat_left);
    GetCellIndices(lon_c, lat_left, z_c, ilyz, ily, ilyx);
    T puff_right = y_c + (coefficient_y * sigma_y) / 2;
    T lon_right, lat_right;
    CartesianToLatLon(x_c, puff_right, lon_right, lat_right);
    GetCellIndices(lon_c, lat_right, z_c, iryz, iry, iryx);

    // Computes puff horizontal extent.
    int ilxz, ilxy, ilx, irxz, irxy, irx;
    puff_left = x_c - (coefficient_y * sigma_x) / 2;
    CartesianToLatLon(puff_left, y_c, lon_left, lat_left);
    GetCellIndices(lon_left, lat_c, z_c, ilxz, ilxy, ilx);
    puff_right = x_c + (coefficient_y * sigma_x) / 2;
    CartesianToLatLon(puff_right, y_c, lon_right, lat_right);
    GetCellIndices(lon_right, lat_c, z_c, irxz, irxy, irx);

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
            max(Concentration_out(s, icz, icy, icx) + GaussianModel
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
                    quantity(s) = GaussianModel
                      ->GetPuffQuantity(puff_index, s);
                    if (quantity(s) != 0.)
                      {
                        concentration_erf =
                          GaussianModel->
                          ComputePuffIntegral(puff_index, s, x_c, y_c, z_c,
                                              cell_width_x, cell_width_y,
                                              cell_width_z);
                        PuffConcentration(s, k, j, i) += concentration_erf;

                        total_mass(s) += concentration_erf * cell_volume;
                      }
                  }
              }

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


  /*! Computes the Gaussian puff chemistry and the subsequent feedback to
    background cells concentrations (in case there is feedback with
    chemistry).*/
  /*!
    \param PuffCellList list of coordinates of the cells containing puffs.
    \param PuffConcentration Concentration perturbation to be added
    to each cell after puff chemistry.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::ComputePuffChemistry(Array<int, 2> PuffCellList,
                         Array<T, 4>& PuffConcentration)
  {
    int l, k, i, s, puff_index;
    int index_x = 0;
    int index_y = 0;
    int index_z = 0;
    int Npuff = GaussianModel->GetGaussianSourceCount();
    // Loop on cells to take them into account for background chemistry.
    bool is_in_list;
    list<int> pufflist_tmp;
    vector<list<int> > PuffList;
    vector<T> PuffCellVolume;
    Array<int, 1> Coord(3);
    vector<Array<int, 1> > PuffCellCoordinates;
    Array<vector<T>, 1 > PuffCellConcentration;
    PuffCellConcentration.resize(this->Ns);
    int Ncell = 0;
    for (l = 0; l < this->Nz; l++)
      for (k = 0; k < this->Ny; k++)
        for (i = 0; i < this->Nx; i++)
          {
            is_in_list = false;
            pufflist_tmp.clear();
            // Is there at least one puff in the cell?
            for (puff_index = 0; puff_index < Npuff; puff_index++)
              if (PuffCellList(puff_index, 0) == l
                  && PuffCellList(puff_index, 1) == k
                  && PuffCellList(puff_index, 2) == i)
                {
                  index_z = l;
                  index_y = k;
                  index_x = i;
                  is_in_list = true;
                  pufflist_tmp.push_back(puff_index);
                }
            if (is_in_list)
              {
                Ncell++;
                // Computing cell volume.
                T cell_width_z, cell_width_y, cell_width_x, cell_volume;
                ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
                                 cell_width_y, cell_width_x, cell_volume);
                PuffCellVolume.push_back(cell_volume);
                Coord(0) = index_z;
                Coord(1) = index_y;
                Coord(2) = index_x;
                PuffCellCoordinates.push_back(Coord);
                for (s = 0; s < this->Ns; s++)
                  PuffCellConcentration(s).push_back
                    (Model.GetConcentration()(s, index_z, index_y, index_x));
                PuffList.push_back(pufflist_tmp);
              }
          }
    GaussianModel->ComputeChemistry(PuffList, PuffCellVolume,
                                    PuffCellConcentration);
    // Adding the remaining concentrations in a buffer array.
    for (i = 0; i < Ncell; i++)
      {
        index_z = PuffCellCoordinates[i](0);
        index_y = PuffCellCoordinates[i](1);
        index_x = PuffCellCoordinates[i](2);

        for (s = 0; s < this->Ns; s++)
          {
            PuffConcentration(s, index_z, index_y, index_x)
              += PuffCellConcentration(s)[i];
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
  void PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
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
  Data<T, 4>& PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::GetConcentration()
  {
    return this->Concentration;
  }


  //! Computes concentration at a given point.
  /*!
    \return The concentration at the point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  T  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::GetConcentration(int species, T z, T y, T x)
  {
    T concentration = 0.;

    // Conversion to cartesian coordinates.
    T abscissa, ordinate, z_grid = z;
    int index_x, index_y, index_z;
    GetCellIndices(x, y, z, index_z, index_y, index_x);
    LatLonToCartesian(x, y, abscissa, ordinate);
    if (z < this->GridZ4D(0))
      z_grid = this->GridZ4D(0);
    concentration = Model.GetConcentration(species, z_grid, y, x);

    concentration +=
      GaussianModel->GetConcentration(species, z, ordinate , abscissa);
    return concentration;
  }


  //! Computes the mean concentration over a given volume.
  /*!
    \return The concentration over the given volume.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  T  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::GetIntegratedConcentration(int species, T z, T y, T x, T lz, T ly, T lx)
  {
    T concentration = 0.;

    // Conversion to cartesian coordinates.
    T abscissa, ordinate;
    LatLonToCartesian(x, y, abscissa, ordinate);

    if (gaussian_type == "puff")
      concentration += Model.GetConcentration(species, z, y, x);
    else if (gaussian_type == "plume-line")
      concentration += EulerianConcentration(species, z, y, x);

    concentration += GaussianModel->GetIntegratedConcentration(species, z,
                                                               ordinate,
                                                               abscissa,
                                                               lz, ly, lx);
    return concentration;
  }


  //! Initialization of the Gaussian Model.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::InitGaussianModel(const vector<string>& photolysis_reaction_list)
  {
    if (gaussian_type == "puff")
      {
        // Gaussian model initialization.
        GaussianModel = new ClassLocalModel(gaussian_file);
        GaussianModel->Init();
        GaussianModel->SetDateMin(this->Date_min);
        delta_t_local = min(GaussianModel->GetDelta_t(), delta_t_eulerian);
        GaussianModel->SetTimeStep(delta_t_local);
        Nt_local = int(delta_t_eulerian / delta_t_local);
        GaussianModel->SetModelParameter("Nt", this->Nt * Nt_local);

        // Is there puff deposition?
        if (!GaussianModel->WithDeposition())
          option_deposition = false;
        else if (!option_deposition)
          throw "Deposition cannot be activated in the GaussianPuff model "
            "without deposition in Eulerian model.";

        // Is there puff scavenging?
        if (!GaussianModel->WithScavenging())
          option_scavenging = false;
        else if (!option_scavenging)
          throw "Scavenging cannot be activated in the GaussianPuff model "
            "without scavenging in the Eulerian model.";

        // Is there puff chemistry?
        if (!GaussianModel->WithChemistry())
          option_chemistry = false;
        else if (!option_chemistry)
          throw "Chemistry cannot be activated in the GaussianPuff model "
            "without chemistry and photolysis in the Eulerian model.";

        if (option_chemistry)
          GaussianModel->InitPhotolysis(Nr_photolysis,
                                        photolysis_reaction_list);

        // Sources initialization in Gaussian puff model.
        GaussianModel->InitSource();

        // Conversion of sources coordinates from (latitude,longitude)
        // to Cartesian system.
        Array<T, 2> source_coordinates;
        GaussianModel->GetSourceCoordinates(source_coordinates);
        int Nemis = source_coordinates.extent(0);

        // Array containing the sources Cartesian coordinates.
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
        GaussianModel->SetSourceCoordinates(source_coordinates_tmp);

        // Clears the list of puffs and initializes positional variables.
        GaussianModel->InitPosition();
      }
    else if (gaussian_type == "plume-line")
      {
        // Gaussian Model creation.
        GaussianModel = new ClassLocalModel(gaussian_file);

        // Set option_cartesian.
        GaussianModel->SetOption("cartesian", option_cartesian);

        // Use finite length plumes.
        GaussianModel->SetOption("infinite_plume", false);

        // Gaussian model initialization.
        GaussianModel->Init();
        delta_t_local = int(GaussianModel->GetTimeStep());
        if (delta_t_local < delta_t_eulerian)
          delta_t_local = delta_t_eulerian;
        else
          delta_t_local = int(delta_t_local / delta_t_eulerian + 0.5)
            * int(delta_t_eulerian);

        GaussianModel->SetDateMin(this->Date_min);
        // Unlike with puff model, 'delta_t_eulerian' is the time step of the
        // Gaussian model. delta_t_local is the time step at which gaussian
        // dispersion is updated.
        GaussianModel->SetTimeStep(delta_t_eulerian);

        // If option_proportional_emission is set to true, only one species
        // is specified in the emission file.
        // Aditional species must be specified with a proportional coefficient.
        if (option_proportional_emission)
          {
            int s = GaussianModel->GetSpecies(0);
            int N = GaussianModel->GetGaussianSourceCount();

            // Loops over all sources.
            for (int plume_index = 1; plume_index < N; ++plume_index)
              if (s != GaussianModel->GetSpecies(plume_index))
                throw string("When \"option_proportional_emission\" is set to ")
                  + "yes, only one species is specified in the emission "
                  + "file. Aditional species must be specified with a proportional "
                  + "coefficient.";
          }

        // Is there plume deposition?
        if (!GaussianModel->WithDeposition())
          option_deposition = false;
        else if (!option_deposition)
          throw string("Depostion cannot be activated in GaussianPlume ")
            + "model without deposition "
            + "in Eulerian model.";

        // Is there plume scavenging?
        if (!GaussianModel->WithScavenging())
          option_scavenging = false;
        else if (!option_scavenging)
          throw string("Scavenging cannot be activated in GaussianPlume ")
            + "model without scavenging "
            + "in Eulerian model.";

        // Is there plume chemistry?
        if (!GaussianModel->WithChemistry())
          option_chemistry = false;
        else if (!option_chemistry)
          throw string("Chemistry cannot be activated in GaussianPuff ")
            + "model without chemistry and photolysis "
            + "in Eulerian model.";

        // Spliting of sources if there are across several grid cell.
        SplitGaussianSource();

        // Conversion of sources coordinates from lat/lon to cartesian.
        Array<T, 2> source_coordinates;
        GaussianModel->GetSourceCoordinates(source_coordinates);
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

        GaussianModel->SetSourceCoordinates(source_coordinates_tmp);
        GaussianModel->EmptyPointSourceList();

        if (option_temporal_profile)
          read_temporal_profile();

        rank = 0;
        Nproc = 1;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &Nproc);

        GaussianModel->SplitSourceList(rank, Nproc);
#endif
      }
  }


  //! Split all line sources that is in several grid cell.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::SplitGaussianSource()
  {
    int index_z_1, index_y_1, index_x_1, index_z_2,
      index_y_2, index_x_2, N, species;
    T x_1, x_2, x_3 = 0., y_1, y_2, y_3 = 0., z_1, z_2, z_3, rate, length,
      VehicleVelocity, Area, Density, length_1 = 0, length_2 = 0, width,
      eps = 0.000001;
    vector<int> index_list;
    bool Split = false, Split_i;

    // Loop over all sources until there is no source left to split.
    while (!Split)
      {
        Split = true;
        N = GaussianModel->GetGaussianSourceCount();

        // Loop over all sources.
        for (int plume_index = 0; plume_index < N; plume_index++)
          {
            Split_i = true;
            if (GaussianModel->GetEmissionType(plume_index)
                == "continuous_line")
              {
                // Read source data.
                x_1 = GaussianModel->GetSourceParameter("X", plume_index);
                y_1 = GaussianModel->GetSourceParameter("Y", plume_index);
                z_1 = GaussianModel->GetSourceParameter("Z", plume_index);
                x_2 = GaussianModel->GetSourceParameter("X2", plume_index);
                y_2 = GaussianModel->GetSourceParameter("Y2", plume_index);
                z_2 = GaussianModel->GetSourceParameter("Z2", plume_index);
                species = GaussianModel->GetSpecies(plume_index);
                length = sqrt((x_1 - x_2) * (x_1 - x_2)
                              + (y_1 - y_2) * (y_1 - y_2));

                // Gets corresponding cell in Eulerian grid.
                GetCellIndices(x_1, y_1, z_1, index_z_1,
                               index_y_1, index_x_1);
                GetCellIndices(x_2, y_2, z_2, index_z_2,
                               index_y_2, index_x_2);

                // Split source into two parts.
                if (index_y_1 != index_y_2)
                  {
                    if (y_1 > y_2)
                      if (abs(this->y_min + this->Delta_y
                              * (index_y_1 - 0.5 - eps) - y_1) > eps)
                        y_3 = this->y_min + this->Delta_y
                          * (index_y_1 - 0.5 - eps);
                      else
                        y_3 = this->y_min + this->Delta_y
                          * (index_y_1 - 1.5 - eps);
                    else if (y_1 < y_2)
                      if (abs(this->y_min + this->Delta_y
                              * (index_y_1 + 0.5 + eps) - y_1) > eps)
                        y_3 = this->y_min + this->Delta_y
                          * (index_y_1 + 0.5 + eps);
                      else
                        y_3 = this->y_min + this->Delta_y
                          * (index_y_1 + 1.5 + eps);
                    else
                      "Error in SplitSource()";

                    x_3 = (x_2 - x_1) * (y_3 - y_1) / (y_2 - y_1) + x_1;
                    length_1 = sqrt((x_1 - x_3) * (x_1 - x_3) + (y_1 - y_3)
                                    * (y_1 - y_3));
                    length_2 = sqrt((x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2)
                                    * (y_3 - y_2));
                    Split_i = length_1 < eps || eps > length_2  ||
                      length < length_1 || length_2 > length;
                  }
                if (index_x_1 != index_x_2 && Split_i)
                  {
                    if (x_1 > x_2)
                      if (abs(this->x_min + this->Delta_x
                              * (index_x_1 - 0.5 - eps) - x_1) > eps)
                        x_3 = this->x_min + this->Delta_x
                          * (index_x_1 - 0.5 - eps);
                      else
                        x_3 = this->x_min + this->Delta_x
                          * (index_x_1 - 1.5 - eps);
                    else if (x_1 < x_2)
                      if (abs(this->x_min + this->Delta_x
                              * (index_x_1 + 0.5 + eps) - x_1) > eps)
                        x_3 = this->x_min + this->Delta_x
                          * (index_x_1 + 0.5 + eps);
                      else
                        x_3 = this->x_min + this->Delta_x
                          * (index_x_1 + 1.5 + eps);
                    y_3 = (y_2 - y_1) * (x_3 - x_1) / (x_2 - x_1) + y_1;
                    length_1 = sqrt((x_1 - x_3) * (x_1 - x_3) + (y_1 - y_3)
                                    * (y_1 - y_3));
                    length_2 = sqrt((x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2)
                                    * (y_3 - y_2));
                    Split_i = length_1 < eps || eps > length_2  ||
                      length < length_1 || length_2 > length;
                  }

                // If a source must be split.
                if (!Split_i)
                  {
                    // Add two sources.
                    z_3 = z_1;

                    rate = GaussianModel->GetSourceParameter("rate",
                                                             plume_index);
                    width = GaussianModel->GetSourceParameter("width",
                                                              plume_index);
                    VehicleVelocity  = GaussianModel->
                      GetSourceParameter("VehicleVelocity", plume_index);
                    Area = GaussianModel->GetSourceParameter("Area",
                                                             plume_index);
                    Density = GaussianModel->GetSourceParameter("Density",
                                                                plume_index);
                    if (length_1 != 0.)
                      GaussianModel->AddLineSource(plume_index, x_1, y_1, z_1,
                                                   x_3, y_3, z_3, rate,
                                                   width, VehicleVelocity,
                                                   Area, Density, species);
                    if (length_2 != 0.)
                      GaussianModel->AddLineSource(plume_index, x_3, y_3, z_3,
                                                   x_2, y_2, z_2, rate,
                                                   width, VehicleVelocity,
                                                   Area, Density, species);
                    index_list.push_back(plume_index);
                  }
              }
            Split = Split && Split_i;
          }

        // Delete duplicates sources.
        if (!Split)
          GaussianModel->EraseSource(index_list);
        index_list.clear();
      }

    // Remove sources outside of the domain.
    int index_x, index_y, index_z;
    for (int plume_index = 0; plume_index < N; plume_index++)
      if (GaussianModel->GetEmissionType(plume_index) == "continuous_line")
        {
          x_1 = GaussianModel->GetSourceParameter("X", plume_index);
          y_1 = GaussianModel->GetSourceParameter("Y", plume_index);
          z_1 = GaussianModel->GetSourceParameter("Z", plume_index);

          // Gets corresponding cell in Eulerian grid.
          GetCellIndices(x_1, y_1, z_1, index_z, index_y, index_x);

          // Tests if puff has reached the end of the domain.
          if (index_x == 0 || index_x == this->Nx - 1
              || index_y == 0 || index_y == this->Ny - 1
              || index_z == this->Nz)
            index_list.push_back(plume_index);
        }
    GaussianModel->EraseSource(index_list);
    index_list.clear();
  }


  //! Compute the contribution of a given source to each grid cell.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeGaussianSourceContribution(int id_source,
                                      Data<T, 4>& PlumeConcentration)
  {
    const T pi = 3.14159265358979323846264;
    int species = GaussianModel->GetSpecies(id_source);
    vector<int> species_list;
    bool rate_NO = false, rate_O3 = false, rate_NO2 = false;

    species_list.push_back(species);

    if (gaussian_type == "plume-line" && species_list.size() > 0)
      {
        T release_distance = 200.;
        T plume_length, alpha;
        if (GaussianModel->GetModelParameter("rural"))
          {
            if (GaussianModel->GetModelParameter("stability") == 0)
              {
                plume_length = 67.;
                alpha = 2 * atan(10. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 1)
              {
                plume_length = 107.;
                alpha = 2 * atan(7. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 2)
              {
                plume_length = 161.;
                alpha = 2 * atan(5. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 3)
              {
                plume_length = 240.;
                alpha = 2 * atan(4. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 4)
              {
                plume_length = 409.;
                alpha = 2 * atan(3. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 5)
              {
                plume_length = 744.;
                alpha = 2 * atan(2. / 15.);
              }
            else
              throw "Stability class is not valid.";
          }
        else
          {
            if (GaussianModel->GetModelParameter("stability") == 0)
              {
                plume_length = 53.;
                alpha = 2 * atan(11. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 1)
              {
                plume_length = 53.;
                alpha = 2 * atan(11. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 2)
              {
                plume_length = 67.;
                alpha = 2 * atan(10. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 3)
              {
                plume_length = 96.;
                alpha = 2 * atan(7. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 4)
              {
                plume_length = 175.;
                alpha = 2 * atan(5. / 15.);
              }
            else if (GaussianModel->GetModelParameter("stability") == 5)
              {
                plume_length = 175.;
                alpha = 2 * atan(5. / 15.);
              }
            else
              throw "Stability class is not valid.";
          }

        Array<T, 1> Concentration, lon_, lat_;
        T cos_angle = GaussianModel->GetModelParameter("cos_angle");
        T sin_angle = GaussianModel->GetModelParameter("sin_angle");
        T x_, y_, sigma_z, W = 1., length = 1., Q;
        int index_x, index_y, index_z, Np;
        T z = GaussianModel->GetSourceParameter("Z", id_source);
        Array<T, 2> Quantity;

        if (GaussianModel->GetEmissionType(id_source) == "continuous_line")
          {
            // Initialize Gaussian parameters.
            T width = GaussianModel->GetSourceParameter("width", id_source);
            T x1 = GaussianModel->GetSourceParameter("X", id_source);
            T y1 = GaussianModel->GetSourceParameter("Y", id_source);
            T x2 = GaussianModel->GetSourceParameter("X2", id_source);
            T y2 = GaussianModel->GetSourceParameter("Y2", id_source);

            length = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            T wind_angle = GaussianModel->GetModelParameter("wind_angle");

            if (width > 0.)
              W = width;

            // Number of points to be computed per sources.
            T Nx_ = 10;
            T Ny_ = 10;
            Np = Nx_ * Ny_;
            lon_.resize(Np * Np);
            lat_.resize(Np * Np);
            Concentration.resize(Np * Np);
            Concentration = 0.;

            // Discretization parameters.
            if (abs(cos_angle) < 0.00001)
              cos_angle = 0.;
            if (abs(sin_angle) < 0.00001)
              sin_angle = 0.;

            T cos_beta = (y2 - y1) / length;
            T sin_beta = (x2 - x1) / length;
            T cos_gamma = sin_angle * sin_beta - cos_angle * cos_beta;
            T sin_gamma = sin_beta * cos_angle + cos_beta * sin_angle;
            if (abs(cos_gamma) < 0.00001)
              cos_gamma = 0.;
            if (abs(sin_gamma) < 0.00001)
              sin_gamma = 0.;

            T Dr1 = release_distance + length * abs(sin_gamma) / 2.;
            T Dr2 = abs(release_distance - length * abs(sin_gamma) / 2.);
            T D;
            if (cos_gamma * sin_gamma >= 0.)
              D = Dr1 / cos(alpha / 2.);
            else
              D = Dr2 / cos(alpha / 2.);
            T D3 = abs(cos_gamma) * length + (Dr1 + Dr2) * tan(alpha / 2.);
            T D4 = plume_length / cos(alpha / 2.);

            T x0, y0;
            if (cos_gamma > 0. && sin_gamma <= 0. ||
                cos_gamma >= 0. && sin_gamma > 0.)
              {
                x0 = x1 + cos(alpha / 2. + wind_angle) * D;
                y0 = y1 + sin(alpha / 2. + wind_angle) * D;
              }
            else
              {
                x0 = x2 + cos(alpha / 2. + wind_angle) * D;
                y0 = y2 + sin(alpha / 2. + wind_angle) * D;
              }

            T dx_1 = D3 / Nx_ * sin_angle;
            T dy_1 = -D3 / Nx_ * cos_angle;
            T dx_2 = cos(alpha / 2. + wind_angle) * D4  / Ny_;
            T dy_2 = sin(alpha / 2. + wind_angle) * D4  / Ny_;

            // Loop on x and y axis.
            int id = 0;
            T D5;
            for (int j = 0; j < Ny_; j++)
              {
                if (j > 0)
                  {
                    D3 = abs(cos_gamma) * length + 2. *
                      (release_distance + (j - 1) * plume_length / Ny_) *
                      tan(alpha / 2.);
                    D5 = abs(cos_gamma) * length + 2. *
                      (release_distance + j * plume_length / Ny_) *
                      tan(alpha / 2.);
                    Nx_ = Nx_  + (D5 - D3) /
                      (sqrt(dx_1 * dx_1 + dy_1 * dy_1));
                  }
                for (int i = 0; i < Nx_; i++)
                  {
                    // Compute coordinates.
                    x_ = x0 + i * dx_1 + j * dx_2;
                    y_ = y0 + i * dy_1 + j * dy_2;

                    CartesianToLatLon(x_, y_, lon_(id), lat_(id));

                    // Compute concentration.
                    Concentration(id) = GaussianModel->
                      ComputeGaussianConcentration(species, z, y_, x_);
                    id++;
                  }
              }
            Np = id;
          }
        else if (GaussianModel->GetEmissionType(id_source) == "continuous")
          {
            // Compute realase distance coordinates.
            T x_source = GaussianModel->GetSourceParameter("X", id_source);
            T y_source = GaussianModel->GetSourceParameter("Y", id_source);
            T wind_angle = GaussianModel->GetModelParameter("wind_angle");

            // Number of points to be computed per sources.
            T Nx_ = 10;
            T Ny_ = 10;
            Np = Nx_ * Ny_;
            lon_.resize(Np * Np);
            lat_.resize(Np * Np);
            Concentration.resize(Np * Np);
            Concentration = 0.;

            // Discretization parameters.
            if (abs(cos_angle) < 0.00001)
              cos_angle = 0.;
            if (abs(sin_angle) < 0.00001)
              sin_angle = 0.;

            T D = release_distance / cos(alpha / 2.);
            T D3 = 2 * release_distance * tan(alpha / 2.);
            T D4 = plume_length / cos(alpha / 2.);
            T x0 = x_source + cos(alpha / 2. + wind_angle) * D;
            T y0 = y_source + sin(alpha / 2. + wind_angle) * D;

            T dx_1 = D3 / Nx_ * sin_angle;
            T dy_1 = -D3 / Nx_ * cos_angle;
            T dx_2 = cos(alpha / 2. + wind_angle) * D4  / Ny_;
            T dy_2 = sin(alpha / 2. + wind_angle) * D4  / Ny_;

            // Loop on x and y axis.
            int id = 0;
            T D5;
            for (int j = 0; j < Ny_; j++)
              {
                if (j > 0)
                  {
                    D3 = 2. * (release_distance + (j - 1) *
                               plume_length / Ny_) *
                      tan(alpha / 2.);
                    D5 = 2. * (release_distance + j *
                               plume_length / Ny_) *
                      tan(alpha / 2.);
                    Nx_ = Nx_  + (D5 - D3) /
                      (sqrt(dx_1 * dx_1 + dy_1 * dy_1));
                  }
                for (int i = 0; i < Nx_; i++)
                  {
                    // Compute coordinates.
                    x_ = x0 + i * dx_1 + j * dx_2;
                    y_ = y0 + i * dy_1 + j * dy_2;

                    CartesianToLatLon(x_, y_, lon_(id), lat_(id));

                    // Compute concentration.
                    Concentration(id) = GaussianModel->
                      ComputeGaussianConcentration(species, z, y_, x_);
                    id++;
                  }
              }
            Np = id;
          }
        else
          throw string(GaussianModel->GetEmissionType(id_source))
            + " is not a valid emission type";

        T S = sum(Concentration);
        if (S > 0.)
          {
            Quantity.resize(int(species_list.size()), Np);
            Quantity(0, Range::all()) = GaussianModel->
              GetSourceParameter("rate", id_source) * W * length
              * delta_t_eulerian * Concentration(Range::all()) / S;

            bool isday;
            for (int i = 0; i < Np; i++)
              {
                LatLonToCartesian(lon_(i), lat_(i), x_, y_);
                sigma_z = GaussianModel->ComputeSigmaZ(id_source, x_, y_);

                isday = IsDay(lon_(i), lat_(i), GaussianModel->
                              GetCurrentDate());

                GetCellIndices(lon_(i), lat_(i), z, index_z,
                               index_y, index_x);

                // Tests if plume has reached the end of the domain.
                if (index_x != 0 && index_x != this->Nx - 1
                    && index_y != 0 && index_y != this->Ny - 1
                    && index_z != this->Nz)
                  // Loop on all required species.
                  for (int sp = 0; sp < int(species_list.size()); sp++)
                    GaussianTransfer(Quantity(sp, i), sigma_z,
                                     species_list[sp], z, lat_(i), lon_(i),
                                     isday, PlumeConcentration);
              }
          }
      }
  }


  //! Reads the temporal profile coefficients.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::read_temporal_profile()
  {
    string data_file, temporal_profile_file;
    Date date_min_;
    T delta_t_;
    this->config.PeekValue("Data_description", data_file);
    ConfigStream data(data_file);
    data.SetSection("[temporal_profile]");
    data.PeekValue("Gaussian_temporal_profile_file", temporal_profile_file);
    date_min_ = data.PeekValue("Date_min");
    data.PeekValue("Delta_t", delta_t_);

    ConfigStream temporal_profile_stream(temporal_profile_file);
    vector<string> tmp;
    while (!temporal_profile_stream.IsEmpty())
      tmp.push_back(temporal_profile_stream.GetElement());

    int N = tmp.size();
    Date date_max_ref, date_max_, date_min_ref;
    date_max_ref = this->Date_min;
    date_max_ = date_min_;
    date_max_ref.AddSeconds(this->Delta_t * this->Nt);
    date_max_.AddSeconds(delta_t_ * tmp.size());
    if (this->Date_min < date_min_ || date_max_ref > date_max_)
      throw string("'PlumeInGrid::read_temporal_profile()': There are "
                   "missing values in the temporal profile file \"")
        + temporal_profile_file + "\".";

    temporal_profile_coefficients.resize(this->Nt);
    if (delta_t_ == this->Delta_t && this->Date_min == date_min_)
      for (int i = 0; i < this->Nt; i++)
        temporal_profile_coefficients(i) = to_num<T>(tmp[i]);
    else
      {
        // Interpolation needed.
        int date_min_sec = 0.;
        int current_t = this->Date_min.GetSecondsFrom(date_min_);
        int current_t_tmp, id;
        T y0, y1;
        for (int i = 0; i < this->Nt; i++)
          {
            current_t = this->Date_min.GetSecondsFrom(date_min_) +
              i * this->Delta_t;
            current_t_tmp = date_min_sec;
            id = 0;
            y0 = to_num<T>(tmp[0]);
            y1 = to_num<T>(tmp[1]);
            while (current_t_tmp > current_t
                   || current_t_tmp + delta_t_ < current_t)
              {
                current_t_tmp = current_t_tmp + delta_t_;
                id = id + 1;
                y0 = to_num<T>(tmp[id]);
                y1 = to_num<T>(tmp[id + 1]);
              }

            if (current_t_tmp == current_t)
              temporal_profile_coefficients(i) = y0;
            else if (current_t_tmp + delta_t_ == current_t)
              temporal_profile_coefficients(i) = y1;
            else if (current_t_tmp < current_t &&
                     current_t_tmp + delta_t_ > current_t)
              temporal_profile_coefficients(i) = (y1 - y0) / delta_t_
                * current_t + y0 - (y1 - y0) / delta_t_ * current_t_tmp;
            else
              throw "'PlumeInGrid::read_temporal_profile()': Error in "
                "temporal profile interpolation";
          }
      }
  }


  // ! Get gaussian meteorological parameters.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void  PlumeInGrid<T, ClassEulerianModel, ClassLocalModel>
  ::GetGaussianMeteo(T z_c, T lat_c, T lon_c, T &SpecificHumidity,
                     T &Temperature, T &Pressure, int &Stability, int &Rural,
                     T& Wind, T& Attenuation)
  {
    // Night or day.
    bool isday = IsDay(lon_c, lat_c, GaussianModel->GetCurrentDate());

    // Extracting meteorological data.
    bool rural;
    string stability;
    map<string, T> meteorological_data;
    bool option_similarity = GaussianModel->WithSimilarity();

    ExtractMeteo(z_c, lat_c, lon_c, isday, option_similarity,
                 meteorological_data, stability, rural);

    SpecificHumidity = meteorological_data["specific_humidity"];
    Temperature = meteorological_data["temperature"];
    Pressure = meteorological_data["pressure"];

    if (stability == "A")
      Stability = 0;
    else if (stability == "B")
      Stability = 1;
    else if (stability == "C")
      Stability = 2;
    else if (stability == "D")
      Stability = 3;
    else if (stability == "E")
      Stability = 4;
    else if (stability == "F")
      Stability = 5;

    Rural = int(rural);
    Wind = meteorological_data["wind"];
    Attenuation = meteorological_data["attenuation"];
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PLUMEINGRID_CXX
#endif
