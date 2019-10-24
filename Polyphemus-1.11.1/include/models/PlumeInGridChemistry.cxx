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


#ifndef POLYPHEMUS_FILE_MODELS_PLUMEINGRIDCHEMISTRY_CXX


#include "PlumeInGridChemistry.hxx"


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
  PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::PlumeInGridChemistry(string config_file):
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>(config_file)
  {
  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::~PlumeInGridChemistry()
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
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::ReadConfiguration();

    this->config.SetSection("[options]");
    bool option_photolysis;
    this->config.PeekValue("With_chemistry", option_chemistry);
    this->config.PeekValue("With_photolysis", option_photolysis);

    option_chemistry = (option_chemistry && option_photolysis);

    this->config.SetSection("[data]");

    string data_file, emission_file;
    this->config.PeekValue("Data_description", data_file);

    // Plume-in-grid general data.
    ConfigStream data(data_file);

    data.SetSection("[plume-in-grid]");
    if (option_chemistry)
      data.PeekValue("With_chemistry_feedback", option_chemical_feedback);

  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::CheckConfiguration()
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>::CheckConfiguration();
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.

  }


  //! Initialization.
  /*! Initializes the Eulerian model and the meteorological conditions. Reads
    the sources and creates the corresponding Gaussian model instances.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::Init();

    /*** Initializations ***/

    Nr_photolysis = 0;
    vector<string> photolysis_reaction_list;

    if (option_chemistry)
      {
        // Extracts information useful for chemistry from Eulerian model.
        Nr_photolysis = this->Model.GetNr_photolysis();
        photolysis_reaction_list = this->Model.GetPhotolysisReactionList();

        // Photolysis rates.
        if (this->Model.HasField("PhotolysisRate"))
          PhotolysisRate_i.Copy(this->Model.D4("PhotolysisRate_i"));
        else
          throw string("Eulerian model doesn't have a field ")
            + "\"PhotolysisRate\".";

        // Meteorological data.
        Attenuation_i.Copy(this->Model.D3("Attenuation_i"));
        Pressure_i.Copy(this->Model.D3("Pressure_i"));
        SpecificHumidity_i.Copy(this->Model.D3("SpecificHumidity_i"));
      }

    // Is there puff chemistry?
    if (!this->GaussianPuffModel->WithChemistry())
      option_chemistry = false;
    else if (!option_chemistry)
      throw string("Chemistry cannot be activated in GaussianPuff ")
        + "model without chemistry and photolysis "
        + "in Eulerian model.";

    if (option_chemistry)
      this->GaussianPuffModel->InitPhotolysis(Nr_photolysis,
                                              photolysis_reaction_list);
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */  template<class T, class ClassEulerianModel, class ClassLocalModel>
   void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
   ::InitStep()
   {
     PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
       ::InitStep();

     if (option_chemistry)
       {
         Attenuation_i.Copy(this->Model.D3("Attenuation_i"));
         Pressure_i.Copy(this->Model.D3("Pressure_i"));
         SpecificHumidity_i.Copy(this->Model.D3("SpecificHumidity_i"));

         // Photolysis rates.
         PhotolysisRate_i.Copy(this->Model.D4("PhotolysisRate_i"));
       }

   }


  //! Performs one step forward.
  /*! Calls the Eulerian model. For each major point emission, the
    corresponding Gaussian Model instance is called.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    this->Concentration.Copy(this->Model.GetConcentration());

    int Npuff;
    int index_x, index_y, index_z, puff_index;
    bool puff_transfer;
    bool isday;

    /*** Tests if puff has to be transfered to Eulerian model. ***/

    Npuff = this->GaussianPuffModel->GetPuffNumber();

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
          }

        // Tests if puff size has reached the cell width.
        else if (this->coefficient_y * sigma_y >= cell_width_y)
          puff_transfer = 1;

        // Puff transfer.
        if (puff_transfer)
          {
            if (this->injection_integrated)
              this->PuffIntegratedTransfer(puff_index, this->Model.GetConcentration());
            else
              for (int s = 0; s < this->Ns; s++)
                this->PuffTransfer(this->GaussianPuffModel->
                                   GetPuffQuantity(puff_index, s),
                                   sigma_z, s, z_c, lat_c, lon_c,
                                   isday, this->Model.GetConcentration());
            this->GaussianPuffModel->ErasePuff(puff_index);
          }
      }

    /*** Inner time-loop for Gaussian models. ***/

    Array<T, 4> PuffConcentration(this->Ns, this->Nz, this->Ny, this->Nx);
    PuffConcentration = 0.;
    int l, k, j, i, s;

    for (j = 0; j < this->Nt_local; j++)
      {
        this->GaussianPuffModel->InitStep();
        Npuff = this->GaussianPuffModel->GetPuffNumber();
        Array<int, 2> PuffCellList(Npuff, 3);

        // Loop on puffs to get the meteorological data.
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
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
              this->UpdateMeteo(puff_index);
          }

        // Puff effective height.
        this->GaussianPuffModel->ComputePlumeRise();

        // Advection and diffusion for puffs.
        this->GaussianPuffModel->Advection();
        this->GaussianPuffModel->Diffusion();
        this->GaussianPuffModel->ComputeLossFactor();

        // Loop on puffs to update data if necessary.
        for (puff_index = 0; puff_index < Npuff; puff_index++)
          {
            bool has_changed = false;
            // Puff center coordinates.
            T x_c, y_c, z_c, lon_c, lat_c, puff_distance, puff_time;
            this->GaussianPuffModel->
              GetPuffPosition(puff_index, x_c, y_c, z_c, puff_distance,
                              puff_time);

            // Conversion to longitude/latitude (degrees).
            this->CartesianToLatLon(x_c, y_c, lon_c, lat_c);

            // Gets corresponding cell in eulerian grid.
            this->GetCellIndices(lon_c, lat_c, z_c, index_z, index_y, index_x);

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
            ComputeChemistry(PuffCellList, PuffConcentration);
          else
            this->GaussianPuffModel->Chemistry();
        this->GaussianPuffModel->AddTime(this->delta_t_local);
      }

    // Eulerian model.
    this->Model.Forward();

    if (option_chemistry)
      if (option_chemical_feedback)
        for (l = 0; l < this->Nz; l++)
          for (k = 0; k < this->Ny; k++)
            for (i = 0; i < this->Nx; i++)
              for (s = 0; s < this->Ns; s++)
                this->Model.GetConcentration()(s, l, k, i)
                  += PuffConcentration(s, l, k, i);

    this->Concentration.Copy(this->Model.GetConcentration());

    /* Adding all puffs to concentration Data in case concentrations are saved
       on the domain.*/

    this->GaussianPuffModel->SubtractTime(this->delta_t_local);
    Npuff = this->GaussianPuffModel->GetPuffNumber();
    for (puff_index = 0; puff_index < Npuff; puff_index++)
      {
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
            T sigma_x, sigma_y, sigma_z;
            this->GaussianPuffModel->GetPuffSigma(puff_index, sigma_x,
                                                  sigma_y, sigma_z);

            isday = IsDay(lon_c, lat_c, this->GaussianPuffModel->
                          GetCurrentDate());
            // Adding puff concentration to the Concentration Data.
            if (this->injection_integrated)
              this->PuffIntegratedTransfer(puff_index, this->Concentration);
            else
              for (int s = 0; s < this->Ns; s++)
                this->PuffTransfer(this->GaussianPuffModel->
                                   GetPuffQuantity(puff_index, s),
                                   sigma_z, s, z_c, lat_c, lon_c,
                                   isday, this->Concentration);
          }
      }
    this->GaussianPuffModel->AddTime(this->delta_t_local);
    this->AddTime(this->delta_t_eulerian);
    this->step++;
  }



  //! Updates the puff meteorological data.
  /*! Updates the puff meteorological data taking them from the center cell.
    \param puff_index index of the puff.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
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

    if (option_chemistry)
      this->GaussianPuffModel->SetPuffAdditionalMeteo
        (puff_index,
         meteorological_data["attenuation"],
         meteorological_data["pressure"],
         meteorological_data["specific_humidity"]);

    // Data depending on the puff species.
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

    if (this->option_deposition)
      this->GaussianPuffModel
        ->SetPuffDepositionVelocity(puff_index,
                                    deposition_velocity);

    if (this->option_scavenging)
      this->GaussianPuffModel
        ->SetPuffScavengingCoefficient(puff_index,
                                       scavenging_coefficient);

    if (option_chemistry)
      {
        this->GaussianPuffModel
          ->SetPuffPhotolysisRate(puff_index,
                                  photolysis_rate);
        this->GaussianPuffModel->SetPuffBackgroundConcentration
          (puff_index, background_concentration);
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
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractMeteo(T height, T lat, T lon, bool isday,
                 bool option_similarity, map<string, T> &met_data,
                 string &stability, bool& rural)
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::ExtractMeteo(height, lat, lon, isday,
                     option_similarity, met_data,
                     stability, rural);
    // const T pi = 3.14159265358979323846264;
    int index_x, index_y, index_z;
    this->GetCellIndices(lon, lat, height, index_z, index_y, index_x);

    Array<T, 1> Coord3D(3);
    Coord3D(0) = height;
    Coord3D(1) = lat;
    Coord3D(2) = lon;

    if (option_chemistry)
      {
        T attenuation, pressure, specific_humidity;
        if (!this->interpolated)
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
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractSpeciesData(T height, T lat, T lon,
                       Array<T, 1>& deposition_velocity,
                       Array<T, 1>& scavenging_coefficient,
                       Array<T, 1>& photolysis_rate,
                       Array<T, 1>& background_concentration)
  {
    PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::ExtractSpeciesData(height, lat, lon,
                           deposition_velocity,
                           scavenging_coefficient);

    int r, s, index_x, index_y, index_z;
    this->GetCellIndices(lon, lat, height, index_z, index_y, index_x);

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
            if (!this->interpolated)
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
            if (!this->interpolated)
              background_concentration(s) =
                this->Model.GetConcentration()(s, index_z, index_y, index_x);
            else
              LinearInterpolationPoint(this->Model.GetConcentration(),
                                       Coord4D, background_concentration(s));
            // To correct extrapolation that could give concentrations < 0.
            background_concentration(s) =
              max(background_concentration(s), 0.);
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
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ComputeChemistry(Array<int, 2> PuffCellList,
                     Array<T, 4>& PuffConcentration)
  {
    int l, k, i, s, puff_index;
    int index_x = 0;
    int index_y = 0;
    int index_z = 0;
    int Npuff = this->GaussianPuffModel->GetPuffNumber();
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
                this->ComputeCellWidth(index_z, index_y, index_x, cell_width_z,
                                       cell_width_y, cell_width_x, cell_volume);
                PuffCellVolume.push_back(cell_volume);
                Coord(0) = index_z;
                Coord(1) = index_y;
                Coord(2) = index_x;
                PuffCellCoordinates.push_back(Coord);
                for (s = 0; s < this->Ns; s++)
                  PuffCellConcentration(s).push_back
                    (this->Model.GetConcentration()(s, index_z, index_y, index_x));
                PuffList.push_back(pufflist_tmp);
              }
          }
#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
    this->GaussianPuffModel->Chemistry(PuffList, PuffCellVolume,
                                       PuffCellConcentration);
#endif
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

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDCHEMISTRY_CXX
#endif
