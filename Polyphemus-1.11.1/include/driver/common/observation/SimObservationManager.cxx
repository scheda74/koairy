// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_OBSERVATION_SIMOBSERVATIONMANAGER_CXX


#include "SimObservationManager.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template <class T>
  SimObservationManager<T>
  ::SimObservationManager(string config_file)
    : BaseObservationManager<T>(config_file)
  {
  }


  //! Destructor.
  template <class T>
  SimObservationManager<T>::~SimObservationManager()
  {
  }


  ////////////////////
  // CONFIGURATIONS //
  ////////////////////


  //! Reads configurations.
  /*! It reads the option that specifies how observations are provided, and
    the attributes for the simulation data files. For station option, it reads
    the numbers of stations for all the observed species, the names of the
    files that describe the stations, and the paths that contain the data
    files.
  */
  template <class T>
  template <class ClassModel>
  void SimObservationManager<T>::ReadConfiguration(ClassModel& Model)
  {
    int i;

    string obs_manager_file;
    this->config.SetSection("[observation_management]");
    this->config.PeekValue("Configuration_file", obs_manager_file);
    ConfigStream config_obs(obs_manager_file);

    this->Nobs_array.resize(this->Ns_obs);

    /*** Simulation manager section ***/

    config_obs.SetSection("[simulation_manager]");
    config_obs.PeekValue("Simulation_option", simulation_option);
    sim_date_min = config_obs.PeekValue("Date_min");
    config_obs.PeekValue("Delta_t", sim_delta_t);

    // Settings for simulation data files.
    string generic_name, file_name;
    config_obs.PeekValue("Input_file", generic_name);
    sim_input_file_list.clear();
    for (i = 0; i < this->Ns_obs; i++)
      {
        file_name
          = find_replace(generic_name, "&f", this->obs_species_list[i]);
        file_name = find_replace(file_name, "&s", this->obs_species_list[i]);
        sim_input_file_list.push_back(file_name);
      }

    config_obs.Find("Levels");
    vector<int> sim_levels_vector;
    split(config_obs.GetLine(), sim_levels_vector);
    Nz_sim = int(sim_levels_vector.size());
    sim_levels.resize(Nz_sim);
    for (i = 0; i < Nz_sim; i++)
      sim_levels(i) = sim_levels_vector[i];

    /*** Treatment of simulation options ***/

    if (simulation_option == "with_station")
      // Simulated concentrations are interpolated at station locations.
      {
        config_obs.SetSection("[stations]");

        // Reads array of station numbers.
        config_obs.Find("Nstation");
        vector<int> Nstation_array_vector;
        split(config_obs.GetLine(), Nstation_array_vector);
        int num_species = int(Nstation_array_vector.size());
        if (num_species != this->Ns_obs)
          throw string("Found ") + to_str(this->Ns_obs)
            + string(" observed species, but ") + to_str(num_species)
            + string(" number(s) of station number were provided ")
            + string(" in entry \"Nstation\".");

        Nstation = 0;
        Nstation_array.resize(this->Ns_obs);
        for (i = 0; i < this->Ns_obs; i++)
          {
            Nstation_array(i) = Nstation_array_vector[i];
            Nstation += Nstation_array_vector[i];
          }

        // Reads station files.
        string generic_name, file_name;
        config_obs.PeekValue("Stations_file", generic_name);
        station_file_list.clear();
        for (i = 0; i < this->Ns_obs; i++)
          {
            file_name
              = find_replace(generic_name, "&f", this->obs_species_list[i]);
            file_name
              = find_replace(file_name, "&s", this->obs_species_list[i]);
            station_file_list.push_back(file_name);
          }

        // Reads input directories for data files.
        config_obs.PeekValue("Input_directory", generic_name);
        input_dir_list.clear();
        for (i = 0; i < this->Ns_obs; i++)
          {
            file_name
              = find_replace(generic_name, "&f", this->obs_species_list[i]);
            file_name
              = find_replace(file_name, "&s", this->obs_species_list[i]);
            input_dir_list.push_back(file_name);
          }
      }
    else if (simulation_option == "with_level")
      // Includes 2D layers of simulated data.
      {
        config_obs.SetSection("[simulation_manager_level]");
        config_obs.Find("Observation_levels");
        vector<int> obs_levels_vector;
        split(config_obs.GetLine(), obs_levels_vector);
        Nz_obs = int(obs_levels_vector.size());
        obs_levels.resize(Nz_obs);
        for (i = 0; i < Nz_obs; i++)
          obs_levels(i) = obs_levels_vector[i];

        // Sets the number of observations for the species.
        this->Nobs_array = Nz_obs * this->model_Ny * this->model_Nx;
        Nobs_total = this->Ns_obs * Nz_obs * this->model_Ny * this->model_Nx;

        this->ObservationStateIndex_s.resize(Nobs_total);
        this->ObservationIndex_s.resize(Nobs_total);
        this->ObservationIndex_z.resize(Nobs_total);
        this->ObservationIndex_y.resize(Nobs_total);
        this->ObservationIndex_x.resize(Nobs_total);
        ObservationIndex_in_state.resize(Nobs_total);

        this->ObservationStateSequenceIndex_z.resize(Nobs_total);
        this->ObservationStateSequenceIndex_z
          = Model.GetStateLevelSequenceIndex(0);

        int pos_s, k, j;
        int r = 0;
        for (int s = 0; s < this->Ns_obs; s++)
          {
            pos_s = Model.GetStateSpeciesIndex(this->obs_species_list[s]);
            pos_s *= Model.GetNz() * this->model_Ny * this->model_Nx;
            for (k = 0; k < Nz_obs; k++)
              for (j = 0; j < this->model_Ny; j++)
                for (i = 0; i < this->model_Nx; i++)
                  {
                    ObservationIndex_in_state(r) = pos_s
                      + obs_levels(k) * this->model_Ny * this->model_Nx
                      + j * this->model_Nx + i;
                    this->ObservationIndex_s(r) = this->obs_species_index(s);
                    this->ObservationStateIndex_s(r)
                      = Model.GetStateSpeciesIndex(this->obs_species_list[s]);
                    this->ObservationIndex_z(r) = obs_levels(k);
                    this->ObservationIndex_y(r) = j;
                    this->ObservationIndex_x(r) = i;
                    r++;
                  }
          }
      }
    else if (simulation_option == "with_location")
      throw string("Simulation option \"") + simulation_option +
        string("\" not supported.");
    else
      throw string("Simulation option \"") + simulation_option +
        string("\" unknown.");
  }


  ////////////////////
  // INITIALIZATION //
  ////////////////////


  //! Reads station descriptions and removes useless stations.
  /*! It reads the station files to get the coordinates of the stations. The
    stations out of the model domain are filtered out.
  */
  template <class T>
  template <class ClassModel>
  void SimObservationManager<T>::ReadStation(ClassModel& Model)
  {
    // Are stations in the model domain?
    Array<bool, 1> StationAvailability(Nstation);
    // Index of the left corner along x of the enclosing cell in the
    // model, for all stations.
    Array<int, 1> StationIndex_x(Nstation);
    // Index of the left corner along y of the enclosing cell in the
    // model, for all stations.
    Array<int, 1> StationIndex_y(Nstation);

    Array<T, 1> station_longitude, station_latitude;
    station_longitude.resize(Nstation);
    station_latitude.resize(Nstation);

    // Reads coordinates.
    int s, i;
    int Nstation_accumulate = 0;
    int Nobs_accumulate = 0;

    for (s = 0; s < this->Ns_obs; s++)
      {

        /*** For species 's' ***/

        // Reads station latitude and longitude.
        FormatFormattedText Emep(string("<c 0 4><c 5 44><c 49 2><c 52 2>")
                                 + string("<c 55 2><c 63 2><c 66 2><c 69 2>")
                                 + string("<c 72 1><e>"));

        ExtStream station_list_stream(station_file_list[s]);
        streampos position = station_list_stream.tellg();

        Array<float, 1> station_tmp0(Nstation_array(s));
        Array<float, 1> station_tmp1(Nstation_array(s));
        Array<float, 1> station_tmp2(Nstation_array(s));
        Array<char, 1> station_longitude_EW(Nstation_array(s));

        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "2", station_tmp0);
        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "3", station_tmp1);
        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "4", station_tmp2);

        for (i = 0; i < Nstation_array(s); i++)
          station_latitude(i + Nstation_accumulate) = station_tmp0(i) +
            station_tmp1(i) / 60. + station_tmp2(i) / 3600.;

        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "5", station_tmp0);
        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "6", station_tmp1);
        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "7", station_tmp2);

        for (i = 0; i < Nstation_array(s); i++)
          station_longitude(i + Nstation_accumulate) = station_tmp0(i) +
            station_tmp1(i) / 60. + station_tmp2(i) / 3600.;

        station_list_stream.seekg(position);
        Emep.Read(station_list_stream, "8", station_longitude_EW);

        for (i = 0; i < Nstation_array(s); i++)
          if (station_longitude_EW(i) == 'w'
              || station_longitude_EW(i) == 'W')
            station_longitude(i + Nstation_accumulate) =
              -station_longitude(i + Nstation_accumulate);

        // Filters out the stations out of model domain.
        for (i = 0; i < Nstation_array(s); i++)
          {
            StationIndex_x(i + Nstation_accumulate) =
              int((station_longitude(i + Nstation_accumulate)
                   - this->model_x_min) / this->model_Delta_x + 0.5);

            StationIndex_y(i + Nstation_accumulate) =
              int((station_latitude(i + Nstation_accumulate)
                   - this->model_y_min) / this->model_Delta_y + 0.5);

            // Filters out the stations out of model domain.
            StationAvailability(i + Nstation_accumulate) =
              StationIndex_x(i + Nstation_accumulate) >= 0
              && StationIndex_x(i + Nstation_accumulate) < this->model_Nx
              && StationIndex_y(i + Nstation_accumulate) >= 0
              && StationIndex_y(i + Nstation_accumulate) < this->model_Ny;
          }

        // Counts the available stations for species 's'.
        this->Nobs_array(s) = 0;
        for (i = 0; i < Nstation_array(s); i++)
          this->Nobs_array(s) += StationAvailability(i + Nstation_accumulate);

        // Sets the observation index array.
        this->ObservationIndex_x.resizeAndPreserve(Nobs_accumulate +
                                                   this->Nobs_array(s));
        this->ObservationIndex_y.resizeAndPreserve(Nobs_accumulate +
                                                   this->Nobs_array(s));
        this->ObservationStateIndex_s.resizeAndPreserve(Nobs_accumulate +
                                                        this->Nobs_array(s));

        int obs_state_index_s
          = Model.GetStateSpeciesIndex(this->obs_species_list[s]);
        int r = 0;
        for (i = 0; i < Nstation_array(s); i++)
          {
            if (StationAvailability(i + Nstation_accumulate))
              {
                this->ObservationIndex_x(r + Nobs_accumulate)
                  = StationIndex_x(i + Nstation_accumulate);
                this->ObservationIndex_y(r + Nobs_accumulate)
                  = StationIndex_y(i + Nstation_accumulate);
                this->ObservationStateIndex_s(r + Nobs_accumulate)
                  = obs_state_index_s;
                r++;
              }
          }

        // Updates indices.
        Nobs_accumulate += this->Nobs_array(s);
        Nstation_accumulate += Nstation_array(s);
      }

    // Computes the total number of observations.
    Nobs_total = 0;
    for (s = 0; s < this->Ns_obs; s++)
      Nobs_total += this->Nobs_array(s);

    // Sets station levels to zero.
    this->ObservationIndex_z.resize(Nobs_total);
    this->ObservationIndex_z = 0;
    this->ObservationStateSequenceIndex_z.resize(Nobs_total);
    this->ObservationStateSequenceIndex_z
      = Model.GetStateLevelSequenceIndex(0);

    // Computes the corresponding indices in state vector for all available
    // observations.
    ObservationIndex_in_state.resize(Nobs_total);
    for (int r = 0; r < Nobs_total; r++)
      ObservationIndex_in_state(r) = this->ObservationStateIndex_s(r) +
        this->ObservationStateSequenceIndex_z(r) * this->model_Ny
        * this->model_Nx + this->ObservationIndex_y(r) * this->model_Nx +
        this->ObservationIndex_x(r);
  }


  //! Initialization of the simulation observation manager.
  /*! It read the configuration, species indices in model state vector, and
    filters out stations out of model domain.
    \param Model model with the following interface:
    <ul>
    <li> GetNx()
    <li> GetNy()
    <li> GetX_min()
    <li> GetY_min()
    <li> GetDelta_x()
    <li> GetDelta_y()
    <li> GetStateNz()
    <li> GetStateLevels()
    <li> GetStateNs()
    <li> GetStateSpeciesList()
    <li> GetStateSpeciesIndex(string)
    <li> GetNstate()
    <li> GetStateLevelSequenceIndex(int)
    </ul>
  */
  template <class T>
  template <class ClassModel>
  void SimObservationManager<T>::Init(ClassModel& Model)
  {
    BaseObservationManager<T>::Init(Model);
    ReadConfiguration(Model);

    if (simulation_option == "with_station")
      ReadStation(Model);
  }


  /////////////////////////////
  // OBSERVATIONS MANAGEMENT //
  /////////////////////////////


  //! Sets the simulation observations at a given date.
  /*! It retrieves the observation data from simulation results, computes the
    number of observations, sets the observation error covariance, and
    perturbs observations.
    \param date the given date.
  */
  template <class T>
  void SimObservationManager<T>::SetDate(Date date)
  {
    double seconds = date.GetSecondsFrom(sim_date_min);
    int istep = int(seconds / sim_delta_t + .5);

    if (seconds < 0 || !is_multiple(seconds, sim_delta_t))
      {
        this->Nobs = 0;
        this->error_covariance.resize(0, 0);
        this->Observation.resize(0);
        this->availability = false;
        return;
      }

    // Reads data from simulated data file for all species.
    int s, r;
    int Nobs_accumulate = 0;
    for (s = 0; s < this->Ns_obs; s++)
      {
        FormatBinary<float> Input;
        Data<T, 3> data(Nz_sim, this->model_Ny, this->model_Nx);
        Input.ReadRecord(sim_input_file_list[s], istep, data);

        this->Observation.resizeAndPreserve(Nobs_accumulate
                                            + this->Nobs_array(s));

        for (r = 0; r < this->Nobs_array(s); r++)
          {
            int z_index = 0;
            while (z_index < Nz_sim &&
                   this->ObservationIndex_z(Nobs_accumulate + r)
                   != sim_levels(z_index))
              z_index++;

            if (z_index == Nz_sim)
              throw string("Level \"")
                + to_str(this->ObservationIndex_z(Nobs_accumulate + r))
                + "\" not in the data file.";

            this->Observation(Nobs_accumulate + r)
              = data(z_index,
                     this->ObservationIndex_y(Nobs_accumulate + r),
                     this->ObservationIndex_x(Nobs_accumulate + r));
          }

        Nobs_accumulate += this->Nobs_array(s);
      }

    // Successfully read all data.
    this->availability = true;
    this->Nobs = Nobs_total;

    this->error_covariance.resize(this->Nobs, this->Nobs);
    SetCovariance();
  }


  //! Operator that maps model state to one observation.
  /*! Linear pointwise observation operator.
    \param r index of the observation.
    \param state the model state vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T SimObservationManager<T>::MltObsOperator(int r, const Array<T, 1>& State)
  {
    return State(ObservationIndex_in_state(r));
  }


  //! Tangent linear operator that maps model state to one observation.
  /*! Linear pointwise observation operator, the same as MltObsOperator.
    \param r index of the observation.
    \param state the model state vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T SimObservationManager<T>::TLMMltObsOperator(int r,
                                                const Array<T, 1>& State)
  {
    return State(ObservationIndex_in_state(r));
  }


  //! Adjoint operator that maps departure vector to one state component.
  /*! Adjoint of linear pointwise observation operator.
    \param i index of the corresponding component in state vector.
    \param Departure the departure vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T SimObservationManager<T>
  ::ADJMltObsOperator(int i, const Array<T, 1>& Departure)
  {
    T sum = T(0.);
    for (int r = 0; r < this->Nobs; r++)
      if (i == ObservationIndex_in_state(r))
        sum += Departure(r);
    return sum;
  }


  //! Returns one entry of the linearized operator.
  /*!
    \param i row index.
    \param j column index.
    \return Element (i, j) of the linearized operator.
  */
  template <class T>
  T SimObservationManager<T>::Linearized(int i, int j)
  {
    if (j == ObservationIndex_in_state(i))
      return T(1);
    else
      return T(0);
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Sets error covariance.
  template <class T>
  void SimObservationManager<T>::SetCovariance()
  {
    this->error_covariance.resize(this->Nobs, this->Nobs);

    for (int i = 0; i < this->Nobs; i++)
      for (int j = 0; j < this->Nobs; j++)
        this->error_covariance(i, j) = 0.;

    int Nobs_accumulate = 0;
    for (int s = 0; s < this->Ns_obs; s++)
      {
        for (int r = 0; r < this->Nobs_array(s); r++)
          this->error_covariance(Nobs_accumulate + r, Nobs_accumulate + r)
            = this->error_variance_array(s);
        Nobs_accumulate += this->Nobs_array(s);
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OBSERVATION_SIMOBSERVATIONMANAGER_CXX
#endif
