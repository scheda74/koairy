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


#ifndef POLYPHEMUS_FILE_OBSERVATION_GROUNDOBSERVATIONMANAGER_CXX


#include "GroundObservationManager.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template <class T>
  GroundObservationManager<T>
  ::GroundObservationManager(string config_file)
    : BaseObservationManager<T>(config_file)
  {
  }


  //! Destructor.
  template <class T>
  GroundObservationManager<T>::~GroundObservationManager()
  {
  }


  ////////////////////
  // CONFIGURATIONS //
  ////////////////////


  //! Read configurations.
  /*! It reads the observation species, all the available stations, and the
    name and coordinates of these stations, as well as the names of data files
    associated with these stations.
  */
  template <class T>
  void GroundObservationManager<T>::ReadConfiguration()
  {
    BaseObservationManager<T>::ReadConfiguration();

    // File path to the configuration of the observation manager.
    string obs_manager_file;
    this->config.SetSection("[observation_management]");
    this->config.PeekValue("Configuration_file", obs_manager_file);
    ConfigStream config_obs(obs_manager_file);

    this->Nobs_array.resize(this->Ns_obs);

    /*** Main configuration ***/

    config_obs.SetSection("[general]");

    config_obs.PeekValue("With_spatial_interpolation",
                         with_spatial_interpolation);

    /*** Stations ***/

    config_obs.SetSection("[stations]");

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
    for (int i = 0; i < this->Ns_obs; i++)
      {
        Nstation_array(i) = Nstation_array_vector[i];
        Nstation += Nstation_array_vector[i];
      }

    // Reads station files.
    string generic_name, file_name;
    config_obs.PeekValue("Stations_file", generic_name);
    station_file_list.clear();
    for (int i = 0; i < this->Ns_obs; i++)
      {
        file_name
          = find_replace(generic_name, "&f", this->obs_species_list[i]);
        file_name = find_replace(file_name, "&s", this->obs_species_list[i]);
        station_file_list.push_back(file_name);
      }

    // Reads input directories for data files.
    config_obs.PeekValue("Input_directory", generic_name);
    input_dir_list.clear();
    for (int i = 0; i < this->Ns_obs; i++)
      {
        file_name
          = find_replace(generic_name, "&f", this->obs_species_list[i]);
        file_name = find_replace(file_name, "&s", this->obs_species_list[i]);
        input_dir_list.push_back(file_name);
      }
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
  void GroundObservationManager<T>::ReadStation(ClassModel& Model)
  {
    Availability.resize(Nstation);
    StationAvailability.resize(Nstation);
    StationIndex_x.resize(Nstation);
    StationIndex_y.resize(Nstation);
    StationWeight_x.resize(Nstation);
    StationWeight_y.resize(Nstation);

    station_raw_data.resize(Nstation);
    station_latitude.resize(Nstation);
    station_longitude.resize(Nstation);
    station_altitude.resize(Nstation);

    station_name.resize(Nstation);
    station_code.resize(Nstation);
    station_file.resize(Nstation);

    station_country.resize(Nstation);
    station_type.resize(Nstation);
    station_network.resize(Nstation);

    state_index_s.resize(this->Ns_obs);

    int s, i;
    int Nstation_accumulate = 0;

    for (s = 0; s < this->Ns_obs; s++)
      {

        /*** For species s ***/

        Array<string, 1> station_code_s(Nstation_array(s));
        Array<string, 1> station_name_s(Nstation_array(s));
        Array<T, 1> station_latitude_s(Nstation_array(s));
        Array<T, 1> station_longitude_s(Nstation_array(s));
        Array<T, 1> station_altitude_s(Nstation_array(s));
        Array<string, 1> station_country_s(Nstation_array(s));
        Array<string, 1> station_type_s(Nstation_array(s));
        Array<string, 1> station_network_s(Nstation_array(s));

        FormatFormattedText reader(string("<e><e><e><e><e><e><e><e>"));
        reader.SetDelimiters(";");
        reader.Read(station_file_list[s], "0", station_code_s);
        reader.Read(station_file_list[s], "1", station_name_s);
        reader.Read(station_file_list[s], "2", station_latitude_s);
        reader.Read(station_file_list[s], "3", station_longitude_s);
        reader.Read(station_file_list[s], "4", station_altitude_s);
        reader.Read(station_file_list[s], "5", station_country_s);
        reader.Read(station_file_list[s], "6", station_type_s);
        reader.Read(station_file_list[s], "7", station_network_s);

        for (i = 0; i < Nstation_array(s); i++)
          {
            station_code(i + Nstation_accumulate) = trim(station_code_s(i));
            station_name(i + Nstation_accumulate) = trim(station_name_s(i));
            station_file(i + Nstation_accumulate) = input_dir_list[s] + "/" +
              trim(station_code_s(i));
            station_latitude(i + Nstation_accumulate) = station_latitude_s(i);
            station_longitude(i + Nstation_accumulate)
              = station_longitude_s(i);
            station_altitude(i + Nstation_accumulate) = station_altitude_s(i);

            station_country(i + Nstation_accumulate)
              = trim(station_country_s(i));
            station_type(i + Nstation_accumulate) = trim(station_type_s(i));
            station_network(i + Nstation_accumulate)
              = trim(station_network_s(i));
          }


        // Filters out the stations out of model domain.
        for (i = 0; i < Nstation_array(s); i++)
          if (with_spatial_interpolation)
            {
              StationIndex_x(i + Nstation_accumulate) =
                int((station_longitude(i + Nstation_accumulate)
                     - this->model_x_min) / this->model_Delta_x);
              StationWeight_x(i + Nstation_accumulate) =
                (station_longitude(i + Nstation_accumulate) -
                 StationIndex_x(i + Nstation_accumulate)
                 * this->model_Delta_x - this->model_x_min)
                / this->model_Delta_x;

              StationIndex_y(i + Nstation_accumulate) =
                int((station_latitude(i + Nstation_accumulate)
                     - this->model_y_min) / this->model_Delta_y);
              StationWeight_y(i + Nstation_accumulate) =
                (station_latitude(i + Nstation_accumulate) -
                 StationIndex_y(i + Nstation_accumulate)
                 * this->model_Delta_y - this->model_y_min)
                / this->model_Delta_y;


              // Filters out the stations out of model domain.
              StationAvailability(i + Nstation_accumulate) =
                StationIndex_x(i + Nstation_accumulate) >= 0
                && StationIndex_x(i + Nstation_accumulate)
                < this->model_Nx - 1
                && StationIndex_y(i + Nstation_accumulate) >= 0
                && StationIndex_y(i + Nstation_accumulate)
                < this->model_Ny - 1;
            }
          else
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


        state_index_s(s)
          = Model.GetStateSpeciesIndex(this->obs_species_list[s]);

        // Updates indices.
        Nstation_accumulate += Nstation_array(s);
      }

    // Zero level for station observations.
    state_sequence_index_z = Model.GetStateLevelSequenceIndex(0);

    // Observation availability array initialized with station availability
    // array.
    for (i = 0; i < Nstation; i++)
      Availability(i) = StationAvailability(i);
  }


  //! Initialization of the ground observation manager.
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
  void GroundObservationManager<T>::Init(ClassModel& Model)
  {
    BaseObservationManager<T>::Init(Model);

    ReadStation(Model);
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the total number of stations (that is, for all species).
  /*! Note that the total number of stations is counted with respect to
    species, thus one station may be counted several times if it provides
    observations for several species.
    \return The total number of stations (that is, for all species).
  */
  template <class T>
  int GroundObservationManager<T>::GetNstation() const
  {
    return Nstation;
  }


  //! Returns the array of the number of stations for all species.
  /*!
    \return The array of the number of stations for all species.
  */
  template <class T>
  Array<int, 1>& GroundObservationManager<T>::GetNstationArray()
  {
    return Nstation_array;
  }


  //! Returns the array of station latitude for all species.
  /*!
    \return The array of station latitude for all species.
  */
  template <class T>
  Array<T, 1>& GroundObservationManager<T>::GetStationLatitude()
  {
    return station_latitude;
  }


  //! Returns the array of station longitudes for all species.
  /*!
    \return The array of station longitudes for all species.
  */
  template <class T>
  Array<T, 1>& GroundObservationManager<T>::GetStationLongitude()
  {
    return station_longitude;
  }


  //! Returns the array of station names.
  /*!
    \return The array of station names.
  */
  template <class T>
  Array<string, 1>& GroundObservationManager<T>::GetStationName()
  {
    return station_name;
  }


  //! Returns the flags indicating whether stations are in the model domain.
  /*!
    \return The flags indicating whether stations are in the model domain.
  */
  template <class T>
  Array<bool, 1>& GroundObservationManager<T>::GetStationAvailability()
  {
    return StationAvailability;
  }


  //! Returns the availabilities of stations at given date.
  /*!
    \return The availabilities of stations at given date.
  */
  template <class T>
  Array<bool, 1>& GroundObservationManager<T>::GetAvailability()
  {
    return Availability;
  }


  //! Returns the raw observations from all stations for all species.
  /*!
    \return The raw observations from all stations  for all species.
  */
  template <class T>
  Array<T, 1>& GroundObservationManager<T>::GetStationRawData()
  {
    return station_raw_data;
  }


  //! Returns the array of station types.
  /*!
    \return The array of station types.
  */
  template <class T>
  Array<string, 1>& GroundObservationManager<T>::GetStationType()
  {
    return station_type;
  }


  /////////////////////////////
  // OBSERVATIONS MANAGEMENT //
  /////////////////////////////


  //! Sets the ground observations at a given date.
  /*! It retrieves the observation data from station data files, computes the
    number of observations, and sets the observation error covariance.
    \param date the given date.
  */
  template <class T>
  void GroundObservationManager<T>::SetDate(Date date)
  {

    /*** Hourly observations only ***/

    if (date.GetMinutes() != 0 || date.GetSeconds() != 0)
      {
        this->Nobs = 0;
        this->availability = false;
        this->error_covariance.resize(0, 0);
        return;
      }

    /*** Reads observations from data file ***/

    string str_date, file_date;
    str_date = date.GetDate("%y%m%d%h");
    int int_date = to_num<int>(str_date);
    T value;

    station_raw_data = -999.;

    int s, i;
    int Nstation_accumulate = 0;
    int Nobs_accumulate = 0;

    for (s = 0; s < this->Ns_obs; s++)
      {
        for (i = 0; i < Nstation_array(s); i++)
          {
            if (!StationAvailability(i + Nstation_accumulate))
              continue;

            if (exists(station_file(i + Nstation_accumulate)))
              {
                ExtStream data(station_file(i + Nstation_accumulate));
                do
                  data >> file_date >> value;
                while (to_num<int>(file_date) < int_date && data);

                if (str_date == file_date)
                  station_raw_data(i + Nstation_accumulate) = value;
              }

            Availability(i + Nstation_accumulate)
              = station_raw_data(i + Nstation_accumulate) >= 0;
          }

        // Computes the number of observations.
        this->Nobs_array(s) = 0;
        for (i = 0; i < Nstation_array(s); i++)
          this->Nobs_array(s) += Availability(i + Nstation_accumulate);

        // Sets observation array and observation index array.
        if (this->Nobs_array(s) > 0)
          {
            if (with_spatial_interpolation)
              {
                ObservationWeight_x.resizeAndPreserve(Nobs_accumulate +
                                                      this->Nobs_array(s));
                ObservationWeight_y.resizeAndPreserve(Nobs_accumulate +
                                                      this->Nobs_array(s));
              }

            this->ObservationIndex_x.resizeAndPreserve(Nobs_accumulate +
                                                       this->Nobs_array(s));
            this->ObservationIndex_y.resizeAndPreserve(Nobs_accumulate +
                                                       this->Nobs_array(s));
            this->ObservationIndex_s.resizeAndPreserve(Nobs_accumulate +
                                                       this->Nobs_array(s));
            this->ObservationStateIndex_s
              .resizeAndPreserve(Nobs_accumulate + this->Nobs_array(s));
            this->Observation.resizeAndPreserve(Nobs_accumulate
                                                + this->Nobs_array(s));

            int r = 0;
            for (i = 0; i < Nstation_array(s); i++)
              {
                if (Availability(i + Nstation_accumulate))
                  {
                    if (with_spatial_interpolation)
                      {
                        ObservationWeight_x(r + Nobs_accumulate)
                          = StationWeight_x(i + Nstation_accumulate);
                        ObservationWeight_y(r + Nobs_accumulate)
                          = StationWeight_y(i + Nstation_accumulate);
                      }

                    this->ObservationIndex_x(r + Nobs_accumulate)
                      = StationIndex_x(i + Nstation_accumulate);
                    this->ObservationIndex_y(r + Nobs_accumulate)
                      = StationIndex_y(i + Nstation_accumulate);
                    this->ObservationIndex_s(r + Nobs_accumulate)
                      = this->obs_species_index(s);
                    this->ObservationStateIndex_s(r + Nobs_accumulate)
                      = state_index_s(s);
                    this->Observation(r + Nobs_accumulate)
                      = station_raw_data(i + Nstation_accumulate);
                    r++;
                  }
              }
          }

        // Update indices.
        Nobs_accumulate += this->Nobs_array(s);
        Nstation_accumulate += Nstation_array(s);
      } // Loop over species.


    // Computes the number of all observations.
    this->Nobs = 0;
    for (s = 0; s < this->Ns_obs; s++)
      this->Nobs += this->Nobs_array(s);

    if (this->Nobs > 0)
      {
        this->availability = true;

        this->ObservationIndex_z.resize(this->Nobs);
        this->ObservationIndex_z = 0;

        this->ObservationStateSequenceIndex_z.resize(this->Nobs);
        this->ObservationStateSequenceIndex_z = state_sequence_index_z;

        // Sets the error covariance matrix.
        this->error_covariance.resize(this->Nobs, this->Nobs);
        SetCovariance();

        // Calculates indices in state vector for all available observations.
        if (with_spatial_interpolation)
          {
            ObservationIndex_in_state.resize(this->Nobs, 4);
            for (int r = 0; r < this->Nobs; r++)
              {
                ObservationIndex_in_state(r, 0) =
                  this->ObservationStateIndex_s(r) +
                  this->ObservationStateSequenceIndex_z(r) * this->model_Ny *
                  this->model_Nx + this->ObservationIndex_y(r) *
                  this->model_Nx + this->ObservationIndex_x(r);
                ObservationIndex_in_state(r, 1) =
                  ObservationIndex_in_state(r, 0) + 1;

                ObservationIndex_in_state(r, 2) =
                  this->ObservationStateIndex_s(r) +
                  this->ObservationStateSequenceIndex_z(r) * this->model_Ny *
                  this->model_Nx + (this->ObservationIndex_y(r) + 1) *
                  this->model_Nx + this->ObservationIndex_x(r);
                ObservationIndex_in_state(r, 3) =
                  ObservationIndex_in_state(r, 2) + 1;
              }
            SpatialInterpolation_coefficient.resize(this->Nobs, 4);
            for (int r = 0; r < this->Nobs; r++)
              {
                SpatialInterpolation_coefficient(r, 0) =
                  (1 - ObservationWeight_x(r)) * (1 - ObservationWeight_y(r));
                SpatialInterpolation_coefficient(r, 1) =
                  ObservationWeight_x(r) * (1 - ObservationWeight_y(r));
                SpatialInterpolation_coefficient(r, 2) =
                  (1 - ObservationWeight_x(r)) * ObservationWeight_y(r);
                SpatialInterpolation_coefficient(r, 3) =
                  ObservationWeight_x(r) * ObservationWeight_y(r);
              }
          }
        else
          {
            ObservationIndex_in_state.resize(this->Nobs, 1);
            for (int r = 0; r < this->Nobs; r++)
              ObservationIndex_in_state(r, 0) =
                this->ObservationStateIndex_s(r) +
                this->ObservationStateSequenceIndex_z(r) * this->model_Ny *
                this->model_Nx + this->ObservationIndex_y(r) *
                this->model_Nx + this->ObservationIndex_x(r);
          }
      }
  }


  //! Operator that maps model state to one observation.
  /*! Linear pointwise observation operator.
    \param r index of the observation.
    \param state the model state vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T GroundObservationManager<T>
  ::MltObsOperator(int r, const Array<T, 1>& State)
  {
    if (with_spatial_interpolation)
      return State(ObservationIndex_in_state(r, 0))
        * SpatialInterpolation_coefficient(r, 0) +
        State(ObservationIndex_in_state(r, 1))
        * SpatialInterpolation_coefficient(r, 1) +
        State(ObservationIndex_in_state(r, 2))
        * SpatialInterpolation_coefficient(r, 2) +
        State(ObservationIndex_in_state(r, 3))
        * SpatialInterpolation_coefficient(r, 3);
    else
      return State(ObservationIndex_in_state(r, 0));
  }


  //! Tangent linear operator that maps model state to one observation.
  /*! Linear pointwise observation operator, the same as MltObsOperator.
    \param r index of the observation.
    \param state the model state vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T GroundObservationManager<T>
  ::TLMMltObsOperator(int r, const Array<T, 1>& State)
  {
    if (with_spatial_interpolation)
      return State(ObservationIndex_in_state(r, 0))
        * SpatialInterpolation_coefficient(r, 0) +
        State(ObservationIndex_in_state(r, 1))
        * SpatialInterpolation_coefficient(r, 1) +
        State(ObservationIndex_in_state(r, 2))
        * SpatialInterpolation_coefficient(r, 2) +
        State(ObservationIndex_in_state(r, 3))
        * SpatialInterpolation_coefficient(r, 3);
    else
      return State(ObservationIndex_in_state(r, 0));
  }


  //! Adjoint operator that maps departure vector to one state component.
  /*! Adjoint of linear pointwise observation operator.
    \param i index of the corresponding component in state vector.
    \param Departure the departure vector.
    \return The mapping result (scalar).
  */
  template <class T>
  T GroundObservationManager<T>
  ::ADJMltObsOperator(int i, const Array<T, 1>& Departure)
  {
    T sum = T(0.);
    if (with_spatial_interpolation)
      {
        for (int r = 0; r < this->Nobs; r++)
          for (int j = 0; j < 4; j++)
            if (i == ObservationIndex_in_state(r, j))
              sum += SpatialInterpolation_coefficient(r, j) * Departure(r);
      }
    else
      {
        for (int r = 0; r < this->Nobs; r++)
          if (i == ObservationIndex_in_state(r, 0))
            sum += Departure(r);
      }
    return sum;
  }


  //! Returns one entry of the linearized operator.
  /*!
    \param i row index.
    \param j column index.
    \return Element (i, j) of the linearized operator.
  */
  template <class T>
  T GroundObservationManager<T>::Linearized(int i, int j)
  {
    if (with_spatial_interpolation)
      if (j == ObservationIndex_in_state(i, 0))
        return SpatialInterpolation_coefficient(i, 0);
      else if (j == ObservationIndex_in_state(i, 1))
        return SpatialInterpolation_coefficient(i, 1);
      else if (j == ObservationIndex_in_state(i, 2))
        return SpatialInterpolation_coefficient(i, 2);
      else if (j == ObservationIndex_in_state(i, 3))
        return SpatialInterpolation_coefficient(i, 3);
      else
        return T(0);
    else if (j == ObservationIndex_in_state(i, 0))
      return T(1);
    else
      return T(0);
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Sets error covariance.
  template <class T>
  void GroundObservationManager<T>::SetCovariance()
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


#define POLYPHEMUS_FILE_OBSERVATION_GROUNDOBSERVATIONMANAGER_CXX
#endif
