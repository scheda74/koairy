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


#ifndef POLYPHEMUS_FILE_OBSERVATION_BASEOBSERVATIONMANAGER_CXX


#include "BaseObservationManager.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template <class T>
  BaseObservationManager<T>::BaseObservationManager(string config_file):
    config(config_file), availability(false)
  {
  }


  //! Destructor.
  template <class T>
  BaseObservationManager<T>::~BaseObservationManager()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads observed species list and observation error parameters.
   */
  template <class T>
  void BaseObservationManager<T>::ReadConfiguration()
  {
    string obs_manager_file;

    this->config.SetSection("[observation_management]");

    this->config.PeekValue("Configuration_file", obs_manager_file);

    /*** Observation configuration ***/

    ConfigStream config_obs(obs_manager_file);

    config_obs.SetSection("[general]");

    config_obs.Find("Species");
    obs_species_list = split(config_obs.GetLine());
    Ns_obs = int(obs_species_list.size());

    vector<T> error_variance_vector;
    config_obs.Find("Error_variance");
    split(config_obs.GetLine(), error_variance_vector);
    int Ns_err = int(error_variance_vector.size());
    if (Ns_obs != Ns_err)
      throw string("In configuration file \"") + obs_manager_file +
        string("\", the number of observed species (") + to_str(Ns_obs) +
        string(") is not equal to the number of error variances.");
    error_variance_array.resize(Ns_err);
    for (int i = 0; i < Ns_err; i++)
      error_variance_array(i) = error_variance_vector[i];
  }


  ////////////////////
  // INITIALIZATION //
  ////////////////////


  //! Initialization for observation manager.
  /*!
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
    <li> GetNstate()
    </ul>
  */
  template <class T>
  template <class ClassModel>
  void BaseObservationManager<T>::Init(ClassModel& Model)
  {
    // Model domain.
    model_x_min = Model.GetX_min();
    model_y_min = Model.GetY_min();
    model_Delta_x = Model.GetDelta_x();
    model_Delta_y = Model.GetDelta_y();
    model_Nx = Model.GetNx();
    model_Ny = Model.GetNy();

    // Information on model state variable.
    Nz_state = Model.GetStateNz();
    state_levels.resize(Nz_state);
    state_levels = Model.GetStateLevels();

    Ns_state = Model.GetStateNs();
    state_species_list = Model.GetStateSpeciesList();

    Nstate = Model.GetNstate();

    // Observation manager configuration.
    this->ReadConfiguration();

    // Observation settings.
    Nobs = 0;
    obs_species_index.resize(Ns_obs);
    for (int s = 0; s < Ns_obs; s++)
      obs_species_index(s) = Model.GetSpeciesIndex(obs_species_list[s]);
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the observation number.
  /*!
    \return Observation dimension.
  */
  template <class T>
  int BaseObservationManager<T>::GetNobs() const
  {
    return Nobs;
  }


  //! Returns the number of observed species.
  /*!
    \return Number of observed species.
  */
  template <class T>
  int BaseObservationManager<T>::GetNs() const
  {
    return Ns_obs;
  }


  //! Checks whether observations are available.
  /*!
    \return True if observations are available, false otherwise.
  */
  template <class T>
  bool BaseObservationManager<T>::IsAvailable() const
  {
    return availability;
  }


  //! Returns the values of observations.
  /*!
    \return The values of observations.
  */
  template <class T>
  const Array<T, 1>& BaseObservationManager<T>::GetObservation() const
  {
    return Observation;
  }


  //! Returns the values of observations.
  /*!
    \return The values of observations.
  */
  template <class T>
  Array<T, 1>& BaseObservationManager<T>::GetObservation()
  {
    return Observation;
  }


  //! Returns the values of observation error covariance.
  /*!
    \return Observation error covariance.
  */
  template <class T>
  const Array<T, 2>& BaseObservationManager<T>::GetCovariance() const
  {
    return error_covariance;
  }


  //! Returns the values of observation error covariance.
  /*!
    \return Observation error covariance.
  */
  template <class T>
  Array<T, 2>& BaseObservationManager<T>::GetCovariance()
  {
    return error_covariance;
  }


  /*! \brief Returns indices (along x) of the closest model-grid-points for
    observations (at current date).
  */
  /*!
    \return The indices (along x) of the closest model-grid-points for
    observations (at current date).
  */
  template <class T>
  Array<int, 1>& BaseObservationManager<T>::GetObservationIndex_x()
  {
    return ObservationIndex_x;
  }


  /*! \brief Returns indices (along y) of the closest model-grid-points for
    observations (at current date).
  */
  /*!
    \return The indices (along y) of the closest model-grid-points for
    observations (at current date).
  */
  template <class T>
  Array<int, 1>& BaseObservationManager<T>::GetObservationIndex_y()
  {
    return ObservationIndex_y;
  }


  /*! \brief Returns indices (along z) of the closest model-grid-points for
    observations (at current date).
  */
  /*!
    \return The indices (along z) of the closest model-grid-points for
    observations (at current date).
  */
  template <class T>
  Array<int, 1>& BaseObservationManager<T>::GetObservationIndex_z()
  {
    return ObservationIndex_z;
  }


  /*! \brief Returns species indices of the observations (at current date).
   */
  /*!
    \return The species indices of the observations (at current date).
  */
  template <class T>
  Array<int, 1>& BaseObservationManager<T>::GetObservationIndex_s()
  {
    return ObservationIndex_s;
  }


  /*! \brief Returns indices of species in state array for all observations
    (at current date).
  */
  /*!
    \return The indices of species in state array for all observations
    (at current date)
  */
  template <class T>
  Array<int, 1>& BaseObservationManager<T>::GetObservationStateIndex_s()
  {
    return ObservationStateIndex_s;
  }


  /*! \brief Returns the array of level storage sequence indices in state
    array for all observations (at current date).
  */
  /*!
    \return The array of level storage sequence indices in state array for
    all observations (at current date).
  */
  template <class T>
  Array<int, 1>&
  BaseObservationManager<T>::GetObservationStateSequenceIndex_z()
  {
    return ObservationStateSequenceIndex_z;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OBSERVATION_BASEOBSERVATIONMANAGER_CXX
#endif
