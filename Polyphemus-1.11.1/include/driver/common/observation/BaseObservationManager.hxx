// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_OBSERVATION_BASEOBSERVATIONMANAGER_HXX


#include <iostream>
#include "AtmoData.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////
  // BASEOBSERVATIONMANAGER //
  ////////////////////////////


  /*! \brief This class is the base class for observation managers. It defines
    interfaces for all observation managers, and reads model
    information. Manager configurations are read in derived classes.
  */
  template<class T>
  class BaseObservationManager
  {

  protected:

    //! Configuration stream.
    ConfigStream config;

    /*** Observations ***/

    //! Observation data.
    Array<T, 1> Observation;
    //! Index along x in the model grid domain for all available observations.
    Array<int, 1> ObservationIndex_x;
    //! Index along y in the model grid domain for all available observations.
    Array<int, 1> ObservationIndex_y;
    //! Level index in the model grid domain for all available observations.
    Array<int, 1> ObservationIndex_z;
    //! Species index in the model for all available observations.
    Array<int, 1> ObservationIndex_s;

    /*! \brief Array of level storage sequence indices in model state array
      for all observations.
    */
    Array<int, 1> ObservationStateSequenceIndex_z;
    //! The species index in model state array for all available observations.
    Array<int, 1> ObservationStateIndex_s;

    /*! \brief Number of total observations (that is, with all species) at
      current date.
    */
    int Nobs;
    //! Array of the number of observations for all species.
    Array<int, 1> Nobs_array;
    //! Availability of observations at current date.
    bool availability;

    //! Array of observation error variance for all species.
    Array<T, 1> error_variance_array;
    //! Observation error covariance matrix.
    Array<T, 2> error_covariance;

    //! Number of observed species.
    int Ns_obs;
    //! Observed species list.
    vector<string> obs_species_list;
    //! Global index array for observed species.
    Array<int, 1> obs_species_index;


    /*** Model domain ***/

    //! Model origin along x.
    T model_x_min;
    //! Model space increment along x.
    T model_Delta_x;
    //! Number of points along x in the model.
    int model_Nx;
    //! Model origin along y.
    T model_y_min;
    //! Model space increment along y.
    T model_Delta_y;
    //! Number of points along y in the model.
    int model_Ny;
    //! Number of levels in the model state.
    int Nz_state;
    //! Levels list of model state.
    Array<int, 1> state_levels;

    //! State species list.
    vector<string> state_species_list;
    //! Number of state species.
    int Ns_state;

    //! Dimension of state variable.
    int Nstate;

  public:

    /*** Constructor ***/

    BaseObservationManager(string config_file);
    virtual ~BaseObservationManager();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    /*** Access methods ***/

    int GetNobs() const;
    int GetNs() const;
    bool IsAvailable() const;

    const Array<T, 1>& GetObservation() const;
    Array<T, 1>& GetObservation();
    const Array<T, 2>& GetCovariance() const;
    Array<T, 2>& GetCovariance();

    Array<int, 1>& GetObservationIndex_x();
    Array<int, 1>& GetObservationIndex_y();
    Array<int, 1>& GetObservationIndex_z();
    Array<int, 1>& GetObservationIndex_s();
    Array<int, 1>& GetObservationStateIndex_s();
    Array<int, 1>& GetObservationStateSequenceIndex_z();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OBSERVATION_BASEOBSERVATIONMANAGER_HXX
#endif
