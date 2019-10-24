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


#ifndef POLYPHEMUS_FILE_OBSERVATION_SIMOBSERVATIONMANAGER_HXX


#include <iostream>
#include "BaseObservationManager.cxx"


namespace Polyphemus
{


  using namespace std;


  ///////////////////////////
  // SIMOBSERVATIONMANAGER //
  ///////////////////////////


  //! This class manages synthetic observations.
  template<class T>
  class SimObservationManager: public BaseObservationManager<T>
  {

  protected:

    /*! \brief Specifies whether observations are provided over 2D-layers, at
      given stations or at given grid points.
    */
    string simulation_option;

    /*** Simulation data file attributes ***/

    //! Name of the data files that store the model simulation results.
    vector<string> sim_input_file_list;
    //! Starting date for the simulation results in data files.
    Date sim_date_min;
    //! Time step in seconds for the simulation results in data files.
    T sim_delta_t;
    //! Levels list for the simulated data in files.
    Array<int, 1> sim_levels;
    //! Number of levels for the simulated data in files.
    int Nz_sim;

    /*** Levels ***/

    //! Number of levels for the simulated observations.
    int Nz_obs;
    //! Levels list for the simulated observations.
    Array<int, 1> obs_levels;

    /*** Stations ***/

    //! Total number of stations (that is, for all species).
    int Nstation;
    //! Array of the number of stations for all species.
    Array<int, 1> Nstation_array;
    //! For all species, names of the files that describe the stations.
    vector<string> station_file_list;
    //! For all species, paths to the station data files.
    vector<string> input_dir_list;

    //! Total number of observations read from data files.
    int Nobs_total;

    /*! Array that stores the corresponding indices in state vector for all
      available observations. */
    Array<int, 1> ObservationIndex_in_state;

  public:

    /*** Constructor and destructor ***/

    SimObservationManager(string config_file);
    virtual ~SimObservationManager();

    /*** Configuration ***/

    template <class ClassModel>
    void ReadConfiguration(ClassModel& Model);

    /*** Initialization ***/

    template <class ClassModel>
    void ReadStation(ClassModel& Model);
    template <class ClassModel>
    void Init(ClassModel& Model);

    /*** Observations management ***/

    void SetDate(Date date);

    T MltObsOperator(int r, const Array<T, 1>& State);
    T TLMMltObsOperator(int r, const Array<T, 1>& State);
    T ADJMltObsOperator(int i, const Array<T, 1>& Departure);
    T Linearized(int i, int j);

  private:

    void SetCovariance();

  };


} // namespace Polyphemus


#define POLYPHEMUS_FILE_OBSERVATION_SIMOBSERVATIONMANAGER_HXX
#endif
