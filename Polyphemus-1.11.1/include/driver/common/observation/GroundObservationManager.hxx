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


#ifndef POLYPHEMUS_FILE_OBSERVATION_GROUNDOBSERVATIONMANAGER_HXX


#include "BaseObservationManager.cxx"
#include <iostream>


namespace Polyphemus
{


  using namespace std;


  //////////////////////////////
  // GROUNDOBSERVATIONMANAGER //
  //////////////////////////////


  //! This class manages ground observations for a single species.
  template<class T>
  class GroundObservationManager: public BaseObservationManager<T>
  {

  protected:

    /*** Stations ***/

    //! Total number of stations (that is, for all species).
    int Nstation;
    //! Array of the number of stations for all species.
    Array<int, 1> Nstation_array;
    //! For all species, names of the files that describe the stations.
    vector<string> station_file_list;
    //! For all species, paths to the station data file.
    vector<string> input_dir_list;

    //! Files associated with the stations.
    Array<string, 1> station_file;
    //! Names of the stations.
    Array<string, 1> station_name;
    //! Codes of the stations.
    Array<string, 1> station_code;

    //! Raw observations from stations, without any filtering.
    Array<T, 1> station_raw_data;
    //! Stations latitudes.
    Array<T, 1> station_latitude;
    //! Stations longitudes.
    Array<T, 1> station_longitude;
    //! Stations altitudes.
    Array<T, 1> station_altitude;
    //! Stations countries.
    Array<string, 1> station_country;
    //! Stations types.
    Array<string, 1> station_type;
    //! Stations networks.
    Array<string, 1> station_network;

    //! Are stations in the model domain?
    Array<bool, 1> StationAvailability;
    //! Availabilities of stations at current date.
    Array<bool, 1> Availability;

    /*! \brief Index of the left corner along x of the enclosing cell in the
      model, for all stations.
    */
    Array<int, 1> StationIndex_x;
    /*! \brief Index of the left corner along y of the enclosing cell in the
      model, for all stations.
    */
    Array<int, 1> StationIndex_y;

    //! Array of indices in model state array for all observed species.
    Array<int, 1> state_index_s;
    /*! \brief Level storage sequence index in model state array for all
      observed species.
    */
    int state_sequence_index_z;

    /*** Spatial interpolation ***/

    //! Flag for spatial interpolation.
    bool with_spatial_interpolation;
    //! Spatial interpolation weights along x.
    Array<T, 1> StationWeight_x;
    //! Spatial interpolation weights along y.
    Array<T, 1> StationWeight_y;
    //! Spatial interpolation weights along x for observations.
    Array<T, 1> ObservationWeight_x;
    //! Spatial interpolation weights along y for observations.
    Array<T, 1> ObservationWeight_y;

    /*! Array that stores the corresponding indices in state vector for all
      available observations. */
    Array<int, 2> ObservationIndex_in_state;
    /*! Array that stores the spatial interpolation coefficients for all
      available observations. */
    Array<T, 2> SpatialInterpolation_coefficient;

  public:

    /*** Constructor and destructor ***/

    GroundObservationManager(string config_file);
    virtual ~GroundObservationManager();

    /*** Configurations ***/

    virtual void ReadConfiguration();

    /*** Initialization ***/

    template <class ClassModel>
    void ReadStation(ClassModel& Model);

    template <class ClassModel>
    void Init(ClassModel& Model);

    /*** Access methods ***/

    int GetNstation() const;
    Array<int, 1>& GetNstationArray();
    Array<T, 1>& GetStationLatitude();
    Array<T, 1>& GetStationLongitude();
    Array<string, 1>& GetStationName();
    Array<bool, 1>& GetStationAvailability();
    Array<bool, 1>& GetAvailability();
    Array<T, 1>& GetStationRawData();
    Array<string, 1>& GetStationType();

    /*** Observations managements ***/

    void SetDate(Date date);
    T MltObsOperator(int r, const Array<T, 1>& State);
    T TLMMltObsOperator(int r, const Array<T, 1>& State);
    T ADJMltObsOperator(int i, const Array<T, 1>& Departure);
    T Linearized(int i, int j);

  protected:

    void SetCovariance();

  };


} // namespace Polyphemus


#define POLYPHEMUS_FILE_OBSERVATION_GROUNDOBSERVATIONMANAGER_HXX
#endif
