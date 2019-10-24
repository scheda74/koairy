// Copyright (C) 2010-2011, INRIA
// Author(s): Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_OBSERVATION_GROUNDNETWORKOBSERVATIONMANAGER_HXX


#include "TalosHeader.hxx"
using namespace Talos;

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#include "VerdandiHeader.hxx"

#include "seldon/vector/Vector2.hxx"
#include "seldon/matrix_sparse/Matrix_ArraySparse.hxx"


namespace Polyphemus
{

  using namespace Verdandi;
  using Seldon::Vector;
  using Verdandi::is_multiple;


  //! This class manages observations provided by an air quality network.
  /*!
    \tparam T the type of floating-point numbers.
  */
  template <class T>
  class GroundNetworkObservationManager: public Verdandi::VerdandiBase
  {
  public:
    //! Type of the tangent linear operator.
    typedef Matrix<T> tangent_linear_operator;
    //! Type of the observation error covariance matrix.
    typedef Matrix<T> error_variance;
    //! Type of a row of the tangent linear operator.
    typedef Vector<T> tangent_linear_operator_row;
    //! Type of the observation vector.
    typedef Vector<T> observation;

  protected:
    //! First available date.
    Date date_begin_;
    //! Last available date plus one time step.
    Date date_end_;
    //! Time step in seconds.
    T Delta_t_;
    //! Number of time steps.
    int Nt_;

    /*! \brief Is a structured mesh used to define the observation operator
      and to filter out the stations outside the corresponding domain? */
    bool with_domain_;
    //! Abscissa of the domain lower-left corner.
    T x_min_;
    //! Ordinate of the domain lower-left corner.
    T y_min_;
    T Delta_x_;
    T Delta_y_;
    T Nx_;
    T Ny_;
    //! Abscissa of the domain upper-right corner.
    T x_max_;
    //! Ordinate of the domain upper-right corner.
    T y_max_;

    //! Network type.
    string network_type_;

    //! Path to the file that describes the stations.
    string station_description_file_;
    /*! \brief Path to the station data file(s). The special markup
      "&l" is replaced with the short name of the station. */
    string station_data_file_;

    //! List of stations.
    // To be implemented.
    //! Station complete names.
    vector<string> station_name_;
    //! Station short names, used in data file names.
    vector<string> station_short_name_;
    //! Abscissa of stations.
    Vector<T> station_x_;
    //! Ordinate of stations.
    Vector<T> station_y_;
    //! Total number of stations.
    int Nstation_;

    /*! Mode: "none", "hourly", "daily" or "daily_peak".
      - "none": keeps the data as read.
      - "daily": compute daily averages.
      - "hourly": retrieves values at a given hour ('hourly_selected_hour_').
      - "daily_peak": compute daily peaks with data in ['hour_begin',
      'hour_end'[ if at least 'Ndaily_peak_hour_' values are available.
    */
    string mode_;
    //! In case of hourly data, the hour to select.
    int hourly_selected_hour_;
    /*! \brief Minimum number of available values per day required to compute
      a meaningful daily average. */
    int Ndaily_average_;
    //! First hour at which a daily peak is sought.
    int daily_peak_hour_begin_;
    //! Last hour, plus one, at which a daily peak is sought.
    int daily_peak_hour_end_;
    /*! \brief Minimum number of available data in [daily_peak_hour_begin_,
      daily_peak_hour_end_[ for the daily peak to be reliable. */
    int Ndaily_peak_hour_;
    /*! \brief First hour in the day (in [0, 23]). It has impact on the daily
      averages and on the daily peaks since 'daily_peak_hour_begin_' and
      'daily_peak_hour_end_' are relative to them. */
    int first_hour_in_day_;
    //! Values put in place of unknown data.
    Vector<T> unknown_value_;

    /*** For data assimilation ***/

    //! Variance of observation errors.
    T error_variance_value_;

    /*** Observations ***/

    //! In order to deactivate the observations.
    bool with_observation_;

    //! Date at which the observations may be requested.
    Date date_;

    //! Time index of the current date.
    int time_index_;

    //! Station indexes of all observations.
    Vector2<int> all_station_index_;

    //! All observations.
    Vector2<T> all_observation_;

    //! Interpolation indices for all stations.
    Matrix<int> interpolation_index_;

    //! Interpolation weights for all stations.
    Matrix<T> interpolation_weight_;

    //! Size of the model state (to define the observation operator).
    int model_Nstate_;

  public:
    // Constructors.
    GroundNetworkObservationManager();
    GroundNetworkObservationManager(string configuration_filename);
    void Initialize(string configuration_filename);
    void Initialize(Verdandi::VerdandiOps& configuration);
    template <class Model>
    void Initialize(Model& model, string configuration_filename);

    void LoadStationInformation(string description_file);

    // Loading the data.
    void LoadStationData(vector<string>& file_data);
    void LoadStationData(vector<string>& file_data,
                         Vector2<int>& station_index_, Vector2<T>& data_);
    // Loads all observations in memory.
    void LoadAllObservation();
    // Prepares to load observations at a given date.
    template <class Model>
    void SetDate(Model& model, string date);
    void SetDate(string date);
    template <class Model>
    void SetTime(Model& model, double time);
    void SetTime(string date);

    // Provides the observations at the selected date.
    bool HasObservation() const;
    int GetNobservation() const;
    void GetObservation(Vector<T>& data) const;
    void GetObservation(Vector<int>& station_index, Vector<T>& data) const;
    void GetObservation(Vector<T>& abscissa, Vector<T>& ordinate,
                        Vector<T>& data) const;
    template <class state>
    void GetInnovation(const state& x, observation& innovation);
    Vector2<int>& GetAllStationIndex();
#ifndef SWIG
    const Vector2<int>& GetAllStationIndex() const;
#endif
    Vector2<T>& GetAllObservation();
#ifndef SWIG
    const Vector2<T>& GetAllObservation() const;
#endif

    template <class state>
    void ApplyOperator(const state& x, observation& y) const;
#ifndef SWIG
    Matrix<T, General, RowSparse> GetOperator() const;
#endif
    T GetTangentLinearOperator(int i, int j) const;
    void GetTangentLinearOperatorRow(int row, tangent_linear_operator_row&
                                     tangent_operator_row) const;
#ifndef SWIG
    tangent_linear_operator GetTangentLinearOperator() const;
#endif

    T GetErrorVariance(int i, int j) const;
#ifndef SWIG
    error_variance GetErrorVariance() const;
#endif

    // Processing.
    void FilterTimeSequence(Vector<int>& time_index, Vector<T>& data);
    void KeepDailyPeak(Vector<int>& time_index, Vector<T>& data);
    void KeepHour(int hour, Vector<int>& time_index, Vector<T>& data);
    void KeepHour(int hour, Vector2<int>& station_index_, Vector2<T>& data_);
    void ComputeDailyAverage(Vector<int>& time_index, Vector<T>& data);
    void FromHourlyToDaily(Vector2<int>& station_index_, Vector2<T>& data_,
                           int selected_hour = -1);
    void FromDailyToHourly(int hour,
                           Vector2<int>& station_index_, Vector2<T>& data_);

    // Access methods.
    int GetNstation() const;
    T GetXMin() const;
    T GetYMin() const;
    T GetXMax() const;
    T GetYMax() const;
    void SetNstate(int Nstate);
    Date GetDateBegin() const;
    Date GetDateEnd() const;
    T GetTimeStep() const;
    const vector<string>& GetStationName() const;
    string GetStationName(int i) const;
    const vector<string>& GetStationShortName() const;
    string GetStationShortName(int i) const;
    int GetStationIndex(string name, bool is_short_name = true) const;
    void GetCoordinate(int station, Vector<T>& coordinate) const;
    Vector<T>& GetStationX();
    Vector<T>& GetStationY();
  };


}


#define POLYPHEMUS_FILE_OBSERVATION_GROUNDNETWORKOBSERVATIONMANAGER_HXX
#endif
