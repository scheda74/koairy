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


#ifndef POLYPHEMUS_FILE_OBSERVATION_GROUNDNETWORKOBSERVATIONMANAGER_CXX


#include "GroundNetworkObservationManager.hxx"

#include "Talos.hxx"
#include "Verdandi.hxx"
#include "seldon/vector/Vector2.cxx"
#include "seldon/matrix_sparse/Matrix_ArraySparse.cxx"


namespace Polyphemus
{
  using Talos::to_str;
  using Talos::split;
  using Talos::is_num;
  using Talos::to_num;
  using Talos::trim;
  using Talos::find_replace;
  using Seldon::Error;


  //////////////////
  // CONSTRUCTORS //
  //////////////////


  //! Default constructor.
  template <class T>
  GroundNetworkObservationManager<T>::GroundNetworkObservationManager()
  {
  }


  //! Main constructor.
  /*!
    \param[in] configuration_filename path to the configuration file.
    \param[in] type type of network data.
    \param[in] index index associated with the network data.
  */
  template <class T>
  GroundNetworkObservationManager<T>
  ::GroundNetworkObservationManager(string configuration_filename)
  {
    Initialize(configuration_filename);
  }


  //! Main constructor.
  /*! It loads the main information about the network.
    \param[in] configuration_filename path to the configuration file.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::Initialize(string configuration_filename)
  {
    Verdandi::VerdandiOps configuration(configuration_filename);
    configuration.SetPrefix("ground_network_observation.");
    Initialize(configuration);
  }


  //! Main constructor.
  /*! It loads the main information about the network.
    \param[in] configuration the configuration available in an Ops object.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::Initialize(Verdandi::VerdandiOps& configuration)
  {

    /*** Mandatory fields ***/

    configuration.Set("option.with_observation", with_observation_);
    configuration.Set("error.variance", error_variance_value_);

    if (!with_observation_)
      return;

    /*** Input ***/

    configuration.Set("input.network_type",
                      "ops_in(v, {'default', 'none'})", network_type_);

    if (network_type_ == "default")
      {
        // Stations description and data files.
        configuration.Set("input.station_description_file",
                          station_description_file_);
        configuration.Set("input.station_data_file", station_data_file_);
      }

    // Time step of the input data.
    configuration.Set("input.Delta_t", Delta_t_);

    configuration.Set("input.unknown_values", unknown_value_);

    /*** Output ***/

    // Data can be filtered in time.
    if (configuration.Is<string>("output.mode"))
      {
        configuration.Set("output.mode",
                          "ops_in(v, {'none', 'daily', 'daily_peak'})",
                          mode_);
        if (mode_ == "daily")
          Ndaily_average_ = 24;
        else if (mode_ == "daily_peak")
          {
            daily_peak_hour_begin_ = 0;
            daily_peak_hour_end_ = 25;
            Ndaily_peak_hour_ = 24;
          }
      }
    else // it is assumed the entry is a table, and the first element is a
      // string.
      {
        configuration.Set("output.mode[1]",
                          "ops_in(v, {'hourly', 'daily', 'daily_peak',})",
                          mode_);
        if (mode_ == "hourly")
          configuration.Set("output.mode[2]", "v >= 0 and v < 24",
                            hourly_selected_hour_);
        else if (mode_ == "daily")
          configuration.Set("output.mode[2]", "v >= 0 and v <= 24",
                            Ndaily_average_);
        else if (mode_ == "daily_peak")
          {
            configuration.Set("output.mode[2]", "v >= 0 and v < 24",
                              daily_peak_hour_begin_);
            configuration.Set("output.mode[3]",  "v > "
                              + to_str(daily_peak_hour_begin_)
                              + " and v < 25",
                              daily_peak_hour_end_);
            configuration.Set("output.mode[4]", "v > 0 and v <= "
                              + to_str(daily_peak_hour_end_
                                       - daily_peak_hour_begin_),
                              daily_peak_hour_end_ - daily_peak_hour_begin_,
                              Ndaily_peak_hour_);
          }
      }

    if (Delta_t_ != T(3600) && mode_ != "none")
      throw Error("GroundNetworkObservationManager::Init",
                  "Only hourly input data can be handled, "
                  "unless the output mode is not \"none\".");

    if (mode_ == "daily" || mode_ == "daily_peak")
      configuration.Set("output.first_hour_in_day", "v >= -24 and v <= 24",
                        first_hour_in_day_);
    else
      first_hour_in_day_ = 0.;

    // Time period.
    string tmp_date;
    date_begin_ = configuration.Get<string>("output.date_begin");
    date_end_ = configuration.Get<string>("output.date_end");
    Nt_ = int(T(date_end_.GetSecondsFrom(date_begin_))
              / Delta_t_ + 0.5);

    // Horizontal output.
    configuration.Set("output.with_domain", with_domain_);

    if (with_domain_)
      {
        configuration.Set("output.x_min", x_min_);
        configuration.Set("output.y_min", y_min_);
        configuration.Set("output.Delta_x", Delta_x_);
        configuration.Set("output.Delta_y", Delta_y_);
        configuration.Set("output.Nx", Nx_);
        configuration.Set("output.Ny", Ny_);
        x_max_ = x_min_ + (Nx_ - 1.) * Delta_x_;
        y_max_ = y_min_ + (Ny_ - 1.) * Delta_y_;
      }

    /*** Load stations information and observations ***/

    if (network_type_ == "default")
      {
        LoadStationInformation(station_description_file_);
        vector<string> file_list(Nstation_);
        for (int s = 0; s < Nstation_; s++)
          file_list[s] = find_replace(station_data_file_, "&l",
                                      station_short_name_[s]);
        LoadStationData(file_list);
      }
  }


  //! Constructor.
  /*! It loads the main information about the network.
    \param[in] configuration_filename path to the configuration file.
  */
  template <class T>
  template <class Model>
  void GroundNetworkObservationManager<T>
  ::Initialize(Model& model, string configuration_filename)
  {
    model_Nstate_ = model.GetNstate();
    Initialize(configuration_filename);
  }


  //! Loads the stations information.
  /*! It loads the main information about the network. It filters out the
    stations that are outside the domain (if a domain has been defined). It
    also pre-computes the observation operator.
    \param[in] description_file the file that describes the network.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::LoadStationInformation(string description_file)
  {

    /*** Station locations ***/

    if (!exists(description_file))
      throw Error("GroundNetworkObservationManager::LoadStationInformation",
                  "Unable to open file \"" + description_file + "\".");
    ExtStream description_stream(description_file);
    string station_line, short_name, name;
    vector<string> split_line;
    string x_string, y_string;
    T x, y;
    Nstation_ = 0;
    while (description_stream.GetLine(station_line))
      {
        split_line = split(station_line, ";");
        short_name = trim(split_line[0]);
        name = trim(split_line[1]);
        y_string = trim(split_line[2]);
        x_string = trim(split_line[3]);
        if (!is_num(x_string) || !is_num(y_string))
          continue;
        to_num(y_string, y);
        to_num(x_string, x);

        // Filters out the stations outside the selected domain.
        if (!with_domain_ || x >= x_min_ && x <= x_max_
            && y >= y_min_ && y <= y_max_)
          {
            Nstation_++;
            station_x_.Append(x);
            station_y_.Append(y);
            station_short_name_.push_back(short_name);
            station_name_.push_back(name);
          }
      }

    /*** For observation operator ***/

    if (with_domain_ && with_observation_)
      {
        interpolation_index_.Reallocate(Nstation_, 4);
        interpolation_weight_.Reallocate(Nstation_, 4);

        int i, j;
        T weight_x, weight_y;
        for (int station = 0; station < Nstation_; station++)
          {
            i = int((station_x_(station) - x_min_) / Delta_x_);
            j = int((station_y_(station) - y_min_) / Delta_y_);
            interpolation_index_(station, 0) = j * int(Nx_) + i;
            interpolation_index_(station, 1) = j * int(Nx_) + i + 1;
            interpolation_index_(station, 2) = (j + 1) * int(Nx_) + i;
            interpolation_index_(station, 3) = (j + 1) * int(Nx_) + i + 1;

            weight_x = (station_x_(station) - x_min_
                        - T(i) * Delta_x_) / Delta_x_;
            weight_y = (station_y_(station) - y_min_
                        - T(j) * Delta_y_) / Delta_y_;
            interpolation_weight_(station, 0) = (T(1) - weight_x)
              * (T(1) - weight_y);
            interpolation_weight_(station, 1) = (T(1) - weight_x) * weight_y;
            interpolation_weight_(station, 2) = weight_x * (T(1) - weight_y);
            interpolation_weight_(station, 3) = weight_x * weight_y;
          }
      }
  }


  //////////////////////
  // LOADING THE DATA //
  //////////////////////


  //! Loads station data from data files.
  /*! It reads the data available in a list of files.
    \param[in] file_list the paths to the data files.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::LoadStationData(vector<string>& file_list)
  {
    LoadStationData(file_list, all_station_index_, all_observation_);
  }


  //! Loads station data from data files.
  /*! It reads the data available in a list of files.
    \param[in] file_list the paths to the data files.
    \param[out] station_index the station indexes, indexed by the time and the
    station.
    \param[out] data the values, indexed by the time and the station.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::LoadStationData(vector<string>& file_list,
                    Vector2<int>& station_index, Vector2<T>& data)
  {
    Nstation_ = int(file_list.size());

    station_index.Clear();
    station_index.Reallocate(Nt_);
    data.Clear();
    data.Reallocate(Nt_);

    int date_begin_int = to_num<int>(date_begin_.GetDate("%y%m%d%h"));
    Date last_date = date_end_;
    last_date.AddSeconds(-Delta_t_);
    int last_date_int = to_num<int>(last_date.GetDate("%y%m%d%h"));

    ifstream data_stream;
    T tmp_value;
    int tmp_date;
    // Temporary vector that stores the data at one station, and the
    // associated date indices.
    Vector<T> data_station(Nt_);
    Vector<int> time_index(Nt_);
    int t;
    for (int l = 0; l < Nstation_; l++)
      {
        data_station.Reallocate(Nt_);
        time_index.Reallocate(Nt_);
        // All dates are supposed to be there in first place.
        time_index.Fill();

        if (!exists(file_list[l]))
          continue;
        data_stream.open(file_list[l].c_str());
        // Searches for 'date_begin_' in the dates.
        while (data_stream >> tmp_date && tmp_date < date_begin_int)
          data_stream >> tmp_value;
        // Now starting reading the data of interest.
        t = 0;
        do
          {
            data_stream >> data_station(t);
            t++;
          }
        while (t < Nt_ && data_stream >> tmp_date
               && data_stream.good());
        if (t != Nt_)
          throw Error("GroundNetworkObservationManager::LoadStationData",
                      "There must be missing dates in \"" + file_list[l]
                      + "\": only " + to_str(t) + " date(s) was/were found, "
                      "instead of " + to_str(Nt_) + ".");
        if (tmp_date != last_date_int)
          throw Error("GroundNetworkObservationManager::LoadStationData",
                      "There must be missing dates in \"" + file_list[l]
                      + "\": the " + to_str(Nt_) + "th date "
                      + "that was read is " + to_str(tmp_date) + " instead of"
                      + last_date.GetDate(" %y-%m-%d %h:%i:%s"));
        data_stream.close();
        // Copies the values to 'data', possibly applying with some filter.
        FilterTimeSequence(time_index, data_station);
        for (t = 0; t < time_index.GetLength(); t++)
          {
            data(time_index(t)).Append(data_station(t));
            station_index(time_index(t)).Append(l);
          }
      }

    if (mode_ == "daily_peak" || mode_ == "daily" || mode_ == "hourly")
      FromHourlyToDaily(station_index, data);
  }


  //! Loads all observations in memory.
  /*! The observations of the whole time period are loaded. */
  template <class T>
  void GroundNetworkObservationManager<T>::LoadAllObservation()
  {
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] model this argument is not used.
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  template <class Model>
  void GroundNetworkObservationManager<T>
  ::SetDate(Model& model, string date)
  {
    model_Nstate_ = model.GetNstate();
    SetDate(date);
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  void GroundNetworkObservationManager<T>::SetDate(string date)
  {
    Date date_obj(date);
    date_ = date_obj;

    if (with_observation_
        && is_multiple(date_.GetSecondsFrom(date_begin_), Delta_t_))
      time_index_ = int(date_.GetSecondsFrom(date_begin_) / Delta_t_);
    else
      time_index_ = -1; // No observations.
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] model this argument is not used.
    \param[in] time the time at which the observations may be requested.
  */
  template <class T>
  template <class Model>
  void GroundNetworkObservationManager<T>
  ::SetTime(Model& model, double time)
  {
    model_Nstate_ = model.GetNstate();
    Date date = model.GetDate(time);
    SetDate(date.GetDate("%y%m%d%h%i"));
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  void GroundNetworkObservationManager<T>::SetTime(string date)
  {
    SetDate(date);
  }


  //! Checks whether observations are available.
  /*!
    \return True if observations are available, false otherwise.
  */
  template <class T>
  bool GroundNetworkObservationManager<T>::HasObservation() const
  {
    return GetNobservation() != 0;
  }


  //! Returns the number of available observations.
  /*!
    \return The number of available observations.
  */
  template <class T>
  int GroundNetworkObservationManager<T>::GetNobservation() const
  {
    if (time_index_ >= Nt_ || time_index_ < 0)
      return 0;
    else
      return all_observation_(time_index_).GetLength();
  }


  //! Retrieves the observations at current date.
  /*!
    \param[out] data all observations at current date.
  */
  template <class T>
  void GroundNetworkObservationManager<T>::GetObservation(Vector<T>& data)
    const
  {
    if (time_index_ >= Nt_ || time_index_ < 0)
      data.Clear();
    else
      data = all_observation_(time_index_);
  }


  //! Retrieves the observations and their stations at current date.
  /*!
    \param[out] station_index station of the observations.
    \param[out] data all observations at current date.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::GetObservation(Vector<int>& station_index, Vector<T>& data)
    const
  {
    if (time_index_ >= Nt_ || time_index_ < 0)
      {
        station_index.Clear();
        data.Clear();
      }
    else
      {
        station_index = all_station_index_(time_index_);
        data = all_observation_(time_index_);
      }
  }


  //! Retrieves the observations and their stations at current date.
  /*!
    \param[out] abscissa position along x of the observations.
    \param[out] ordinate position along y of the observations.
    \param[out] data all observations at current date.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::GetObservation(Vector<T>& abscissa, Vector<T>& ordinate, Vector<T>& data)
    const
  {
    if (time_index_ >= Nt_ || time_index_ < 0)
      {
        abscissa.Clear();
        ordinate.Clear();
        data.Clear();
      }
    else
      {
        data = all_observation_(time_index_);
        abscissa.Reallocate(data.GetLength());
        ordinate.Reallocate(data.GetLength());
        for (int i = 0; i < data.GetLength(); i++)
          abscissa(i) = station_x_(all_station_index_(time_index_)(i));
        for (int i = 0; i < data.GetLength(); i++)
          ordinate(i) = station_y_(all_station_index_(time_index_)(i));
      }
  }


  //! Computes the innovation.
  /*!
    \param[in] state state vector.
    \param[out] innovation innovation vector.
  */
  template <class T>
  template <class state>
  void GroundNetworkObservationManager<T>
  ::GetInnovation(const state& x, observation& innovation)
  {
    if (time_index_ >= Nt_ || time_index_ < 0)
      {
        innovation.Clear();
        return;
      }

    observation observation = all_observation_(time_index_);
    int Nobservation = observation.GetLength();

    innovation.Reallocate(Nobservation);
    GetObservation(observation);
    ApplyOperator(x, innovation);
    Mlt(T(-1), innovation);
    Seldon::Add(T(1), observation, innovation);
  }


  //! Returns all the available stations.
  /*!
    \return All available stations, that is, the attribute
    'all_station_index_'.
  */
  template <class T>
  Vector2<int>& GroundNetworkObservationManager<T>::GetAllStationIndex()
  {
    return all_station_index_;
  }


  //! Returns all the available stations.
  /*!
    \return All available stations, that is, the attribute
    'all_station_index_'.
  */
  template <class T>
  const Vector2<int>& GroundNetworkObservationManager<T>
  ::GetAllStationIndex() const
  {
    return all_station_index_;
  }


  //! Returns all the available observations.
  /*!
    \return All available observations, that is, the attribute
    'all_observation_'.
  */
  template <class T>
  Vector2<T>& GroundNetworkObservationManager<T>::GetAllObservation()
  {
    return all_observation_;
  }


  //! Returns all the available observations.
  /*!
    \return All available observations, that is, the attribute
    'all_observation_'.
  */
  template <class T>
  const Vector2<T>& GroundNetworkObservationManager<T>::GetAllObservation()
    const
  {
    return all_observation_;
  }


  //! Applies the observation operator.
  /*!
    \param[in] x state vector.
    \param[out] y output of the observation operator applied to \a x.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  template <class state>
  void GroundNetworkObservationManager<T>
  ::ApplyOperator(const state& x, observation& y) const
  {
    if (!with_domain_)
      throw Error("GroundNetworkObservationManager::ApplyOperator(x, y)",
                  "The domain is missing so that "
                  "the observation operator cannot be defined.");

    if (time_index_ >= Nt_ || time_index_ < 0)
      {
        y.Clear();
        return;
      }

    int Nobservation = all_observation_(time_index_).GetLength();
    y.Reallocate(Nobservation);

    for (int l = 0; l < Nobservation; l++)
      {
        int index = all_station_index_(time_index_, l);
        y(l) = interpolation_weight_(index, 0)
          * x(interpolation_index_(index, 0))
          + interpolation_weight_(index, 1)
          * x(interpolation_index_(index, 1))
          + interpolation_weight_(index, 2)
          * x(interpolation_index_(index, 2))
          + interpolation_weight_(index, 3)
          * x(interpolation_index_(index, 3));
      }
  }


  //! Returns the observation operator.
  /*!
    \return The observation operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  Matrix<T, General, RowSparse> GroundNetworkObservationManager<T>
  ::GetOperator() const
  {
    if (!with_domain_)
      throw Error("GroundNetworkObservationManager::GetOperator()",
                  "The domain is missing so that "
                  "the observation operator cannot be defined.");

    if (time_index_ >= Nt_ || time_index_ < 0)
      return Matrix<T, General, RowSparse>();

    Matrix<T, General, ArrayRowSparse>
      H(all_station_index_(time_index_).GetLength(), model_Nstate_);

    for (int l = 0; l < all_station_index_(time_index_).GetLength(); l++)
      {
        int index = all_station_index_(time_index_, l);
        H.Set(l, interpolation_index_(index, 0),
              interpolation_weight_(index, 0));
        H.Set(l, interpolation_index_(index, 1),
              interpolation_weight_(index, 1));
        H.Set(l, interpolation_index_(index, 2),
              interpolation_weight_(index, 2));
        H.Set(l, interpolation_index_(index, 3),
              interpolation_weight_(index, 3));
      }

    Matrix<T, General, RowSparse> output;
    Copy(H, output);

    return output;
  }


  //! Returns one element of the tangent linear operator.
  /*!
    \param[in] i row index.
    \param[in] j column index.
    \return The element (\a i, \a j) of the linearized operator.
  */
  template <class T>
  T GroundNetworkObservationManager<T>
  ::GetTangentLinearOperator(int i, int j) const
  {
    if (!with_domain_)
      throw Error("GroundNetworkObservationManager::"
                  "GetTangentLinearOperator(i, j)",
                  "The domain is missing so that "
                  "the observation operator cannot be defined.");

    if (time_index_ >= Nt_ || time_index_ < 0)
      throw Error("GroundNetworkObservationManager::"
                  "GetTangentLinearOperator(i, j)",
                  "The observation operator is empty: element ("
                  + to_str(i) + ", " + to_str(j) + ") does not exist.");

    int index = all_station_index_(time_index_, i);

    if (j == interpolation_index_(index, 0))
      return interpolation_weight_(index, 0);
    else if (j == interpolation_index_(index, 1))
      return interpolation_weight_(index, 1);
    else if (j == interpolation_index_(index, 2))
      return interpolation_weight_(index, 2);
    else if (j == interpolation_index_(index, 3))
      return interpolation_weight_(index, 3);
    else
      return T(0);
  }


  //! Returns the tangent linear operator.
  /*!
    \param[in] row row index.
    \param[out] tangent_operator_row the row \a row of the linearized
    operator.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::GetTangentLinearOperatorRow(int row,
                                tangent_linear_operator_row&
                                tangent_operator_row)
    const
  {
    int index = all_station_index_(time_index_, row);
    tangent_operator_row.Reallocate(model_Nstate_);
    tangent_operator_row.Zero();
    tangent_operator_row(interpolation_index_(index, 0))
      = interpolation_weight_(index, 0);
    tangent_operator_row(interpolation_index_(index, 1))
      = interpolation_weight_(index, 1);
    tangent_operator_row(interpolation_index_(index, 2))
      = interpolation_weight_(index, 2);
    tangent_operator_row(interpolation_index_(index, 3))
      = interpolation_weight_(index, 3);
  }


  //! Returns the tangent linear operator.
  /*!
    \return The tangent linear operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  typename GroundNetworkObservationManager<T>::tangent_linear_operator
  GroundNetworkObservationManager<T>::GetTangentLinearOperator() const
  {
    Matrix<T, General, RowSparse> sparse = GetOperator();
    Matrix<T> output(sparse.GetM(), sparse.GetN());
    for (int i = 0; i < sparse.GetM(); i++)
      for (int j = 0; j < sparse.GetN(); j++)
        output(i, j) = sparse(i, j);
    return output;
  }


  //! Observation error covariance.
  /*!
    \param[in] i row index.
    \param[in] j column index.
    \return The element (\a i, \a j) of the observation error covariance.
  */
  template <class T>
  T GroundNetworkObservationManager<T>
  ::GetErrorVariance(int i, int j) const
  {
    if (i == j)
      return error_variance_value_;
    else
      return T(0);
  }


  //! Returns the error variance.
  /*! \return The error variance. */
  template <class T>
  typename GroundNetworkObservationManager<T>::error_variance
  GroundNetworkObservationManager<T>::GetErrorVariance() const
  {
    int Nobservation = GetNobservation();
    Matrix<T> R(Nobservation, Nobservation);
    R.SetIdentity();
    Mlt(error_variance_value_, R);
    return R;
  }


  ////////////////
  // PROCESSING //
  ////////////////


  //! Filters a time sequence.
  /*! It removes the values considered as "unknown values", and it may extract
    daily peaks or a given hour from the time sequence.
    \param[in] time_index time indices associated with the data.
    \param[in,out] data on entry, the sequence of data to be filtered; on
    exit, the filtered sequence.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::FilterTimeSequence(Vector<int>& time_index, Vector<T>& data)
  {
    // First removes useless data.
    int t = 0, i;
    bool remove;
    for (int h = 0; h < Nt_; h++)
      {
        remove = false;
        for (i = 0; i < unknown_value_.GetLength(); i++)
          remove = remove || unknown_value_(i) == data(h);
        if (!remove)
          {
            data(t) = data(h);
            time_index(t) = time_index(h);
            t++;
          }
      }
    time_index.Reallocate(t);
    data.Reallocate(t);

    // Second keeps only the peaks if required.
    if (mode_ == "daily_peak")
      KeepDailyPeak(time_index, data);
    else if (mode_ == "hourly")
      KeepHour(hourly_selected_hour_, time_index, data);
    else if (mode_ == "daily")
      ComputeDailyAverage(time_index, data);
  }


  //! Filters data to keep the daily peaks only.
  /*!
    \param[in] time_index time indices associated with the data.
    \param[in,out] data on entry, the sequence of data to be filtered; on
    exit, the sequence of daily peaks.
    \warning It is assumed that this function takes hourly data on input.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::KeepDailyPeak(Vector<int>& time_index, Vector<T>& data)
  {
    // Old and new indices in 'time_index'.
    int t(0), new_t(0);
    int hour, day, next_day_index;
    T data_max;
    // Number of data available in [daily_peak_hour_begin_,
    // daily_peak_hour_end_[.
    int number_data;
    while (t < time_index.GetLength())
      {
        day = (time_index(t) + date_begin_.GetHour()
               - first_hour_in_day_ + 24) / 24 - 1;
        next_day_index = (day + 1) * 24 - date_begin_.GetHour()
          + first_hour_in_day_;
        data_max = -999.;
        number_data = 0;
        while (t < time_index.GetLength() && time_index(t) < next_day_index)
          {
            hour = (time_index(t) + date_begin_.GetHour()
                    - first_hour_in_day_ + 24) % 24;
            if (hour >= daily_peak_hour_begin_ && hour < daily_peak_hour_end_)
              {
                data_max = max(data_max, data(t));
                number_data++;
              }
            t++;
          }
        if (number_data != 0 && number_data >= Ndaily_peak_hour_)
          {
            data(new_t) = data_max;
            // The peak is moved to midnight, in order to ease latter
            // cleaning.
            time_index(new_t) = max(0, next_day_index - 24);
            new_t++;
          }
      }
    time_index.Reallocate(new_t);
    data.Reallocate(new_t);
  }


  //! Filters data to keep the data at a given hour only.
  /*! Only the data at hour \a hour are kept. The rest is removed.
    \param[in] hour the target hour.
    \param[in] time_index time indices associated with the data.
    \param[in,out] data on entry, the sequence of data to be filtered; on
    exit, the sequence of data at hour \a hour.
    \warning It is assumed that the input data is hourly data.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::KeepHour(int hour, Vector<int>& time_index, Vector<T>& data)
  {
    if (hour < 0 || hour > 23)
      throw Error("GroundNetworkObservationManager::KeepHour",
                  "The hour should be in [0, 23], but " + to_str(hour)
                  + " was provided.");
    // New index in 'time_index'.
    int new_t(0);
    int data_hour;
    if (!is_multiple(Delta_t_, 3600.))
      throw Error("GroundNetworkObservationManager::KeepHour",
                  "The input time step size is not multiple of 3600.");
    int delta = int(Delta_t_) / 3600;
    for (int t = 0; t < time_index.GetLength(); t++)
      {
        data_hour = (time_index(t) * delta
                     + date_begin_.GetHour()) % 24;
        if (data_hour == hour)
          {
            data(new_t) = data(t);
            // The value is moved to midnight, in order to ease latter
            // cleaning.
            time_index(new_t) = max(0, time_index(t) - hour);
            new_t++;
          }
      }
    time_index.Reallocate(new_t);
    data.Reallocate(new_t);
  }


  //! Filters data to keep the container data at a given hour only.
  /*! Only the data at hour \a hour are kept. The rest is removed.
    \param[in] hour the target hour.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::KeepHour(int hour, Vector2<int>& station_index, Vector2<T>& data)
  {
    if (hour < 0 || hour > 23)
      throw Error("GroundNetworkObservationManager::KeepHour",
                  "The hour should be in [0, 23], but " + to_str(hour)
                  + " was provided.");
    Date current_date(date_begin_);
    for (int t = 0; t < Nt_; t++)
      {
        if (current_date.GetHour() != hour)
          {
            data(t).Clear();
            station_index(t).Clear();
          }
        current_date.AddSeconds(double(Delta_t_));
      }
  }


  //! Computes the daily average.
  /*!
    \param[in] time_index time indices associated with the data.
    \param[in,out] data on entry, the sequence of data to be processed; on
    exit, the sequence of daily averages.
    \warning It is assumed that the input data is hourly data.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::ComputeDailyAverage(Vector<int>& time_index, Vector<T>& data)
  {
    // Old and new indices in 'time_index'.
    int t(0), new_t(0);
    int hour, day, next_day_index;
    T data_mean;
    // Number of data available in the day.
    int number_data;
    while (t < time_index.GetLength())
      {
        day = (time_index(t) + date_begin_.GetHour()
               - first_hour_in_day_ + 24) / 24 - 1;
        next_day_index = (day + 1) * 24 - date_begin_.GetHour()
          + first_hour_in_day_;
        data_mean = 0.;
        number_data = 0;
        while (t < time_index.GetLength() && time_index(t) < next_day_index)
          {
            data_mean += data(t);
            number_data++;
            t++;
          }
        if (number_data != 0 && number_data >= Ndaily_average_)
          {
            data(new_t) = data_mean / T(number_data);
            // The daily mean is moved to midnight, in order to ease latter
            // cleaning.
            time_index(new_t) = max(0, next_day_index - 24);
            new_t++;
          }
      }
    time_index.Reallocate(new_t);
    data.Reallocate(new_t);
  }


  //! Switches from a per-hour indexation to a daily indexation.
  /*! Checks that, at any station, there is no more than one value per day (at
    \a selected_hour if applicable) and changes the time indices to switch to
    an indexation per day. On exit, the data is available at midnight each
    day.
    \param[in] selected_hour the hour at which the daily data may be
    available. If it is set to -1, the data may be at any hour in the day.
    \warning In general, no data should be loaded after this method was
    called.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::FromHourlyToDaily(Vector2<int>& station_index_, Vector2<T>& data_,
                      int selected_hour)
  {
    if (Delta_t_ != 3600.)
      throw Error("GroundNetworkObservationManager::"
                  "FromHourlyToDaily",
                  "The current data is not hourly data.");

    Vector2<int> station_index_copy;
    station_index_copy.Copy(station_index_);
    Vector2<T> data_copy;
    data_copy.Copy(data_);

    Date last_date = date_end_;
    last_date.AddSeconds(-Delta_t_);
    int Nday = last_date.GetDaysFrom(date_begin_) + 1;

    data_.Clear();
    data_.Reallocate(Nday);
    station_index_.Clear();
    station_index_.Reallocate(Nday);

    int hour, previous_day_index(-1), day_index;
    // Was data already selected for the current day?
    bool selected_data = false;
    for (int t = 0; t < Nt_; t++)
      {
        hour = (t + date_begin_.GetHour()) % 24;
        day_index = (t + date_begin_.GetHour()) / 24;
        selected_data = selected_data && previous_day_index == day_index;
        if ((selected_hour == -1 && station_index_copy(t).GetLength() != 0)
            || hour == selected_hour)
          {
            if (selected_data)
              throw Error("GroundNetworkObservationManager::"
                          "FromHourlyToDaily",
                          "Found data twice in day "
                          + to_str(day_index) + ".");
            data_(day_index) = data_copy(t);
            station_index_(day_index) = station_index_copy(t);
            selected_data = true;
          }
        else if (station_index_copy(t).GetLength() != 0)
          throw Error("GroundNetworkObservationManager::"
                      "FromHourlyToDaily",
                      "Found data at another hour than "
                      + to_str(selected_hour) + ".");
        previous_day_index = day_index;
      }

    date_begin_.AddHours(-first_hour_in_day_);
    date_begin_.SetHour(0);
    last_date.AddHours(-first_hour_in_day_);
    last_date.SetHour(0);
    date_end_ = last_date;
    date_end_.AddDays(1);
    Delta_t_ = 86400.;
    Nt_ = Nday;
  }


  //! Switches from a daily indexation to an hourly indexation.
  /*! On exit, the container indexes the data per hour. The daily values have
    been moved at \a hour. Other hours have no data. The last date stored is
    at 23:00 hours in the last day.
    \param[in] hour the hour at which all daily values are to be placed.
    \warning In general, no data should be loaded after this method was
    called.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::FromDailyToHourly(int hour,
                      Vector2<int>& station_index_, Vector2<T>& data_)
  {
    if (Delta_t_ != 86400.)
      throw Error("GroundNetworkObservationManager::"
                  "FromDailyToHourly",
                  "The current data is not daily data.");

    Delta_t_ = 3600.;
    Nt_ *= 24;
    // 'date_begin_' is supposed to be at midnight in the first day with data,
    // and the 'date_end_' is supposed to be at midnight in the first day
    // without data. So, no changes to these dates should be necessary.

    Vector2<int> station_index_copy;
    station_index_copy.Copy(station_index_);
    Vector2<T> data_copy;
    data_copy.Copy(data_);

    data_.Clear();
    data_.Reallocate(Nt_);
    station_index_.Clear();
    station_index_.Reallocate(Nt_);

    for (int t = hour, i = 0; t < Nt_; t += 24, i++)
      {
        data_(t) = data_copy(i);
        station_index_(t) = station_index_copy(i);
      }
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the total number of stations.
  /*!
    \return The total number of stations.
  */
  template <class T>
  inline int GroundNetworkObservationManager<T>::GetNstation() const
  {
    return Nstation_;
  }


  //! Returns the abscissa of the domain lower-left corner.
  /*!
    \return The abscissa of the domain lower-left corner.
  */
  template <class T>
  inline T GroundNetworkObservationManager<T>::GetXMin() const
  {
    return x_min_;
  }


  //! Returns the ordinate of the domain lower-left corner.
  /*!
    \return The ordinate of the domain lower-left corner.
  */
  template <class T>
  inline T GroundNetworkObservationManager<T>::GetYMin() const
  {
    return y_min_;
  }


  //! Returns the abscissa of the domain upper-right corner.
  /*!
    \return The abscissa of the domain upper-right corner.
  */
  template <class T>
  inline T GroundNetworkObservationManager<T>::GetXMax() const
  {
    return x_max_;
  }


  //! Returns the ordinate of the domain upper-right corner.
  /*!
    \return The ordinate of the domain upper-right corner.
  */
  template <class T>
  inline T GroundNetworkObservationManager<T>::GetYMax() const
  {
    return y_max_;
  }


  //! Sets the dimension of the model's state vector.
  /*!
    \param[in] Nstate the size of the model's state vector.
  */
  template <class T>
  inline void GroundNetworkObservationManager<T>::SetNstate(int Nstate)
  {
    model_Nstate_ = Nstate;
  }


  //! Returns the initial date.
  /*! There can be no observations before that date.
    \return The first considered date.
  */
  template <class T>
  inline Date GroundNetworkObservationManager<T>::GetDateBegin() const
  {
    return date_begin_;
  }


  //! Returns the last date.
  /*! There can be no observations after that date.
    \return The last considered date.
  */
  template <class T>
  inline Date GroundNetworkObservationManager<T>::GetDateEnd() const
  {
    return date_end_;
  }


  //! Returns the time step.
  /*!
    \return The time step.
  */
  template <class T>
  inline T GroundNetworkObservationManager<T>::GetTimeStep() const
  {
    return Delta_t_;
  }


  //! Returns the list of station complete names.
  /*!
    \return The station complete names.
  */
  template <class T>
  inline const vector<string>&
  GroundNetworkObservationManager<T>::GetStationName() const
  {
    return station_name_;
  }


  //! Returns a station complete name.
  /*!
    \param[in] i the station index.
    \return The complete name of station #\a i.
  */
  template <class T>
  inline string GroundNetworkObservationManager<T>
  ::GetStationName(int i) const
  {
    return station_name_.at(i);
  }


  //! Returns the list of station short names.
  /*!
    \return The station short names.
  */
  template <class T>
  inline const vector<string>&
  GroundNetworkObservationManager<T>::GetStationShortName() const
  {
    return station_short_name_;
  }


  //! Returns a station short name.
  /*!
    \param[in] i the station index.
    \return The short name of station #\a i.
  */
  template <class T>
  inline string GroundNetworkObservationManager<T>
  ::GetStationShortName(int i)
    const
  {
    return station_short_name_.at(i);
  }


  //! Returns the index of a station based on its name.
  /*! An exception is thrown if the station is not found.
    \param[in] name short or complete name of the station, depending on \a
    is_short_name.
    \param[in] is_short_name true if \a name is the short name, false if \a
    name is the complete name.
    \return The index of the station named \a name.
  */
  template <class T>
  inline int GroundNetworkObservationManager<T>
  ::GetStationIndex(string name, bool is_short_name) const
  {
    std::vector<string>::const_iterator position;
    if (is_short_name)
      {
        position = find(station_short_name_.begin(),
                        station_short_name_.end(), name);
        if (position == station_short_name_.end())
          throw Error("GroundNetworkObservationManager::GetStationIndex",
                      "Unable to find a station with short name \""
                      + name + "\".");
        return int(position - station_short_name_.begin());
      }
    else
      {
        position = find(station_name_.begin(),
                        station_name_.end(), name);
        if (position == station_name_.end())
          throw Error("GroundNetworkObservationManager::GetStationIndex",
                      "Unable to find a station with complete name \""
                      + name + "\".");
        return int(position - station_name_.begin());
      }
  }


  //! Returns the coordinates of a station.
  /*!
    \param[in] station the station index.
    \param[out] coordinate the coordinates of the station.
  */
  template <class T>
  void GroundNetworkObservationManager<T>
  ::GetCoordinate(int station, Vector<T>& coordinate) const
  {
    coordinate.Reallocate(2);
    coordinate(0) = station_x_(station);
    coordinate(1) = station_y_(station);
  }


  //! Returns the abscissa of stations.
  /*!
    \return The abscissa of stations.
  */
  template <class T>
  Vector<T>& GroundNetworkObservationManager<T>::GetStationX()
  {
    return station_x_;
  }


  //! Returns the ordinates of stations.
  /*!
    \return The ordinates of stations.
  */
  template <class T>
  Vector<T>& GroundNetworkObservationManager<T>::GetStationY()
  {
    return station_y_;
  }


}


#define POLYPHEMUS_FILE_OBSERVATION_GROUNDNETWORKOBSERVATIONMANAGER_CXX
#endif
