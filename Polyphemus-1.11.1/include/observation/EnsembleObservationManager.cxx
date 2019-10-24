// Copyright (C) 2011-2012, INRIA
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


#ifndef POLYPHEMUS_FILE_OBSERVATION_ENSEMBLEOBSERVATIONMANAGER_CXX


#include "EnsembleObservationManager.hxx"

#include "TalosHeader.hxx"
#include "Verdandi.hxx"
#include "seldon/vector/Vector2.cxx"

#include "GroundNetworkObservationManager.cxx"


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
  EnsembleObservationManager<T>::EnsembleObservationManager()
  {
  }


  //! Main constructor.
  /*!
    \param[in] configuration_filename path to the configuration file.
  */
  template <class T>
  EnsembleObservationManager<T>
  ::EnsembleObservationManager(string configuration_filename)
  {
    Initialize(configuration_filename);
  }


  //! Main constructor.
  /*! It loads the main information about the network.
    \param[in] configuration_filename path to the configuration file.
  */
  template <class T>
  void EnsembleObservationManager<T>
  ::Initialize(string configuration_filename)
  {
    Verdandi::VerdandiOps configuration(configuration_filename);
    Initialize(configuration);
  }


  //! Main constructor.
  /*! It loads the main information about the network.
    \param[in] configuration the configuration available in a VerdandiOps
    object.
  */
  template <class T>
  void EnsembleObservationManager<T>
  ::Initialize(Verdandi::VerdandiOps& configuration)
  {

    /*** Configuration ***/

    configuration.SetPrefix("ensemble_observation.");

    configuration.Set("simulation_path", simulation_path_);
    Nmember_ = int(simulation_path_.size());
    simulation_.resize(Nmember_);

    string tmp_date;
    configuration.Set("date_begin", tmp_date);
    date_begin_ = Date(tmp_date);
    configuration.Set("date_end", tmp_date);
    date_end_ = Date(tmp_date);
    configuration.Set("Delta_t", Delta_t_);
    Nt_ = int(T(date_end_.GetSecondsFrom(date_begin_)) / Delta_t_ + 0.5);

    if (configuration.Is<string>("mode"))
      configuration.Set("mode", "ops_in(v, {'none'})", mode_);
    else // it is assumed the entry is a table, and the first element is a
      // string.
      {
        configuration.Set("mode[1]", "ops_in(v, {'daily_peak'})", mode_);
        if (mode_ == "daily_peak")
          configuration.Set("mode[2]", "v >= 0 and v < 24",
                            peak_selected_hour_);
      }

    if (mode_ == "daily_peak")
      {
        Delta_t_ = 3600. * 24.;
        simulation_first_hour_ = date_begin_.GetHour();
        date_begin_.SetHour(0);
        date_end_.SetHour(0);
        date_end_.AddDays(1);
        Nt_ = int(T(date_end_.GetSecondsFrom(date_begin_))
                  / Delta_t_ + 0.5);
      }

    configuration.Set("Nx", Nx_);
    configuration.Set("Ny", Ny_);
    configuration.Set("Nz", Nz_);
    Nstate_ = Nx_ * Ny_ * Nz_;

    configuration.Set("error_variance", error_variance_);

    configuration
      .SetPrefix("ensemble_observation.ground_network_observation.");
    network_observation_manager_.Initialize(configuration);
    network_observation_manager_.SetNstate(Nstate_);
  }


  //! Constructor.
  /*! It loads the main information about the observation manager.
    \param[in] configuration_filename path to the configuration file.
  */
  template <class T>
  template <class Model>
  void EnsembleObservationManager<T>
  ::Initialize(Model& model, string configuration_filename)
  {
    Verdandi::VerdandiOps configuration(configuration_filename);
    configuration.SetPrefix("ensemble_observation.network_observation");
    Initialize(configuration_filename);
  }


  //! Constructor.
  /*! It loads the main information about the observation manager.
    \param[in] configuration the configuration available in a VerdandiOps
    object.
  */
  template <class T>
  template <class Model>
  void EnsembleObservationManager<T>
  ::Initialize(Model& model, Verdandi::VerdandiOps& configuration)
  {
    configuration.SetPrefix("ensemble_observation.network_observation");
    Initialize(configuration);
  }


  ////////////////////
  // DATE SELECTION //
  ////////////////////


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] model this argument is not used.
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  template <class Model>
  void EnsembleObservationManager<T>
  ::SetDate(Model& model, string date)
  {
    SetDate(date);
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  void EnsembleObservationManager<T>::SetDate(string date)
  {
    Date date_obj(date);
    if ((date_obj > date_begin_
         || date_obj == date_begin_
         && simulation_first_hour_ <= peak_selected_hour_)
        && date_obj < date_end_
        && is_multiple(date_obj.GetSecondsFrom(date_begin_), Delta_t_))
      time_index_ = int(date_obj.GetSecondsFrom(date_begin_) / Delta_t_);
    else
      time_index_ = -1; // No simulation data.

    LoadEnsemble();

    network_observation_manager_.SetDate(date);
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] model this argument is not used.
    \param[in] time the time at which the observations may be requested.
  */
  template <class T>
  template <class Model>
  void EnsembleObservationManager<T>
  ::SetTime(Model& model, double time)
  {
    Date date = date_begin_;
    date.AddSeconds(time);
    SetDate(date.GetDate("%y%m%d%h"));
  }


  //! Sets the date at which the observations should be later requested.
  /*!
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  void EnsembleObservationManager<T>::SetTime(string date)
  {
    SetDate(date);
  }


  ////////////////////////
  // OBSERVATION ACCESS //
  ////////////////////////


  //! Checks whether observations are available.
  /*!
    \return True if observations are available, false otherwise.
  */
  template <class T>
  bool EnsembleObservationManager<T>::HasObservation() const
  {
    return network_observation_manager_.GetNobservation() != 0
      && time_index_ >= 0;
  }


  //! Returns the number of available observations.
  /*!
    \return The number of available observations.
  */
  template <class T>
  int EnsembleObservationManager<T>::GetNobservation() const
  {
    if (time_index_ < 0)
      return 0;
    else
      return network_observation_manager_.GetNobservation();
  }


  //! Retrieves the observations at current date.
  /*!
    \param[out] data all observations at current date.
  */
  template <class T>
  void EnsembleObservationManager<T>::GetObservation(Vector<T>& data)
    const
  {
    network_observation_manager_.GetObservation(data);
  }


  //! Retrieves the observations and their stations at current date.
  /*!
    \param[out] station_index station of the observations.
    \param[out] data all observations at current date.
  */
  template <class T>
  void EnsembleObservationManager<T>
  ::GetObservation(Vector<int>& station_index, Vector<T>& data) const
  {
    network_observation_manager_.GetObservation(station_index, data);
  }


  //! Returns all the available observations.
  /*!
    \return All available observations.
  */
  template <class T>
  Vector2<T>& EnsembleObservationManager<T>::GetAllObservation()
  {
    return network_observation_manager_.GetAllObservation();
  }


  //! Returns all the available observations.
  /*!
    \return All available observations.
  */
  template <class T>
  const Vector2<T>& EnsembleObservationManager<T>::GetAllObservation() const
  {
    return network_observation_manager_.GetAllObservation();
  }


  //! Applies the observation operator.
  /*!
    \param[in] x weight vector.
    \param[out] y output of the observation operator applied to \a x.
  */
  template <class T>
  template <class state>
  void EnsembleObservationManager<T>
  ::ApplyOperator(const state& x, observation& y) const
  {
    observation tmp(GetNobservation());
    y.Reallocate(GetNobservation());
    y.Zero();
    for (int s = 0; s < Nmember_; s++)
      {
        network_observation_manager_.ApplyOperator(simulation_[s], tmp);
        Add(T(x(s)), tmp, y);
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
  Matrix<T> EnsembleObservationManager<T>::GetOperator() const
  {
    Matrix<T> H(GetNobservation(), Nmember_);

    observation H_col(GetNobservation());
    for (int s = 0; s < Nmember_; s++)
      {
        network_observation_manager_.ApplyOperator(simulation_[s], H_col);
        SetCol(H_col, s, H);
      }

    return H;
  }


  //! Returns the observation operator.
  /*!
    \param[out] H the observation operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  void EnsembleObservationManager<T>::GetOperator(Matrix<T>& H) const
  {
    H.Reallocate(GetNobservation(), Nmember_);

    observation H_col(GetNobservation());
    for (int s = 0; s < Nmember_; s++)
      {
        network_observation_manager_.ApplyOperator(simulation_[s], H_col);
        SetCol(H_col, s, H);
      }
  }


  //! Returns a column of the observation operator.
  /*!
    \param[in] c index of the column to be returned.
    \return The column \a c of the observation operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  typename EnsembleObservationManager<T>::observation
  EnsembleObservationManager<T>::GetOperatorCol(int c) const
  {
    if (c < 0 || c >= Nmember_)
      throw Error("EnsembleObservationManager::GetOperatorCol(int)",
                  "Column index is " + to_str(c)
                  + ", but it should be in [0, " + to_str(Nmember_ - 1)
                  + "].");

    observation H_col(GetNobservation());
    network_observation_manager_.ApplyOperator(simulation_[c], H_col);

    return H_col;
  }


  //! Returns a column of the observation operator.
  /*!
    \param[in] c index of the column to be returned.
    \param[out] H_col on exit, the column \a c of the observation operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  void EnsembleObservationManager<T>
  ::GetOperatorCol(int c, observation& H_col) const
  {
    if (c < 0 || c >= Nmember_)
      throw
        Error("EnsembleObservationManager::GetOperatorCol(int, observation)",
              "Column index is " + to_str(c)
              + ", but it should be in [0, " + to_str(Nmember_ - 1)
              + "].");

    H_col.Reallocate(GetNobservation());
    network_observation_manager_.ApplyOperator(simulation_[c], H_col);
  }


  //! Returns one element of the tangent linear operator.
  /*!
    \param[in] i row index.
    \param[in] j column index.
    \return The element (\a i, \a j) of the linearized operator.
  */
  template <class T>
  T EnsembleObservationManager<T>
  ::GetTangentLinearOperator(int i, int j) const
  {
    observation H_col(GetNobservation());
    network_observation_manager_.ApplyOperator(simulation_[j], H_col);
    return H_col(i);
  }


  //! Returns the tangent linear operator.
  /*!
    \param[in] row row index.
    \param[out] tangent_operator_row the row \a row of the linearized
    operator.
  */
  template <class T>
  void EnsembleObservationManager<T>
  ::GetTangentLinearOperatorRow(int row,
                                tangent_linear_operator_row&
                                tangent_operator_row)
    const
  {
    Matrix<T> H = GetOperator();
    GetRow(H, row, tangent_operator_row);
  }


  //! Returns the tangent linear operator.
  /*!
    \return The tangent linear operator.
    \note The observation operator is defined only if a domain was provided in
    the configuration file. If the observation operator is not defined, an
    exception is thrown.
  */
  template <class T>
  typename EnsembleObservationManager<T>::tangent_linear_operator
  EnsembleObservationManager<T>::GetTangentLinearOperator() const
  {
    return GetOperator();
  }


  //! Observation error covariance.
  /*!
    \param[in] i row index.
    \param[in] j column index.
    \return The element (\a i, \a j) of the observation error covariance.
  */
  template <class T>
  T EnsembleObservationManager<T>
  ::GetErrorVariance(int i, int j) const
  {
    if (i == j)
      return error_variance_;
    else
      return T(0);
  }


  //! Returns the error variance.
  /*! \return The error variance. */
  template <class T>
  typename EnsembleObservationManager<T>::error_variance
  EnsembleObservationManager<T>::GetErrorVariance() const
  {
    int Nobservation = GetNobservation();
    Matrix<T> R(Nobservation, Nobservation);
    R.SetIdentity();
    Mlt(error_variance_, R);
    return R;
  }


  //! Returns the underlying ground network observation manager.
  /*! \return The underlying ground network observation manager. */
  template <class T>
  GroundNetworkObservationManager<T>&
  EnsembleObservationManager<T>::GetGroundNetworkObservationManager()
  {
    return network_observation_manager_;
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Loads the ensemble of simulations at current time.
  /*!
    \param[in] date the date at which the observations may be requested.
  */
  template <class T>
  void EnsembleObservationManager<T>::LoadEnsemble()
  {
    if (time_index_ == -1)
      {
        for (int s = 0; s < Nmember_; s++)
          simulation_[s].Clear();
        return;
      }
    int record_size = Nstate_ * 4;
    int position = record_size * time_index_;
    if (mode_ == "daily_peak")
      {
        position *= 24;
        position += (peak_selected_hour_ - simulation_first_hour_)
          * record_size;
      }
    Vector<float> tmp(Nstate_);
    for (int s = 0; s < Nmember_; s++)
      {
        ifstream input(simulation_path_[s].c_str());
        if (!input.good())
          throw ErrorConfiguration("EnsembleObservationManager::LoadEnsemble",
                                   "Unable to open simulation file \""
                                   + simulation_path_[s] + "\".");
        input.seekg(position);
        tmp.Read(input, false);
        input.close();

        simulation_[s].Reallocate(Nstate_);
        for (int i = 0; i < Nstate_; i++)
          simulation_[s](i) = tmp(i);
      }
  }


}


#define POLYPHEMUS_FILE_OBSERVATION_ENSEMBLEOBSERVATIONMANAGER_CXX
#endif
