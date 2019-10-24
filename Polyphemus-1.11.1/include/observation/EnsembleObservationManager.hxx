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


#ifndef POLYPHEMUS_FILE_OBSERVATION_ENSEMBLEOBSERVATIONMANAGER_HXX


#include "TalosHeader.hxx"
using namespace Talos;

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#include "VerdandiHeader.hxx"
#include "seldon/vector/Vector2.hxx"

#include "GroundNetworkObservationManager.hxx"


namespace Polyphemus
{

  using namespace Verdandi;
  using Seldon::Vector;
  using Verdandi::is_multiple;


  //! This class implements an ensemble-based observation manager.
  /*!
    \tparam T the type of floating-point numbers.
  */
  template <class T>
  class EnsembleObservationManager: public Verdandi::VerdandiBase
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
    //! Underlying observation manager.
    GroundNetworkObservationManager<T> network_observation_manager_;

    //! Ensemble simulations.
    vector<Vector<T> > simulation_;

    //! Paths to the ensemble simulations.
    vector<string> simulation_path_;
    //! Number of members in the ensemble;
    int Nmember_;

    /*! Mode: "none" or "daily_peak".
      - "none": keeps the data as read.
      - "daily_peak": returns the peak value (supposed to be at
      'peak_selected_hour_').
    */
    string mode_;
    //! In case of hourly data, the hour to select.
    int peak_selected_hour_;
    //! First hour of simulations.
    int simulation_first_hour_;

    //! First available date in the simulations.
    Date date_begin_;
    //! Last available date plus one time step in the simulations.
    Date date_end_;
    //! Time step in seconds of the simulations.
    T Delta_t_;
    //! Number of time steps in the simulations.
    int Nt_;

    //! Dimension along x of the simulations data.
    int Nx_;
    //! Dimension along y of the simulations data.
    int Ny_;
    //! Dimension along z of the simulations data.
    int Nz_;
    //! Sizes of the state.
    int Nstate_;

    //! Variance of the observational error.
    T error_variance_;

    //! Current time index.
    int time_index_;

  public:
    // Constructors.
    EnsembleObservationManager();
    EnsembleObservationManager(string configuration_filename);
    void Initialize(string configuration_filename);
    void Initialize(Verdandi::VerdandiOps& configuration);
    template <class Model>
    void Initialize(Model& model, string configuration_filename);
    template <class Model>
    void Initialize(Model& model, Verdandi::VerdandiOps& configuration);

    // Date selection.
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
    Vector2<T>& GetAllObservation();
#ifndef SWIG
    const Vector2<T>& GetAllObservation() const;
#endif

    template <class state>
    void ApplyOperator(const state& x, observation& y) const;
#ifndef SWIG
    Matrix<T> GetOperator() const;
#endif
    void GetOperator(Matrix<T>& H) const;
    observation GetOperatorCol(int c) const;
    void GetOperatorCol(int c, observation& H_col) const;
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

    GroundNetworkObservationManager<T>& GetGroundNetworkObservationManager();

  protected:
    void LoadEnsemble();
  };


}


#define POLYPHEMUS_FILE_OBSERVATION_ENSEMBLEOBSERVATIONMANAGER_HXX
#endif
