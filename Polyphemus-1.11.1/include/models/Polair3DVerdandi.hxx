// Copyright (C) 2010, INRIA
// Author(s): Anne Tilloy, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_INCLUDE_MODEL_POLAIR3DVERDANDI_HXX

#include <vector>
#include <string>
#include <map>
#include "AtmoDataHeader.hxx"

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
#include <mpi.h>
#endif

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#include "seldon/SeldonHeader.hxx"
#include "seldon/vector/VectorCollection.hxx"


namespace Polyphemus
{

  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // POLAIR3DVERDANDI //
  //////////////////////


  //! This class provides the Polair3D interface to Verdandi.
  template<class T, class ClassModel, class ClassOutputSaver>
  class Polair3DVerdandi: public Verdandi::VerdandiBase
  {
  public:

    typedef Seldon::Vector<T> state;
    typedef Seldon::Vector<Seldon::Vector<T>, Seldon::Collection>
    uncertain_variable;
    typedef Seldon::Vector<T> parameter;
    typedef Seldon::Matrix<T> tangent_linear_operator;
    typedef Seldon::Matrix<T> state_error_variance;
    typedef Seldon::Vector<T> state_error_variance_row;
    typedef Seldon::Matrix<T> matrix_state_observation;
    typedef Seldon::Matrix<T> error_variance;

  protected:

    /*** Main components ***/

    //! Underlying model.
    ClassModel Model;
    //! Output saver.
    ClassOutputSaver OutputSaver;

    //! Pertubation manager.
    Verdandi::NewranPerturbationManager perturbation_manager;

    /*** Configuration ***/

    //! Display options.
    map<string, bool> option_display;

    //! MPI rank.
    int rank;

    /*** Assimilation-related variables ***/

    //! State vector.
    state state_;

    /*! \brief The tangent linear model TLM is approximated with finite
      differences. TLM(z) is approximated by [M(x + tlm_delta z) - M(x)] /
      tlm_delta_ at point x. */
    T delta_tlm;

    //! Square root of the state error variance.
    state_error_variance P_sqrt;

    //! Square root of the error variance.
    state_error_variance Q_sqrt;

    //! Number of columns in Q_sqrt.
    int Nq_sqrt;

    //! Perturbations vectors for the computation of Q_sqrt.
    vector<vector<uncertain_variable> > perturbation;

    //! Type of parameterized-covariance matrix.
    string error_covariance_model_;
    //! Background error variance.
    T background_error_variance_;
    //! Balgovind horizontal scale for the background error variance.
    T balgovind_horizontal_scale_background_;
    //! Balgovind vertical scale for the background error variance.
    T balgovind_vertical_scale_background_;

    //! Buffer for the row of the background error covariance matrix.
    state_error_variance_row buffer_background_row_;
    //! Row index associated with 'buffer_background_row_'.
    int buffer_background_row_index_;

    //! Species included in the state.
    vector<string> species_state_;
    //! Number of species in the state.
    int Nspecies_state_;
    //! Global index of the species that are in the state.
    Seldon::Vector<int> species_index_state_;

    //! Vertical layers included in the state.
    Seldon::Vector<int> level_state_;
    //! Number of vertical levels in the state.
    int Nlevel_state_;

    /*** Parameters ***/

    /*! \brief Name of input parameters (note: state vector is not considered
      as a parameter). */
    vector<string> parameter_name;
    //! Dimensions of input parameters.
    vector<Seldon::Vector<int> > parameter_dimension;
    //! Value of input parameters.
    vector<parameter> parameter_value;

    /*** Maps ***/

    //! Map from field names to their data vectors.
    map<string, Seldon::Vector<T>*> vector_map;

    //! Map from field names to their data vectors.
    map<string, uncertain_variable*> uncertain_variable_list;
    //! Keeps track of the perturbations on the boundary conditions.
    vector<Seldon::Vector<T>* > bc_perturbation;

    //! Map from field names to their collection vectors.
    map<string, Seldon::Vector<T>*> uncertain_variable_correlation;

    //! Map from perturbed field names to their distribution.
    map<string, string> uncertain_variable_pdf;

    //! Map from perturbed field names to their covariance matrix.
    map<string, Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>*>
    uncertain_variable_variance;

    //! Map from perturbed field names to their distribution parameters.
    map<string, Seldon::Vector<T>* > uncertain_variable_parameter;

    //! Map from perturbed field names to their perturbation options.
    map<string, string> uncertain_variable_option;

    //! Vector of uncertain species for boundary conditions.
    vector<string> bc_species_list;

    //! Matrix saving unperturbed east boundary conditions.
    map<string, Seldon::Matrix<T> > bc_backup_east;
    //! Matrix saving unperturbed west boundary conditions.
    map<string, Seldon::Matrix<T> > bc_backup_west;
    //! Matrix saving unperturbed north boundary conditions.
    map<string, Seldon::Matrix<T> > bc_backup_north;
    //! Matrix saving unperturbed south boundary conditions.
    map<string, Seldon::Matrix<T> > bc_backup_south;

  public:

    /*** Constructors and destructor ***/

    Polair3DVerdandi();
    Polair3DVerdandi(string config_file);
    virtual ~Polair3DVerdandi();

    void Initialize(string config_file);
    void InitializeStep();
    void Forward();
    void ApplyOperator(state& x,
                       bool forward = false, bool preserve_state = true);
    void ApplyTangentLinearOperator(state& z);
    void FinalizeStep();

    bool HasFinished() const;
    void Finalize();

    /*** Access methods ***/

    ClassModel& GetModel();
    ClassOutputSaver& GetOutputSaver();

    double GetTime() const;
    void SetTime(double time);
    void SetTime(string time);
    void SetTime(Date time);
    Date GetDate(double time) const;

    int GetNstate() const;
    state& GetState();
    void StateUpdated();
    void GetState(state& state) const;
    void SetState(const Seldon::Vector<T>& state);

    int GetNfull_state() const;
    void GetFullState(state& state) const;
    void SetFullState(const Seldon::Vector<T>& state);

    map<string, Seldon::Vector<T>*>& GetVectorMap();
    int GetNuncertain();
    uncertain_variable& GetUncertainVariable(int i);
    Seldon::Vector<T>& GetPDFCorrelation(int i);
    string GetPDF(int i);
    Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>&
    GetPDFVariance(int i);
    Seldon::Vector<T>& GetPDFParameter(int i);
    string GetPerturbationOption(int i);

    int GetNparameter() const;
    string GetParameterName(int i) const;
    int GetParameterIndex(string name) const;
    Seldon::Vector<int> GetParameterDimension(int i) const;
    parameter& GetParameterValue(int i);
    void SetParameterValue(int i, parameter& value);

    state_error_variance& GetErrorVarianceSqrt();
    state_error_variance& GetStateErrorVarianceSqrt();

    void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                  state_error_variance_row);

  protected:

    /*** Methods to fill the maps ***/

    void FillVarianceXYMap(string map_name, T std, T spatial_decorrelation);
    void FillVarianceXMap(string map_name, T std, T spatial_decorrelation);
    void FillVarianceYMap(string map_name, T std, T spatial_decorrelation);
    void FillParameterMap(string map_name, double clip);
    void FillPDFMap(string map_name, string PDF);
    void FillUncertainVariableD3Map(string name, string field_name,
                                    Seldon::Vector<T>* correlation);
    void FillUncertainVariableD3Map(string name, string field_name,
                                    string map_name);
    void FillUncertainVariableD4Map(string name, string field_name,
                                    string map_name,
                                    Seldon::Vector<T>* correlation);
    void FillUncertainVariableBCMap(string name, string field_name,
                                    string map_name_0, string map_name_1,
                                    Seldon::Vector<T>* correlation_0,
                                    Seldon::Vector<T>* correlation_1);
    void SaveBC(string name);
    void RestoreBC(string name);
    void PerturbBC(string name);

    template <class T0, class Storage0, class Allocator0>
    void SetDimension(Seldon::Vector<T0, Storage0, Allocator0>& in,
                      Seldon::Vector<T0, Storage0, Allocator0>& out);
    template <class T0, class Allocator0>
    void
    SetDimension(Seldon::Vector<T0, Seldon::Collection, Allocator0>& in,
                 Seldon::Vector<T0, Seldon::Collection, Allocator0>& out);

    template <class T0, class Allocator0>
    void Fill(Seldon::Vector<T0, Seldon::Collection, Allocator0>& in,
              string pdf);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_INCLUDE_MODEL_POLAIR3DVERDANDI_HXX
#endif
