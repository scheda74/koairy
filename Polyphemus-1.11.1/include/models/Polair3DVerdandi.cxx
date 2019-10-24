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


#ifndef POLYPHEMUS_FILE_INCLUDE_MODELS_POLAIR3DVERDANDI_CXX

#include "Polair3DVerdandi.hxx"

#include "AtmoData.hxx"
#include "seldon/Seldon.hxx"
#include "seldon/vector/VectorCollection.cxx"

namespace Polyphemus
{


  /////////////////////////////////
  // CONSTRUCTORS AND DESTRUCTOR //
  /////////////////////////////////


  //! Default constructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::Polair3DVerdandi(): Model(), buffer_background_row_index_(-1)
  {
  }


  //! Main constructor.
  /*! Builds the underlying model.
    \param config_file configuration file of the underlying model.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::Polair3DVerdandi(string config_file): Model(config_file),
                                          buffer_background_row_index_(-1)
  {
  }


  //! Destructor.
  /*! Deletes the maps used for the communication with Verdandi.
   */
  template<class T, class ClassModel, class ClassOutputSaver>
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::~Polair3DVerdandi()
  {
    for (size_t i = 0; i < parameter_value.size(); i++)
      parameter_value[i].Nullify();

    typename map<string, Seldon::Vector<T>*>::iterator i;
    for (i = vector_map.begin(); i != vector_map.end(); i++)
      {
        i->second->Nullify();
        delete i->second;
      }

    typename map<string, uncertain_variable*>::iterator j;
    for (j = uncertain_variable_list.begin();
         j != uncertain_variable_list.end(); j++)
      {
        j->second->Nullify();
        delete j->second;
      }

    typename map < string,
                   typename Seldon::Matrix < T, Seldon::Symmetric,
                                             Seldon::RowSymPacked > * >
      ::iterator k;
    for (k = uncertain_variable_variance.begin();
         k != uncertain_variable_variance.end(); k++)
      delete k->second;

    typename map<string, typename Seldon::Vector<T>*>::iterator l;
    for (l = uncertain_variable_parameter.begin();
         l != uncertain_variable_parameter.end(); l++)
      delete l->second;

    typename map<string, typename Seldon::Vector<T>*>::iterator m;
    for (m = uncertain_variable_correlation.begin();
         m != uncertain_variable_correlation.end(); m++)
      delete m->second;

    for (int n = 0; n < int(bc_perturbation.size()); n++)
      delete bc_perturbation[n];

    for (int n = 0; n < int(perturbation.size()); n++)
      for (int o = 0; o < int(perturbation[n].size()); o++)
        perturbation[n][o].Deallocate();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Finalize();
#endif
  }


  //! Performs the initialization.
  /*! Initializes the underlying model and fills the uncertain-variables maps
    for Verdandi.
    \param config_file configuration file of the underlying model.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::Initialize(string config_file)
  {
    Model.Construct(config_file);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Init();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    Model.Init();
    if (rank == 0)
      OutputSaver.Init(Model);

    /*** Display options ***/

    // Configuration file.
    ConfigStream configuration(config_file);

    configuration.SetSection("[display]");
    // Should iterations be displayed on screen?
    configuration.PeekValue("Show_iterations",
                            option_display["show_iterations"]);
    // Should current date be displayed on screen?
    configuration.PeekValue("Show_date", option_display["show_date"]);

    /*** Assimilation-related options ***/

    configuration.SetSection("[data_assimilation]");

    // Scaling factor for the tangent linear model.
    configuration.PeekValue("delta_tlm", delta_tlm);

    // Number of columns in the square root of the error variance.
    configuration.PeekValue("Nq_sqrt", Nq_sqrt);

    // Parameterized background error covariance matrix.
    configuration.PeekValue("Error_covariance_model",
                            error_covariance_model_);
    configuration.PeekValue("Background_error_variance",
                            background_error_variance_);
    configuration.PeekValue("Balgovind_horizontal_scale_background",
                            balgovind_horizontal_scale_background_);
    configuration.PeekValue("Balgovind_vertical_scale_background",
                            balgovind_vertical_scale_background_);

    configuration.SetSection("[state]");

    configuration.Find("Species");
    species_state_ = split(configuration.GetLine());
    Nspecies_state_ = int(species_state_.size());
    species_index_state_.Reallocate(Nspecies_state_);
    for (int i = 0; i < Nspecies_state_; i++)
      species_index_state_(i) = Model.GetSpeciesIndex(species_state_[i]);

    configuration.Find("Levels");
    vector<string> level_state = split(configuration.GetLine());
    Nlevel_state_ = int(level_state.size());
    level_state_.Reallocate(Nlevel_state_);
    for (int i = 0; i < Nlevel_state_; i++)
      level_state_(i) = to_num<int>(level_state[i]);

    /*** Input parameters ***/

    parameter_value.resize(Model.GetNparameter());
    for (int i = 0; i < Model.GetNparameter(); i++)
      {
        string name = Model.GetParameterName(i);
        parameter_name.push_back(name);
        Seldon::Vector<int> new_parameter_dimension;
        int dimension = Model.GetFieldDimension(name);
        new_parameter_dimension.Reallocate(dimension);
        if (dimension == 2)
          {
            parameter_value[i].SetData(Model.A2(name).numElements(),
                                       Model.A2(name).data());
            new_parameter_dimension(0) = Model.A2(name).extent(0);
            new_parameter_dimension(1) = Model.A2(name).extent(1);
          }
        if (dimension == 3)
          {
            parameter_value[i].SetData(Model.A3(name).numElements(),
                                       Model.A3(name).data());
            new_parameter_dimension(0) = Model.A3(name).extent(0);
            new_parameter_dimension(1) = Model.A3(name).extent(1);
            new_parameter_dimension(2) = Model.A3(name).extent(2);
          }
        if (dimension == 4)
          {
            parameter_value[i].SetData(Model.A4(name).numElements(),
                                       Model.A4(name).data());
            new_parameter_dimension(0) = Model.A4(name).extent(0);
            new_parameter_dimension(1) = Model.A4(name).extent(1);
            new_parameter_dimension(2) = Model.A4(name).extent(2);
            new_parameter_dimension(3) = Model.A4(name).extent(3);
          }
        if (dimension == 5)
          {
            parameter_value[i].SetData(Model.A5(name).numElements(),
                                       Model.A5(name).data());
            new_parameter_dimension(0) = Model.A5(name).extent(0);
            new_parameter_dimension(1) = Model.A5(name).extent(1);
            new_parameter_dimension(2) = Model.A5(name).extent(2);
            new_parameter_dimension(3) = Model.A4(name).extent(3);
            new_parameter_dimension(4) = Model.A4(name).extent(4);
          }
        parameter_dimension.push_back(new_parameter_dimension);
      }

    /*** Reading the perturbation file ***/

    // Perturbation file.
    string pert_file;
    configuration.SetSection("[perturbation]");
    configuration.PeekValue("Configuration_file", pert_file);
    ConfigStream perturbation_file(pert_file);

    string line, name;
    vector<string> word_vector;

    // List of fields to be perturbed.
    string field_name;
    vector<string> uncertain_field_list;

    // Names of uncertain variables in maps.
    string uncertain_variable_name;
    string uncertain_variable_name_z;
    string uncertain_variable_name_y_south;
    string uncertain_variable_name_y_north;
    string uncertain_variable_name_x_east;
    string uncertain_variable_name_x_west;

    // Probability distribution and parameters.
    string PDF;
    T std, spatial_decorrelation;

    // Pointers on data, used to fill maps.
    Seldon::Vector<T>* tmp;
    Seldon::Vector<T>* correlation = NULL;
    Seldon::Vector<T>* correlation_z = NULL;
    Seldon::Vector<T>* correlation_y_south = NULL;
    Seldon::Vector<T>* correlation_y_north = NULL;
    Seldon::Vector<T>* correlation_x_east = NULL;
    Seldon::Vector<T>* correlation_x_west = NULL;

    /*** Perturbed fields ***/

    // Filling the uncertain-fields list.
    perturbation_file.SetSection("[general]");
    perturbation_file.Find("Fields");
    uncertain_field_list = split(perturbation_file.GetLine());
    if (uncertain_field_list.size() == 1
        && uncertain_field_list[0] == "---")
      uncertain_field_list.clear();
    uncertain_field_list.push_back("AdditionalField");

    // Perturbation loop on uncertain-fields list.
    for (int it = 0; it < int(uncertain_field_list.size()); it++)
      {
        field_name = uncertain_field_list[it];
        perturbation_file.SetSection("[" + field_name + "]");

        while (!perturbation_file.IsEmpty() &&
               perturbation_file.PeekElement()[0] != '[')
          {
            line = perturbation_file.GetFullLine();
            while (line[0] == '#')
              line = perturbation_file.GetFullLine();
            split(line, word_vector);

            /*** New variable ***/

            if (word_vector.size() == 4)
              {
                // Short name of the variable.
                name = word_vector[0];

                // Parameters registering.
                PDF = word_vector[1];
                std = convert<T>(word_vector[2]);
                if (PDF == "LN" || PDF == "LNH")
                  std = log(std);
                spatial_decorrelation = convert<T>(word_vector[3]);

                // Boundary conditions.
                if (field_name == "BoundaryCondition")
                  {
                    if (PDF != "LN" && PDF != "LNH")
                      throw "Only log-normal distributions are supported for "
                        "boundary conditions in Polair3DVerdandi.";

                    this->bc_species_list.push_back(name);

                    // Five correlation vectors.
                    correlation_z = new Seldon::Vector<T>();
                    correlation_x_east = new Seldon::Vector<T>();
                    correlation_x_west = new Seldon::Vector<T>();
                    correlation_y_south = new Seldon::Vector<T>();
                    correlation_y_north = new Seldon::Vector<T>();

                    // Full names of variables. Key names in maps.
                    uncertain_variable_name_z
                      = field_name + "_z_" + name;
                    uncertain_variable_name_x_east
                      = field_name + "_x_east_" + name;
                    uncertain_variable_name_x_west
                      = field_name + "_x_west_" + name;
                    uncertain_variable_name_y_south
                      = field_name + "_y_south_" + name;
                    uncertain_variable_name_y_north
                      = field_name + "_y_north_" + name;

                    // Five Vector collections.
                    this->uncertain_variable_list
                      [uncertain_variable_name_z]
                      = new Seldon
                      ::Vector<Seldon::Vector<T>, Seldon::Collection>;

                    this->uncertain_variable_list
                      [uncertain_variable_name_x_east]
                      = new Seldon
                      ::Vector<Seldon::Vector<T>, Seldon::Collection>;

                    this->uncertain_variable_list
                      [uncertain_variable_name_x_west]
                      = new Seldon
                      ::Vector<Seldon::Vector<T>, Seldon::Collection>;

                    this->uncertain_variable_list
                      [uncertain_variable_name_y_south]
                      = new Seldon
                      ::Vector<Seldon::Vector<T>, Seldon::Collection>;

                    this->uncertain_variable_list
                      [uncertain_variable_name_y_north]
                      = new Seldon
                      ::Vector<Seldon::Vector<T>, Seldon::Collection>;

                    // Filling of maps.
                    FillUncertainVariableD3Map(name,
                                               "BoundaryCondition_z",
                                               uncertain_variable_name_z);
                    FillUncertainVariableBCMap
                      (name,
                       "BoundaryCondition_x",
                       uncertain_variable_name_x_west,
                       uncertain_variable_name_x_east,
                       correlation_x_west,
                       correlation_x_east);
                    FillUncertainVariableBCMap
                      (name, "BoundaryCondition_y",
                       uncertain_variable_name_y_south,
                       uncertain_variable_name_y_north,
                       correlation_y_south,
                       correlation_y_north);

                    FillVarianceXYMap(uncertain_variable_name_z,
                                      std / 2., spatial_decorrelation);
                    FillVarianceYMap(uncertain_variable_name_x_west,
                                     std / 2., spatial_decorrelation);
                    FillVarianceYMap(uncertain_variable_name_x_east,
                                     std / 2., spatial_decorrelation);
                    FillVarianceXMap(uncertain_variable_name_y_south,
                                     std / 2., spatial_decorrelation);
                    FillVarianceXMap(uncertain_variable_name_y_north,
                                     std / 2., spatial_decorrelation);

                    FillParameterMap(uncertain_variable_name_z, 2.);
                    FillParameterMap(uncertain_variable_name_y_south, 2.);
                    FillParameterMap(uncertain_variable_name_y_north, 2.);
                    FillParameterMap(uncertain_variable_name_x_west, 2.);
                    FillParameterMap(uncertain_variable_name_x_east, 2.);

                    // Default option for the moment.
                    this->uncertain_variable_option
                      [uncertain_variable_name_z] = "every_step";
                    this->uncertain_variable_option
                      [uncertain_variable_name_x_west] = "init_step";
                    this->uncertain_variable_option
                      [uncertain_variable_name_x_east] = "init_step";
                    this->uncertain_variable_option
                      [uncertain_variable_name_y_south] = "init_step";
                    this->uncertain_variable_option
                      [uncertain_variable_name_y_north] = "init_step";
                  }
                // Other fields.
                else
                  {
                    // One correlation vector.
                    correlation = new Seldon::Vector<T>();

                    // Full name of the variable.
                    uncertain_variable_name = field_name + "_" + name;

                    // One vector collection.
                    this->uncertain_variable_list
                      [uncertain_variable_name]
                      = new Seldon::Vector < Seldon::Vector<T>,
                                             Seldon::Collection >;

                    // Fields independent of chemical species.
                    if (field_name == "AdditionalField")
                      FillUncertainVariableD3Map(name, field_name,
                                                 correlation);
                    // Fields dependent of chemical species.
                    else
                      {
                        if (field_name == "SurfaceEmission" ||
                            field_name == "DepositionVelocity")
                          FillUncertainVariableD3Map
                            (name, field_name,
                             uncertain_variable_name);
                        else if (field_name == "PhotolysisRate")
                          FillUncertainVariableD4Map
                            (name, field_name,
                             uncertain_variable_name,
                             correlation);
                      }

                    // Variance and parameters maps filling.
                    FillVarianceXYMap(uncertain_variable_name, std / 2.,
                                      spatial_decorrelation);
                    FillParameterMap(uncertain_variable_name, 2.);
                    this->uncertain_variable_option
                      [uncertain_variable_name] = "every_step";
                  }

                line = perturbation_file.PeekFullLine();
                while (line[0] == '#')
                  {
                    line = perturbation_file.GetFullLine();
                    line = perturbation_file.PeekFullLine();
                  }
                split(line, word_vector);

                /*** Correlated variable ***/

                while (word_vector.size() == 2)
                  {
                    line = perturbation_file.GetFullLine();
                    split(line, word_vector);
                    name = word_vector[0];

                    // Boundary conditions.
                    if (field_name == "BoundaryCondition")
                      {
                        this->bc_species_list.push_back(name);

                        FillUncertainVariableD3Map
                          (name,
                           "BoundaryCondition_z",
                           uncertain_variable_name_z);

                        // Filling of correlation vectors.  Only
                        // correlations to 1 are supported for x and y.
                        correlation_z->
                          PushBack(convert<T>(word_vector[1]));
                        correlation_x_west->PushBack
                          (convert<T>(word_vector[1]));
                        correlation_x_east->PushBack
                          (convert<T>(word_vector[1]));
                        correlation_y_south->PushBack
                          (convert<T>(word_vector[1]));
                        correlation_y_north->PushBack
                          (convert<T>(word_vector[1]));

                        // Filling of maps.
                        FillUncertainVariableBCMap
                          (name, "BoundaryCondition_x",
                           uncertain_variable_name_x_west,
                           uncertain_variable_name_x_east,
                           correlation_x_west, correlation_x_east);
                        FillUncertainVariableBCMap
                          (name, "BoundaryCondition_y",
                           uncertain_variable_name_y_south,
                           uncertain_variable_name_y_north,
                           correlation_y_south, correlation_y_north);
                      }
                    // Other fields.
                    else
                      {
                        // Fields independent of chemical species.
                        if (field_name == "AdditionalField")
                          {
                            correlation->
                              PushBack(convert<T>(word_vector[1]));
                            FillUncertainVariableD3Map(name, field_name,
                                                       correlation);
                          }
                        // Fields dependent of chemical species.
                        else
                          {
                            if (field_name == "SurfaceEmission" ||
                                field_name == "DepositionVelocity")
                              {
                                correlation->PushBack
                                  (convert<T>(word_vector[1]));
                                FillUncertainVariableD3Map
                                  (name, field_name,
                                   uncertain_variable_name);
                              }
                            else if (field_name == "PhotolysisRate")
                              {
                                correlation->PushBack
                                  (convert<T>(word_vector[1]));
                                FillUncertainVariableD4Map
                                  (name, field_name,
                                   uncertain_variable_name,
                                   correlation);
                              }
                          }
                      }

                    line = perturbation_file.PeekFullLine();
                    split(line, word_vector);

                  } // End of correlated variable loop.

                // Filling probability distributions maps and correlation
                // maps, for boundary conditions.
                if (field_name == "BoundaryCondition")
                  {
                    FillPDFMap(uncertain_variable_name_z, PDF);
                    FillPDFMap(uncertain_variable_name_x_west, PDF);
                    FillPDFMap(uncertain_variable_name_x_east, PDF);
                    FillPDFMap(uncertain_variable_name_y_south, PDF);
                    FillPDFMap(uncertain_variable_name_y_north, PDF);

                    this->uncertain_variable_correlation
                      [uncertain_variable_name_z] = correlation_z;
                    this->uncertain_variable_correlation
                      [uncertain_variable_name_x_west]
                      = correlation_x_west;
                    this->uncertain_variable_correlation
                      [uncertain_variable_name_x_east]
                      = correlation_x_east;
                    this->uncertain_variable_correlation
                      [uncertain_variable_name_y_south]
                      = correlation_y_south;
                    this->uncertain_variable_correlation
                      [uncertain_variable_name_y_north]
                      = correlation_y_north;
                  }
                // Other fields.
                else
                  {
                    this->uncertain_variable_correlation
                      [uncertain_variable_name] = correlation;
                    FillPDFMap(uncertain_variable_name, PDF);
                  }
              } // End of reference variable loop.
          }
      }

    /*** Initializes the perturbations for Q_sqrt ***/

    string perturbation_manager_configuration_file;
    configuration.SetSection("[perturbation_manager]");
    configuration.PeekValue("Configuration_file",
                            perturbation_manager_configuration_file);

    perturbation_manager.Initialize(perturbation_manager_configuration_file);

    perturbation.resize(Nq_sqrt);

    for (int q = 0; q < Nq_sqrt; q++)
      for (int i = 0; i < GetNuncertain(); i++)
        {
          uncertain_variable output;
          bool allocate;

          if (GetPDF(i) == "Normal"
              || GetPDF(i) == "LogNormal"
              || GetPDF(i) == "BlockNormal"
              || GetPDF(i) == "BlockLogNormal")
            {
              SetDimension(GetUncertainVariable(i), output);
              Fill(output, GetPDF(i));
              perturbation_manager.Sample(GetPDF(i),
                                          GetPDFVariance(i),
                                          GetPDFParameter(i),
                                          GetPDFCorrelation(i),
                                          output);
            }
          else if (GetPDF(i) == "NormalHomogeneous"
                   || GetPDF(i) == "LogNormalHomogeneous"
                   || GetPDF(i) == "BlockNormalHomogeneous"
                   || GetPDF(i) == "BlockLogNormalHomogeneous")
            {
              SetDimension(GetUncertainVariable(i), output);
              Fill(output, GetPDF(i));
              perturbation_manager.Sample(GetPDF(i),
                                          GetPDFVariance(i)(0, 0),
                                          GetPDFParameter(i),
                                          GetPDFCorrelation(i),
                                          output);
            }
          else
            throw
              Verdandi::ErrorConfiguration("Polair3DVerdandi"
                                           "::Initialize(string)",
                                           "The probability distribution \""
                                           + GetPDF(i)
                                           + "\" is not supported.");
          perturbation[q].push_back(output);
        }

    perturbation_manager.Finalize();
  }


  //! Performs the initialization of a time step.
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::InitializeStep()
  {
    if (rank == 0 && option_display["show_iterations"])
      cout << "Performing iteration #" << " " << endl;

    if (rank == 0 && option_display["show_date"])
      cout << "Current date: "
           << Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

    Model.InitStep();
    if (rank == 0)
      OutputSaver.InitStep(Model);
  }


  //! Performs a model time step.
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::Forward()
  {
    // Perturbations of lateral boundary conditions.
    for (int i = 0; i < int((this->bc_species_list).size()); i++)
      PerturbBC(this->bc_species_list[i]);

    Model.Forward();
    if (rank == 0)
      OutputSaver.Save(Model);
  }


  //! Applies the model to a given state vector.
  /*!
    \param[in,out] c on entry, the state vector to which the model is
    applied; on exit, the state vector after the model is applied.
    \param[in] forward Boolean to indicate if the model has to go on to the
    next step.
    \param[in] preserve_state Boolean to indicate if the model state has to
    be preserved.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::ApplyOperator(state& x, bool forward, bool preserve_state)
  {
    state current_state;
    double current_time = GetTime();
    SetTime(current_time);
    InitializeStep();
    if (preserve_state)
      GetFullState(current_state);
    SetState(x);
    Model.Forward();
    if (!forward)
      SetTime(current_time);
    GetState(x);
    if (preserve_state)
      SetFullState(current_state);
  }


  //! Applies the tangent linear model to a given vector.
  /*!
    \param[in,out] x on entry, a vector to which the tangent linear model
    should be applied; on exit, the result.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::ApplyTangentLinearOperator(state& z)
  {
    state current_state;
    GetState(current_state);

    // Computes M(x).
    state M_x = current_state;
    ApplyOperator(M_x, false, true);
    // Computes M(x + delta_tlm * z).
    T delta = Norm1(current_state) / Norm1(z) * delta_tlm;
    Seldon::Add(delta, z, current_state);
    z = current_state;
    for (int j = 0; j < z.GetLength(); j++)
      if (z(j) < T(0))
        z(j) = T(0);
    ApplyOperator(z, false, true);

    // Computes the finite difference.
    Seldon::Add(T(-1), M_x, z);
    Seldon::Mlt(T(1) / delta, z);
  }


  //! Performs operations at the end of one time step.
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::FinalizeStep()
  {
  }


  //! Asks the model whether all time steps have been performed.
  template<class T, class ClassModel, class ClassOutputSaver>
  bool
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::HasFinished() const
  {
    return Model.HasFinished();
  }


  //! Performs operations at the end of the simulation.
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::Finalize()
  {
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the underlying model.
  /*!
    \return The underlying model.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  ClassModel& Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetModel()
  {
    return Model;
  }


  //! Returns the underlying output saver.
  /*!
    \return The underlying output saver.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  ClassOutputSaver& Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetOutputSaver()
  {
    return OutputSaver;
  }


  /*! \brief Returns the number of seconds between the simulation beginning
    and the current date. */
  /*!
    \return The number of seconds between the simulation beginning and the
    current date.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  double Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetTime() const
  {
    return Model.GetCurrentTime();
  }


  //! Moves the model to a given date.
  /*!
    \param[in] time time at which the model should be moved. It should be a
    number of seconds between the simulation beginning and the target date.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetTime(double time)
  {
    Date new_date = Model.GetDate_min();
    new_date.AddSeconds(time);
    Model.SetDate(new_date);
  }


  //! Moves the model to a given date.
  /*!
    \param[in] time time at which the model should be moved.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetTime(string time)
  {
    Date new_date(time);
    Model.SetDate(new_date);
  }


  //! Moves the model to a given date.
  /*!
    \param[in] time time at which the model should be moved.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetTime(Date time)
  {
    Model.SetDate(time);
  }


  //! Converts a time in seconds to a date.
  /*!
    \param[in] time time to be converted. It should be a number of seconds
    between the simulation beginning and the target date.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Date Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetDate(double time)
    const
  {
    Date output = Model.GetDate_min();
    output.AddSeconds(time);
    return output;
  }


  //! Returns the dimension of the state.
  /*!
    \return The dimension of the state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  int Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetNstate() const
  {
    return Nspecies_state_ * Nlevel_state_ * Model.GetNx() * Model.GetNy();
  }


  //! Returns the state vector.
  /*!
    \param[out] state the state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::state&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetState()
  {
    int position = 0;
    state_.Reallocate(GetNstate());
    for (int s = 0; s < Nspecies_state_; s++)
      for (int k = 0; k < Nlevel_state_; k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            state_(position++)
              = Model.GetConcentration()(species_index_state_(s),
                                         level_state_(k), j, i);
    return state_;
  }


  //! Copies an updated state vector to the concentration array.
  template<class T, class ClassModel, class ClassOutputSaver>
  void
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::StateUpdated()
  {
    int position = 0;
    for (int s = 0; s < Nspecies_state_; s++)
      for (int k = 0; k < Nlevel_state_; k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            Model.GetConcentration()(species_index_state_(s),
                                     level_state_(k), j, i)
              = state_(position++);
  }


  //! Provides the state vector.
  /*!
    \param[out] state the state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetState(typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
             ::state& state) const
  {
    int position = 0;
    state.Reallocate(GetNstate());
    for (int s = 0; s < Nspecies_state_; s++)
      for (int k = 0; k < Nlevel_state_; k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            state(position++)
              = Model.GetConcentration()(species_index_state_(s),
                                         level_state_(k), j, i);
  }


  //! Sets the state vector.
  /*!
    \param[in] state the new state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetState(const typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
             ::state& state)
  {
    if (state.GetLength() != GetNstate())
      throw Verdandi::ErrorArgument("Polair3DVerdandi::SetState(state&)",
                                    "The input vector has "
                                    + to_str(state.GetLength())
                                    + " elements, but the state has "
                                    + to_str(GetNstate()) + " elements.");
    int position = 0;
    for (int s = 0; s < Nspecies_state_; s++)
      for (int k = 0; k < Nlevel_state_; k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            Model.GetConcentration()(species_index_state_(s),
                                     level_state_(k), j, i)
              = max(0., state(position++));
  }


  //! Returns the dimension of the full state.
  /*!
    \return The dimension of the full state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  int Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetNfull_state() const
  {
    return Model.GetNs() * Model.GetNx() * Model.GetNy() * Model.GetNz();
  }


  //! Provides the full state vector.
  /*!
    \param[out] state the full state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetFullState(typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
                 ::state& state) const
  {
    int position = 0;
    state.Reallocate(Model.GetNs() * Model.GetNx() * Model.GetNy()
                     * Model.GetNz());
    for (int n = 0; n < Model.GetNs(); n++)
      for (int k = 0; k < Model.GetNz(); k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            state(position++) = Model.GetConcentration()(n, k, j, i);
  }


  //! Sets the full state vector.
  /*!
    \param[in] state the new full state vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetFullState(const typename
                 Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
                 ::state& state)
  {
    if (state.GetLength() != GetNfull_state())
      throw Verdandi::ErrorArgument("Polair3DVerdandi::SetFullState(state&)",
                                    "The input vector has "
                                    + to_str(state.GetLength())
                                    + " elements, but the full state has "
                                    + to_str(GetNfull_state())
                                    + " elements.");
    int position = 0;
    for (int n = 0; n < Model.GetNs(); n++)
      for (int k = 0; k < Model.GetNz(); k++)
        for (int j = 0; j < Model.GetNy(); j++)
          for (int i = 0; i < Model.GetNx(); i++)
            Model.GetConcentration()(n, k, j, i) = state(position++);
  }


  //! Returns the map from field names to their data vectors.
  /*!
    \return The map from field names to their data vectors.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  map<string, Seldon::Vector<T>*>&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetVectorMap()
  {
    return vector_map;
  }


  /*! Returns the number of vectors to be perturbed.
    \return The number of vectors to be perturbed.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  int Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetNuncertain()
  {
    return uncertain_variable_list.size();
  }


  //! Returns the i-th uncertain variable.
  /*!
    \param i variable index in the map of uncertain variables.
    \return The vector of the i-th uncertain variable.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::uncertain_variable&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetUncertainVariable(int i)
  {
    typename map < string, typename Polair3DVerdandi < T, ClassModel,
                                                       ClassOutputSaver >::uncertain_variable* >::iterator
      it = uncertain_variable_list.begin();
    for (int j = 0; j < i; ++j, ++it);
    return *it->second;
  }


  //! Returns the correlation between the vectors of the i-th collection.
  /*! An uncertain variable is associated with a vector collection. This
    method returns the correlation between the vectors of the collection.
    \param i uncertain-variable index.
    \return Correlation vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Seldon::Vector<T>& Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetPDFCorrelation(int i)
  {
    typename map<string, Seldon::Vector<T>* >
      ::iterator it = uncertain_variable_correlation.begin();
    for (int j = 0; j < i; ++j, ++it);
    return *it->second;
  }


  //! Returns the PDF of the i-th uncertain variable.
  /*!
    \param i uncertain-variable index.
    \return The PDF of the i-th uncertain variable.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  string Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetPDF(int i)
  {
    map<string, string>::iterator it = uncertain_variable_pdf.begin();
    for (int j = 0; j < i; ++j, ++it);
    return it->second;
  }


  /*! \brief Returns the covariance matrix associated with the i-th uncertain
    variable. */
  /*!
    \param i uncertain-variable index.
    \return The covariance matrix associated with the i-th variable.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetPDFVariance(int i)
  {
    typename map < string, Seldon::Matrix < T, Seldon::Symmetric,
                                            Seldon::RowSymPacked > * >
      ::iterator it = uncertain_variable_variance.begin();
    for (int j = 0; j < i; ++j, ++it);
    return *it->second;
  }


  //! Returns parameters associated with the PDF of some uncertain variable.
  /*! In case of normal or log-normal distribution, the parameters are
    clipping parameters.
    \param i uncertain-variable index.
    \return The parameters associated with the i-th variable.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Seldon::Vector<T>& Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetPDFParameter(int i)
  {
    typename map<string, Seldon::Vector<T>* >::iterator
      it = uncertain_variable_parameter.begin();
    for (int j = 0; j < i; ++j, ++it);
    return *it->second;
  }


  //! Returns the perturbation option of the i-th uncertain variable.
  /*!
    \param i uncertain-variable index.
    \return The perturbation option of the i-th uncertain variable.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  string Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetPerturbationOption(int i)
  {
    map<string, string>::iterator it = uncertain_variable_option.begin();
    for (int j = 0; j < i; ++j, ++it);
    return it->second;
  }


  //! Returns the total number of parameters.
  /*!
    \return The total number of parameters.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  int Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetNparameter() const
  {
    return int(parameter_name.size());
  }


  //! Returns the name of a parameter.
  /*!
    \param[in] i index of the parameter.
    \return The name of the parameter \a i.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  string Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetParameterName(int i) const
  {
    if (i < 0 || i >= GetNparameter())
      throw Verdandi::ErrorArgument("Polair3DVerdandi::GetParameterName(int)",
                                    "Parameter index should be in [0, "
                                    + to_str(GetNparameter()) + "], but "
                                    + to_str(i) + " was provided.");
    return parameter_name[i];
  }


  //! Returns the index of a parameter.
  /*!
    \param[in] name the name of the parameter.
    \return The index of the parameter named \a name.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  int Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetParameterIndex(string name) const
  {
    for (int i = 0; i < GetNparameter(); i++)
      if (GetParameterName(i) == name)
        return i;
    throw Verdandi::ErrorArgument("Polair3DVerdandi"
                                  "::GetParameterNameIndex(string)",
                                  "Parameter named \"" + name + "\" unknown.");
  }


  //! Returns the dimensions of a parameter.
  /*!
    \param[in] i index of the parameter.
    \return The dimensions of the parameter \a i.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  Seldon::Vector<int> Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetParameterDimension(int i) const
  {
    if (i < 0 || i >= GetNparameter())
      throw Verdandi::ErrorArgument("Polair3DVerdandi::"
                                    "GetParameterDimension(int)",
                                    "Parameter index should be in [0, "
                                    + to_str(GetNparameter()) + "], but "
                                    + to_str(i) + " was provided.");
    return parameter_dimension[i];
  }


  //! Returns the value of a parameter.
  /*!
    \param[in] i index of the parameter.
    \return The value of the parameter \a i.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::parameter&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>::GetParameterValue(int i)
  {
    if (i < 0 || i >= GetNparameter())
      throw Verdandi::ErrorArgument("Polair3DVerdandi"
                                    "::GetParameterValue(int)",
                                    "Parameter index should be in [0, "
                                    + to_str(GetNparameter()) + "], but "
                                    + to_str(i) + " was provided.");
    return parameter_value[i];
  }


  //! Sets the value of a parameter.
  /*!
    \param[in] i index of the parameter.
    \param[in] value the new value of the parameter \a i.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetParameterValue(int i, parameter& value)
  {
    if (i < 0 || i >= GetNparameter())
      throw Verdandi::ErrorArgument("Polair3DVerdandi"
                                    "::SetParameterValue(int, parameter&)",
                                    "Parameter index should be in [0, "
                                    + to_str(GetNparameter()) + "], but "
                                    + to_str(i) + " was provided.");
    parameter_value[i].Copy(value);
  }


  //! Returns the state error variance.
  /*!
    \return The state error variance.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::error_variance&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetErrorVarianceSqrt()
  {
    // Current state.
    state current_state;
    GetState(current_state);
    state current_full_state;
    GetFullState(current_full_state);
    double current_time = GetTime();

    SetTime(current_time);
    InitializeStep();

    // Saves the lateral boundary conditions.
    for (int i = 0; i < int((this->bc_species_list).size()); i++)
      SaveBC(this->bc_species_list[i]);

    // Computes M(x).
    Model.Forward();
    state M_x;
    GetState(M_x);

    // Back at the current time.
    SetFullState(current_full_state);
    SetTime(current_time);

    Q_sqrt.Reallocate(GetNstate(), Nq_sqrt);

    // Makes a copy of the perturbed variables.
    map < string, typename Polair3DVerdandi < T, ClassModel,
                                              ClassOutputSaver >::uncertain_variable > uncertain_variable_copy;
    typename map < string, typename Polair3DVerdandi < T, ClassModel,
                                                       ClassOutputSaver >::uncertain_variable* >::iterator it;
    for (it = uncertain_variable_list.begin();
         it != uncertain_variable_list.end(); ++it)
      {
        typename Polair3DVerdandi < T, ClassModel,
                                    ClassOutputSaver >::uncertain_variable collection_duplicate;
        Seldon::Vector<T>* vector_duplicate;
        for (int i = 0; i < it->second->GetNvector(); i++)
          {
            vector_duplicate
              = new Seldon::Vector<T>(it->second->GetVector(i));
            uncertain_variable_copy[it->first].AddVector(*vector_duplicate);
            vector_duplicate->Nullify();
            delete vector_duplicate;
          }
      }

    for (int q = 0; q < Nq_sqrt; q++)
      {
        InitializeStep();

        for (int i = 0; i < GetNuncertain(); i++)
          if (GetPerturbationOption(i) == "every_step")
            if (GetPDF(i) == "Normal"
                || GetPDF(i) == "BlockNormal"
                || GetPDF(i) == "NormalHomogeneous"
                || GetPDF(i) == "BlockNormalHomogeneous")
              Seldon::Add(1., perturbation[q][i], GetUncertainVariable(i));
            else if (GetPDF(i) == "LogNormal"
                     || GetPDF(i) == "BlockLogNormal"
                     || GetPDF(i) == "LogNormalHomogeneous"
                     || GetPDF(i) == "BlockLogNormalHomogeneous")
              for (int k = 0; k < perturbation[q][i].GetM(); k++)
                GetUncertainVariable(i)(k) *= perturbation[q][i](k);

        // Computes M(x) with perturbed parameters.
        Forward();
        state M_x_perturbed;
        GetState(M_x_perturbed);

        Seldon::Add(T(-1), M_x, M_x_perturbed);
        Seldon::SetCol(M_x_perturbed, q, Q_sqrt);

        // Back at the current time.
        SetFullState(current_full_state);
        SetTime(current_time);

        // Restores unperturbed data.
        for (it = uncertain_variable_list.begin();
             it != uncertain_variable_list.end(); ++it)
          for (int i = 0; i < it->second->GetNvector(); i++)
            it->second->GetVector(i)
              = uncertain_variable_copy[it->first].GetVector(i);
        // Restores the lateral boundary conditions.
        for (int i = 0; i < int((this->bc_species_list).size()); i++)
          RestoreBC(this->bc_species_list[i]);
      }

    // Deallocations.
    for (it = uncertain_variable_list.begin();
         it != uncertain_variable_list.end(); ++it)
      for (int i = 0; i < it->second->GetNvector(); i++)
        uncertain_variable_copy[it->first].Deallocate();

    return Q_sqrt;
  }


  //! Returns the square root of the state error variance.
  /*!
    \return The square root of the state error variance.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  typename Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::state_error_variance&
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetStateErrorVarianceSqrt()
  {
    P_sqrt.Reallocate(GetNstate(), 1);
    P_sqrt.Fill(T(0));
    return P_sqrt;
  }


  //! Returns a row of the state error variance.
  /*!
    \param[in] row row index.
    \param[out] P_row the row with index \a row in the state error variance.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void
  Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::GetStateErrorVarianceRow(int row, state_error_variance_row& P_row)
  {
    if (row == -1)
      {
        P_row = buffer_background_row_;
        return;
      }

    P_row.Reallocate(GetNstate());

    int s_ref = row / (Model.GetNx() * Model.GetNy() * Nlevel_state_);
    row -= s_ref * Model.GetNx() * Model.GetNy() * Nlevel_state_;
    int k_ref = row / (Model.GetNx() * Model.GetNy());
    row -= k_ref * Model.GetNx() * Model.GetNy();
    int j_ref = row / Model.GetNx();
    int i_ref = row - j_ref * Model.GetNx();
    T altitude_ref = Model.GetGridZArray1D()(k_ref);
    T abscissa = Model.GetX_min() + T(i_ref) * Model.GetDelta_x();
    T ordinate = Model.GetY_min() + T(j_ref) * Model.GetDelta_y();

    T abscissa_tmp, ordinate_tmp;
    T distance_h, distance_v;
    // We assume no correlation between the species, so everything is zero but
    // for species 's_ref'.
    P_row.Zero();
    int index = s_ref * Model.GetNx() * Model.GetNy() * Nlevel_state_;
    for (int k = 0; k < Nlevel_state_; k++)
      for (int j = 0; j < Model.GetNy(); j++)
        for (int i = 0; i < Model.GetNx(); i++)
          {
            abscissa_tmp = Model.GetX_min() + T(i) * Model.GetDelta_x();
            ordinate_tmp = Model.GetY_min() + T(j) * Model.GetDelta_y();

            distance_h = sqrt((ordinate_tmp - ordinate)
                              * (ordinate_tmp - ordinate) +
                              (abscissa_tmp - abscissa)
                              * (abscissa_tmp - abscissa));

            P_row(index) = background_error_variance_
              * (1. + distance_h / balgovind_horizontal_scale_background_)
              * exp(-distance_h / balgovind_horizontal_scale_background_);

            if (k != k_ref)
              {
                distance_v = abs(Model.GetGridZArray1D()(k) - altitude_ref);
                P_row(index)
                  *= (1. + distance_v / balgovind_vertical_scale_background_)
                  * exp(-distance_v / balgovind_vertical_scale_background_);
              }

            index++;
          }
  }


  //////////////////////////////
  // METHODS TO FILL THE MAPS //
  //////////////////////////////


  //! Builds the covariance matrix of a horizontal field (Ny by Nx).
  /* The covariance matrix is in Balgovind form.
     \param map_name key name of the perturbed variable in the map.
     \param std standard deviation.
     \param spatial_decorrelation spatial decorrelation scale in degree.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillVarianceXYMap(string map_name, T std, T spatial_decorrelation)
  {
    Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>* variance
      = new Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>
      (Model.GetNx() * Model.GetNy(), Model.GetNx() * Model.GetNy());

    Seldon::Vector<T, Seldon::VectFull>
      abscissa(Model.GetNx() * Model.GetNy());
    Seldon::Vector<T, Seldon::VectFull>
      ordinate(Model.GetNx() * Model.GetNy());

    for (int j = 0; j < Model.GetNy(); j++)
      for (int i = 0; i < Model.GetNx(); i++)
        {
          abscissa(j * Model.GetNx() + i) = Model.GetX_min()
            + T(i) * Model.GetDelta_x();
          ordinate(j * Model.GetNx() + i) = Model.GetY_min()
            + T(j) * Model.GetDelta_y();
        }

    T distance;
    for (int i = 0; i < Model.GetNx() * Model.GetNy(); i++)
      for (int j = 0; j <= i; j++)
        {
          distance =
            sqrt((ordinate(j) - ordinate(i))
                 * (ordinate(j) - ordinate(i)) +
                 (abscissa(j) - abscissa(i))
                 * (abscissa(j) - abscissa(i)));
          (*variance)(i, j) = std * std
            * (1. + distance / spatial_decorrelation)
            * exp(-distance / spatial_decorrelation);
        }

    this->uncertain_variable_variance[map_name] = variance;
  }


  //! Builds the covariance matrix of a 1D field along Nx.
  /* The covariance matrix is in Balgovind form.
     \param map_name key name of the perturbed variable in the map.
     \param std standard deviation.
     \param spatial_decorrelation spatial decorrelation scale in degree.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillVarianceXMap(string map_name, T std, T spatial_decorrelation)
  {
    Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>* variance
      = new Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>
      (Model.GetNx(), Model.GetNx());

    Seldon::Vector<T, Seldon::VectFull> abscissa(Model.GetNx());
    for (int i = 0; i < Model.GetNx(); i++)
      abscissa(i) = Model.GetX_min() + T(i) * Model.GetDelta_x();

    T distance;
    for (int i = 0; i < Model.GetNx(); i++)
      for (int j = 0; j <= i; j++)
        {
          distance = abs(abscissa(j) - abscissa(i));
          (*variance)(i, j) = std * std
            * (1. + distance / spatial_decorrelation)
            * exp(-distance / spatial_decorrelation);
        }

    this->uncertain_variable_variance[map_name] = variance;
  }


  //! Builds the covariance matrix of a 1D field along Ny.
  /* The covariance matrix is in Balgovind form.
     \param map_name key name of the perturbed variable in the map.
     \param std standard deviation.
     \param spatial_decorrelation spatial decorrelation scale in degree.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillVarianceYMap(string map_name, T std, T spatial_decorrelation)
  {
    Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>* variance
      = new Seldon::Matrix<T, Seldon::Symmetric, Seldon::RowSymPacked>
      (Model.GetNy(), Model.GetNy());

    Seldon::Vector<T, Seldon::VectFull> ordinate(Model.GetNy());
    for (int j = 0; j < Model.GetNy(); j++)
      ordinate(j) = Model.GetY_min() + T(j) * Model.GetDelta_y();

    T distance;
    for (int i = 0; i < Model.GetNy(); i++)
      for (int j = 0; j <= i; j++)
        {
          distance = abs(ordinate(j) - ordinate(i));
          (*variance)(i, j) = std * std
            * (1. + distance / spatial_decorrelation)
            * exp(-distance / spatial_decorrelation);
        }

    this->uncertain_variable_variance[map_name] = variance;
  }


  //! Fills the parameter map.
  /*! This method puts the clipping parameters in the parameter map. It means
    that the distribution should be clipped outside [-clip, +clip].
    \param[in] map_name key name of the perturbed field in the parameter
    map.
    \param[in] clip clipping value.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillParameterMap(string map_name, double clip)
  {
    Seldon::Vector<double>* parameter = new Seldon::Vector<double>(2);
    (*parameter)(0) = -clip;
    (*parameter)(1) = clip;
    this->uncertain_variable_parameter[map_name] = parameter;
  }


  //! Fills the probability distribution function map.
  /*!
    \param[in] map_name key name of the perturbed field in the parameter
    map.
    \param[in] PDF name of the probability density function, that is, "N",
    "NH", "LN" or "LNH".
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillPDFMap(string map_name, string PDF)
  {
    if (PDF == "N")
      this->uncertain_variable_pdf[map_name] = "Normal";
    else if (PDF == "LN")
      this->uncertain_variable_pdf[map_name] = "LogNormal";
    else if (PDF == "NH")
      this->uncertain_variable_pdf[map_name] = "NormalHomogeneous";
    else if (PDF == "LNH")
      this->uncertain_variable_pdf[map_name] = "LogNormalHomogeneous";
  }


  /*! \brief Fills the uncertain variable map, from a 3D field independent
    of species. */
  /*! The 3D field is of shape (Nz, Ny, Nx). All vertical levels are
    correlated to the first level. The correlation value is set to one.
    \param[in] name variable name.
    \param[in] field_name field name.
    \param[out] correlation vector of correlation values between vectors of
    the collection.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillUncertainVariableD3Map(string name, string field_name,
                               Seldon::Vector<T>* correlation)
  {
    Seldon::Vector<T>* tmp;
    string map_name = field_name + "_" + name;

    // All vertical levels are put in the same collection.
    for (int k = 0; k <  Model.D3(name).GetLength(0); k++)
      {
        tmp = new Seldon::Vector<T>
          (Model.D3(name).GetNbElements() / Model.D3(name).GetLength(0),
           &Model.D3(name).Value(k, 0, 0));
        this->uncertain_variable_list[map_name]->AddVector(*tmp);
        tmp->Nullify();
        delete tmp;
        if (k > 0)
          correlation->PushBack(1.);
      }
  }


  /*! \brief Fills the uncertain variable map, from a 3D field dependent of
    species. */
  /*! The 3D field is of shape (Ns, Ny, Nx).
    \param[in] name species name.
    \param[in] field_name field name.
    \param[in] map_name is the key name of the perturbed field in the
    map of uncertain variables.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillUncertainVariableD3Map(string name, string field_name,
                               string map_name)
  {
    Seldon::Vector<T>* tmp;
    int species_index;
    species_index = Model.GetSpeciesIndex(field_name, name);
    tmp = new Seldon::Vector<T>
      (Model.D3(field_name).GetNbElements()
       / Model.D3(field_name).GetLength(0),
       &Model.D3(field_name).Value(species_index, 0, 0));
    this->uncertain_variable_list[map_name]->AddVector(*tmp);
    tmp->Nullify();
    delete tmp;
  }


  /*! \brief Fills the uncertain variable map, from a 4D field dependent of
    species. */
  /*! The 4D field is of shape (Ns, Nz, Ny, Nx). All vertical levels are
    correlated to the first level. The correlation value is set to one.
    \param[in] name species name.
    \param[in] field_name field name.
    \param[in] map_name key name of the perturbed field in the map of
    uncertain variables.
    \param[out] correlation vector of correlation values between vectors of
    the collection.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillUncertainVariableD4Map(string name, string field_name,
                               string map_name,
                               Seldon::Vector<T>* correlation)
  {
    Seldon::Vector<T>* tmp;
    int species_index;

    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      {
        species_index = Model
          .GetSpeciesIndex(field_name, name);
        tmp = new Seldon::Vector<T>
          (Model.D4(field_name).GetNbElements()
           / (Model.D4(field_name).GetLength(0)
              * Model.D4(field_name).GetLength(1)),
           &Model.D4(field_name).Value(species_index, k, 0, 0));
        this->uncertain_variable_list[map_name]->AddVector(*tmp);
        tmp->Nullify();
        delete tmp;
        if (k > 0)
          correlation->PushBack(1.);
      }
  }


  /*! \brief Fills the uncertain variable map, from a 4D field dependent of
    species. */
  /*! Function used for lateral boundary conditions only. Two maps are
    filled. The size of the two reference fields is Ny in case of the field
    name "BoundaryCondition_x" and Nx in case of the field name
    "BoundaryCondition_y".  All vertical levels are correlated to the first
    level. The correlation value is set to one.
    \param[in] name species name.
    \param[in] field_name field name.
    \param[in] map_name_0 key name of the perturbed variable in the map of
    uncertain variables for the specific lateral fields: west or south.
    \param[in] map_name_1 key name of the perturbed variable in the map of
    uncertain variables for the specific lateral fields: east or north.
    \param[out] correlation_0 vector of correlation values between vectors of
    the collection for the specific lateral fields: west or south.
    \param[out] correlation_1 vector of correlation values between vectors of
    the collection for the specific lateral fields: east or north.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::FillUncertainVariableBCMap(string name, string field_name,
                               string map_name_0, string map_name_1,
                               Seldon::Vector<T>* correlation_0,
                               Seldon::Vector<T>* correlation_1)
  {
    Seldon::Vector<T>* tmp;
    int species_index;

    if (field_name == "BoundaryCondition_x")
      {
        species_index = Model.GetSpeciesIndex(field_name, name);

        // Levels loop.
        for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
          {
            tmp = new
              Seldon::Vector<T>(Model.D4(field_name).GetLength(2));
            tmp->Fill(T(1.));
            this->uncertain_variable_list[map_name_0]->AddVector(*tmp);
            bc_perturbation.push_back(tmp);
            if (k > 0)
              correlation_0->PushBack(1.);
          }
        // Levels loop.
        for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
          {
            tmp = new
              Seldon::Vector<T>(Model.D4(field_name).GetLength(2));
            tmp->Fill(T(1.));
            this->uncertain_variable_list[map_name_1]->AddVector(*tmp);
            bc_perturbation.push_back(tmp);
            if (k > 0)
              correlation_1->PushBack(1.);
          }
      }

    else if (field_name == "BoundaryCondition_y")
      {
        species_index = Model.GetSpeciesIndex(field_name, name);

        // Levels loop.
        for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
          {
            tmp = new
              Seldon::Vector<T>(Model.D4(field_name).GetLength(3));
            tmp->Fill(T(1.));
            this->uncertain_variable_list[map_name_0]->AddVector(*tmp);
            bc_perturbation.push_back(tmp);
            if (k > 0)
              correlation_0->PushBack(1.);
          }
        // Levels loop.
        for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
          {
            tmp = new
              Seldon::Vector<T>(Model.D4(field_name).GetLength(3));
            tmp->Fill(T(1.));
            this->uncertain_variable_list[map_name_1]->AddVector(*tmp);
            bc_perturbation.push_back(tmp);
            if (k > 0)
              correlation_1->PushBack(1.);
          }
      }
  }


  //! Saves the lateral boundary conditions of some species.
  /*!
    \param[in] name name of the species.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SaveBC(string name)
  {
    bc_backup_east[name]
      .Reallocate(Model.D4("BoundaryCondition_x").GetLength(1),
                  Model.D4("BoundaryCondition_x").GetLength(2));
    bc_backup_west[name]
      .Reallocate(Model.D4("BoundaryCondition_x").GetLength(1),
                  Model.D4("BoundaryCondition_x").GetLength(2));
    bc_backup_south[name]
      .Reallocate(Model.D4("BoundaryCondition_y").GetLength(1),
                  Model.D4("BoundaryCondition_y").GetLength(3));
    bc_backup_north[name]
      .Reallocate(Model.D4("BoundaryCondition_y").GetLength(1),
                  Model.D4("BoundaryCondition_y").GetLength(3));

    int species_index;
    string field_name, map_name_0, map_name_1;

    field_name = "BoundaryCondition_x";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_east_" + name;
    map_name_1 = field_name + "_west_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        bc_backup_east[name](k, j)
          = Model.D4(field_name).Value(species_index, k, j, 0);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        bc_backup_west[name](k, j)
          = Model.D4(field_name).Value(species_index, k, j, 1);

    field_name = "BoundaryCondition_y";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_south_" + name;
    map_name_1 = field_name + "_north_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        bc_backup_south[name](k, i)
          = Model.D4(field_name).Value(species_index, k, 0, i);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        bc_backup_north[name](k, i)
          = Model.D4(field_name).Value(species_index, k, 1, i);
  }


  //! Restores the lateral boundary conditions of some species.
  /*! This method can be used after SaveBC and PerturbBC so as to retrieve the
    unperturbed boundary conditions.
    \param[in] name name of the species.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::RestoreBC(string name)
  {
    int species_index;
    string field_name, map_name_0, map_name_1;

    field_name = "BoundaryCondition_x";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_east_" + name;
    map_name_1 = field_name + "_west_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        Model.D4(field_name).Value(species_index, k, j, 0)
          = bc_backup_east[name](k, j);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        Model.D4(field_name).Value(species_index, k, j, 1)
          = bc_backup_west[name](k, j);

    field_name = "BoundaryCondition_y";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_south_" + name;
    map_name_1 = field_name + "_north_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        Model.D4(field_name).Value(species_index, k, 0, i)
          = bc_backup_south[name](k, i);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        Model.D4(field_name).Value(species_index, k, 1, i)
          = bc_backup_north[name](k, i);
  }


  //! Perturbs the lateral boundary conditions of some species.
  /*! The values of the boundary conditions in the model are not structured in
    memory so that it would be possible for the Verdandi Monte Carlo driver to
    apply correlated perturbations. Consequently, intermediate vectors
    initially filled with ones are perturbed by the Verdandi Monte Carlo
    driver, and the model boundary conditions are multiplied by these
    perturbed vectors.
    \param[in] name name of the species.
    \warning Only the log-normal distribution is supported.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::PerturbBC(string name)
  {
    int species_index;
    string field_name, map_name_0, map_name_1;

    field_name = "BoundaryCondition_x";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_east_" + name;
    map_name_1 = field_name + "_west_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        Model.D4(field_name).Value(species_index, k, j, 0) =
          Model.D4(field_name).Value(species_index, k, j, 0) *
          (*this->uncertain_variable_list[map_name_0])
          .GetVector(k)(j);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int j = 0; j < Model.D4(field_name).GetLength(2); j++)
        Model.D4(field_name).Value(species_index, k, j, 1) =
          Model.D4(field_name).Value(species_index, k, j, 1) *
          (*this->uncertain_variable_list[map_name_1])
          .GetVector(k)(j);

    field_name = "BoundaryCondition_y";
    species_index = Model.GetSpeciesIndex(field_name, name);
    map_name_0 = field_name + "_south_" + name;
    map_name_1 = field_name + "_north_" + name;
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        Model.D4(field_name).Value(species_index, k, 0, i) =
          Model.D4(field_name).Value(species_index, k, 0, i) *
          (*this->uncertain_variable_list[map_name_0])
          .GetVector(k)(i);
    // Levels loop.
    for (int k = 0; k < Model.D4(field_name).GetLength(1); k++)
      for (int i = 0; i < Model.D4(field_name).GetLength(3); i++)
        Model.D4(field_name).Value(species_index, k, 1, i) =
          Model.D4(field_name).Value(species_index, k, 1, i) *
          (*this->uncertain_variable_list[map_name_1])
          .GetVector(k)(i);
  }


  //! Allocates an output vector to the dimension of the input vector.
  /*!
    \param[in] in input vector.
    \param[out] out output vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  template <class T0, class Storage0, class Allocator0>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetDimension(Seldon::Vector<T0, Storage0, Allocator0>& in,
                 Seldon::Vector<T0, Storage0, Allocator0>& out)
  {
    out.Reallocate(in.GetLength());
  }


  /*! \brief Allocates an output vector collection to the dimension of the
    input vector collection. */
  /*!
    \param[in] in input collection vector.
    \param[out] out output collection vector.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  template <class T0, class Allocator0>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::SetDimension(Seldon::Vector<T0, Seldon::Collection, Allocator0>& in,
                 Seldon::Vector<T0, Seldon::Collection, Allocator0>& out)
  {
    T0 suboutput;
    for (int i = 0; i < in.GetNvector(); i++)
      {
        suboutput.Reallocate(in.GetVector(i).GetLength());
        out.AddVector(suboutput);
        suboutput.Nullify();
      }
  }


  /*! Fills an input vector collection according to its probability
    distribution. */
  /*!
    \param[in,out] in input collection vector.
    \param[in] pdf probability density function: Normal, NormalHomogeneous,
    LogNormal or LogNormalHomogeneous.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  template <class T0, class Allocator0>
  void Polair3DVerdandi<T, ClassModel, ClassOutputSaver>
  ::Fill(Seldon::Vector<T0, Seldon::Collection, Allocator0>& in, string pdf)
  {
    if (pdf == "Normal" || pdf == "NormalHomogeneous")
      for (int i = 0; i < in.GetNvector(); i++)
        in.GetVector(i).Fill(typename T0::value_type(0));
    else if (pdf == "LogNormal" || pdf == "LogNormalHomogeneous")
      for (int i = 0; i < in.GetNvector(); i++)
        in.GetVector(i).Fill(typename T0::value_type(1));
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_INCLUDE_MODELS_POLAIR3DVERDANDI_CXX
#endif
