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

// This file is part of the Eulerian model Polair3D.


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRYASSIMCONC_CXX


#include "Polair3DChemistryAssimConc.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::Polair3DChemistryAssimConc(string config_file):
    Polair3DChemistry < T, ClassAdvection,
    ClassDiffusion, ClassChemistry > (config_file),
    current_row(-1), current_column(-1)
  {
  }


  //! Destructor.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::~Polair3DChemistryAssimConc()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >
  ::ReadConfiguration()
  {
    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::ReadConfiguration();

    /*** Species ***/

    this->config.SetSection("[state]");
    this->config.Find("Species");
    state_species_list = split(this->config.GetLine());
    if (state_species_list[0] == "all")
      state_species_list = this->species_list;
    Ns_state = int(state_species_list.size());

    /*** Levels ***/

    this->config.SetSection("[state]");
    this->config.Find("Levels");
    vector<int> state_levels_vector;
    split(this->config.GetLine(), state_levels_vector);
    Nz_state = int(state_levels_vector.size());
    state_levels.resize(Nz_state);
    for (int i = 0; i < Nz_state; i++)
      state_levels(i) = state_levels_vector[i];

    /*** Covariance model ***/

    this->config.SetSection("[data_assimilation]");
    this->config.PeekValue("Error_covariance_model",
                           "Balgovind | diagonal_constant",
                           error_covariance_model);

    int length;

    this->config.SetSection("[data_assimilation]");
    this->config.Find("Background_error_variance");
    vector<T> background_error_variance_vector;
    split(this->config.GetLine(), background_error_variance_vector);
    if (Ns_state != (length = int(background_error_variance_vector.size())))
      throw string("Invalid length (") + to_str(length) +
        string(") of \"Background_error_variance\", should be the number (")
        + to_str(Ns_state) + string(") of species")
        + " in state vector.";
    background_error_variance.resize(Ns_state);
    for (int i = 0; i < Ns_state; i++)
      background_error_variance(i) = background_error_variance_vector[i];

    this->config.SetSection("[data_assimilation]");
    this->config.Find("Model_error_variance");
    vector<T> model_error_variance_vector;
    split(this->config.GetLine(), model_error_variance_vector);
    if (Ns_state != (length = int(model_error_variance_vector.size())))
      throw string("Invalid length (") + to_str(length) +
        string(") of \"Model_error_variance\", should be the number (")
        + to_str(Ns_state) + string(") of species")
        + " in state vector.";
    model_error_variance.resize(Ns_state);
    for (int i = 0; i < Ns_state; i++)
      model_error_variance(i) = model_error_variance_vector[i];

    if (error_covariance_model == "Balgovind")
      {
        this->config.SetSection("[data_assimilation]");
        this->config.Find("Balgovind_scale_background");
        vector<T> Balgovind_scale_background_vector;
        split(this->config.GetLine(), Balgovind_scale_background_vector);
        if (Ns_state !=
            (length = int(Balgovind_scale_background_vector.size())))
          throw string("Invalid length (") + to_str(length) +
            string(") of \"Balgovind_scale_background\", should be")
            + string(" the number (") + to_str(Ns_state)
            + string(") of species")
            + " in state vector.";
        Balgovind_scale_background.resize(Ns_state);
        for (int i = 0; i < Ns_state; i++)
          Balgovind_scale_background(i)
            = Balgovind_scale_background_vector[i];

        this->config.SetSection("[data_assimilation]");
        this->config.Find("Balgovind_vertical_scale_background");
        vector<T> Balgovind_vertical_scale_background_vector;
        split(this->config.GetLine(),
              Balgovind_vertical_scale_background_vector);
        if (Ns_state !=
            (length = int(Balgovind_vertical_scale_background_vector.size())))
          throw string("Invalid length (") + to_str(length) +
            string(") of \"Balgovind_vertical_scale_background\", should be")
            + string(" the number (") + to_str(Ns_state)
            + string(") of species")
            + " in state vector.";
        Balgovind_vertical_scale_background.resize(Ns_state);
        for (int i = 0; i < Ns_state; i++)
          Balgovind_vertical_scale_background(i)
            = Balgovind_vertical_scale_background_vector[i];

        this->config.SetSection("[data_assimilation]");
        this->config.Find("Balgovind_scale_model");
        vector<T> Balgovind_scale_model_vector;
        split(this->config.GetLine(), Balgovind_scale_model_vector);
        if (Ns_state != (length = int(Balgovind_scale_model_vector.size())))
          throw string("Invalid length (") + to_str(length) +
            string(") of \"Balgovind_scale_model\", should be")
            + string(" the number (") + to_str(Ns_state)
            + string(") of species")
            + " in state vector.";
        Balgovind_scale_model.resize(Ns_state);
        for (int i = 0; i < Ns_state; i++)
          Balgovind_scale_model(i) = Balgovind_scale_model_vector[i];

        this->config.SetSection("[data_assimilation]");
        this->config.Find("Balgovind_vertical_scale_model");
        vector<T> Balgovind_vertical_scale_model_vector;
        split(this->config.GetLine(), Balgovind_vertical_scale_model_vector);
        if (Ns_state !=
            (length = int(Balgovind_vertical_scale_model_vector.size())))
          throw string("Invalid length (") + to_str(length) +
            string(") of \"Balgovind_vertical_scale_model\", should be")
            + string(" the number (") + to_str(Ns_state)
            + string(") of species")
            + " in state vector.";
        Balgovind_vertical_scale_model.resize(Ns_state);
        for (int i = 0; i < Ns_state; i++)
          Balgovind_vertical_scale_model(i)
            = Balgovind_vertical_scale_model_vector[i];
      }
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates the array that stores a row of B.
   */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >::Allocate()
  {
    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::Allocate();

    state_species_index.resize(Ns_state);

    Nstate = Ns_state * Nz_state * this->Ny * this->Nx;
    error_covariance_vector.resize(Nstate);
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >::Init()
  {
    Polair3DChemistry<T, ClassAdvection, ClassDiffusion, ClassChemistry>
      ::Init();

    for (int s = 0; s < Ns_state; s++)
      state_species_index(s) = this->GetSpeciesIndex(state_species_list[s]);
  }


  //////////////////////
  // GET OR SET STATE //
  //////////////////////


  //! Returns the state.
  /*!
    \param state_vector (output) state vector.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >
  ::GetState(Array<T, 1>& state_vector)
  {
    state_vector.resize(Nstate);

    int s, k, j, i;
    int index = 0;
    for (s = 0; s < Ns_state; s++)
      for (k = 0; k < Nz_state; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              state_vector(index) =
                this->Concentration(state_species_index(s),
                                    state_levels(k), j, i);
              index++;
            }
  }


  //! Sets the state.
  /*!
    \param state_vector new model state.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >
  ::SetState(const Array<T, 1>& state_vector)
  {
    if ((int) state_vector.size() != Nstate)
      throw Error("Polair3DChemistryAssimConc::SetState",
                  "Index error.");

    int s, k, j, i;
    int index = 0;
    for (s = 0; s < Ns_state; s++)
      for (k = 0; k < Nz_state; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              this->Concentration(state_species_index(s),
                                  state_levels(k), j, i) = state_vector(index);
              index++;
            }
  }


  //! Returns the adjoint state vector.
  /*!
    \param state_vector (output) adjoint state vector.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >
  ::GetState_ccl(Array<T, 1>& state_vector)
  {
    state_vector.resize(Nstate);

    int s, k, j, i;
    int index = 0;
    for (s = 0; s < Ns_state; s++)
      for (k = 0; k < Nz_state; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              state_vector(index) =
                this->Concentration_ccl(state_species_index(s),
                                        state_levels(k), j, i);
              index++;
            }
  }


  //! Sets the adjoint state vector.
  /*!
    \param state_vector new adjoint state vector.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  void Polair3DChemistryAssimConc < T, ClassAdvection,
                                    ClassDiffusion, ClassChemistry >
  ::SetState_ccl(const Array<T, 1>& state_vector)
  {
    if ((int) state_vector.size() != Nstate)
      throw Error("Polair3DChemistryAssimConc::SetState",
                  "Index error.");

    int s, k, j, i;
    int index = 0;
    for (s = 0; s < Ns_state; s++)
      for (k = 0; k < Nz_state; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              this->Concentration_ccl(state_species_index(s),
                                      state_levels(k), j, i)
                = state_vector(index);
              index++;
            }
  }


  ////////////
  // ACCESS //
  ////////////


  //! Returns the model for concentrations error covariances.
  /*!
    \return The model for concentrations error covariances.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  string Polair3DChemistryAssimConc < T, ClassAdvection,
                                      ClassDiffusion, ClassChemistry >
  ::GetErrorCovarianceModel() const
  {
    return error_covariance_model;
  }


  //! Returns the number of elements in the state.
  /*!
    \return The state dimension.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DChemistryAssimConc < T, ClassAdvection,
                                   ClassDiffusion, ClassChemistry >
  ::GetNstate() const
  {
    return Nstate;
  }


  //! Returns the number of levels in the state.
  /*!
    \return The state levels number.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DChemistryAssimConc < T, ClassAdvection,
                                   ClassDiffusion, ClassChemistry >
  ::GetStateNz() const
  {
    return Nz_state;
  }


  //! Returns the vertical levels included in the state.
  /*!
    \return The indices of the vertical levels included in the state.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  const Array<int, 1>&
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::GetStateLevels() const
  {
    return state_levels;
  }


  //! Returns the number of species in the state.
  /*!
    \return The state species number.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DChemistryAssimConc < T, ClassAdvection,
                                   ClassDiffusion, ClassChemistry >
  ::GetStateNs() const
  {
    return Ns_state;
  }


  //! Returns the species name list of the state.
  /*!
    \return The species name list of the state.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  vector<string>
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::GetStateSpeciesList() const
  {
    return state_species_list;
  }


  //! Returns the index of a species that is in the state array.
  /*!
    \param species species name.
    \return The species index in the state array.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DChemistryAssimConc < T, ClassAdvection,
                                   ClassDiffusion, ClassChemistry >
  ::GetStateSpeciesIndex(string species)
  {
    vector<string>::iterator iter;
    iter = find(state_species_list.begin(), state_species_list.end(),
                species);
    if (iter == state_species_list.end())
      throw string("Species \"") + species + "\" not in the state.";

    return int(iter - state_species_list.begin())
      * Nz_state * this->Ny * this->Nx;
  }


  /*! \brief Returns the storage sequence number of a level that is in the
    state array. */
  /*!
    \param level level number.
    \return The storage sequence index of the given level in the state array.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  int Polair3DChemistryAssimConc < T, ClassAdvection,
                                   ClassDiffusion, ClassChemistry >
  ::GetStateLevelSequenceIndex(int level)
  {
    for (int index = 0; index < Nz_state; index++)
      if (level == state_levels(index))
        return index;

    throw "Observation level not included in the state.";
  }


  //////////////////////
  // COVARIANCE MODEL //
  //////////////////////


  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  const Array<T, 1>&
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::GetBackgroundErrVar() const
  {
    return background_error_variance;
  }


  //! Computes a row of the background error covariance matrix B.
  /*! The row is computed and stored. If the next call to this method is
    performed with the same row number, computations are then saved.
    \param row row number.
    \return The value of row number \a row.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  const Array<T, 1>&
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::BackgroundErrorCovariance(int row)
  {
    // The row has already been computed.
    if (row == current_row)
      return error_covariance_vector;
    current_row = row;
    current_column = -1;

    if (error_covariance_model == "Balgovind")
      {
        /*** Only for one species ***/

        for (int i = 0; i < Nstate; i++)
          error_covariance_vector(i) = 0;

        int species_point_number = Nz_state * this->Nx * this->Ny;
        int species_index =  row / species_point_number;
        int species_alignment = species_index * species_point_number;

        int level_point_number = this->Nx * this->Ny;
        int row_alignment = species_alignment +
          (row - species_alignment) / level_point_number * level_point_number;
        // Index along z.
        int k_row = (row - species_alignment) / level_point_number;

        // Indices along x and y in the model domain.
        int j_row = (row - row_alignment) / this->Nx;
        int i_row = (row - row_alignment) - j_row * this->Nx;

        T distance, distance_x, distance_y, distance_z, balgovind_z;
        for (int k = 0; k < Nz_state; k++)
          {
            if (k == k_row)
              balgovind_z = 1.;
            else if (Balgovind_vertical_scale_background(species_index) == 0.)
              balgovind_z = 0.;
            else
              {
                distance_z = abs(this->GridZ3D(k) - this->GridZ3D(k_row));
                balgovind_z
                  = (1. + distance_z
                     / Balgovind_vertical_scale_background(species_index))
                  * exp(-distance_z /
                        Balgovind_vertical_scale_background(species_index));
              }
            for (int j = 0; j < this->Ny; j++)
              for (int i = 0; i < this->Nx; i++)
                {
                  distance_x = T(i - i_row) * this->Delta_x;
                  distance_y = T(j - j_row) * this->Delta_y;
                  distance = sqrt(distance_x * distance_x
                                  + distance_y * distance_y);
                  error_covariance_vector(species_alignment
                                          + k * this->Ny * this->Nx
                                          + j * this->Nx + i)
                    = (1. + distance
                       / Balgovind_scale_background(species_index))
                    * exp(- distance
                          / Balgovind_scale_background(species_index))
                    * balgovind_z
                    * background_error_variance(species_index);
                }
          }
        return error_covariance_vector;
      }
    else if (error_covariance_model == "diagonal_constant")
      {

        /*** Only for one species ***/

        int species_point_number = Nz_state * this->Nx * this->Ny;
        int species_index =  row / species_point_number;

        for (int i = 0; i < Nstate; i++)
          error_covariance_vector(i) = 0.;
        error_covariance_vector(row)
          = background_error_variance(species_index);
        return error_covariance_vector;
      }
    else
      throw string("Error: unknown error covariance model. This")
        + " part of the code should never have been reached.";
  }


  //! Computes a column of the model error covariance matrix Q.
  /*! The column is computed and stored. If the next call to this method is
    performed with the same column number, computations are then saved.
    \param column column number.
    \return The value of column number \a column.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  const Array<T, 1>&
  Polair3DChemistryAssimConc < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  ::ModelErrorCovariance(int column)
  {
    // The column has already been computed.
    if (column == current_column)
      return error_covariance_vector;
    current_column = column;
    current_row = -1;

    if (error_covariance_model == "Balgovind")
      {
        /*** Only for one species ***/

        for (int i = 0; i < Nstate; i++)
          error_covariance_vector(i) = 0;

        int species_point_number = Nz_state * this->Nx * this->Ny;
        int species_index = column / species_point_number;
        int species_alignment = species_index * species_point_number;

        int level_point_number = this->Nx * this->Ny;
        int row_alignment = species_alignment +
          (column - species_alignment) / level_point_number
          * level_point_number;
        // Index along z.
        int k_row = (column - species_alignment) / level_point_number;

        // Indices along x and y in the model domain.
        int j_row = (column - row_alignment) / this->Nx;
        int i_row = (column - row_alignment) - j_row * this->Nx;

        T distance, distance_x, distance_y, distance_z, balgovind_z;
        for (int k = 0; k < Nz_state; k++)
          {
            if (k == k_row)
              balgovind_z = 1.;
            else if (Balgovind_vertical_scale_model(species_index) == 0.)
              balgovind_z = 0.;
            else
              {
                distance_z = abs(this->GridZ3D(k) - this->GridZ3D(k_row));
                balgovind_z
                  = (1. + distance_z
                     / Balgovind_vertical_scale_model(species_index))
                  * exp(-distance_z
                        / Balgovind_vertical_scale_model(species_index));
              }
            for (int j = 0; j < this->Ny; j++)
              for (int i = 0; i < this->Nx; i++)
                {
                  distance_x = T(i - i_row) * this->Delta_x;
                  distance_y = T(j - j_row) * this->Delta_y;
                  distance = sqrt(distance_x * distance_x
                                  + distance_y * distance_y);
                  error_covariance_vector(species_alignment
                                          + k * this->Ny * this->Nx
                                          + j * this->Nx + i)
                    = (1. + distance / Balgovind_scale_model(species_index))
                    * exp(- distance / Balgovind_scale_model(species_index))
                    * balgovind_z * model_error_variance(species_index);
                }
          }
        return error_covariance_vector;
      }
    else if (error_covariance_model == "diagonal_constant")
      {

        /*** Only for one species ***/

        int species_point_number = Nz_state * this->Nx * this->Ny;
        int species_index =  column / species_point_number;

        for (int i = 0; i < Nstate; i++)
          error_covariance_vector(i) = 0.;
        error_covariance_vector(column) = model_error_variance(species_index);
        return error_covariance_vector;
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRYASSIMCONC_CXX
#endif
