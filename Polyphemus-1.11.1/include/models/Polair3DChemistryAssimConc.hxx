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


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRYASSIMCONC_HXX


#include <vector>
#include "AtmoData.hxx"
#include "Polair3DChemistry.cxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////////
  // POLAIR3DCHEMISTRYASSIMCONC //
  ////////////////////////////////


  /*! \brief This class is a solver for an advection-diffusion-reaction
    equation, and it is dedicated to data assimilation.
  */
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  class Polair3DChemistryAssimConc:
    public Polair3DChemistry < T, ClassAdvection,
                               ClassDiffusion, ClassChemistry >
  {

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >::iterator
    input_files_iterator;

  protected:

    /*** Domain ***/

    //! State species list.
    vector<string> state_species_list;
    //! Number of state species.
    int Ns_state;
    //! State species indices.
    Array<int, 1> state_species_index;
    //! Number of levels.
    int Nz_state;
    //! Levels list of state.
    Array<int, 1> state_levels;

    //! Dimension of state vector.
    int Nstate;

    /*** Covariance model ***/

    //! Model for concentrations error covariances.
    string error_covariance_model;
    //! Balgovind scale for background covariance.
    Array<T, 1> Balgovind_scale_background;
    //! Balgovind vertical scale for background covariance.
    Array<T, 1> Balgovind_vertical_scale_background;
    //! Background error variance.
    Array<T, 1> background_error_variance;

    //! Balgovind scale for model covariance.
    Array<T, 1> Balgovind_scale_model;
    //! Balgovind vertical scale for model covariance.
    Array<T, 1> Balgovind_vertical_scale_model;
    //! Model error variance.
    Array<T, 1> model_error_variance;

    //! Number of the row of B currently stored.
    int current_row;
    //! Number of the column of Q currently stored.
    int current_column;
    //! Value of the row of B or of the column of Q currently stored.
    Array<T, 1> error_covariance_vector;

  public:

    /*** Constructor and destructor ***/

    Polair3DChemistryAssimConc(string config_file);
    virtual ~Polair3DChemistryAssimConc();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Initializations ***/

    void Allocate();
    void Init();

    /*** Set or get state ***/

    void GetState(Array<T, 1>& state_vector);
    void SetState(const Array<T, 1>& state_vector);

    void GetState_ccl(Array<T, 1>& state_vector);
    void SetState_ccl(const Array<T, 1>& state_vector);

    /*** Access ***/

    string GetErrorCovarianceModel() const;

    int GetNstate() const;
    int GetStateNz() const;
    const Array<int, 1>& GetStateLevels() const;
    int GetStateNs() const;
    vector<string> GetStateSpeciesList() const;

    int GetStateSpeciesIndex(string species);
    int GetStateLevelSequenceIndex(int level);

    /*** Covariance model ***/

    const Array<T, 1>& BackgroundErrorCovariance(int row);
    const Array<T, 1>& ModelErrorCovariance(int column);

    const Array<T, 1>& GetBackgroundErrVar() const;

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRYASSIMCONC_HXX
#endif
