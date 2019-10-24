// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi, Vivien Mallet, Denis Quélo
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_DECAY_HXX


#include <vector>
#include <map>
#include <math.h>
#include "AtmoData.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ///////////
  // DECAY //
  ///////////


  //! This class is a solver for decay (radioactive or so).
  template<class T>
  class Decay: public BaseModule
  {

  protected:

    //! Number of species.
    int Ns;
    //! Species list.
    vector<string> species_list;
    //! Number of aerosol species.
    int Ns_aer;
    //! Aerosol species list.
    vector<string> species_list_aer;
    //! Number of aerosol bins.
    int Nbin_aer;

    //! Is a filiation matrix used?
    bool with_filiation_matrix;
    //! Are two half lives used?
    bool with_time_dependence;

    //! Species half-lives (in days).
    Array<T, 1> half_life;
    //! Species half-lives during the day (in days).
    Array<T, 1> half_life_day;
    //! Species half-lives during the night (in days).
    Array<T, 1> half_life_night;

    //! Aerosols half-lives (in days).
    Array<T, 1> half_life_aer;
    //! Aerosols half-lives during the day (in days).
    Array<T, 1> half_life_day_aer;
    //! Aerosols half-lives during the night (in days).
    Array<T, 1> half_life_night_aer;

    //! Filiation matrix: decay constants.
    Array<T, 2> filiation_matrix;
    //! Filiation matrix for aerosols: decay constants.
    Array<T, 2> filiation_matrix_aer;

  public:

    /*** Constructor ***/

    Decay();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);

    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_DECAY_HXX
#endif
