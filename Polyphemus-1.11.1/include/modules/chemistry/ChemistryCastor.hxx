// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYCASTOR_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModule.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // CHEMISTRYCASTOR //
  /////////////////////


  //! This class is the default chemical module for Castor.
  template<class T>
  class ChemistryCastor: public BaseModule
  {

  protected:

    //! Number of vertical levels in the model.
    int Nz;
    //! Number of cells along y in the model.
    int Ny;
    //! Number of cells along x in the model.
    int Nx;

    //! Number of species.
    int Ns;
    //! List of species names.
    vector<string> species_list;

    //! Path to the file with chemical reactions.
    string reaction_file;
    //! Path to stoichiometry description.
    string stoichiometry_file;
    //! Path to tabulated photolysis rates.
    string photolysis_file;
    //! Path to the file containing reaction rates.
    string rates_file;

    //! Number of reactions.
    int Nr;
    //! Number of external species.
    int Ns_ext;
    int Ntemps;

    Data<int, 1> kreacp;
    Data<int, 1> kreacl;
    Data<int, 2> ireacp;
    Data<int, 2> ireacl;
    Data<int, 1> nreactants;
    Data<int, 1> nprods;
    Data<int, 2> irctt;
    Data<T, 3> stoi;

    int Ntab_phot;
    int Nphot;
    int Nwave;
    int Nz_phot;

    Data<T, 1> zenang;
    Data<T, 1> zetaref;
    Data<T, 1> altiphot;
    Data<T, 3> photoj;
    Data<int, 1> iphoto;

    Data<int, 1> ltabrate;
    Data<T, 2> tabrate;
    Data<int, 1> ityperate;

    Data<int, 3> istoit;
    Data<T, 3> wgstl;
    Data<T, 3> wgsth;

    //! Chemical rates.
    Data<T, 4> rate;
    //! Photolysis rates.
    Data<T, 4> photorate;

  public:

    /*** Constructor ***/

    ChemistryCastor();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void InitStep(ClassModel& Model);

    template<class ClassModel>
    void ComputePhotorate(ClassModel& Model);
    template<class ClassModel>
    void ComputeRate(ClassModel& Model);

    template<class ClassModel>
    void LossProduction(ClassModel& Model, int s, int k, int j, int i,
                        T& loss, T& production);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYCASTOR_HXX
#endif
