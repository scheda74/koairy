// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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

// This file is part of the Eulerian model Polair3D.


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRY_HXX


#include <vector>

#include "AtmoDataHeader.hxx"
#include "Polair3DTransport.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ///////////////////////
  // POLAIR3DCHEMISTRY //
  ///////////////////////


  //! This class is a solver for an advection-diffusion-reaction equation.
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  class Polair3DChemistry:
    public Polair3DTransport<T, ClassAdvection, ClassDiffusion>
  {

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >::iterator
    input_files_iterator;

  public:

    /*** Configuration ***/

    //! List of species with forced concentrations.
    vector<string> species_list_forced;
    //! List of species with photolysis reactions.
    vector<string> photolysis_reaction_list;
    //! List of altitudes at which photolysis rates are available.
    vector<string> altitudes_photolysis;
    //! With source splitting?
    bool source_splitting;

    /*** Photolysis rates and attenuation ***/

    //! Number of photolysis reactions.
    int Nr_photolysis;
    //! Grid for photolysis reactions.
    RegularGrid<T> GridR_photolysis;

    //! Starting date of photolysis-rates files.
    Date photolysis_date_min;
    //! Time step of photolysis-rates files in seconds.
    T photolysis_delta_t;
    //! Number of days in photolysis-rates files.
    int Nphotolysis_days;

    //! First time angle in photolysis-rates data.
    T photolysis_time_angle_min;
    //! Time-angle step in photolysis-rates data.
    T photolysis_delta_time_angle;
    //! Number of time angles in photolysis-rates data.
    int Nphotolysis_time_angle;
    //! Grid for time angles of photolysis rates.
    RegularGrid<T> Grid_time_angle_photolysis;

    //! First latitude in photolysis-rates data.
    T photolysis_latitude_min;
    //! Latitude step in photolysis-rates data.
    T photolysis_delta_latitude;
    //! Number of latitudes in photolysis-rates data.
    int Nphotolysis_latitude;
    //! Grid for latitudes of photolysis rates.
    RegularGrid<T> Grid_latitude_photolysis;

    //! Number of vertical levels in photolysis-rates data.
    int Nphotolysis_z;
    //! Grid for altitudes of photolysis rates.
    RegularGrid<T> GridZ_photolysis;

    //! Photolysis rates at current date.
    Data<T, 4> PhotolysisRate_i;
    //! Photolysis rates at next date.
    Data<T, 4> PhotolysisRate_f;

    //! Attenuation coefficients at current date.
    Data<T, 3> Attenuation_i;
    //! Attenuation coefficients at next date.
    Data<T, 3> Attenuation_f;
    //! Attenuation buffer.
    Data<T, 3> FileAttenuation_i;
    //! Attenuation buffer.
    Data<T, 3> FileAttenuation_f;
    //! Photolysis buffer (used when photolysis is computed).
    Data<T, 4> FilePhotolysis_i;
    //! Photolysis buffer (used when photolysis is computed).
    Data<T, 4> FilePhotolysis_f;

    // Option to determine how photolysis rates are computed.
    string computed_photolysis;

    /*** Forced concentrations ***/

    //! Number of species with forced concentrations.
    int Ns_forced;
    //! Grid for species with forced concentrations.
    RegularGrid<T> GridS_forced;
    //! Forced concentrations at current date.
    Data<T, 4> ForcedConcentration_i;
    //! Forced concentrations at next date.
    Data<T, 4> ForcedConcentration_f;
    //! Forced concentrations buffer.
    Data<T, 4> FileForcedConcentration_i;
    //! Forced concentrations buffer.
    Data<T, 4> FileForcedConcentration_f;

    /*** Sources (source splitting) ***/

    //! Sources at current date (for source splitting).
    Data<T, 4> Source_i;
    //! Sources at next date (for source splitting).
    Data<T, 4> Source_f;

    /*** Chemical mechanism ***/

    //! Chemical mechanism and chemical numerical scheme.
    ClassChemistry Chemistry_;

    /*** Adjoint data ***/

    //! Adjoint data for sources at current date (for source splitting).
    Data<T, 4> Source_i_ccl;
    //! Adjoint data for sources at next date (for source splitting).
    Data<T, 4> Source_f_ccl;

  public:

    /*** Constructor and destructor ***/

    Polair3DChemistry();
    Polair3DChemistry(string config_file);
    void Construct(string config_file);
    virtual ~Polair3DChemistry();

    /*** Configuration ***/

    void ReadConfiguration();
    virtual void CheckConfiguration();

    bool HasForcedConcentration(int s) const;
    bool HasForcedConcentration(string name) const;
    int ForcedConcentrationIndex(int s) const;
    int ForcedConcentrationIndex(string name) const;
    string ForcedConcentrationName(int s) const;
    int ForcedConcentrationGlobalIndex(int s) const;

    /*** Initializations ***/

    void Allocate();
    void Init();
    void InitPhotolysis(Date date, Data<T, 4>& Rates);
    void InitStep();

    virtual void SetDate(Date date);

    /*** Integration ***/

    void Chemistry();
    void Chemistry_b();

    void Forward();
    void SetBackward(bool flag);
    void Backward();

    /*** Access methods ***/

    vector<string> GetPhotolysisReactionList() const;
    int GetNr_photolysis() const;

    int GetNs_source();
    int GetNz_source();
    int SourceGlobalIndex(int s) const;
    Data<T, 4>& GetSource_i();
    Data<T, 4>& GetSource_f();

  protected:

    virtual void InitAllData();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DCHEMISTRY_HXX
#endif
