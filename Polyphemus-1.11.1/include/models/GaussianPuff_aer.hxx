// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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

// This file is part of a Gaussian puff model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_AER_HXX


#include "GaussianPuff.cxx"
#include <list>


namespace Polyphemus
{


  using namespace std;


  ///////////////////
  // GAUSSIAN PUFF //
  ///////////////////


  //! Gaussian puff model for aerosol species.
  template<class T>
  class GaussianPuff_aer: public GaussianPuff<T>
  {

  protected:

    //! File containing the aerosol diameters.
    string file_diameter;
    //! List of diameters.
    vector<T> diameter_list;

    //! File containing the puff data for aerosol species.
    string file_puff_aer;
    //! Aerosol puff list.
    vector< vector< Puff<T>* > > PuffList_aer;
    //! Number of aerosol puffs.
    int Npuff_aer;

    //! Array of aerosol half-lives (s).
    Array<T, 1> half_life_time_aer;

    //! Array of aerosol biological half-lives (s).
    Array<T, 2> biological_half_life_time_aer;

    //! Array of aerosol scavenging coefficients ( s^(-1) ).
    Array<T, 2> scavenging_coefficient_aer;

    //! Array of aerosol deposition velocities (m/s).
    Array<T, 2> deposition_velocity_aer;

    //! 5D grid along x.
    RegularGrid<T> GridX5D;
    //! 5D grid along y.
    RegularGrid<T> GridY5D;
    //! 5D grid along z.
    RegularGrid<T> GridZ5D;

    //! 5D grid for aerosol species.
    RegularGrid<T> GridS5D_aer;

    //! 5D grid for diameters.
    RegularGrid<T> GridB5D_aer;

    //! Concentrations of aerosol species.
    Data<T, 5> Concentration_aer;

  public:

    /*** Main constructor and destructor ***/

    GaussianPuff_aer(string config_file);
    ~GaussianPuff_aer();

    /*** Main methods ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    void ClearPuffList_aer();

    virtual void Init();
    virtual void InitStep();
    virtual void InitPuffSource_aer(string file_puff_aer);
    virtual void InitMeteo(ConfigStream& meteo, bool show_meteo);

    void Advection_aer();
    void Diffusion_aer();
    void Forward();

    virtual T GetConcentration_aer(int species, int diameter, T z, T y, T x);
    Data<T, 5>& GetConcentration_aer();

    void ComputeLossFactor_aer(int index, int diam);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_AER_HXX
#endif
