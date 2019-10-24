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

// This file is part of a Gaussian plume model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_AER_HXX


#include "GaussianPlume.cxx"
#include <list>
#include <vector>


namespace Polyphemus
{


  using namespace std;


  ////////////////////
  // GAUSSIAN PLUME //
  ////////////////////


  //! Plume Gaussian model for aerosols.
  template<class T>
  class GaussianPlume_aer: public GaussianPlume<T>
  {

  protected:

    //! File containing the aerosol diameters.
    string file_diameter;
    //! List of aerosol diameters.
    vector<T> diameter_list;

    //! File containing the source data for aerosol species.
    string file_source_aer;
    //! List of aerosol sources.
    list<PlumeSource<T> > SourceList_aer;
    //! Number of aerosol sources.
    int Nsource_aer;

    //! List of aerosol half-lives (s).
    Array<T, 1> half_life_time_aer;

    //! List of aerosol biological half-lives (s).
    Array<T, 1> biological_half_life_time_aer;

    //! List of aerosol scavenging coefficients ( s^(-1) ).
    Array<T, 2> scavenging_coefficient_aer;

    //! List of aerosol deposition velocities ( m/s ).
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

    /*** Main constructor ***/

    GaussianPlume_aer(string config_file);

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    virtual void Init(bool read_all_input = false);
    virtual void InitSource_aer(string source_file_aer);
    virtual void InitMeteo(ConfigStream& meteo, bool show_meteo);

    /*** Computations ***/

    virtual void InitCompute();
    virtual void ComputePlumeRise_aer();
    virtual void Compute();
    virtual T GetConcentration_aer(int species, int diameter, T z, T y, T x);
    Data<T, 5>& GetConcentration_aer();

    void ComputeLossFactor_aer(T distance, T z,
                               int species, int diam,
                               T& loss_factor, T& overcamp_factor);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_AER_HXX
#endif
