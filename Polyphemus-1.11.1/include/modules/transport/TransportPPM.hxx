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


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_TRANSPORTPPM_HXX


#include <vector>
#include "AtmoData.hxx"
#include "BaseModule.cxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////
  // TRANSPORTPPM //
  //////////////////


  //! This class is a numerical solver for transport.
  /*! It uses the PPM (piecewise parabolic method) third order scheme for
    advection. It is also possible to use the upwind algorithm for given
    species.
  */
  template<class T>
  class TransportPPM: public BaseModule
  {

  protected:

    //! List of species with PPM advection.
    vector<string> species_list_ppm;
    //! Number of species with PPM.
    int Ns_ppm;
    //! Array of indices of species with PPM.
    Array<int, 1> species_list_ppm_index;

    //! Number of vertical levels in the model.
    int Nz;
    //! Number of cells along y in the model.
    int Ny;
    //! Number of cells along x in the model.
    int Nx;

    //! Fluxes on cells West sides.
    Data<T, 3> WestFlux;
    //! Fluxes on cells East sides.
    Data<T, 3> EastFlux;
    //! Fluxes on cells South sides.
    Data<T, 3> SouthFlux;
    //! Fluxes on cells North sides.
    Data<T, 3> NorthFlux;

    //! Incoming vertical fluxes.
    Data<T, 4> VerticalFlux_in;
    //! Outgoing vertical fluxes.
    Data<T, 3> VerticalFlux_out;

    //! Wind on cells West sides.
    Data<T, 3> WestWind;
    //! Wind on cells East sides.
    Data<T, 3> EastWind;
    //! Wind on cells South sides.
    Data<T, 3> SouthWind;
    //! Wind on cells North sides.
    Data<T, 3> NorthWind;
    //! Vertical wind.
    Data<T, 3> VerticalWind;

    //! Extended air density.
    Data<T, 3> AirDensity_ext;

  public:

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void InitStep(ClassModel& Model);

    template<class ClassModel>
    void LossProduction(ClassModel& Model, int s, int k, int j, int i,
                        T& loss, T& production);
    template<class ClassModel>
    void LossProductionUpwind(ClassModel& Model, int s, int k, int j, int i,
                              T& loss, T& production);
    template<class ClassModel>
    void LossProductionPPM(ClassModel& Model, int s, int k, int j, int i,
                           T& loss, T& production);

    template<class ClassModel>
    void Reconstruct(ClassModel& Model, int s, int k, int j,
                     int i, T y, T& conc_left, T& conc_right,
                     T& delta_conc, T& parab, T& x, bool zonal);

    T PPMInterpolation(T c1, T c2, T c3, T c4);
    T PPMDelta(T c1, T c2, T c3);
    template<class ClassModel>
    T MeanConcLeft(ClassModel& Model, int s, int k, int j, int i,
                   T y, bool zonal);
    template<class ClassModel>
    T MeanConcRight(ClassModel& Model, int s, int k, int j, int i,
                    T y, bool zonal);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_TRANSPORTPPM_HXX
#endif
