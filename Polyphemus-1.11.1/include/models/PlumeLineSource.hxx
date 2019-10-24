// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): R�gis Briant
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

#ifndef POLYPHEMUS_FILE_MODELS_PLUMELINESOURCE_HXX


#include <list>
#include "PlumeSource.cxx"

namespace Polyphemus
{


  using namespace std;


  ///////////////////////
  // PLUME LINE SOURCE //
  ///////////////////////


  //! This class stores a line source description.
  template<class T>
  class PlumeLineSource: public PlumeSource<T>
  {
  protected:

    //! Abscissa of the second extremity of line source.
    T x2_;

    //! Ordinate of the second extremity of line source.
    T y2_;

    //! Height of the second extremity of line source.
    T z2_;

    // Coefficient for the combination between point sources and line sources.
    T combination_factor_;

    //! Width of line source.
    T width_;

    //! Velocity of line source.
    T VehicleVelocity_;

    //! Area of line source.
    T Area_;

    //! Density of line source.
    T Density_;

    //! Id of the source section.
    T id_section_;


  public:

    static T  vertical_spread(T z);
    static T  vertical_spread(T z, T width, T area, T density,
                              T vehicle_velocity);

  public:

    /*** Constructor and destructor ***/
    PlumeLineSource(T rate, T x, T y, T z, T x2, T y2, T z2, T width,
                    T VehicleVelocity, T Area, T Density, int species_index,
                    int id_source = 0, int id_section = 0);
    virtual ~PlumeLineSource();

    /*** Methods ***/

    T GetX2() const;
    void SetX2(T x2);
    T GetY2() const;
    void SetY2(T y2);
    T GetZ2() const;
    void SetZ2(T z2);
    T GetCombinationFactor() const;
    void SetCombinationFactor(T combination_factor);
    T GetWidth() const;
    virtual int GetIdSection() const;
    T GetVehicleVelocity() const;
    T GetArea() const;
    T GetDensity() const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PLUMELINESOURCE_HXX
#endif
