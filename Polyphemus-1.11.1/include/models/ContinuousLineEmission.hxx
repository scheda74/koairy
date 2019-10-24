// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_HXX



//////////////
// INCLUDES //
//////////////


#include "PointEmissionUnit.hxx"

namespace Polyphemus
{
  using namespace AtmoData;


  ////////////////////
  // LINE EMISSIONS //
  ////////////////////


  //! This class manages point emissions.
  template<class T>
  class ContinuousLineEmission: public PointEmissionUnit<T>
  {

  protected:

    //! Ending date.
    Date date_end;

    //! Emission rate (mass/s).
    vector<vector<T> > rate;

    //! List of Coordinate.
    list<Array<T, 1> > coordinate_list_;

    //! Line sources widths (m).
    vector<T> Width;

    //! Line sources Velocitys (m/s).
    vector<T> VehicleVelocity;

    //! Line sources Area (mÂ²).
    vector<T> Area;

    //! Line sources Density (veh/m).
    vector<T> Density;

  public:

    ContinuousLineEmission();
    virtual ~ContinuousLineEmission();
    void Init(ConfigStream& config, vector<string> species_list);
    bool HasEnded(Date date);
    Date GetEndDate() const;
    void GetCoordinates(list<Array<T, 1> >& coordinate_list) const;
    T GetRate(int species, int id_section) const;
    T GetDiameter() const;
    T GetWidth(int id_section) const;
    void SetCoordinates(const list<Array<T, 1> >& coordinate_list);
    void GetPlumeRiseParam(T& velocity, T& temperature, T& diameter);
    T GetVehicleVelocity(int id_section) const;
    T GetArea(int id_section) const;
    T GetDensity(int id_section) const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_HXX
#endif
