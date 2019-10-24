// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_PUFFEMISSION_HXX


//////////////
// INCLUDES //
//////////////


#include "PointEmissionUnit.hxx"

namespace Polyphemus
{
  using namespace AtmoData;


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages point emissions.
  template<class T>
  class PuffEmission: public PointEmissionUnit<T>
  {

  protected:

    //! Ending date.
    Date date_end;

    //! Emitted mass (mass).
    vector<T> quantity;

    //! Source abscissa.
    T x;

    //! Source ordinate.
    T y;

    //! Source altitude.
    T z;

  public:

    PuffEmission();
    void Init(ConfigStream& config, vector<string> species_list);
    void GetCoordinates(T& abscissa, T& ordinate, T& height) const;
    void SetCoordinates(T abscissa, T ordinate, T height);
    bool HasEnded(Date date);
    void GetEmission(Date date_beg, Date date_end, int s,
                     Array<T, 2>& point_emission);
    void GetEmission(int s, Array<T, 2>& point_emission);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFEMISSION_HXX
#endif
