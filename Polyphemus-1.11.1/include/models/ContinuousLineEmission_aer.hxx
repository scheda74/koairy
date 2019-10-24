// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_AER_HXX



//////////////
// INCLUDES //
//////////////


#include "ContinuousLineEmission.cxx"

namespace Polyphemus
{
  using namespace AtmoData;


  ////////////////////
  // LINE EMISSIONS //
  ////////////////////


  //! This class manages point emissions.
  template<class T>
  class ContinuousLineEmission_aer: public ContinuousLineEmission<T>
  {

  protected:

    // //! Ending date.
    // Date date_end;

    //! Aerosol emission rate (mass/s).
    vector<T> rate_aer;

    // //! List of Coordinate.
    // list<Array<T, 1> > coordinate_list_;

  public:

    ContinuousLineEmission_aer();
    void Init(ConfigStream& config, vector<string> species_list,
              vector<string> species_list_aer, int Nbin_aer);
    // bool HasEnded(Date date);
    // void GetCoordinates(list<Array<T, 1> >& coordinate_list) const;
    // T GetRate(int species) const;
    // T GetDiameter() const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CONTINUOUSLINEEMISSION_AER_HXX
#endif
