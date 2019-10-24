// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_AER_HXX



//////////////
// INCLUDES //
//////////////


#include "PointEmissionUnit.cxx"
#include "ContinuousEmission.cxx"

namespace Polyphemus
{
  using namespace AtmoData;


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages point emissions.
  template<class T>
  class ContinuousEmission_aer: public ContinuousEmission<T>
  {

  protected:

    //! Aerosol emission rate (mass/s).
    vector<T> rate_aer;

  public:

    ContinuousEmission_aer();
    void Init(ConfigStream& config, vector<string> species_list,
              vector<string> species_list_aer, int Nbin_aer);
    void GetEmission_aer(Date current_date,
                         Date next_date, int s, T& quantity);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CONTINUOUSEMISSION_AER_HXX
#endif
