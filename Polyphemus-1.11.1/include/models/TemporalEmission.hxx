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


#ifndef POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_HXX



//////////////
// INCLUDES //
//////////////


#include "PointEmissionUnit.hxx"
#include "ContinuousEmission.hxx"

namespace Polyphemus
{
  using namespace AtmoData;


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages point emissions.
  template<class T>
  class TemporalEmission: public ContinuousEmission<T>
  {

  protected:

    //! Temporal factors.
    Array<T, 1> TemporalFactor;

    //! File containing temporal factors.
    string input_file;

    //! File time step.
    T Delta_t_file;

    //! File beginning date.
    Date date_min_file;

    //! Current time step in the file.
    int current_record;

  public:

    virtual void Init(ConfigStream& config, vector<string> species_list);
    void GetEmission(Date date_beg, Date date_end, int s,
                     Array<T, 2>& point_emission);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_TEMPORALEMISSION_HXX
#endif
