// Copyright (C) 2009, ENPC - INRIA - EDF R&D
// Author(s): Pierre Tran
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

// This file is part of a Lagrangian model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_FOKKERPLANCK_HXX


#include "ParticleDIFPAR_Horker.cxx"


namespace Polyphemus
{


  ////////////////////////////////////////////////////////
  // DIFPAR PARTICLE WITH THE FOKKER-PLANCK FORMULATION //
  ////////////////////////////////////////////////////////


  //! Class for DIFPAR particles with the Fokker-Planck formulation.
  template<class T>
  class ParticleDIFPAR_FokkerPlanck: public ParticleDIFPAR_Horker<T>
  {

  public:

    /*** Constructor and destructor ***/
    ParticleDIFPAR_FokkerPlanck(T z, T lat, T lon,
                                const vector<int>& species_index,
                                const vector<T>& quantity, T age);
    ~ParticleDIFPAR_FokkerPlanck() {};

    /*** Methods ***/

    template<class ClassModel>
    void Transport(ClassModel* Model);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_FOKKERPLANCK_HXX
#endif
