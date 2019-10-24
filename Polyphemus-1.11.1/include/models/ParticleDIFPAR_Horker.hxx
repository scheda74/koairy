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


#ifndef POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_HORKER_HXX


#include "BaseParticle.hxx"


namespace Polyphemus
{


  ////////////////////////////////////////////////////////
  // DIFPAR PARTICLE WITH THE FOKKER-PLANCK FORMULATION //
  ////////////////////////////////////////////////////////


  //! Class for DIFPAR particles with the Horker formulation.
  template<class T>
  class ParticleDIFPAR_Horker: public BaseParticle<T>
  {

  protected:

    /*** Data ***/

    //! Altitude in the curvilinear coordinate system (meters).
    T z_;

    //! Ordinate in the curvilinear coordinate system (meters).
    T y_;

    //! Abscissa in the curvilinear coordinate system (meters).
    T x_;

    //! Indices of species.
    vector<int> species_index_;

    //! Total mass per species (mass).
    vector<T> quantity_;

    //! Time since release (seconds).
    T age_;

    //! Horizontal standard deviation (meters).
    T sigma_h2_;

    //! Vertical standard deviation (meters).
    T sigma_v2_;



  public:

    /*** Constructor and destructor ***/

    ParticleDIFPAR_Horker(T z, T lat, T lon,
                          const vector<int>& species_index,
                          const vector<T>& quantity, T age);
    virtual ~ParticleDIFPAR_Horker() {};

    /*** Methods ***/

    bool IsOutside(T x_min, T x_max, T y_min, T y_max) const;
    virtual T GetConcentrationContributionOnPoint(int s, T z, T lat, T lon)
      const;

    template<class ClassModel>
    void Transport(ClassModel* Model);

    virtual void Update(T time_elapsed);

  protected:

    T GetLongitude() const;
    T GetLatitude() const;
    virtual T GetInterpolation(Data<T, 3>& D3) const;
    virtual T GetInterpolationDerivativeFirstDim(Data<T, 3>& D3) const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_HORKER_HXX
#endif
