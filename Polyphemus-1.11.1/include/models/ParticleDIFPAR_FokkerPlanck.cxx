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



#ifndef POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_FOKKERPLANCK_CXX


#include "ParticleDIFPAR_FokkerPlanck.hxx"


namespace Polyphemus
{

  /////////////////
  // CONSTRUCTOR //
  /////////////////


  //! Main constructor.
  /*!
    \param z Altitude of the particle (meter).
    \param lat Latitude of the particle (degree).
    \param lon Longitude of the particle (degree).
    \species_index Indices of the species carried by the particle.
    \quantity Quantities of the species carried by the particle.
    \age Age of the particle (s).
  */
  template<class T>
  ParticleDIFPAR_FokkerPlanck<T>
  ::ParticleDIFPAR_FokkerPlanck(T z, T lat, T lon,
                                const vector<int>& species_index,
                                const vector<T>& quantity, T age):
    ParticleDIFPAR_Horker<T>(z, lat, lon, species_index, quantity, age)
  {
  }


  //! Performs the stochastic transport of the particle given the related
  // model.
  /*! param Model Model object.
   */
  template<class T>
  template<class ClassModel>
  void ParticleDIFPAR_FokkerPlanck<T>
  ::Transport(ClassModel* Model)
  {
    T U(GetInterpolation(Model->D3("ZonalWind_i")));
    T V(GetInterpolation(Model->D3("MeridionalWind_i")));
    T W(GetInterpolation(Model->D3("VerticalWind_i")));
    T dwx(Model->GetRandomNumber());
    T dwy(Model->GetRandomNumber());
    T dwz(Model->GetRandomNumber());
    T delta_t(min(Model->GetDelta_t(), this->age_));

    T rho_vertical_diffusion
      (GetInterpolation(Model->D3("VerticalDiffusionCoefficient_i")));
    rho_vertical_diffusion = (rho_vertical_diffusion < 0.2) ?
      0.2 : rho_vertical_diffusion;

    T gradient_vertical_diffusion
      (GetInterpolationDerivativeFirstDim
       (Model->D3("VerticalDiffusionCoefficient_i")));
    T horizontal_diffusion(Model->GetHorizontalDiffusion());
    T gaussian_kernel_horizontal_diffusion
      (Model->GetGaussianKernelHorizontalDiffusion());
    T air_density(GetInterpolation(Model->D3("AirDensity_i")));

    T pi(3.1416);
    T earth_radius(6371229.);
    T lat = this->GetLatitude() * pi / 180.;

    // Computes the displacement.
    this->x_ += U * delta_t
      + sqrt(2. * horizontal_diffusion * delta_t) / cos(lat) * dwx;
    this->y_ += V * delta_t
      - 2. * horizontal_diffusion / earth_radius * sin(lat) * delta_t
      + cos(lat) * sqrt(2. * horizontal_diffusion * delta_t) * dwy;
    this->z_ += (W + gradient_vertical_diffusion / air_density) * delta_t
      + sqrt(2. * rho_vertical_diffusion * delta_t / air_density) * dwz;

    // Computes the vertical standard deviation.
    this->sigma_h2_ +=
      delta_t * max(gaussian_kernel_horizontal_diffusion, 0.2);
    this->sigma_v2_ += delta_t * max(rho_vertical_diffusion / air_density, 0);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_FOKKERPLANCK_CXX
#endif
