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



#ifndef POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_HORKER_CXX


#include "ParticleDIFPAR_Horker.hxx"


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
  ParticleDIFPAR_Horker<T>
  ::ParticleDIFPAR_Horker(T z, T lat, T lon,
                          const vector<int>& species_index,
                          const vector<T>& quantity, T age):
    z_(z), species_index_(species_index), quantity_(quantity),
    age_(age), sigma_h2_(0.), sigma_v2_(0.)
  {
    const T earth_radius(6371229.);
    const T pi(3.1416);

    x_ = earth_radius * lon * pi / 180.;
    y_ = earth_radius * sin(lat * pi / 180.);
  }


  /////////////
  // METHODS //
  /////////////


  //! Checks whether the particle is outside the simulation domain.
  /*!
    \param x_min Abscissa in the curvilinear coordinate system of the lower
    domain boundary.
    \param x_max Abscissa in the curvilinear coordinate system of the upper
    domain boundary.
    \param y_min Ordinate in the curvilinear coordinate system of the lower
    domain boundary.
    \param y_max Ordinate in the curvilinear coordinate system of the upper
    domain boundary.
  */
  template<class T>
  bool ParticleDIFPAR_Horker<T>
  ::IsOutside(T x_min, T x_max, T y_min, T y_max) const
  {
    return (x_ < x_min || x_ > x_max || y_ < y_min || y_ > y_max);
  }


  //! Computes the particle contribution to concentration on a given point.
  /*!
    \param s Index of the specie whose concentration is to be computed.
    \param z Altitude of the point where the contribution is computed.
    \param lat Latitude of the point where the contribution is computed.
    \param lon Longitude of the point where the contribution is computed.
  */
  template<class T>
  T ParticleDIFPAR_Horker<T>
  ::GetConcentrationContributionOnPoint(int s, T z, T lat, T lon)
    const
  {
    const T earth_radius(6371229.);
    const T pi(3.1416);

    T lon_distance = earth_radius * lon * pi / 180. - x_;
    T lat_distance = earth_radius * sin(lat * pi / 180.) - y_;

    T horizontal_distance_square = lon_distance * lon_distance
      + lat_distance * lat_distance;

    horizontal_distance_square /= (2 * sigma_h2_);
    T vertical_distance_square = (z - z_) * (z - z_);
    vertical_distance_square /= (2 * sigma_v2_);

    T impact_distance = horizontal_distance_square + vertical_distance_square;
    T coeff = pow(2 * pi, 1.5);
    T result =  quantity_[s] / (coeff * sigma_h2_ * sqrt(sigma_v2_))
      * exp(-impact_distance);

    return result;
  }


  //! Performs the stochastic transport of the particle given the related
  // model.
  /*! param Model Model object.
   */
  template<class T>
  template<class ClassModel>
  void ParticleDIFPAR_Horker<T>
  ::Transport(ClassModel* Model)
  {
    T U(GetInterpolation(Model->D3("ZonalWind_i")));
    T V(GetInterpolation(Model->D3("MeridionalWind_i")));
    T W(GetInterpolation(Model->D3("VerticalWind_i")));
    T dwx(Model->GetRandomNumber());
    T dwy(Model->GetRandomNumber());
    T dwz(Model->GetRandomNumber());
    T delta_t(min(Model->GetDelta_t(), age_));

    T rho_vertical_diffusion
      (GetInterpolation(Model->D3("VerticalDiffusionCoefficient_i")));
    rho_vertical_diffusion = (rho_vertical_diffusion < 0.2) ?
      0.2 : rho_vertical_diffusion;

    T gradient_vertical_diffusion
      (GetInterpolationDerivativeFirstDim
       (Model->D3("VerticalDiffusionCoefficient_i")));
    T air_density(GetInterpolation(Model->D3("AirDensity_i")));

    T pi(3.1416);
    T earth_radius(6371229.);
    T lat = GetLatitude() * pi / 180.;

    // Computes the displacement.
    x_ += U * delta_t;
    y_ += V * delta_t;
    z_ += (W + gradient_vertical_diffusion / air_density) * delta_t
      + sqrt(2. * rho_vertical_diffusion * delta_t / air_density) * dwz;

    // Computes the vertical standard deviation.
    sigma_v2_ += delta_t * max(rho_vertical_diffusion / air_density, 0);
    sigma_h2_ = 0.25 * age_ * age_;
  }


  //! Updates the state of the particle.
  /*! param time_elapsed Time elapsed since the previous update.
   */
  template<class T>
  void ParticleDIFPAR_Horker<T>
  ::Update(T time_elapsed)
  {
    age_ += time_elapsed;

    // Reflection at ground level.
    z_ = (z_ >= 0) ? z_ : -z_;
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Returns the longitude of the particle.
  template<class T>
  T ParticleDIFPAR_Horker<T>
  ::GetLongitude() const
  {
    const T earth_radius(6371229.);
    const T pi(3.1416);
    return x_ / earth_radius * 180. / pi;
  }


  //! Returns the latitude of the particle.
  template<class T>
  T ParticleDIFPAR_Horker<T>
  ::GetLatitude() const
  {
    const T earth_radius(6371229.);
    const T pi(3.1416);
    return asin(y_ / earth_radius) * 180. / pi;
  }


  //! Interpolates the values of a Data object at the particle location.
  /*!
    \param D3 3D Data object.
  */
  template<class T>
  T ParticleDIFPAR_Horker<T>
  ::GetInterpolation(Data<T, 3>& D3) const
  {
    Array<T, 1> Coord(3);
    Coord(0) = z_;
    Coord(1) = GetLatitude();
    Coord(2) = GetLongitude();

    T value;
    LinearInterpolationPoint(D3, Coord, value);
    return value;
  }


  //! Computes at the particle location the derivative of a Data object with
  // respect to its first dimension.
  /*!
    \param D3 3D-Data object.
  */
  template<class T>
  T ParticleDIFPAR_Horker<T>
  ::GetInterpolationDerivativeFirstDim(Data<T, 3>& D3) const
  {
    Array<T, 1> Coord(3);
    Coord(1) = GetLatitude();
    Coord(2) = GetLongitude();

    int upper_index(0), lower_index;

    while (upper_index < D3.GetLength(0) && D3[0](upper_index) < z_)
      upper_index++;

    if (upper_index == D3.GetLength(0))
      upper_index--;
    else if (upper_index == 0)
      upper_index = 1;

    lower_index = upper_index - 1;

    T upper_value, lower_value;
    Coord(0) = D3[0](upper_index);
    LinearInterpolationPoint(D3, Coord, upper_value);
    Coord(0) = D3[0](lower_index);
    LinearInterpolationPoint(D3, Coord, lower_value);

    T gradient_value = (upper_value - lower_value) /
                         (D3[0](upper_index) - D3[0](lower_index));
    return gradient_value;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PARTICLEDIFPAR_HORKER_CXX
#endif
