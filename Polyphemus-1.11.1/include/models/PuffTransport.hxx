// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Ir?e Korsakissok, Vivien Mallet, Youngseob Kim
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

// This file is part of a Gaussian puff model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_PUFFTRANSPORT_HXX


namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


#include <list>
  using namespace std;


  //////////////
  // GAUSSIAN //
  //////////////


  //! Class that stores all data about a puff.
  template<class T>
  class Puff
  {

  protected:

    //! Time at which the puff is released (seconds).
    T release_time_;

    //! Time since puff release (seconds).
    T puff_time_;

    //! Total mass released (mass).
    T quantity_;

    /*** Data about the source that released the puff. ***/

    //! Efflux speed of gases.
    T source_velocity_;
    //! Temperature.
    T source_temperature_;
    //! Diameter.
    T source_diameter_;
    //! Source abscissa (meters).
    T source_abscissa_;
    //! Source ordinate (meters).
    T source_ordinate_;
    //! Source altitude (meters).
    T source_height_;
    //! Source water (kg).
    T source_water_;
    //! Volume at previous timestep
    T volume_prev_;
    //! Source id
    string source_id_;

    //! Height of the fraction of puff above BL.
    T z_above_;
    //! Fraction of the puff above BL.
    T penetration_factor_;
    //! Current abscissa (meters).
    T x_;
    //! Current ordinate (meters).
    T y_;
    //! Current altitude (meters).
    T z_;

    //! Downwind distance to the source (meters).
    T distance_;

    // Diffusion coefficients.
    //! Horizontal diffusion coefficient in the downwind direction (meters).
    T sigma_x_;
    //! Horizontal diffusion coefficient in the crosswind direction (meters).
    T sigma_y_;
    //! Vertical diffusion coefficient (meters).
    T sigma_z_;
    //! Square of the initial horizontal plume spread.
    T sigma_y_2_;
    T sigma_x_2_;
    //! Square of the initial vertical plume spread.
    T sigma_z_2_;

    //! Number of species.
    int Ns;

    //! Species index.
    int species_index_;

    //! Temperature of ambient air (Celsius degrees).
    T temperature_;
    //! Angle between positive x-axis and wind, counterclockwise.
    T wind_angle_;
    //! Wind (m/s).
    T wind_;

    //! Friction velocity (m/s).
    T friction_velocity_;
    //! Convective velocity (m/s).
    T convective_velocity_;
    //! Boundary layer height (m).
    T boundary_height_;
    //! Monin Obukhov length (m).
    T lmo_;
    //! Coriolis parameter.
    T coriolis_;

    //! Cosine of wind angle.
    T cos_angle_;
    //! Sine of wind angle.
    T sin_angle_;

    //! Stability class (A, B, C, D, E or F).
    string stability_class_;
    //! Stability (corresponding integer: from 0 for A to 5 for F).
    int stability_;

    //! longitude.
    T longitude_;

    //! latitude.
    T latitude_;

    //! Night or day.
    bool isday_;

    //! Rural or urban.
    bool rural_;

    //! Scavenging coefficient.
    T scavenging_coefficient_;

    //! Deposition velocity.
    T deposition_velocity_;

    //! Fraction of the puff reflected by the ground.
    T reflection_factor_;

    //! Are meteorological data provided to the puff?
    bool has_meteo;

    //! Is the puff from a volume source?
    bool is_volume_source_;

  public:

    /*** Constructor and destructor ***/

    Puff(T time_puff, T velocity, T temperature,
	 T diameter, T quantity, T abscissa, T ordinate, T height,
	 T source_water,T volume_prev,
	 int species_index, string source_id);

    Puff(T time_puff, T velocity, T temperature, T width, T length,
	 T quantity, T abscissa, T ordinate, T height,
	 T source_water,T volume_prev,
	 int species_index, string source_id); 

    Puff(T time_puff, T velocity, T temperature,
         T diameter, T width, T length,
	 T quantity, T abscissa, T ordinate, T height,
	 T source_water, T volume_prev,
	 int species_index, bool is_volume_source, string source_id); 

    Puff(T time_puff, T velocity, T temperature,
         T diameter, T quantity, T abscissa, T ordinate, T height,
         int species_index, string source_id);

    virtual ~Puff();

    /*** Methods ***/

    virtual void InitPuff();

    string GetSourceId() const;
    T GetReleaseTime() const;
    T GetPuffVolume() const;
    T GetSourceVelocity() const;
    T GetSourceTemperature() const;
    T GetSourceDiameter() const;
    T GetSourceHeight() const;
    T GetPuffTime() const;
    int GetSpeciesIndex() const;
    virtual int GetNs();
    virtual int GetNs_aer();
    T GetDistance() const;
    virtual T GetQuantity(int s) const;
    virtual T GetNumberQuantity(int b) const;
    virtual T GetQuantity(int s, int b) const;
    T GetSigma_x() const;
    T GetSigma_y() const;
    T GetSigma_z() const;
    T GetInitialSigma_x() const;
    T GetInitialSigma_y() const;
    T GetInitialSigma_z() const;
    T GetX() const;
    T GetY() const;
    T GetZ() const;
    T GetHeightAboveBL() const;
    T GetHeight() const;
    T GetPenetrationFactor() const;
    T GetEmittedWater() const;
    void SetHeight(T z);
    void SetPuffVolume(T volume_prev);
    void SetHeightAboveBL(T z_above);
    void SetPenetrationFactor(T P);

    bool HasMeteo() const;
    void GetMeteo(T& temperature, T& wind, T& cos_angle, T& sin_angle,
                  int& stability, T& boundary_height,
                  T& longitude, T& latitude, bool& isday, bool& rural);
    void GetMeteo(T& temperature, T& wind, T& cos_angle, T& sin_angle,
                  T& friction_velocity, T& convective_velocity,
                  T& boundary_height, T& lmo, T& coriolis,
                  T& longitude, T& latitude, bool& isday, bool& rural);
    virtual void GetAdditionalMeteo(T& attenuation, T& pressure,
                                    T& specific_humidity) const;
    virtual void GetAdditionalMeteo(T& liquid_water_content) const;
    virtual T GetScavengingCoefficient(int s) const;
    virtual T GetScavengingCoefficient(int s, int b) const;
    virtual T GetDepositionVelocity(int s) const;
    virtual T GetDepositionVelocity(int s, int b) const;
    virtual T GetReflectionFactor(int s) const;
    virtual T GetReflectionFactor(int s, int b) const;
    virtual T GetPhotolysisRate(int s) const;
    virtual T GetBackgroundConcentration(int s) const;
    virtual T GetPreviousBackgroundConcentration(int s) const;
    virtual T GetBackgroundNumberConcentration(int b) const;
    virtual T GetPreviousBackgroundNumberConcentration(int b) const;
    virtual T GetBackgroundConcentration(int s, int b) const;
    virtual T GetPreviousBackgroundConcentration(int s, int b) const;
    virtual T GetSpecificHumidity() const;
    // virtual T GetLiquidWaterContent() const;
    
    void SetMeteo(T temperature, T wind_angle, T wind,
                  string stability_class, T longitude, T latitude,
                  bool isday, bool rural, T boundary_height);
    void SetMeteo(T temperature, T wind_angle, T wind,
                  T friction_velocity, T convective_velocity,
                  T boundary_height, T lmo, T coriolis,
                  T longitude, T latitude, bool isday, bool rural);
    virtual void SetAdditionalMeteo(T attenuation, T pressure,
                                    T specific_humidity);
    virtual void SetAdditionalMeteo(T liquid_water_content);
    virtual void SetScavengingCoefficient(T lambda, int s);
    virtual void SetDepositionVelocity(T vd, int s);
    virtual void SetReflectionFactor(T alpha, int s);
    virtual void SetReflectionFactor(T alpha, int s, int b);
    virtual void SetPhotolysisRate(T rate, int s);
    virtual void SetBackgroundConcentration(T concentration, int s);
    virtual void SetPreviousBackgroundConcentration(T concentration, int s);
    virtual void SetBackgroundNumberConcentration(T concentration, int b);
    virtual void SetPreviousBackgroundNumberConcentration(T concentration, int b);
    virtual void SetBackgroundConcentration(T concentration, int s, int b);
    virtual void SetPreviousBackgroundConcentration(T concentration, int s, int b);

    virtual void SetQuantity(T quantity, int s);
    virtual void SetNumberQuantity(T quantity, int b);
    virtual void SetQuantity(T quantity, int s, int b);
    void SetSigma_x(T sigma_x);
    void SetSigma_y(T sigma_y);
    void SetSigma_z(T sigma_z);
    void SetInitialSigma_x(T sigma_x_2);
    void SetInitialSigma_y(T sigma_y_2);
    void SetInitialSigma_z(T sigma_z_2);
    void SetPuffPosition(T x, T y, T z, T distance, T time);
    bool IsVolumeSource() const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PUFFTRANSPORT_HXX
#endif
