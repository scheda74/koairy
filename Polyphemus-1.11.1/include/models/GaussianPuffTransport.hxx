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

#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFTRANSPORT_HXX


#include "PuffTransport.cxx"
#include <list>
#include "BaseModel.cxx"

#include <list>


namespace Polyphemus
{


  using namespace std;


  //////////////////
  // GAUSSIAN PUFF//
  //////////////////


  //! Gaussian puff model.
  template<class T>
  class GaussianPuffTransport: public BaseModel<T>
  {

  protected:

    //! True if land category is set to "rural".
    bool option_rural;
    //! Is it daytime?.
    bool option_day;

    /********* Standard deviations *********/

    //! True if parameterization to compute standard deviations is "Briggs".
    bool option_briggs;
    //! True if parameterization to compute standard deviations is "Doury".
    bool option_doury;
    /*! \brief True if parameterization for sigma_z above BL is computed with
      Gillani. */
    bool option_gillani;
    //! True if HPDM formulae are used (only useful with similarity theory).
    bool option_hpdm;

    /********* Loss processes *********/

    //! List of species half-lives (s).
    Array<T, 1> half_life_time;

    //! List of species biological half-lives (s).
    Array<T, 1> biological_half_life_time;

    //! List of species scavenging coefficients ( s^(-1) ).
    Array<T, 1> scavenging_coefficient;

    //! List of species deposition velocities (m/s).
    Array<T, 1> deposition_velocity;

    //! Deposition model.
    string deposition_model;

    /********* Plume rise *********/

    //! True if Holland formula for plume rise is used.
    bool option_holland;
    //! True if Concawe formula for plume rise is used.
    bool option_concawe;
    bool compute_trajectory;

    /********* Puff list management *********/

    //! Pointer to the current puff.
    typename list<Puff<T>* >::iterator current_puff;

    //! Index of the current puff.
    int puff_index;

    //! Puff list.
    list<Puff<T>* > PuffList;
    //! Number of puffs.
    int Npuff;
    //! Time step between two puffs.
    T Delta_t_puff;

    /********* Background Meteorological data ************/

    //! Temperature of ambient air (Celsius degrees).
    T background_temperature_;
    //! Angle between positive x-axis and wind, counterclockwise.
    T background_wind_angle_;
    //! Wind (m/s).
    T background_wind_;

    //! Friction velocity (m/s).
    T background_friction_velocity_;
    //! Convective velocity (m/s).
    T background_convective_velocity_;
    //! Boundary layer height (m).
    T background_boundary_height_;
    //! Monin Obukhov length (m).
    T background_lmo_;
    //! Coriolis parameter.
    T background_coriolis_;

    //! Cosine of wind angle.
    T background_cos_angle_;
    //! Sine of wind angle.
    T background_sin_angle_;

    //! Stability class (A, B, C, D, E or F).
    string background_stability_class_;
    //! Stability (corresponding integer: from 0 for A to 5 for F).
    int background_stability_;

    /********* Current Meteorological data ************/

    //! Temperature of ambient air (Celsius degrees).
    T temperature_;
    //! Wind (m/s).
    T wind_;
    //! Inversion height (m).
    T inversion_height_;

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

    //! Night or day.
    bool isday_;

    //! Rural or urban.
    bool rural_;
    //! longitude;
    T longitude_;
    //! latitude;
    T latitude_;

  public:

    /*** Constructors and destructor ***/

    GaussianPuffTransport();
    GaussianPuffTransport(string config_file);
    virtual ~GaussianPuffTransport();

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    void Init();
    void InitStep();
    void InitPuffSource();
    virtual void InitSource(string file_puff, T delta_t_puff);
    void GetSourceCoordinates(Array<T, 2>& source_coordinates);
    void SetSourceCoordinates(Array<T, 2> source_coordinates);
    virtual void InitPuff();
    void InitPosition();

    /*** Access Methods ***/

    void SetDateMin(Date date_min);
    void SetTimeStep(T delta_t);
    void SetNt(int Nt);
    bool WithSimilarity();
    bool WithDeposition();
    bool WithScavenging();

    /*** Puff list management ***/

    int GetPuffNumber() const;
    void SetCurrentPuff(int index);
    void ErasePuff(int index);
    void ClearPuffList();
    string GetPuffSourceId(int index);

    /*** Background Meteorological data ***/

    void InitMeteo(ConfigStream& meteo, bool show_meteo);
    void InitScavenging(ConfigStream& meteo);
    void InitDeposition(ConfigStream& meteo);
    void SetCurrentMeteo(Puff<T>* puff);

    /*** Access methods for Puff attributes ***/

    void SetPuffMeteo(int index, T temperature, T wind_angle, T wind,
                      string stability_class,
                      T longitude, T latitude, bool isday, bool rural,
                      T boundary_height);
    void SetPuffMeteo(int index, T temperature, T wind_angle, T wind,
                      T friction_velocity, T convective_velocity,
                      T boundary_height, T lmo,
                      T coriolis, T longitude,
                      T latitude, bool isday, bool rural);
    void SetPuffDepositionVelocity(int index,
                                   Array<T, 1> deposition_velocity);
    void SetPuffScavengingCoefficient(int index,
                                      Array<T, 1> scavenging_coefficient);

    void GetPuffPosition(int index, T& x, T& y, T& z,
                         T& distance, T& time);
    void GetCurrentPuffTime(int index, T& time);
    void GetPuffSigma(int index, T& sigma_x, T& sigma_y, T& sigma_z);
    T GetPuffReleaseTime(int index);
    T GetPuffQuantity(int index, int s);

    /*** Computational Methods ***/

    void ComputePlumeRise();
    void ComputeEvolutivePlumeRise();
    void ComputePuffEffectiveHeight(Puff<T>* puff);
    void ComputePuffEvolutiveHeight(Puff<T>* puff);
    T GetConcentration(int species, T z, T y, T x);
    T GetIntegratedConcentration(int species, T z, T y, T x,
                                 T lz, T ly, T lx);
    using BaseModel<T>::GetConcentration;

    void Advection();
    void Diffusion();
    void ComputeLossFactor();
    void Forward();

    void Advection(int index);
    void Advection(Puff<T>* puff);
    void Diffusion(int index);
    void Diffusion(Puff<T>* puff);
    T ComputePuffConcentration(int index, int s, T x, T y, T z);
    T ComputePuffConcentration(Puff<T>* puff, int s, T x, T y, T z);
    T ComputePuffIntegral(int index, int s, T x, T y, T z,
                          T lx, T ly, T lz);
    T ComputePuffIntegral(Puff<T>* puff, int s, T x, T y, T z,
                          T lz, T ly, T lx);
    T ComputeVerticalTrajectory();
    void ComputeSigma(T distance, T transfer_time,
                      T z, T& sigma_x, T& sigma_y, T& sigma_z, bool diff = 0);
    void ComputeADMSSigma(T distance, T transfer_time, T plume_rise,
			  T source_diameter, T source_height, T z,
			  T& sigma_x, T& sigma_y, T& sigma_z);
    void ComputeTime(T z, T sigma_y, T& distance, T& transfer_time);

    void ComputeLossFactor(int index);
    void ComputeLossFactor(Puff<T>* puff);
    virtual void ComputeLossFactor(Puff<T>* puff, int s);
    void ComputeLossFactor(T distance, T transfer_time,
                           T z, T rad, T bio, T scav, T dep,
                           T& loss_factor, T& overcamp_factor);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFTRANSPORT_HXX
#endif
