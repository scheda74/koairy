// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Vivien Mallet, Régis Briant
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


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_HXX


#include "PuffChemistry.cxx"
#include <list>
#include "BaseModel.cxx"
#include "Photochemistry.cxx"


namespace Polyphemus
{


  using namespace std;


  //////////////////
  // GAUSSIAN PUFF//
  //////////////////


  //! Gaussian puff model.
  template<class T>
  class GaussianPuff: public BaseModel<T>
  {

  protected:

    //! True if land category is set to "rural".
    bool option_rural;
    //! Is it daytime?.
    bool option_day;

    /********* Standard deviations *********/

    //! Are standard deviations increasing in time?
    bool option_increasing_sigma;
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

    //! Is there radioactive decay?
    bool option_radioactive_decay;
    //! Is there biological decay?
    bool option_biological_decay;
    //! Is there scavenging?
    bool option_scavenging;
    //! Is there dry deposition?
    bool option_dry_deposition;

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

    //! Is there plume rise?
    bool option_plume_rise;
    //! Is there plume rise breakup?
    bool option_breakup;
    //! True if Holland formula for plume rise is used.
    bool option_holland;
    //! True if Concawe formula for plume rise is used.
    bool option_concawe;

    /********* Chemistry *********/

    //! Is there chemistry?
    bool option_chemistry;
    //! Is there chemical interaction between puffs?
    bool option_interaction;
    //! Number of subcycles for chemistry.
    int Ncycle;

    //! Are background concentrations saved along with puff concentrations?
    bool with_background_concentration;

    //! Is the total plume mass to be saved?
    bool option_mass;
    //! File to save the total plume mass.
    string file_mass;

    //! List of photolysis rates.
    Array<T, 1> photolysis_rate_;

    //! List of background concentrations.
    Array<T, 1> background_concentration_;

    //! List of species with photolysis reactions.
    vector<string> photolysis_reaction_list;

    //! Number of photolysis reactions.
    int Nr_photolysis;

    //! Chemical mechanism and chemical numerical scheme.
    Photochemistry<T> Chemistry_;

    /********* Puff list management *********/

    //! Pointer to the current puff.
    typename list<Puff<T>* >::iterator current_puff;

    //! Index of the current puff.
    int puff_index;

    //! Puff list.
    list<Puff<T>* > PuffList;
    //! Number of puffs.
    int Npuff;
    //! Number of sources.
    int Nsource;
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

    //! Attenuation.
    T background_attenuation_;
    //! Pressure.
    T background_pressure_;
    //! Specific humidity.
    T background_specific_humidity_;

    //! longitude;
    T background_longitude;
    //! latitude;
    T background_latitude;

    /********* Current Meteorological data ************/

    //! Temperature of ambient air (Celsius degrees).
    T temperature_;
    //! Angle between positive x-axis and wind, counterclockwise.
    T wind_angle_;
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

    //! Attenuation.
    T attenuation_;
    //! Pressure.
    T pressure_;
    //! Specific humidity.
    T specific_humidity_;

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

    GaussianPuff();
    GaussianPuff(string config_file);
    virtual ~GaussianPuff();

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    void Init();
    void InitStep();
    void InitSource();
    void InitSource(string file_puff, T delta_t_puff);
    void GetSourceCoordinates(Array<T, 2>& source_coordinates);
    void SetSourceCoordinates(Array<T, 2> source_coordinates);
    void InitPuff();
    void InitPosition();

    /*** Access Methods ***/

    int GetNr_photolysis() const
    {
      return Nr_photolysis;
    }
    int GetNs_source() const
    {
      return 0;
    }
    int GetNz_source() const
    {
      return 0;
    }
    int SourceGlobalIndex(int s) const
    {
      return 0;
    }
    vector<string> GetPhotolysisReactionList() const
    {
      return photolysis_reaction_list;
    }
    void SetDateMin(Date date_min);
    void SetTimeStep(T delta_t);
    void SetNt(int Nt);
    bool WithSimilarity();
    bool WithDeposition();
    bool WithScavenging();
    bool WithChemistry();

    /*** Puff list management ***/

    const list<Puff<T>* >& GetPuffList() const;
    int GetInputSourceCount() const;
    int GetGaussianSourceCount() const;
    void SetPuffList(list<Puff<T>* > pufflist);
    void SetCurrentPuff(int index);
    int GetPuffIndex(typename list<Puff<T>* >::iterator puff);
    void ErasePuff(int index);
    void ClearPuffList();

    /*** Background Meteorological data ***/

    void InitMeteo(ConfigStream& meteo, bool show_meteo);
    void InitPhotolysis(int Nr, vector<string> photolysis_list);
    void InitChemistry(ConfigStream& meteo);
    void InitScavenging(ConfigStream& meteo);
    void InitDeposition(ConfigStream& meteo);
    void SetBackgroundMeteo(T temperature, T wind_angle, T wind,
                            string stability_class, T longitude, T latitude,
                            T boundary_height);
    void SetBackgroundMeteo(T temperature, T wind_angle, T wind,
                            T friction_velocity,
                            T convective_velocity,
                            T boundary_height, T lmo,
                            T coriolis, T longitude, T latitude);
    void SetBackgroundAdditionalMeteo(T attenuation, T pressure,
                                      T specific_humidity);
    void SetBackgroundPhotolysisRate(Array<T, 1> photolysis);
    void SetBackgroundConcentration(Array<T, 1> concentration);
    void SetCurrentMeteo(Puff<T>* puff);
    bool IsOption(string option_name);
    T GetModelParameter(string name);
    void SetModelParameter(string name, T value);

    /*** Access methods for Puff attributes ***/

    void SetGaussianMeteo(int index, T temperature, T wind_angle, T wind,
                          string stability_class,
                          T longitude, T latitude, bool isday, bool rural,
                          T boundary_height);
    void SetGaussianMeteo(int index, T temperature, T wind_angle, T wind,
                          T friction_velocity, T convective_velocity,
                          T boundary_height, T lmo,
                          T coriolis, T longitude,
                          T latitude, bool isday, bool rural);
    void SetPuffAdditionalMeteo(int index, T attenuation, T pressure,
                                T specific_humidity);
    void SetGaussianPhotolysisRate(int index,
                                   const Array<T, 1>& photolysis_rate);
    void SetGaussianDepositionVelocity(int index,
                                       const Array<T, 1>& deposition_velocity);
    void SetGaussianScavengingCoefficient(int index,
                                          const Array<T, 1>& scavenging_coefficient);
    void SetGaussianBackgroundConcentration(int index,
                                            const Array<T, 1>& concentration);
    void GetGaussianPosition(int index, T& x, T& y, T& z,
                             T& distance, T& time);
    void GetGaussianSigma(int index, T& sigma_x, T& sigma_y, T& sigma_z);
    T GetPuffReleaseTime(int index);
    T GetPuffQuantity(int index, int s);

    /*** Computational Methods ***/

    void ComputePlumeRise();
    void ComputePuffEffectiveHeight(Puff<T>* puff);
    T GetConcentration(int species, T z, T y, T x);
    T GetIntegratedConcentration(int species, T z, T y, T x,
                                 T lz, T ly, T lx);
    using BaseModel<T>::GetConcentration;

    void Advection();
    void Diffusion();
    void ComputeLossFactor();
    void ComputeChemistry();
    void ComputeChemistry(const vector<list<int> >& CellPuff,
                          const vector<T>& CellVolume,
                          Array<vector<T>, 1 >& CellConcentration);
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
    void ComputeSigma(T distance, T transfer_time,
                      T z, T& sigma_x, T& sigma_y, T& sigma_z, bool diff = 0);
    void ComputeTime(T z, T sigma_y, T& distance, T& transfer_time);

    void ComputeLossFactor(int index);
    void ComputeLossFactor(Puff<T>* puff);
    void ComputeLossFactor(Puff<T>* puff, int s);
    void ComputeLossFactor(T distance, T transfer_time,
                           T z, T rad, T bio, T scav, T dep,
                           T& loss_factor, T& overcamp_factor);
    T ComputeInteractionCoefficient(Puff<T>* puff_alpha, Puff<T>* puff_beta);

    /*** Empty methods required for Continuous-plume-in-grid model ***/

    void Chemistry_Plume_list(Date current_date,
                              Data<T, 2> BackgroundConcentration) {};
    void Compute_list(int index_) {};
    void SetOption(string option, bool value) {};
    void SetZeroConcentrationList() {};
    void RestrictComputationToSource(int computation_index) {};
    void RestoreRate() {};
    void EraseDiscretizedSource() {};
    void EraseSource(const vector<int>& index_list) {};
    void AddLineSource(int index_, T x_1, T y_1, T z_1, T x_2, T y_2, T z_2,
                       T rate_, T width_, T VehicleVelocity, T Area,
                       T Density, int s_) {};
    void Chemistry_Plume(Date current_date, Array<T, 2> latitude,
                         Array<T, 2> longitude,
                         Data<T, 4> BackgroundConcentration,
                         Data<T, 4>& Concentration_out) {};
    int GetSpecies(int index)
    {
      return 0;
    };
    T GetSourceParameter(string name, int index)
    {
      return 0.;
    };
    T ComputeSigmaZ(int index_, T x, T y)
    {
      return 0.;
    };
    T ComputeGaussianConcentration(int species, T z, T y, T x)
    {
      return 0.;
    };
    string GetEmissionType(int index)
    {
      return "";
    };
    Array<T, 2> GetCoordinateList()
    {
      return Array<T, 2>();
    };
    void InitCorrectionCoefficients() {};
    void FillPointSourceList(int begin, int end) {};
    void FillPointSourceList(int index) {};
    void EmptyPointSourceList() {};
    void SplitSourceList(int rank, int Nproc) {};
    Array<T, 2>& ConcentrationList()
    {
      throw "'Array<T, 2>& GaussianPuff::ConcentrationList()' is undefined.";
    }
    const Array<T, 2>& ConcentrationList() const
    {
      throw "'const Array<T, 2>& GaussianPuff::ConcentrationList() const' "
        "is undefined.";
    }
    void MultiplySourcesRate(T coefficient) {};
    void SaveCurrentPlume() {};
    void RestoreCurrentPlume() {};
    T GetTimeStep()
    {
      return T();
    };
    void SetChemistryMeteorologicalParameters(Array<T, 3> Sp, Array<T, 3> Te,
                                              Array<T, 3> P, Array<int, 3> St,
                                              Array<int, 3> R, Array<T, 3> W,
                                              Array<T, 3> A) {};
    void SetChemistryMeteorologicalParameters(Array<T, 1> Sp, Array<T, 1> Te,
                                              Array<T, 1> P, Array<int, 1> St,
                                              Array<int, 1> R, Array<T, 1> W,
                                              Array<T, 1> A) {};
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_HXX
#endif
