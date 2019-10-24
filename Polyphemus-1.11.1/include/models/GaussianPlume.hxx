// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Irène Korsakissok,
// Régis Briant
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

// This file is part of a Gaussian plume model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_HXX


#include "PlumeSource.cxx"
#include "BaseModel.cxx"
#include "PlumeLineSource.cxx"

#include <list>


namespace Polyphemus
{


  using namespace std;


  ////////////////////
  // GAUSSIAN PLUME //
  ////////////////////


  //! Plume Gaussian model.
  template<class T>
  class GaussianPlume: public BaseModel<T>
  {

  protected:

    typedef vector<PlumeSource<T>* > source_list;
    typedef typename source_list::const_iterator source_iterator;

    static const T pi;
    static const T log_2;
    static const T cos_89;
    static const T sin_89;

    //! Read all input data, whatever the parameterizations may be?
    bool read_all;

    //! True if parameterization to compute standard deviations is "Briggs".
    bool option_briggs;

    //! True if parameterization to compute standard deviations is "Doury".
    bool option_doury;

    /*! \brief True if parameterization for sigma_z above BL is computed with
      Gillani. */
    bool option_gillani;

    //! True if HPDM formulae are used (only useful with similarity theory).
    bool option_hpdm;

    //! True if Holland formula for plume rise is used.
    bool option_holland;

    //! True if Concawe formula for plume rise is used.
    bool option_concawe;

    //! True if land category is set to "rural".
    bool option_rural;

    //! Is it daytime?.
    bool option_day;

    //! Is there radioactive decay?
    bool option_radioactive_decay;

    //! Is there biological decay?
    bool option_biological_decay;

    //! Is there scavenging?
    bool option_scavenging;

    //! Is there dry deposition?
    bool option_dry_deposition;

    //! Is there plume rise?
    bool option_plume_rise;

    //! Is there plume rise breakup?
    bool option_breakup;

    /*! \brief If true, increases the precision of the computation of the
      line source width along with the computational time.
    */
    bool option_high_width_precision;

    //! True if the NO2 chemistry should be computed.
    bool option_NO2_chemistry;

    //! True if the OH chemistry should be computed.
    bool option_OH_chemistry;

    //! Kinetic constant of the reaction of NO with O3 (mass / seconds).
    T k1_;

    //! NO2 photolysis rate (seconds^(-1)).
    Array<T, 2> k2_domain;
    Array<T, 1> k2_list;

    //! File storing all source data.
    string file_source;

    //! Source list.
    source_list SourceList;

    //! Point source list.
    source_list PointSourceList;

    //! Number of sources.
    int Nsource;

    //! Temporary number of sources.
    int Nsource_temp;

    //! List of species half-lives (s).
    Array<T, 1> half_life_time;

    //! List of species biological half-lives (s).
    Array<T, 1> biological_half_life_time;

    //! List of species scavenging coefficients ( s^(-1) ).
    Array<T, 1> scavenging_coefficient;

    //! List of species deposition velocities (m/s).
    Array<T, 1> deposition_velocity;
    //! Deposition model ("Chamberlain" or "Overcamp").
    string deposition_model;
    //! Number of points to compute the Chamberlain integral.
    int Nchamberlain;

    //! Temperature of ambient air (Celsius degrees).
    T temperature_;
    //! Angle between positive x-axis and wind, counterclockwise (degrees).
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

    //! Cosinus of wind angle.
    T cos_angle_;
    //! Sinus of wind angle.
    T sin_angle_;

    //! Stability class (A, B, C, D, E or F).
    string stability_class_;

    //! Stability (corresponding integer: from 0 for A to 5 for F).
    int stability_;

    // Correction coefficients for the line source model.
    vector<T> max_1_, mu_1_, sigma_1_,  max_2_, mu_2_, sigma_2_, mu_, a_, b_,
      c_, d_, beta_, lambda_;

    // Correction coefficients path file.
    string correction_coefficient_path_;

    // Threshold for the extremities correction in degrees.
    T threshold_extremities_correction_;

    // Threshold for the downsind correction in degrees.
    int threshold_downwind_correction_1_, threshold_downwind_correction_2_;

    // Is concentrations needs to be computed on the entire domain.
    bool option_compute_domain;

    // Is concentrations needs to be computed on a coordinate list.
    bool option_compute_list;

    // Coordinate list where concentration are computed.
    Array<T, 2> coordinate_list_;

    // Concentration list corresponding to the coordinate_list;
    Data<T, 2> Concentration_list_;

    //! Cartesian or longitude/latitude.
    bool option_cartesian;

    //! True if plumes length are infinite.
    bool option_infinite_plume;

    //! Index of source to be computed. -1 means all sources.
    int computed_source_;

    //! Pointer to the current plume.
    PlumeSource<T>* current_plume;

    //! Pointer to the previously saved plume.
    PlumeSource<T>* previous_plume;

    // Maximum number of point sources per line source for the discretization.
    int Np_max;

    // Discretization step for line sources parallel wind cases.
    int discretization_step;

    //! Pressure.
    T Pressure_;

    //! Specific humidity.
    T RelativeHumidity_;

    // Meteorological parameters for Plume In grid Chemistry.
    Array<T, 3> SpecificHumidity_grid, Temperature_grid,
      Pressure_grid, Wind_grid, Attenuation_grid;
    Array<int, 3> Stability_grid, Rural_grid;

    // Meteorological parameters for Plume In list Chemistry.
    Array<T, 1> SpecificHumidity_list, Temperature_list,
      Pressure_list, Wind_list, Attenuation_list;
    Array<int, 1> Stability_list, Rural_list;

  public:

    /*** Main constructor ***/

    GaussianPlume(string config_file);

    /*** Main destructor ***/

    virtual ~GaussianPlume();

    /*** Access methods ***/

    bool IsOption(string option);
    void SetOption(string option, bool value);
    string GetStringOption(string option);
    void SetStringOption(string option, string value);
    T GetModelParameter(string name);
    void SetModelParameter(string name, T value);
    void GetSourceData(string name, Array<T, 1>& value);
    void SetSourceData(string name, const Array<T, 1>& value);
    void GetSpeciesData(string name, Array<T, 1>& value);
    void SetSpeciesData(string name, const Array<T, 1>& value);
    virtual bool WithSimilarity();
    virtual bool WithDeposition();
    virtual bool WithScavenging();
    virtual bool WithChemistry();
    virtual Data<T, 4>&  GetConcentration();
    virtual void SetCurrentPlume(int index);
    virtual void SetGaussianMeteo(int index, T temperature, T wind_angle,
                                  T wind, string stability_class, T longitude,
                                  T latitude, bool isday, bool rural,
                                  T boundary_height);
    virtual void SetGaussianMeteo(int index, T temperature, T wind_angle,
                                  T wind, T friction_velocity,
                                  T convective_velocity, T boundary_height,
                                  T lmo, T coriolis, T longitude, T latitude,
                                  bool isday, bool rural);
    virtual void
    SetGaussianDepositionVelocity(int index,
                                  const Array<T, 1>& deposition_velocity);
    virtual void
    SetGaussianScavengingCoefficient(int index,
                                     const Array<T, 1>& scavenging_coefficient);
    virtual void GetGaussianPosition(int index, T& x, T& y, T& z,
                                     T& distance, T& time);
    virtual void GetGaussianSigma(int index, T& sigma_x,
                                  T& sigma_y, T& sigma_z);
    virtual string GetEmissionType(int index);
    T GetSourceParameter(string name, int index);
    virtual int GetSpecies(int index);
    virtual void SetDateMin(Date date_min);
    virtual void SetTimeStep(T delta_t);
    int GetInputSourceCount() const;
    int GetGaussianSourceCount() const;
    virtual T GetTimeStep();

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    virtual void InitSource(string source_file);
    virtual void Init(bool read_all_input = false);
    virtual void InitMeteo(ConfigStream& meteo, bool show_meteo);
    virtual void InitScavenging(ConfigStream& meteo);
    virtual void InitDeposition(ConfigStream& meteo);

    /*** Computations ***/

    virtual void InitCompute();
    virtual void InitStep();
    virtual void ComputePlumeRise();
    virtual void ComputeSourcePlumeRise(PlumeSource<T>& source);
    virtual T ComputeGaussianConcentration(int species, T z, T y, T x);
    using BaseModel<T>::GetConcentration;
    virtual T GetConcentration(int species, T z, T y, T x);
    virtual T GetIntegratedConcentration(int species, T z, T y, T x,
                                         T lz, T ly, T lx);

    virtual void Compute();
    virtual void ComputeSigma(T distance, T transfer_time,
                              T z, T& sigma_y, T& sigma_z, bool diff = 0);
    void ComputeLossFactor(T distance, T z, int species,
                           T& loss_factor, T& overcamp_factor);
    void ComputeLossFactor(T distance, T transfer_time,
                           T z, T rad, T bio, T scav, T dep,
                           T& loss_factor, T& overcamp_factor);
    virtual void RestrictComputationToSource(int source_index);
    virtual void RestoreRate();
    virtual void Discretization();
    virtual T ComputeCorrection(T cos_theta, T sin_theta,
                                T tan_theta, T y1_source, T y2_source,
                                T y_source, T  x_source, T distance_x,
                                T effective_height, T sigma_y);
    virtual bool ReceptorDownwindPlume(T x_source, T y_source, T y1_source,
                                       T y2_source, T cos_theta, T sin_theta);
    virtual T ComputeExtremitiesCorrection(T rate, T loss_factor, T wind_,
                                           T inversion_height_,
                                           T effective_height,
                                           T distance_x, T distance_y,
                                           T z, T initial_sigma_z_2,
                                           T overcamp_factor, T minimum_volume,
                                           T cos_theta, int species);
    virtual vector<T> ReadCoefficients(string section, int length);
    virtual void InitCorrectionCoefficients();

    void CheckSourceIndex(int source_index);
    void CartesianToLatLon(T x, T y, T& lon, T& lat);
    void LatLonToCartesian(T lon, T lat, T& x, T& y);

    // Plume In Grid methods
    virtual void InitSource();
    virtual void GetSourceCoordinates(Array<T, 2>& source_coordinates);
    virtual void SetSourceCoordinates(const Array<T, 2>& source_coordinates);
    virtual  void AddLineSource(int index_, T x_1, T y_1, T z_1, T x_2, T y_2,
                                T z_2, T rate_, T width_, T VehicleVelocity_,
                                T Area_, T Density_, int s_);
    virtual void EraseDiscretizedSource();
    virtual void EraseSource(vector<int> index_list);
    virtual T ComputeSigmaZ(int index_, T x, T y);


    // Plume In Grid methods
    virtual void Chemistry_Plume(const Date& current_date,
                                 const Array<T, 2>& latitude,
                                 const Array<T, 2>& longitude,
                                 Data<T, 4> BackgroundConcentration,
                                 Data<T, 4>& Concentration_out);
    virtual void Chemistry_Plume_list(const Date& current_date,
                                      Data<T, 2> BackgroundConcentration);
    virtual Array<T, 2> GetCoordinateList();
    virtual void SetZeroConcentrationList();
    virtual void Compute_list(int index_);
    virtual void FillPointSourceList(int begin, int end);
    virtual void FillPointSourceList(int index);
    virtual void EmptyPointSourceList();
    virtual void SaveCurrentPlume();
    virtual void RestoreCurrentPlume();
    virtual Array<T, 2>& ConcentrationList();
    virtual const Array<T, 2>& ConcentrationList() const;
    virtual void SplitSourceList(int rank, int Nproc);
    virtual void MultiplySourcesRate(T coefficient);
    virtual void ComputeOHChemistry(Data<T, 4>& concentration,
                                    const Array<T, 2>* latitude_matrix = 0,
                                    const Array<T, 2>* longitude_matrix = 0);
    virtual void ComputeOHChemistry(Data<T, 2>& concentration,
                                    bool option_plume_in_grid = false);
    virtual T ComputeO3PhotolysisRate(T zenith_angle, T Attenuation = 1.);
    virtual void
    SetChemistryMeteorologicalParameters(Array<T, 3> SpecificHumidity_grid_,
                                         Array<T, 3> Temperature_grid_,
                                         Array<T, 3> Pressure_grid_,
                                         Array<int, 3> Stability_grid_,
                                         Array<int, 3> Rural_grid_,
                                         Array<T, 3> Wind_grid_,
                                         Array<T, 3> Attenuation_grid_);
    virtual void
    SetChemistryMeteorologicalParameters(Array<T, 1> SpecificHumidity_list_,
                                         Array<T, 1> Temperature_list_,
                                         Array<T, 1> Pressure_list_,
                                         Array<int, 1> Stability_list_,
                                         Array<int, 1> Rural_list_,
                                         Array<T, 1> Wind_list_,
                                         Array<T, 1> Attenuation_list_);

    /*** Empty methods required for Plume-in-grid model ***/

    void InitPosition() {};
    void ComputeLossFactor() {};
    void ErasePuff(int index) {};
    void InitPhotolysis(int Nr, vector<string> photolysis_list) {};
    void SetPuffAdditionalMeteo(int index, T attenuation, T pressure,
                                T specific_humidity) {};
    void SetGaussianPhotolysisRate(int index, Array<T, 1> photolysis_rate) {};
    void SetGaussianBackgroundConcentration(int index,
                                            Array<T, 1> concentration) {};
    void Advection() {};
    void Diffusion() {};
    void ComputeChemistry() {};
    void ComputeChemistry(const vector<list<int> >& PuffCellList,
                          const vector<T>& PuffCellVolume,
                          Array<vector<T>, 1 >& PuffCellConcentration) {};
    T ComputePuffIntegral(int index, int s, T x, T y, T z,
                          T lx, T ly, T lz)
    {
      return 0.;
    };
    T GetPuffReleaseTime(int index)
    {
      return 0.;
    };
    T GetPuffQuantity(int index, int s)
    {
      return 0.;
    };
    // end of empty methods


  protected:

    virtual void ComputeNO2Chemistry(Data<T, 4>& concentration);
    virtual void ComputeNO2Chemistry(Data<T, 2>& concentration);
    virtual T ComputeNO2PhotolysisRate(T zenith_angle, T Attenuation = 1.);
    virtual T ComputeZenithAngle(int nb_sec, T latitude, T longitude);
    virtual T LineSourceConcentration(int species, T z, T y, T x,
                                      PlumeSource<T>& iter, T step = 0.);
    virtual T RombergLineSourceConcentration(T a, T b, int species,
                                             T z, T y, T x, PlumeSource<T>& iter);


  private:
    bool IsPlume(T distance_from_source);
    template<int Dim>
    void NullifyNegativeConcentration(Array<T, Dim>& concentration,
                                      const string& species_name) const;
    template<int Dim>
    void NullifyNegativeConcentration(Array<T, Dim>& concentration) const;
    void DeleteSource();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_HXX
#endif
