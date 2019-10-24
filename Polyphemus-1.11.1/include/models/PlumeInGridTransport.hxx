// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDTRANSPORT_HXX

#include "BaseModel.cxx"
#include <map>
#include <list>
#include <vector>
#include <string>
#include "AtmoData.hxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////////
  // PLUME-IN-GRID DRIVER //
  //////////////////////////


  /*! \brief This class provides a driver to perform a plume-in-grid
    simulation.
  */
  /*! The driver is responsible for the model initialization, for the loop
    over all meteorological conditions and over all time steps, and for the
    calls to the output saver. It also Its reference floating-point precision
    is 'T'. The model is an instance of 'ClassModel' and the output saver is
    an instance of 'ClassOutputSaver'.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  class PlumeInGridTransport: public BaseModel<T>
  {

  protected:

    /*** Main components ***/

    //! Underlying model.
    ClassEulerianModel Model;

    /*** Configuration ***/

    //! Configuration file for gaussian model.
    string gaussian_file;

    //! Cartesian or longitude/latitude.
    bool option_cartesian;

    //! Is there deposition?
    bool option_deposition;

    //! Is there scavenging?
    bool option_scavenging;

    // //! Is there chemistry?
    // bool option_chemistry;

    //! Vertical levels management.
    bool altitude_field;

    /*! \brief Major point sources emissions coordinates
      (longitude, latitude, height).*/

    //! Gaussian puff model.
    ClassLocalModel*  GaussianPuffModel;

    //! Eulerian time step.
    T delta_t_eulerian;

    //! Local time step, number of time steps and time step between puffs.
    T delta_t_local;
    int Nt_local;

    //! Coefficients to multiply standard deviations to obtain puff size.
    T coefficient_y;
    T coefficient_z;

    //! List of species with scavenging.
    vector<string> species_list_scav;
    vector<string> species_list_dep; 
    //! Dry deposition fluxes at current date.
    Data<T, 3> DryDepositionFlux;
  
    //! Wet deposition fluxes at current date.
    Data<T, 3> WetDepositionFlux;
    //! In cloud wet deposition fluxes at current date.
    Data<T, 3> InCloudWetDepositionFlux;
 
    //! Time criteria for reinjection.
    bool option_time;
    //! Reinjection time (in seconds since release).
    T reinjection_time;
    //! Injection method (column or integrated).
    bool injection_integrated;
    //! Merge overlaping puff;
    bool merge_puff;
    bool merge_source_id;
    //! Compute evolutive plume_rise;
    bool evolutive_plume_rise;
    //! Are meteorological data interpolated?
    bool interpolated;

    //! Save gaussian model concentration in the domain
    bool save_gaussian_domain;

    /*** Ground data ***/

    //! Land use coverage.
    Data<T, 3> LUC;

    //! Number of LUC classes.
    int Nc;

    //! LUC file.
    string LUC_file;

    //! Urban index.
    int urban_index;

    //! Option ground ("rural", "urban", or "LUC").
    string land_type;

    //! Proportion of urban LUC in a cell.
    T percent_urban;

    /*** Meteorological data ***/

    //! Zonal wind at current date.
    Data<T, 3> ZonalWind_i;
    //! Meridional wind at current date.
    Data<T, 3> MeridionalWind_i;
    //! Temperature at current date.
    Data<T, 3> Temperature_i;

    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_i;
    //! Zonal wind buffer.
    Data<T, 3> FileZonalWind_f;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_i;
    //! Meridional wind buffer.
    Data<T, 3> FileMeridionalWind_f;
    //! Temperature buffer.
    Data<T, 3> FileTemperature_i;
    //! Temperature buffer.
    Data<T, 3> FileTemperature_f;

    /*** Meteorological data for stability class ***/

    //! First level wind at current date.
    Data<T, 2> FirstLevelWind_i;
    //! Low cloudiness at current date.
    Data<T, 2> LowCloudiness_i;
    //! Medium cloudiness at current date.
    Data<T, 2> MediumCloudiness_i;
    //! High cloudiness at current date.
    Data<T, 2> HighCloudiness_i;
    //! Insolation at current date.
    Data<T, 2> Insolation_i;

    //! First level wind buffer.
    Data<T, 2> FileFirstLevelWind_i;
    //! First level wind buffer.
    Data<T, 2> FileFirstLevelWind_f;
    //! Low cloudiness buffer.
    Data<T, 2> FileLowCloudiness_i;
    //! Low cloudiness buffer.
    Data<T, 2> FileLowCloudiness_f;
    //! Medium cloudiness buffer.
    Data<T, 2> FileMediumCloudiness_i;
    //! Medium cloudiness buffer.
    Data<T, 2> FileMediumCloudiness_f;
    //! High cloudiness buffer.
    Data<T, 2> FileHighCloudiness_i;
    //! High cloudiness buffer.
    Data<T, 2> FileHighCloudiness_f;
    //! Insolation buffer.
    Data<T, 2> FileInsolation_i;
    //! Insolation buffer.
    Data<T, 2> FileInsolation_f;

    /*** Meteorological data for similarity theory ***/

    //! Friction velocity at current date.
    Data<T, 2> FrictionModule_i;
    //! Boundary height at current date.
    Data<T, 2> BoundaryHeight_i;
    //! Monin Obukhov length at current date.
    Data<T, 2> LMO_i;

    //! Friction velocity buffer.
    Data<T, 2> FileFrictionModule_i;
    //! Friction velocity buffer.
    Data<T, 2> FileFrictionModule_f;
    //! Boundary height buffer.
    Data<T, 2> FileBoundaryHeight_i;
    //! Boundary height buffer.
    Data<T, 2> FileBoundaryHeight_f;
    //! Monin Obukhov length buffer.
    Data<T, 2> FileLMO_i;
    //! Monin Obukhov length buffer.
    Data<T, 2> FileLMO_f;

    /*** Other data for loss processes ***/

    //! Deposition velocities.
    Data<T, 3> DepositionVelocity_i;

    //! Scavenging coefficients.
    Data<T, 4> ScavengingCoefficient_i;

    //! Gas compensation, needed for gas/part mass conservation
    Data<T, 4> GasCompensation;

  public:

    /*** Constructor and destructor ***/

    PlumeInGridTransport(string config_file);
    virtual ~PlumeInGridTransport();

    /*** General methods ***/

    void Allocate();
    void ReadConfiguration();
    void CheckConfiguration();
    void Init();
    void InitStep();
    void Forward();

    string GetModelName() const;
    /*** Getting Concentrations ***/

    T GetConcentration(int species, T z, T y, T x);
    T GetIntegratedConcentration(int species, T z, T y, T x,
                                 T lz, T ly, T lx);
    Data<T, 4>& GetConcentration();

    /*** Meteorological and species data ***/

    void UpdateMeteo(int puff_index);
    void ExtractMeteo(T height, T lat, T lon, bool isday,
                      bool option_similarity, map<string, T>& met_data,
                      string& stability, bool& rural);
    void ExtractSpeciesData(T height, T lat, T lon,
                            Array<T, 1>& deposition_velocity,
                            Array<T, 1>& scavenging_coefficient);

    /*** Puff transfer ***/

    void PuffTransfer(T quantity, T sigmaz, int s, T z,
                      T lat, T lon, bool isday,
                      Data<T, 4>& Concentration_out);
    void PuffIntegratedTransfer(int puff_index,
                                Data<T, 4>& Concentration_out);

    /*** Other methods ***/

    void ComputeCellWidth(int k, int j, int i, T& CellWidth_z,
                          T& CellWidth_y, T& CellWidth_x, T& CellVolume);
    void GetCellIndices(T lon, T lat, T height,
                        int& index_z, int& index_y, int& index_x);
    T GetGasCompensation(int species, int z, int y, int x);
    void SetGasCompensation(int species, int z, int y, int x, 
			    T gas_compensation);
    void LatLonToCartesian(T lon, T lat, T& x, T& y);
    void CartesianToLatLon(T x, T y, T& lon, T& lat);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDTRANSPORT_HXX
#endif
