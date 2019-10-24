// Copyright (C) 2016-2007, CEREA- ENPC - EDF R&D
// Author(s): Youngseob Kim
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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


#ifndef POLYPHEMUS_FILE_DRIVER_STREETINGRIDTRANSPORT_HXX

#include "BaseModel.cxx"
#include <map>
#include <list>
#include <vector>
#include <string>
#include "AtmoData.hxx"

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
#include <mpi.h>
#endif

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////////////
  // STREET-IN-GRID TRANSPORT MODEL //
  ////////////////////////////////////


  /*! \brief This class provides an interface to perform a street-in-grid
    simulation.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  class StreetInGridTransport: public BaseModel<T>
  {

  protected:

    static const T pi;

    /*** Main components ***/

    //! Underlying model.
    ClassEulerianModel Model;

    /*** Configuration ***/

    //! Configuration file for a street-network model.
    string street_config;

    //! Cartesian or longitude/latitude.
    bool option_cartesian;

    //! Is there deposition?
    bool option_deposition;

    //! Is there scavenging?
    bool option_scavenging;

    //! Vertical levels management.
    bool altitude_field;

    //! Are meteorological data interpolated?
    bool interpolated;

    //! Eulerian time step.
    T delta_t_eulerian;


    /*** Ground data ***/

    //! Land use coverage.
    Data<T, 3> LUC;

    //! Number of LUC classes.
    int Nc;

    //! LUC file.
    string LUC_file;

    //! Urban index.
    int urban_index;


    /*** Street-network model ***/

    ClassLocalModel*  StreetNetworkModel;

    //! Number of streets.
    int Nstreet;

    //! Number of intersections.
    int Nintersection;

    //! Local time step.
    T delta_t_local;


    /*** Volume correction ***/

    //! corrected volume or not.
    // bool corrected_background;
    
    //! If the grid cell is urban or not.
    Array<bool, 2> is_urban;

    //! Ratio of the urbanized surface. 
    Array<T, 2> urban_fraction;

    //! Length-weighted height of the urban canopy.
    Array<T, 2> weighted_height;

    //! Sum of product height and length for the streets in the urban canopy.
    Array<T, 2> sum_height_length;

    //! Sum of length of the streets in the urban canopy.
    Array<T, 2> sum_length;

    //! Number of the streets in a urban grid cell.
    Array<T, 2> street_number;

    //! Number of the streets in a urban grid cell.
    Array<T, 2> sum_street_volume;

    //! Mass of species in a urban grid cell.
    Array<T, 3> cell_quantity;

    //! Mass flux in a urban grid cell.
    Array<T, 3> mass_flux_grid;

    //! Mass flux in a urban grid cell.
    Array<T, 3> emission_rate;


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

    //! Concentration over canopy.
    Data<T, 3> ConcentrationOverCanopy;

    Data<T, 2> StreetConcentration;

  public:

    /*** Constructor and destructor ***/

    StreetInGridTransport(string config_file);
    virtual ~StreetInGridTransport();

    /*** General methods ***/

    void Allocate();
    void ReadConfiguration();
    void CheckConfiguration();
    void Init();
    void InitStep();
    void Forward();

    /*** Meteorological and species data ***/

    void UpdateBackgroundConcentration();
    void UpdateStreetData(int street_index);
    void UpdateIntersectionData(int intersection_index);
    void ExtractMeteo(T height, T lat, T lon, map<string, T>& met_data);
    void ExtractMeteo(T height, T lat, T lon, map<string, T>& met_data, 
                      bool check);
    void ExtractSpeciesData(T height, T lat, T lon,
			    Array<T, 1>& background_concentration);
    void UpdateMeteo();

    /*** Mass transfer from the streets ***/
 
    void StreetTransfer(Data<T, 4>& Concentration_out);
    void ComputeMassFlux();

    /*** Correctiion of background concentrations ***/

    void AddMassFluxToBackgroundConcentration();
    void ComputeUrbanFraction();
    void CorrectBackgroundConcentration();
    void RestoreBackgroundConcentration();

    /*** Other methods ***/

    void ComputeCellWidth(int k, int j, int i, T& CellWidth_z,
			  T& CellWidth_y, T& CellWidth_x, T& CellVolume);
    void GetCellIndices(T lon, T lat, T height,
			int& index_z, int& index_y, int& index_x);
    int GetNStreet() const;
    Data<T, 2>& GetStreetConcentration();

    void EmissionRateSave();

    /*** Save Concentration over canopy ***/
    
    void ComputeConcentrationOverCanopy();
    Data<T, 3>& GetConcentrationOverCanopy();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETINGRIDTRANSPORT_HXX
#endif
