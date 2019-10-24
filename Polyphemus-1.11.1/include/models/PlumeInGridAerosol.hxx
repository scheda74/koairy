// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDAEROSOL_HXX

#include "PlumeInGridChemistry.cxx"

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
#include <mpi.h>
#endif

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
  class PlumeInGridAerosol:
    public PlumeInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  {

  protected:

    /*** Meteorological data for chemistry ***/

    Data<T, 3> LiquidWaterContent_i;

    Data<T, 4> delta_number_puff;

  public:

    /*** Domain ***/

    //! 5D grid along z.
    RegularGrid<T> GridZ5D;
    //! 5D grid along y.
    RegularGrid<T> GridY5D;
    //! 5D grid along x.
    RegularGrid<T> GridX5D;

    //! 4D grid along z.
    RegularGrid<T> GridZ4D;
    //! 4D grid along y.
    RegularGrid<T> GridY4D;
    //! 4D grid along x.
    RegularGrid<T> GridX4D;

    //! Coordinates of boundary conditions along x.
    RegularGrid<T> GridX5D_interf_bc;
    //! Coordinates of boundary conditions along y.
    RegularGrid<T> GridY5D_interf_bc;

    /*** Species ***/

    vector<string> bin_list;

    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;

    //! 4D grid for aerosol bins.
    RegularGrid<T> GridB4D_aer;
    //! 5D grid for aerosol bins.
    RegularGrid<T> GridB5D_aer;

    //! 4D grid for aerosol species.
    RegularGrid<T> GridS4D_aer;
    //! 5D grid for aerosol species.
    RegularGrid<T> GridS5D_aer;

    //! 4D grid for gas species.
    RegularGrid<T> GridS4D;

    //! 4D grid for number species.
    RegularGrid<T> GridS4D_number;

    //! dry deposition fluxes at current date.
    Data<T, 3> DryDepositionFluxNumber_aer;
    //! Wet deposition fluxes at current date.
    Data<T, 3> WetDepositionFluxNumber_aer;
    //! In cloud wet deposition fluxes at current date.
    Data<T, 3> InCloudWetDepositionFluxNumber_aer;

    //! dry deposition fluxes at current date.
    Data<T, 4> DryDepositionFlux_aer;
    //! Wet deposition fluxes at current date.
    Data<T, 4> WetDepositionFlux_aer;
    //! In cloud wet deposition fluxes at current date.
    Data<T, 4> InCloudWetDepositionFlux_aer;
    //! List of aerosol bins with deposition velocities.
    vector<int> bin_list_dep_aer;
    //! List of aerosol bins with scavenging.
    vector<int> bin_list_scav_aer;


    /*** Constructor and destructor ***/

    PlumeInGridAerosol(string config_file);
    virtual ~PlumeInGridAerosol();

    /*** General methods ***/

    void Allocate();
    void ReadConfiguration();
    void CheckConfiguration();
    void Init();
    void InitStep();
    void Forward();
    void Forward_Gaussian();
    void Forward_Euler();
    void AddGaussianConcentrationsDomain();

    /*** Meteorological and species data ***/

    virtual bool HasNumberConcentration_aer();         
    void UpdateMeteo(int puff_index);         
    void ExtractMeteo(T height, T lat, T lon, bool isday,
                      bool option_similarity, map<string, T>& met_data,
                      string& stability, bool& rural);
    void ExtractSpeciesData(T height, T lat, T lon,
			    Array<T, 1>& deposition_velocity,
			    Array<T, 1>& deposition_velocity_aer,
			    Array<T, 1>& scavenging_coefficient,
			    Array<T, 1>& scavenging_coefficient_aer,
			    Array<T, 1>& photolysis_rate,
			    Array<T, 1>& background_concentration,
			    Array<T, 1>& background_concentration_number,
                            Array<T, 2>& background_concentration_aer);

    /*** Puff transfer ***/

    void PuffTransfer_aer(T quantity, T sigmaz, int s, T z,
                          T lat, T lon, bool isday,
                          Data<T, 5>& Concentration_out_aer, int b);
 
   void PuffTransfer_number(T quantity, T sigmaz, T z,
			    T lat, T lon, bool isday,
			    Data<T, 4>& Concentration_out_number, int b);

    void PuffIntegratedTransfer_aer(int puff_index,
				    Data<T, 5>& Concentration_out_aer,
				    Data<T, 4>& Concentration_out);

    void PuffIntegratedTransfer_number(int puff_index,
				       Data<T, 5> Concentration_out_aer,
				       Data<T, 4>& Concentration_out_number);

    void ComputeChemistry(Array<int, 2> PuffCellList,
                          Array<T, 4>& PuffConcentration,
                          Array<T, 5>& PuffConcentration_aer);

  };

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDAEROSOL_HXX
#endif
