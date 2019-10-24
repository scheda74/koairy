// Copyright (C) 2006-2007,  ENPC - INRIA - EDF R&D
//     Author(s): Edouard Debry
//
// This file is part of the Size Resolved Aerosol Model (SIREAM), a component
// of the air quality modeling system Polyphemus.
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


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SORGAM_HXX

extern "C"
{
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/ipc.h>
}

#include <vector>
#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"

#if !defined(POLYPHEMUS_PARALLEL_WITH_OPENMP)	\
  && !defined(POLYPHEMUS_PARALLEL_WITH_MPI)
// Global variables.
#include "data_shm.h"
#endif

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////

#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _compute_coagulation_coefficient compute_coagulation_coefficient__
#else
#define _compute_coagulation_coefficient compute_coagulation_coefficient_
#endif

#define _chem chem_
#define _aerosol aerosol_
#define _dimensions dimensions_

  extern "C"
  {
    void _chem(int*, int*, int*, int*, int*, int*,
               int*, int*, int*, int*, int*, int*, int*,
               double*, double*, int*,
               double*, double*, double*, double*, double*, double*,
               double*, double*, double*, double*, double*,
               double*, double*, double*, double*, double*, int*,
               double*, double*, double*, int*, int*, int*, int*,
               double*, double*, double*, double*, int*, double*,
               int*, double*, double*, int*);

    void _aerosol(int*, int*, int*, int*,
                  int*, int*, int*, int*, double*, double*, double*,
                  double*, double*, double*, double*, int*, int*,
                  int*, int*, int*, double*, double*, double*, double*,
                  double*, int*, int*, int*, double*, double*, double*,
                  double*, double*, double*);

    void _compute_coagulation_coefficient(int*, double*,
                                          int*, int*, int*, double*);

    void _dimensions(int*, int*, int*);

  }


  ///////////////////////////
  // AEROSOL_SIREAM_SORGAM //
  ///////////////////////////


  //! This class is a numerical solver for the chemical mechanisms
  //  RACM_SIREAM_SORGAM, CB05_SIREAM_SORGAM and RACM2_SIREAM_SORGAM.
  /*! It uses a second-order Rosenbrock method.
   */
  template<class T>
  class Aerosol_SIREAM_SORGAM: public BaseModuleParallel
  {

  protected:

    //! Number of species.
    int Ns;
    //! Number of reactions.
    int Nr;
    //! Number of photolysis reactions.
    int Nr_photolysis;
    //! Number of species with volume sources.
    int Ns_source;
    //! Number of levels with volume sources.
    int Nz_source;

    //! Number of sub-cycles.
    int Ncycle;

    //! Sorted list of species names.
    vector<string> species_list;

    //! Molecular weights of species.
    Array<T, 1> molecular_weight;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 1> ConversionFactor;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 2> ConversionFactorJacobian;

    /*! \brief Map between photolysis reactions names and their indices in
      photolysis reactions.
    */
    map<string, int> photolysis_reaction_name;
    //! Indices of photolysis reactions among other reactions.
    Array<int, 1> photolysis_reaction_index;
    //! Indices of species with volume sources.
    Array<int, 1> source_index;

    //! Number of sub-cycles ofr aerosol processes
    int Ncycle_aer;

    //! Number of aerosol species.
    int Ns_aer;

    //! Number of aerosol bins.
    int Nbin_aer;

    //! Bin index corresponding to fixed cutting diameter.
    int cutting_bin;
    /* bins from 1 to cutting_bin included are at equilibrium
       and are dynamic from cutting_bin+1 to Nbin_aer */

    //! Indices of species with heterogeneous reactions.
    Array<int, 1> heterogeneous_reaction_index;

    //! list of aerosol options.
    int Noptions_aer;
    Array<int, 1> options_aer;

#if !defined(POLYPHEMUS_PARALLEL_WITH_OPENMP)	\
  && !defined(POLYPHEMUS_PARALLEL_WITH_MPI)
    //! Shared memory object.
    SharedMemory shared_memory;
#endif

  public:
    /*** Configuration ***/

    map<string, bool> option_process_aer;

    //! Numerical solver for dynamic bin condensation (etr, ros2 or ebi).
    string dynamic_condensation_solver;

    //! Cutting diameter between equilibrium and dynamic bins.
    double fixed_cutting_diameter;

    //! Sulfate condensation computation method (equilibrium, dynamic).
    string sulfate_computation;

    /*! Redistribution method of lagrangian bins (number-conserving,
      interpolation).
    */
    string redistribution_method;

    //! If nucleation, which nucleation model (binary, ternary).
    string nucleation_model;

    //! Aerosol wet diameter estimation (gerber or isoropia).
    string wet_diameter_estimation;

    //! Name of the aqueous module.
    string aqueous_module;

    //! Thermodynamics module for bulk equilibrium (eqsam or isorropia).
    string thermodynamics_module;

    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;

    //! Aerosol density (kg / m^3).
    double FixedDensity_aer;
    Array<T, 1> Density_aer;

    //! Coagulation coefficients.
    Array<int, 1> CoagulationCouple;
    Array<int, 2> CoagulationFirstIndex;
    Array<int, 2> CoagulationSecondIndex;
    Array<T, 3> CoagulationPartitionCoefficient;

    //! liquid water content threshold for clouds.
    double lwc_cloud_threshold;

    //! Options for adaptive time stepping for gaseous chemistry
    bool with_adaptive_time_step;
    int option_adaptive_time_step;
    double adaptive_time_step_tolerance, min_adaptive_time_step;

    //! Option for photolysis
    int option_photolysis_tabulation;


    /*** Constructor ***/

    Aerosol_SIREAM_SORGAM();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    void DisplayConfiguration();

    template<class ClassModel>
    void Forward(ClassModel& Model);

    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    void Forward(T current_time,
                 Data<T, 3>& Attenuation_i,
                 Data<T, 3>& SpecificHumidity_i,
                 Data<T, 3>& Temperature_i,
                 Data<T, 3>& Pressure_i,
                 Data<T, 4>& VolumeEmission_i,
                 Data<T, 4>& PhotolysisRate_i,
                 T next_time,
                 Data<T, 3>& Attenuation_f,
                 Data<T, 3>& SpecificHumidity_f,
                 Data<T, 3>& Temperature_f,
                 Data<T, 3>& Pressure_f,
                 Data<T, 4>& VolumeEmission_f,
                 Data<T, 4>& PhotolysisRate_f,
                 Array<T, 1>& Longitude,
                 Array<T, 1>& Latitude,
                 Data<T, 4>& Concentration,
                 Data<T, 3>& LiquidWaterContent_i,
                 Data<T, 4>& WetDiameter_aer,
                 Data<T, 5>& Concentration_aer);

    void Forward_aer(T current_time,
                     Data<T, 3>& SpecificHumidity_i,
                     Data<T, 3>& Temperature_i,
                     Data<T, 3>& Pressure_i,
                     T next_time,
                     Data<T, 4>& Concentration,
                     Data<T, 3>& LiquidWaterContent_i,
                     Data<T, 2>& Rain_i,
                     Array<T, 1>& VerticalInterface,
                     Data<T, 5>& Concentration_aer,
                     Data<T, 3>& InCloudWetDepositionFlux,
                     Data<T, 4>& InCloudWetDepositionFlux_aer,
                     Data<T, 3>& pH);

    bool IsRequired(string field);
    bool IsComputed(string field);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SORGAM_HXX
#endif
