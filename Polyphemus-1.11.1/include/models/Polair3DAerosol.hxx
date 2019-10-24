// Copyright (C) 2006-2018, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Shupeng Zhu
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

// This file is part of the Eulerian model Polair3D.


#ifndef POLYPHEMUS_FILE_MODELS_POLAIR3DAEROSOL_HXX


#include <vector>
#include <utility>
#include "AtmoData.hxx"
#include "Polair3DChemistry.cxx"


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
#define _compute_dry_deposition compute_dry_deposition__
#define _compute_scavenging_coefficient_aer	\
  compute_scavenging_coefficient_aer__
#else
#define _compute_dry_deposition compute_dry_deposition_
#define _compute_scavenging_coefficient_aer	\
  compute_scavenging_coefficient_aer_
#endif

  extern "C"
  {
    void _compute_dry_deposition(double*, int*, int*, double*, double*,
                                 double*, double*, double*, double*, double*,
                                 double*, double*, double*, double*, double*,
                                 double*, double*, double*, double*);
    void _compute_scavenging_coefficient_aer(int*, double*, double*, double*,
                                             double*, double*, double*,
                                             double*, double*);
  }


  /////////////////////
  // POLAIR3DAEROSOL //
  /////////////////////


  /*! \brief This class is a solver for an advection-diffusion-reaction
    equation with aerosols.*/
  template < class T, class ClassAdvection,
             class ClassDiffusion, class ClassChemistry >
  class Polair3DAerosol:
    public Polair3DChemistry < T, ClassAdvection, ClassDiffusion,
                               ClassChemistry >
  {

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >::iterator
    input_files_iterator;

  public:

    /*** Configuration ***/

    //! List of aerosol species (with their bins) with initial conditions.
    vector<pair<string, vector<int> > > species_list_ic_aer;
    //! List of aerosol number bins with initial conditions.
    vector<int> ic_bin_list_aer;
    //! List of aerosol species (with their bins) with boundary conditions.
    vector<pair<string, vector<int> > > species_list_bc_aer;
    //! List of aerosol number bins with boundary conditions.
    vector<int> bc_bin_list_aer;
    //! List of aerosol bins with deposition velocities.
    vector<int> bin_list_dep_aer;
    //! List of aerosol bins with scavenging.
    vector<int> bin_list_scav_aer;
    //! List of aerosol species (with their bins) with surface emissions.
    vector<pair<string, vector<int> > > species_list_surf_emis_aer;
	//! List of bin aerosol species with surface emissions
	vector<int> surface_emission_bin_list_aer;
	//! List of aerosol species (with their bins) with volume emissions.
    vector<pair<string, vector<int> > > species_list_vol_emis_aer;
	//! List of bin aerosol species with volume emissions
	vector<int> volume_emission_bin_list_aer;

    /*** Domain ***/

    //! 5D grid along z.
    RegularGrid<T> GridZ5D;
    //! 5D grid along y.
    RegularGrid<T> GridY5D;
    //! 5D grid along x.
    RegularGrid<T> GridX5D;

    //! Coordinates of boundary conditions along x.
    RegularGrid<T> GridX5D_interf_bc;
    //! Coordinates of boundary conditions along x.
    RegularGrid<T> GridX3D_interf_bc;
    //! Coordinates of boundary conditions along y.
    RegularGrid<T> GridY5D_interf_bc;
    //! Coordinates of boundary conditions along y.
    RegularGrid<T> GridY3D_interf_bc;

    /*** Species ***/

    //! List of aerosol bins.
    vector<string> bin_list;

    //! List of aerosol fraction.
    vector<string> fraction_list;
    vector<string> fraction_list_2;
    //! Bins bounds (in m).
    Array<T, 1> Fractionbound_aer;
    Array<T, 1> Fractionbound_aer_2;
    
    //! Number of aerosol fraction sections.
    int Nfraction_aer;
    int Nfraction_aer_2;
    string coagulation_coefficient_file;
    
    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;
    //! relations between aerosol species index and groups index
    Array<int, 1> aerosol_species_group_relation;
    //! List of aerosol compositions
    Data<T, 3> composition_bounds;

    //! 4D grid for aerosol groups.
    RegularGrid<T> GridG3D_aer;
    //! 4D grid for aerosol groups.
    RegularGrid<T> GridG4D_aer;
    //! 5D grid for aerosol groups.
    RegularGrid<T> GridG5D_aer;

    //! 4D grid for aerosol compositions.
    RegularGrid<T> GridC4D_aer;
    //! 5D grid for aerosol compositions.
    RegularGrid<T> GridC5D_aer;
    
    //! 4D grid for aerosol bins.
    RegularGrid<T> GridB4D_aer;
    //! 5D grid for aerosol bins.
    RegularGrid<T> GridB5D_aer;

    //! 4D grid for aerosol bins.(for internal data)
    RegularGrid<T> GridB4D_aer_i;
    //! 5D grid for aerosol bins.(for internal data)
    RegularGrid<T> GridB5D_aer_i;
    
    //! 4D grid for aerosol species.
    RegularGrid<T> GridS4D_aer;
    //! 5D grid for aerosol species.
    RegularGrid<T> GridS5D_aer;
    //! 5D grid for aerosol species without H2O aerosol.
    RegularGrid<T> GridS5D_aer_noH2O;

    //! Fixed density (kg/m^3).
    T fixed_density_aer;
    //! Density per bin (kg/m^3).
    Array<T, 1> Density_aer;

    //!Tag of input data format for initial conditions
    string ic_format;
    //! Number of aerosol species with initial conditions.
    int Ns_ic_aer;
    //! Number of aerosol number bins with initial conditions.
    int Nb_ic_aer;

    //! Dissolution heat (kcal / mol at 298K).
    map<string, T> dissolution_heat;

    /*** Land data ***/

    // Land data are used to compute deposition velocities. They depend on the
    // land use category and sometimes on the season. It is assumed that there
    // are the four seasons winter, spring, summer, fall (in this order), plus
    // a fifth season for snow.

    //! Number of land categories.
    int Nland;
    //! Amount (in [0, 1]) of each land category.
    Data<T, 3> LUC;
    //! Roughness height for all land categories and seasons (m).
    Array<T, 2> LandRoughnessHeight;
    /*! Characteristic radius of small receptors (deposition) for all land
      categories and seasons (m). */
    Array<T, 2> SmallRadius;
    /*! Characteristic radius of large receptors (deposition) for all land
      categories and seasons (m). */
    Array<T, 2> LargeRadius;
    //! Coefficient from Zhang deposition model.
    Array<T, 1> ZhangAlpha;
    //! Coefficient from Zhang deposition model.
    Array<T, 1> ZhangGamma;

    /*** Meteorological fields ***/

    //! Temperature at 2 m at current date (K).
    Data<T, 2> SurfaceTemperature_i;
    //! Temperature at 2 m at next date (K).
    Data<T, 2> SurfaceTemperature_f;
    //! Buffer for temperature at 2 m.
    Data<T, 2> FileSurfaceTemperature_i;
    //! Buffer for temperature at 2 m.
    Data<T, 2> FileSurfaceTemperature_f;

    //! Surface pressure at current date (Pa).
    Data<T, 2> SurfacePressure_i;
    //! Surface pressure at next date (Pa).
    Data<T, 2> SurfacePressure_f;
    //! Surface pressure buffer.
    Data<T, 2> FileSurfacePressure_i;
    //! Surface pressure buffer.
    Data<T, 2> FileSurfacePressure_f;

    //! First level wind module at current date (m / s).
    Data<T, 2> FirstLevelWindModule_i;
    //! First level wind module at next date (m / s).
    Data<T, 2> FirstLevelWindModule_f;

    //! Liquid water content (kg / kg).
    Data<T, 3> LiquidWaterContent_i;
    //! Liquid water content buffer.
    Data<T, 3> FileLiquidWaterContent_i;
    //! Liquid water content buffer.
    Data<T, 3> FileLiquidWaterContent_f;

    //! Relative humidity at current date.
    Data<T, 3> RelativeHumidity_i;
    //! Relative humidity at next date.
    Data<T, 3> RelativeHumidity_f;

    //! Snow height at current date (m).
    Data<T, 2> SnowHeight_i;
    //! Snow height at next date (m).
    Data<T, 2> SnowHeight_f;
    //! First level wind module buffer.
    Data<T, 2> FileSnowHeight_i;
    //! First level wind module buffer.
    Data<T, 2> FileSnowHeight_f;

    //! Cloud Optical Depth.
    Data<T, 3> CloudOpticalDepth_i;
    //!  Cloud Optical Depth buffer.
    Data<T, 3> FileCloudOpticalDepth_i;
    //!  Cloud Optical Depth buffer.
    Data<T, 3> FileCloudOpticalDepth_f;

    //! Ice Optical Depth.
    Data<T, 3> IceOpticalDepth_i;
    //! Ice Optical Depth  buffer.
    Data<T, 3> FileIceOpticalDepth_i;
    //! Ice Optical Depth  buffer.
    Data<T, 3> FileIceOpticalDepth_f;

    /*** Boundary conditions ***/

    //!Tag of input data format for boundary conditions.
    string bc_format;
    //! Number of aerosol species with boundary conditions.
    int Ns_bc_aer;
    //! Grid for aerosol species with boundary conditions.
    RegularGrid<T> GridS_bc_aer;
    //! Number of aerosol number bin (sections) with boundary conditions.
    int Nb_bc_aer;
    int Nbin_bc_aer;
    //! Grid for aerosol bins with boundary conditions.
    RegularGrid<T> GridB_bc_aer;
    //! Grid for aerosol bins with boundary conditions.(internal data)
    RegularGrid<T> GridB_bc_aer_i;
    
    //! Boundary conditions along z at current date for aerosols.
    Data<T, 4> BoundaryCondition_z_aer_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_z_aer_i;
    //! Boundary conditions buffer.
    Data<T, 4> FileBoundaryCondition_z_aer_f;
    //! Number boundary conditions along z at current date for aerosols.
    Data<T, 3> NumberBoundaryCondition_z_aer_i;
    //! Number  boundary conditions buffer.
    Data<T, 3> FileNumberBoundaryCondition_z_aer_i;
    //! Number boundary conditions buffer.
    Data<T, 3> FileNumberBoundaryCondition_z_aer_f;

    //! Boundary conditions along y at current date for aerosols.
    Data<T, 5> BoundaryCondition_y_aer_i;
    //! Boundary conditions buffer.
    Data<T, 5> FileBoundaryCondition_y_aer_i;
    //! Boundary conditions buffer.
    Data<T, 5> FileBoundaryCondition_y_aer_f;
    //! Number boundary conditions along z at current date for aerosols.
    Data<T, 4> NumberBoundaryCondition_y_aer_i;
    //! Number  boundary conditions buffer.
    Data<T, 4> FileNumberBoundaryCondition_y_aer_i;
    //! Number boundary conditions buffer.
    Data<T, 4> FileNumberBoundaryCondition_y_aer_f;


    //! Boundary conditions along x at current date for aerosols.
    Data<T, 5> BoundaryCondition_x_aer_i;
    //! Boundary conditions buffer.
    Data<T, 5> FileBoundaryCondition_x_aer_i;
    //! Boundary conditions buffer.
    Data<T, 5> FileBoundaryCondition_x_aer_f;
    //! Number boundary conditions along z at current date for aerosols.
    Data<T, 4> NumberBoundaryCondition_x_aer_i;
    //! Number  boundary conditions buffer.
    Data<T, 4> FileNumberBoundaryCondition_x_aer_i;
    //! Number boundary conditions buffer.
    Data<T, 4> FileNumberBoundaryCondition_x_aer_f;
    //!temperate storage for internal concentration
    Data<T,4> BoundaryCondition_z_aer_i_tmp;
    Data<T,4> FileBoundaryCondition_z_aer_i_tmp;
    Data<T,4> FileBoundaryCondition_z_aer_f_tmp;
    Data<T,5> BoundaryCondition_y_aer_i_tmp;
    Data<T,5> FileBoundaryCondition_y_aer_i_tmp;
    Data<T,5> FileBoundaryCondition_y_aer_f_tmp;
    Data<T,5> BoundaryCondition_x_aer_i_tmp;
    Data<T,5> FileBoundaryCondition_x_aer_i_tmp;
    Data<T,5> FileBoundaryCondition_x_aer_f_tmp;
	  
    Data<T,3> NumberBoundaryCondition_z_aer_i_tmp;
    Data<T,3> FileNumberBoundaryCondition_z_aer_i_tmp;
    Data<T,3> FileNumberBoundaryCondition_z_aer_f_tmp;
    Data<T,4> NumberBoundaryCondition_y_aer_i_tmp;
    Data<T,4> FileNumberBoundaryCondition_y_aer_i_tmp;
    Data<T,4> FileNumberBoundaryCondition_y_aer_f_tmp;
    Data<T,4> NumberBoundaryCondition_x_aer_i_tmp;
    Data<T,4> FileNumberBoundaryCondition_x_aer_i_tmp;
    Data<T,4> FileNumberBoundaryCondition_x_aer_f_tmp;
    Data<T,5> Concentration_aer_i;
    Data<T,4> NumberConcentration_aer_i;
	
    /*** Loss terms ***/

    //! Number of aerosol bins with deposition velocities.
    int Nbin_dep_aer;
    //! Grid for aerosol bins with deposition velocities.
    RegularGrid<T> GridB_dep_aer;
    //! Deposition velocities at current date for aerosols.
    Data<T, 3> DepositionVelocity_aer_i;
    //! Deposition velocities at next date for aerosols.
    Data<T, 3> DepositionVelocity_aer_f;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_aer_i;
    //! Deposition velocities buffer.
    Data<T, 3> FileDepositionVelocity_aer_f;
    //! Dry deposition fluxes at current date.
    Data<T, 4> DryDepositionFlux_aer;
    //! Dry deposition number fluxes at current date.
    Data<T, 3> DryDepositionFluxNumber_aer;

    //! Number of aerosol bins with scavenging.
    int Nbin_scav_aer;
    //! Grid for aerosol bins with scavenging.
    RegularGrid<T> GridB_scav_aer;
    //! Scavenging coefficients at current date for aerosols.
    Data<T, 4> ScavengingCoefficient_aer_i;
    //! Scavenging coefficients at next date for aerosols.
    Data<T, 4> ScavengingCoefficient_aer_f;
    //! Wet deposition fluxes at current date.
    Data<T, 4> WetDepositionFlux_aer;
    //! In cloud wet deposition fluxes at current date.
    Data<T, 4> InCloudWetDepositionFlux_aer;
    //! Wet deposition number fluxes at current date.
    Data<T, 3> WetDepositionFluxNumber_aer;
    //! In cloud wet deposition number fluxes at current date.
    Data<T, 3> InCloudWetDepositionFluxNumber_aer;
	
    /*** Source terms ***/

    //!Tag of input data format for emissions//SZ
    string point_emis_format;
    string surface_emis_format;
    string volume_emis_format;
    //! List of point emissions.
    vector<map<string, string> > point_emission_list_aer;

    //! Number of aerosol species with surface emissions.
    int Ns_surf_emis_aer;
    //! Grid for aerosol species with surface emissions.
    RegularGrid<T> GridS_surf_emis_aer;
    //! Number of bin aerosol species with surface emissions.
    int Nb_surf_emis_aer;
    //! Grid for aerosol bins with surface emissions.
    RegularGrid<T> GridB_surf_emis_aer;
    //! Surface emissions at current date for aerosols.
    Data<T, 4> SurfaceEmission_aer_i;
    //! Surface emissions at next date for aerosols.
    Data<T, 4> SurfaceEmission_aer_f;
    //! Surface emissions buffer.
    Data<T, 4> FileSurfaceEmission_aer_i;
    //! Surface emissions buffer.
    Data<T, 4> FileSurfaceEmission_aer_f;
    //! Number surface emissions at current date for aerosols.
    Data<T, 3> NumberSurfaceEmission_aer_i;
    //! Number surface emissions at next date for aerosols.
    Data<T, 3> NumberSurfaceEmission_aer_f;
    //! Number surface emissions buffer.
    Data<T, 3> FileNumberSurfaceEmission_aer_i;
    //! Number surface emissions buffer.
    Data<T, 3> FileNumberSurfaceEmission_aer_f;

    //! Number of aerosol species with volume emissions.
    int Ns_vol_emis_aer;
    //! Grid for aerosol species with volume emissions.
    RegularGrid<T> GridS_vol_emis_aer;
    //! Number of aerosol bins with volume emissions.
    int Nb_vol_emis_aer;
    //! Grid for aerosol bins with volume emissions.
    RegularGrid<T> GridB_vol_emis_aer;
    //! Number of aerosol species with volume emissions.
    int Nz_vol_emis_aer;
    //! Grid for altitudes of volume emissions for aerosols.
    RegularGrid<T> GridZ_vol_emis_aer;
    //! Volume emissions at current date for aerosols.
    Data<T, 5> VolumeEmission_aer_i;
    //! Volume emissions at next date for aerosols.
    Data<T, 5> VolumeEmission_aer_f;
    //! Volume emissions buffer.
    Data<T, 5> FileVolumeEmission_aer_i;
    //! Volume emissions buffer.
    Data<T, 5> FileVolumeEmission_aer_f;
    //! Number volume emissions at current date for aerosols.
    Data<T, 4> NumberVolumeEmission_aer_i;
    //! Number volume emissions at next date for aerosols.
    Data<T, 4> NumberVolumeEmission_aer_f;
    //! Number volume emissions buffer.
    Data<T, 4> FileNumberVolumeEmission_aer_i;
    //! Number volume emissions buffer.
    Data<T, 4> FileNumberVolumeEmission_aer_f;

    /*** Aerosol physical properties ***/

    //! Wet aerosol diameter (m).
    Data<T, 4> WetDiameter_aer;

    //! Droplet pH.
    Data<T, 3> pH;

    //! Mass density (\mu g.\mu m^-3).
    vector<T> Mass_Density_aer;


    /*** Sources (source splitting) ***/

    //! Sources for aerosols at current date (for source splitting).
    Data<T, 5> Source_aer_i;
    //! Sources for aerosols at next date (for source splitting).
    Data<T, 5> Source_aer_f;

  public:

    /*** Constructor and destructor ***/

    Polair3DAerosol(string config_file);
    virtual ~Polair3DAerosol();

    /*** Configuration ***/

    //! Threshold on liquid water content for clouds (g / m^3).
    T lwc_cloud_threshold;

    // Options related to radiative computation (photolysis rates).
    int option_well_mixed_index, option_black_carbon_treatment;
    int option_wet_index;
    T time_step_for_computing_photolysis_rates;
    int iteration_radiatif;
    Date next_date_radiatif;
    Date previous_date_radiatif;

    string directory_OPAC, file_water_refractive_index, file_species_polyphemus_OPAC;
    string directory_efficiency_factor, fastJ_parameter_files;

    int black_carbon_index;
    int N_OPAC_wavelength;
    int Nwater_wavelength;
    int tabulation_index_real;
    int tabulation_index_imaginary;
    int index_diameter;
    int Nwavelength;
    int NLegendre;
    vector<T> FastJ_wavelength;

    RegularGrid<T> GridWavelength;
    RegularGrid<T> GridIndexReal;
    RegularGrid<T> GridIndexImaginary;
    RegularGrid<T> GridIndexDiameter;
    RegularGrid<T> GridLegendre;
    Data<T, 2> PureSpeciesIndexReal;
    Data<T, 2> PureSpeciesIndexImaginary;
    Data<T, 1> WaterIndexReal;
    Data<T, 1> WaterIndexImaginary;
    Data<T, 4> AbsorptionEfficiencyFactorTable;
    Data<T, 4> ExtinctionEfficiencyFactorTable;
    Data<T, 5> PhaseFunctionTable;
    Data<T, 4> OpticalDepthAerosol;
    Data<T, 4> SingleScatteringAlbedo;
    Data<T, 4> MeanExtinctionEfficiencyFactor;
    Data<T, 4> MeanAbsorbtionEfficiencyFactor;
    Data<T, 5> PhaseFunction;

    Data<string, 1, T> OPACNames;
    Data<string, 1, T> SpeciesNames_polyphemus_OPAC;
    Data<string, 1, T> OPACNames_tmp;
    Data<string, 1, T> SpeciesNames_polyphemus_OPAC_tmp;

    //! Indices of the photolysis reactions among other reactions.
    Array<int, 1> photolysis_reaction_index;

    //! \brief Maps photolysis reaction names to their reaction indices.
    map<string, int> photolysis_reaction_name;

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    bool HasInitialCondition_aer(int s, int b) const;
    bool HasInitialCondition_aer(string name, int b) const;
    bool HasNumberInitialCondition_aer(int b) const;
    int NumberInitialConditionIndex_aer(int b) const;
    vector<int> InitialConditionBinList_aer(int b);
	
    bool HasBoundaryCondition_aer(int s, int b) const;
    bool HasBoundaryCondition_aer(string name, int b) const;
    bool HasNumberBoundaryCondition_aer(int b) const;
    vector<int> BoundaryConditionIndex_aer(int s, int b) const;
    int NumberBoundaryConditionIndex_aer(int b) const;
    vector<int> BoundaryConditionBinList_aer(int b);

    bool HasDepositionVelocity_aer(int b) const;
    int DepositionVelocityIndex_aer(int b) const;

    bool HasScavenging_aer(int b) const;
    int ScavengingIndex_aer(int b) const;

    bool HasSurfaceEmission_aer(int s, int b) const;
    bool HasSurfaceEmission_aer(string name, int b) const;
    bool HasNumberSurfaceEmission_aer(int b) const;
    vector<int> SurfaceEmissionIndex_aer(int s, int b) const;
    int NumberSurfaceEmissionIndex_aer(int b) const;
    vector<int> SurfaceEmissionBinList_aer(int b);
    bool HasVolumeEmission_aer(int s, int b) const;
    bool HasVolumeEmission_aer(string name, int b) const;
    bool HasNumberVolumeEmission_aer(int b) const;
    int VolumeEmissionIndex_aer(int s) const;
    string VolumeEmissionName_aer(int s) const;
    vector<int> VolumeEmissionBinList_aer(int b);
    int NumberVolumeEmissionIndex_aer(int b) const;
    int VolumeEmissionGlobalIndex_aer(int s) const;
    int Bin_to_size_index_aer(int s) const;
    int Bin_to_composition_index(int s) const;
    int Bin_index_translate_aer(int j, int k) const;

    /*** Initializations ***/

    virtual void Allocate();
    virtual void Init();
    virtual void InitStep();
    void InitWetDiameter_aer(Data<T, 3>& RelativeHumidity_,
			     Data<T, 3>& Temperature_,
			     Data<T, 4>& WetDiameter_aer_);
    void InitWetDiameter_aer(Data<T, 3>& RelativeHumidity_,
			     Data<T, 3>& Temperature_,
			     Data<T, 5>& Concentration_aer_,
			     Data<T, 4>& NumberConcentration_aer_,
			     Data<T, 4>& WetDiameter_aer_,
			     bool with_computation_drydiameter);
    T ComputeDensity(Data<T, 1> Conc_aer_tmp,
		     vector<T> Rho_species, T TotalMass, int Ns);
    void ComputeNumberConcentration_forIC_aer(int b);
    void ComputeNumberBoundaryCondition_z_aer(int b);
    void ComputeNumberBoundaryCondition_y_aer(int b);
    void ComputeNumberBoundaryCondition_x_aer(int b);
    void ComputeNumberSurfaceEmission_aer(int b);
    void ComputeNumberVolumeEmission_aer(int b);

    void InitDepositionVelocity(Data<T, 3>& Temperature_,
				Data<T, 3>& Pressure_,
				Data<T, 2>& SurfaceTemperature_,
				Data<T, 2>& SurfacePressure_,
				Data<T, 2>& FirstLevelWindModule_,
				Data<T, 4>& WetDiameter_aer_,
				Data<T, 2>& SnowHeight_,
				Data<T, 3>& DepositionVelocity_aer_);
    void InitScavengingCoefficient(Data<T, 3>& Temperature_,
                                   Data<T, 3>& Pressure_,
                                   Data<T, 3>& LiquidWaterContent_,
                                   Data<T, 3>& pH_,
                                   Data<T, 2>& CloudBaseHeight_,
                                   Data<T, 2>& Rain_,
                                   Data<T, 4>&
                                   ScavengingCoefficient_);
    void InitScavengingCoefficient_aer(Data<T, 3>& Temperature_,
                                       Data<T, 3>& Pressure_,
                                       Data<T, 4>& WetDiameter_aer_,
                                       Data<T, 2>& CloudBaseHeight_,
                                       Data<T, 2>& Rain_,
                                       Data<T, 4>&
                                       ScavengingCoefficient_aer_);

    virtual void SetDate(Date date);
    
   /*** External composition methods ***///SZ
    int FindExternalCompositionID(string species);
    
    int FindCompositionID(Data<T,2>& total_mass_,
			    Data<T,3>& group_mass_,
			    Data<T,3>& composition_bounds_,
			    int Nc ,int Ng, int y, int x);
			    
    int FindCompositionID(Data<T,3>& total_mass_,
			    Data<T,4>& group_mass_,
			    Data<T,3>& composition_bounds_,
			    int Nc ,int Ng, int z, int y, int x);

    bool BoundaryConditionTransformation();

    /*** Radiative ***/
    void Radiatif(Data<T, 4>& PhotolysisRate, Date date);

    /*** Integration ***/

    virtual void Advection();
    virtual void Diffusion();
    void PointEmission_aer();
    virtual void Chemistry();

    virtual void Forward();

    /*** Access methods ***/

    virtual T GetConcentration_aer(int s, int b, T z, T y, T x);
    using BaseModel<T>::GetConcentration_aer;
    int GetNs_source_aer();
    int GetNz_source_aer();
    int GetNbinMax_source_aer();
    int SourceGlobalBinIndex_aer(int s, int b);
    bool HasSource_aer(int s, int b);
    int SourceGlobalIndex_aer(int s);

    Data<T, 5>& GetSource_aer_i();
    Data<T, 5>& GetSource_aer_f();
	virtual bool HasNumberConcentration_aer();

    Data<T, 3>& GetDryDepositionFluxNumber_aer();
    Data<T, 3>& GetWetDepositionFluxNumber_aer();
    Data<T, 3>& GetInCloudWetDepositionFluxNumber_aer();

    Data<T, 4>& GetDryDepositionFlux_aer();
    Data<T, 4>& GetWetDepositionFlux_aer();
    Data<T, 4>& GetInCloudWetDepositionFlux_aer();

  protected:

    virtual void InitAllData();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POLAIR3DAEROSOL_HXX
#endif
