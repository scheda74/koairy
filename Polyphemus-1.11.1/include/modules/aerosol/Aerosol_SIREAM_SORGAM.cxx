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


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SORGAM_CXX


#include "Aerosol_SIREAM_SORGAM.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  Aerosol_SIREAM_SORGAM<T>::Aerosol_SIREAM_SORGAM()
  {
    _dimensions(&Ns, &Nr, &Nr_photolysis);

    Ncycle = 1;
    molecular_weight.resize(Ns);
    photolysis_reaction_index.resize(Nr_photolysis);
    Ncycle_aer = 1;
    Noptions_aer = 14;
  }


  //! Initialization of the scheme.
  /*! In case the model is incompatible with the chemical mechanism, an
    exception is thrown.
    \param Model model with the following interface:
    <ul>
    <li> GetNs()
    <li> GetSpeciesList()
    <li> GetSpeciesFile()
    <li> GetNr_photolysis()
    <li> GetPhotolysisReactionList()
    <li> GetNs_source()
    <li> GetNz_source()
    <li> SourceGlobalIndex(int)
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SIREAM_SORGAM<T>::Init(ClassModel& Model)
  {
    if (Model.GetNs() != Ns)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNs())
        + string(" species while chemistry has ") + to_str(Ns) + " species.";
    if (Model.GetNr_photolysis() != Nr_photolysis)
      throw string("Incompatibility: model manages ")
        + to_str(Model.GetNr_photolysis())
        + string(" photolysis reactions while chemistry has ")
        + to_str(Nr_photolysis) + " photolysis reactions.";
    species_list = Model.GetSpeciesList();

    Nbin_aer = Model.GetNbin_aer();
    Ns_aer = Model.GetNs_aer();

    // Photolysis reactions.
    const string& config_file = Model.GetConfigurationFile();
    const string& species_file = Model.GetSpeciesFile();
    ConfigStreams config_species(config_file, species_file);

    string molecular_weight_file = species_file;
    config_species.SetSection("[molecular_weight]");

    // Checks if there is an included file.
    if (config_species.GetElement() == "file")
      molecular_weight_file = config_species.GetElement();
    ConfigStream molecular_weight_config(molecular_weight_file);
    molecular_weight_config.SetSection("[molecular_weight]");
    for (int i = 0; i < Ns; i++)
      molecular_weight_config.PeekValue(species_list[i], molecular_weight(i));

    config_species.SetSection("[photolysis_reaction_index]");
    string species;
    for (int i = 0; i < Nr_photolysis; i++)
      {
        species = config_species.GetElement();
        photolysis_reaction_name[species] =
          convert<int>(config_species.GetElement());
      }

    for (int i = 0; i < Nr_photolysis; i++)
      photolysis_reaction_index(i) =
        photolysis_reaction_name[Model.GetPhotolysisReactionList()[i]];

    // Conversions.
    double Navogadro = 6.02213e23;
    ConversionFactor.resize(Ns);
    for (int i = 0; i < Ns; i++)
      ConversionFactor(i) = Navogadro * 1e-12 / molecular_weight(i);
    ConversionFactorJacobian.resize(Ns, Ns);
    for (int i = 0; i < Ns; i++)
      for (int j = 0; j < Ns; j++)
        ConversionFactorJacobian(i, j) = molecular_weight(j)
          / molecular_weight(i);



    // Volumic sources.
    Ns_source = Model.GetNs_source();
    Nz_source = Model.GetNz_source();
    source_index.resize(Ns_source);
    for (int i = 0; i < Ns_source; i++)
      source_index(i) = Model.SourceGlobalIndex(i);

    ConfigStream config(Model.GetConfigurationFile());
    config.SetSection("[options]");

    config.PeekValue("With_coagulation",
                     option_process_aer["with_coagulation"]);
    config.PeekValue("With_condensation",
                     option_process_aer["with_condensation"]);
    config.PeekValue("With_nucleation",
                     option_process_aer["with_nucleation"]);
    config.PeekValue("With_in_cloud_scavenging",
                     option_process_aer["with_in_cloud_scavenging"]);
    config.PeekValue("Collect_wet_flux_aerosol",
                     option_process_aer["collect_wet_flux_aer"]);
    config.PeekValue("Collect_wet_flux",
                     option_process_aer["collect_wet_flux"]);
    config.PeekValue("aqueous_module", aqueous_module);
    aqueous_module = lower_case(aqueous_module);
    config.PeekValue("With_heterogeneous_reactions",
                     option_process_aer["with_heterogeneous_reactions"]);
    config.PeekValue("With_kelvin_effect",
                     option_process_aer["with_kelvin_effect"]);
    config.PeekValue("With_volume_emission_aerosol",
                     option_process_aer["with_volume_emission_aer"]);
    config.PeekValue("Dynamic_condensation_solver",
                     dynamic_condensation_solver);
    dynamic_condensation_solver = lower_case(dynamic_condensation_solver);
    config.PeekValue("Fixed_cutting_diameter", fixed_cutting_diameter);
    fixed_cutting_diameter = 1.e-6 * fixed_cutting_diameter;
    config.PeekValue("Sulfate_computation", sulfate_computation);
    sulfate_computation = lower_case(sulfate_computation);
    config.PeekValue("Redistribution_method", redistribution_method);
    redistribution_method = lower_case(redistribution_method);
    config.PeekValue("Nucleation_model", nucleation_model);
    nucleation_model = lower_case(nucleation_model);
    config.PeekValue("With_fixed_density",
                     option_process_aer["with_fixed_density"]);
    config.PeekValue("Wet_diameter_estimation", wet_diameter_estimation);
    wet_diameter_estimation = lower_case(wet_diameter_estimation);
    config.PeekValue("Thermodynamics_module", thermodynamics_module);
    thermodynamics_module = lower_case(thermodynamics_module);

    config.PeekValue("With_adaptive_time_step_for_gas_chemistry",
                     with_adaptive_time_step);
    if (with_adaptive_time_step)
      {
        option_adaptive_time_step = 1;
        config.PeekValue("Adaptive_time_step_tolerance",
                         adaptive_time_step_tolerance);
        config.PeekValue("Min_adaptive_time_step", min_adaptive_time_step);
      }
    else
      option_adaptive_time_step = 0;

    config.PeekValue("Photolysis_tabulation_option", option_photolysis_tabulation);

    // Density in kg / m^3.
    config.PeekValue("Fixed_aerosol_density", FixedDensity_aer);
    Density_aer.resize(Nbin_aer);
    Density_aer = FixedDensity_aer;

    // Lwc cloud threshold.
    config.PeekValue("Lwc_cloud_threshold", lwc_cloud_threshold);

    // Reads bin bounds.
    vector<string> bin_list;
    config.SetSection("[domain]");
    config.Find("Bin_bounds");
    bin_list = split(config.GetLine());

    // Reads bin bounds in micrometers and converts it to meters.
    BinBound_aer.resize(Nbin_aer + 1);
    for (int i = 0; i < Nbin_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]);

    // Compute coagulation partition coefficient, if coagulation active.
    CoagulationCouple.resize(Nbin_aer);
    CoagulationFirstIndex.resize(4 * Nbin_aer, Nbin_aer);
    CoagulationSecondIndex.resize(4 * Nbin_aer, Nbin_aer);
    CoagulationPartitionCoefficient.resize(Nbin_aer, Nbin_aer,
                                           Nbin_aer);

    if (option_process_aer["with_coagulation"])
      {
        _compute_coagulation_coefficient(&Nbin_aer,
                                         BinBound_aer.data(),
                                         CoagulationCouple.data(),
                                         CoagulationFirstIndex.data(),
                                         CoagulationSecondIndex.data(),
                                         CoagulationPartitionCoefficient.
                                         data());
      }

    // Bin corresponding to cutting diameter (for condensation).
    cutting_bin = 0;
    while (BinBound_aer(cutting_bin) < fixed_cutting_diameter &&
           cutting_bin < Nbin_aer)
      cutting_bin++;

    // List of aerosol options to be passed to fortran routine.
    options_aer.resize(Noptions_aer);

    options_aer(0) = option_process_aer["with_coagulation"] ? 1 : 0;

    options_aer(1) = option_process_aer["with_condensation"] ? 1 : 0;

    options_aer(2) = option_process_aer["with_nucleation"] ? 1 : 0;

    if (aqueous_module == "no")
      options_aer(3) = 0;
    else if (aqueous_module == "vsrm")
      options_aer(3) = 1;
    else if (aqueous_module == "simple") // simple
      options_aer(3) = 2;
    else
      {
        options_aer(3) = 2;
        aqueous_module = "simple";
      }

    options_aer(4) = option_process_aer["with_in_cloud_scavenging"] ? 1 : 0;
    options_aer(5) = option_process_aer["with_kelvin_effect"] ? 1 : 0;
    options_aer(6) = option_process_aer["with_fixed_density"] ? 0 : 1;
    options_aer(7) = cutting_bin;

    if (sulfate_computation == "dynamic")
      options_aer(8) = 1;
    else if (sulfate_computation == "equilibrium")
      options_aer(8) = 0;
    else
      throw string("bad string \"") + sulfate_computation +
                                       string("\" for sulfate computation option,\n") +
                                       string("possibilities are dynamic or equilibrium.");

    if (dynamic_condensation_solver == "etr")
      options_aer(9) = 1;
    else if (dynamic_condensation_solver == "ros2")
      options_aer(9) = 2;
    else if (dynamic_condensation_solver == "ebi")
      options_aer(9) = 3;
    else
      throw string("bad string \"") + dynamic_condensation_solver +
                                       string("\" for dynamic c/e solver option,\n") +
                                       string("possibilities are etr, ros2 or ebi.");

    if (redistribution_method == "number-conserving")
      options_aer(10) = 1;
    else if (redistribution_method == "interpolation")
      options_aer(10) = 2;
    else
      throw string("bad string \"") + redistribution_method +
                                       string("\" for redistribution method option,\n") +
                                       string("possibilities are number-conserving or interpolation.");

    if (nucleation_model == "binary")
      options_aer(11) = 0;
    else if (nucleation_model == "ternary")
      options_aer(11) = 1;
    else
      throw string("bad string \"") + nucleation_model +
                                       string("\" for nucleation model option,\n") +
                                       string("possibilities are binary or ternary.");

    if (wet_diameter_estimation == "isoropia"
        || wet_diameter_estimation == "isorropia")
      options_aer(12) = 0;
    else if (wet_diameter_estimation == "gerber")
      options_aer(12) = 1;
    else
      throw string("bad string \"") + wet_diameter_estimation +
                                       string("\" for wet diameter estimation option,\n") +
                                       string("possibilities are gerber or isorropia.");

    if (thermodynamics_module == "isoropia"
        || thermodynamics_module == "isorropia")
      options_aer(13) = 0;
    else if (thermodynamics_module == "eqsam")
      options_aer(13) = 1;
    else
      throw string("bad string \"") + thermodynamics_module +
                                       string("\" for thermodynamics module option,\n") +
                                       string("possibilities are eqsam or isorropia.");

    // Index of species subject to heterogeneous reactions.
    heterogeneous_reaction_index.resize(4);
    heterogeneous_reaction_index(0) = Model.GetSpeciesIndex("HO2");
    heterogeneous_reaction_index(1) = Model.GetSpeciesIndex("NO2");
    heterogeneous_reaction_index(2) = Model.GetSpeciesIndex("NO3");
    heterogeneous_reaction_index(3) = Model.GetSpeciesIndex("N2O5");

    if (option_process_aer["with_in_cloud_scavenging"])
      {
        if (option_process_aer["collect_wet_flux_aer"])
          {
            if (!Model.HasField("InCloudWetDepositionFlux_aer"))
              throw string("if field \" InCloudWetDepositionFlux_aer \" ") +
                string("does not exist in Model, in cloud wet deposition ") +
                string("fluxes cannot be collected.");
            if (Model.D4("InCloudWetDepositionFlux_aer").GetLength(2)
                != Model.GetNy()
                || Model.D4("InCloudWetDepositionFlux_aer").GetLength(3)
                != Model.GetNx())
              throw string("\" InCloudWetDepositionFlux_aer \" has") +
                string(" wrong dimension.");
          }
        if (option_process_aer["collect_wet_flux"])
          {
            if (!Model.HasField("InCloudWetDepositionFlux"))
              throw string("if field \" InCloudWetDepositionFlux \" ") +
                string("does not exist in Model, in cloud wet deposition ") +
                string("fluxes cannot be collected.");
            if (Model.D3("InCloudWetDepositionFlux").GetLength(1)
                != Model.GetNy()
                || Model.D3("InCloudWetDepositionFlux").GetLength(2)
                != Model.GetNx())
              throw string("\" InCloudWetDepositionFlux \" has") +
                string(" wrong dimension.\n"
                       "(Have you forgotten to set the in cloud scavenging model "
                       "to 'aqueous'?");
          }
      }

    BaseModuleParallel::Init(Model);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Partition along axis x.
    BuildPartition_x();
#endif

    /*** Display options ***/
    config.SetSection("[display]");
    // Should the configuration be displayed on screen?
    config.PeekValue("Show_configuration",
                     option_process_aer["show_configuration"]);

    if (option_process_aer["show_configuration"] && GetRank() == 0)
      DisplayConfiguration();
  }


  //! Display the configuration.
  template<class T>
  void Aerosol_SIREAM_SORGAM<T>::DisplayConfiguration()
  {

    if (with_adaptive_time_step)
      cout << "Module: with_adaptive_time_step_for_gas_chemistry" << endl;

    cout << "Module: photolysis tabulation option: ";
    if (option_photolysis_tabulation == 1)
      cout << "Tabulation by SPACK" << endl;
    else if (option_photolysis_tabulation == 2)
      cout << "Reading the binary files obtained by FastJ" << endl;
    else
      throw string("Wrong number is given. Please put 1 or 2");

    cout << "Module: Fixed_aerosol_density = " << FixedDensity_aer << endl;
    cout << "Module: Lwc_cloud_threshold = " << lwc_cloud_threshold << endl;
    if (option_process_aer["with_coagulation"])
      {
        cout << "Module: computing coagulation coefficients ... ";
        cout << "done" << endl;
      }

    if (option_process_aer["with_coagulation"])
      cout << "Module: with coagulation" << endl;

    if (option_process_aer["with_condensation"])
      cout << "Module: with condensation" << endl;

    if (option_process_aer["with_nucleation"])
      cout << "Module: with nucleation" << endl;

    cout << "Module: cloud chemistry:" << aqueous_module << endl;

    if (option_process_aer["with_in_cloud_scavenging"])
      cout << "Module: with in-cloud scavenging" << endl;

    if (option_process_aer["with_heterogeneous_reactions"])
      cout << "Module: with heterogeneous reactions" << endl;

    if (option_process_aer["with_kelvin_effect"])
      cout << "Module: with kelvin effect" << endl;

    if (option_process_aer["with_fixed_density"])
      cout << "Module: with fixed density" << endl;

    if (cutting_bin == 0)
      cout << "Module: all bins at dynamic." << endl;
    else if (cutting_bin == Nbin_aer)
      cout << "Module: all bins at equilibrium." << endl;
    else
      cout << "Module: bins from 1 to " << cutting_bin <<
        "(included) at equilibrium." << endl;

    cout << "Module: dynamic condensation solver: " <<
      dynamic_condensation_solver << endl;

    cout << "Module: redistribution method: " <<
      redistribution_method << endl;

    if (option_process_aer["with_nucleation"])
      cout << "Module: nucleation model: " <<
        nucleation_model << endl;

    cout << "Module: wet diameter estimation: " <<
      wet_diameter_estimation << endl;

    cout << "Module: thermodynamics module: " <<
      thermodynamics_module << endl;

#ifdef WITHOUT_NACL_IN_THERMODYNAMICS
    cout << "Module: without sea salt in isorropia" << endl;
#else
    cout << "Module: with sea salt in isorropia" << endl;
#endif
  }

  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i")
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetSource_i()
    <li> D4("PhotolysisRate_i")
    <li> D3("Attenuation_f")
    <li> D3("SpecificHumidity_f")
    <li> D3("Temperature_f")
    <li> D3("Pressure_f")
    <li> GetSource_f()
    <li> D4("PhotolysisRate_f")
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D4("WetDiameter_aer")
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SIREAM_SORGAM<T>::Forward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

    Forward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
            Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
            Model.D3("Pressure_i"), Model.GetSource_i(),
            Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
            Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
            Model.D3("Temperature_f"), Model.D3("Pressure_f"),
            Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
            Model.GetGridXArray1D(), Model.GetGridYArray1D(),
            Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
            Model.D4("WetDiameter_aer"), Model.GetConcentration_aer());
  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D2("Rain_i")
    <li> GetLayerInterface()
    <li> GetConcentration_aer()
    <li> D3("pH")
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SIREAM_SORGAM<T>::Forward_aer(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

    Data<T, 3> in_cloud_wet_flux(Ns, Model.GetNy(), Model.GetNx());
    Data<T, 4> in_cloud_wet_flux_aer(Ns_aer, Nbin_aer, Model.GetNy(),
                                     Model.GetNx());

    Forward_aer(T(date_i.GetNumberOfSeconds()),
                Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
                Model.D3("Pressure_i"), T(date_f.GetSecondsFrom(date_i)),
                Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
                Model.D2("Rain_i"), Model.GetLayerInterface(),
                Model.GetConcentration_aer(), in_cloud_wet_flux,
                in_cloud_wet_flux_aer, Model.D3("pH"));

    if (option_process_aer["with_in_cloud_scavenging"])
      {
        if (option_process_aer["collect_wet_flux"])
          for (int s = 0; s < Ns; s++)
            if (Model.HasScavenging(s))
              for (int j = 0; j < Model.GetNy(); j++)
                for (int i = 0; i < Model.GetNx(); i++)
                  Model.D3("InCloudWetDepositionFlux")
                    (Model.ScavengingIndex(s), j, i) =
                    in_cloud_wet_flux(s, j, i);
        if (option_process_aer["collect_wet_flux_aer"])
          for (int b = 0; b < Nbin_aer; b++)
            if (Model.HasScavenging_aer(b))
              for (int s = 0; s < Ns_aer; s++)
                for (int j = 0; j < Model.GetNy(); j++)
                  for (int i = 0; i < Model.GetNx(); i++)
                    Model.D4("InCloudWetDepositionFlux_aer")
                      (s, Model.ScavengingIndex_aer(b), j, i) =
                      in_cloud_wet_flux_aer(s, b, j, i);
      }
  }


  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param next_time time at the end of the time step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param WetDiameter_aer Aerosol wet diameters at the beginning of the time
    step.
    \param Concentration_aer aerosol concentrations.
  */
  template<class T>
  void Aerosol_SIREAM_SORGAM<T>::Forward(T current_time,
                                         Data<T, 3>& Attenuation_i,
                                         Data<T, 3>& SpecificHumidity_i,
                                         Data<T, 3>& Temperature_i,
                                         Data<T, 3>& Pressure_i,
                                         Data<T, 4>& Source_i,
                                         Data<T, 4>& PhotolysisRate_i,
                                         T delta_t,
                                         Data<T, 3>& Attenuation_f,
                                         Data<T, 3>& SpecificHumidity_f,
                                         Data<T, 3>& Temperature_f,
                                         Data<T, 3>& Pressure_f,
                                         Data<T, 4>& Source_f,
                                         Data<T, 4>& PhotolysisRate_f,
                                         Array<T, 1>& Longitude,
                                         Array<T, 1>& Latitude,
                                         Data<T, 4>& Concentration,
                                         Data<T, 3>& LiquidWaterContent_i,
                                         Data<T, 4>& WetDiameter_aer,
                                         Data<T, 5>& Concentration_aer)

  {
    int Nz = Concentration.GetLength(1);
    int Ny = Concentration.GetLength(2);
    int Nx = Concentration.GetLength(3);

    int icld = options_aer(3);
    int iheter = option_process_aer["with_heterogeneous_reactions"] ? 1 : 0;

#if defined(POLYPHEMUS_PARALLEL_WITH_MPI)	\
  || defined(POLYPHEMUS_PARALLEL_WITH_OPENMP)
    int first_index_y(1);
    int first_index_slice_x, last_index_slice_x;
    int Nthreads(1);
    Array<int, 1> first_index_subslice_x(1), last_index_subslice_x(1);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GetEdgePartition_x(first_index_slice_x, last_index_slice_x);
    ScatterSlice_x_MPI(Concentration);
    ScatterSlice_x_MPI(Concentration_aer);
#else
    first_index_slice_x = 0;
    last_index_slice_x = Nx;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    Nthreads = GetNthreads_openmp();
    BuildSubSegment(first_index_slice_x, last_index_slice_x,
                    first_index_subslice_x, last_index_subslice_x);
#pragma omp parallel for num_threads(Nthreads)			\
  shared(Nz, Ny, Nx, icld, iheter, first_index_subslice_x,	\
         last_index_subslice_x, Nthreads)
#else
    first_index_subslice_x(0) = first_index_slice_x + 1;
    last_index_subslice_x(0) = last_index_slice_x;
#endif
    for (int c = 0; c < Nthreads; c++)
      _chem(first_index_subslice_x.data() + c, last_index_subslice_x.data() + c,
            &first_index_y, &Ny,
            &Nx, &Ny, &Nz, &Ns, &Nr,
            &Nr_photolysis, photolysis_reaction_index.data(),
            &Ns_source, source_index.data(),
            ConversionFactor.data(), ConversionFactorJacobian.data(), &Nz_source,
            &lwc_cloud_threshold, molecular_weight.data(),
            &current_time, Attenuation_i.GetData(),
            SpecificHumidity_i.GetData(), Temperature_i.GetData(),
            Pressure_i.GetData(), Source_i.GetData(),
            PhotolysisRate_i.GetData(), &delta_t,
            Attenuation_f.GetData(),
            SpecificHumidity_f.GetData(), Temperature_f.GetData(),
            Pressure_f.GetData(), Source_f.GetData(),
            PhotolysisRate_f.GetData(), &Ncycle, Longitude.data(),
            Latitude.data(), Concentration.GetData(),
            &icld, &iheter, &Ns_aer, &Nbin_aer,
            LiquidWaterContent_i.GetData(),
            BinBound_aer.data(), &FixedDensity_aer,
            WetDiameter_aer.GetData(),
            heterogeneous_reaction_index.data(),
            Concentration_aer.GetData(),
            &option_adaptive_time_step, &adaptive_time_step_tolerance,
            &min_adaptive_time_step, &option_photolysis_tabulation);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(Concentration);
    GatherSlice_x_MPI(Concentration_aer);
#endif

#else
#endif
  }


  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param next_time time at the end of the time step.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param Rain_i rain rate at the beginning of the time step.
    \param VerticalInterface height of vertical interfaces in meter.
    \param Concentration_aer aerosol concentrations.
    \param aerosol pH.
  */
  template<class T>
  void Aerosol_SIREAM_SORGAM<T>::Forward_aer(T current_time,
                                             Data<T, 3>& SpecificHumidity_i,
                                             Data<T, 3>& Temperature_i,
                                             Data<T, 3>& Pressure_i,
                                             T delta_t,
                                             Data<T, 4>& Concentration,
                                             Data<T, 3>& LiquidWaterContent_i,
                                             Data<T, 2>& Rain_i,
                                             Array<T, 1>& VerticalInterface,
                                             Data<T, 5>& Concentration_aer,
                                             Data<T, 3>&
                                             InCloudWetDepositionFlux,
                                             Data<T, 4>&
                                             InCloudWetDepositionFlux_aer,
                                             Data<T, 3>& pH)
  {
    int Nz = Concentration.GetLength(1);
    int Ny = Concentration.GetLength(2);
    int Nx = Concentration.GetLength(3);

#if defined(POLYPHEMUS_PARALLEL_WITH_MPI)	\
  || defined(POLYPHEMUS_PARALLEL_WITH_OPENMP)
    int first_index_y(1);
    int first_index_slice_x, last_index_slice_x;
    int Nthreads(1);
    Array<int, 1> first_index_subslice_x(1), last_index_subslice_x(1);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GetEdgePartition_x(first_index_slice_x, last_index_slice_x);
    ScatterSlice_x_MPI(Concentration);
    ScatterSlice_x_MPI(Concentration_aer);
    ScatterSlice_x_MPI(pH);
    ScatterSlice_x_MPI(InCloudWetDepositionFlux_aer);
    ScatterSlice_x_MPI(InCloudWetDepositionFlux);
#else
    first_index_slice_x = 0;
    last_index_slice_x = Nx;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    Nthreads = GetNthreads_openmp();
    BuildSubSegment(first_index_slice_x, last_index_slice_x,
                    first_index_subslice_x, last_index_subslice_x);
#pragma omp parallel for num_threads(Nthreads)				\
  shared(Nz, Ny, Nx, first_index_subslice_x, last_index_subslice_x, Nthreads)
#else
    first_index_subslice_x(0) = first_index_slice_x + 1;
    last_index_subslice_x(0) = last_index_slice_x;
#endif
    for (int c = 0; c < Nthreads; c++)
      _aerosol(first_index_subslice_x.data() + c,
               last_index_subslice_x.data() + c, &first_index_y, &Ny,
               &Nx, &Ny, &Nz, &Ns,
               &lwc_cloud_threshold, &current_time,
               SpecificHumidity_i.GetData(), Temperature_i.GetData(),
               Pressure_i.GetData(), &delta_t,
               Concentration.GetData(),
               &Noptions_aer, options_aer.data(), &Ns_aer, &Nbin_aer,
               &Ncycle_aer, LiquidWaterContent_i.GetData(),
               Rain_i.GetData(),
               BinBound_aer.data(), &FixedDensity_aer,
               Density_aer.data(),
               CoagulationCouple.data(), CoagulationFirstIndex.data(),
               CoagulationSecondIndex.data(),
               CoagulationPartitionCoefficient.data(),
               VerticalInterface.data(),
               Concentration_aer.GetData(),
               InCloudWetDepositionFlux.GetData(),
               InCloudWetDepositionFlux_aer.GetData(), pH.GetData());

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(Concentration);
    GatherSlice_x_MPI(Concentration_aer);
    GatherSlice_x_MPI(pH);
    GatherSlice_x_MPI(InCloudWetDepositionFlux_aer);
    GatherSlice_x_MPI(InCloudWetDepositionFlux);
#endif

#else
#endif
  }


  //! Checks whether a given field is required by this module.
  /*! Checks whether a given field must be available in the underlying model
    for this module to perform properly.
    \param field the field name.
    \return True is the field is required, false otherwise.
  */
  template<class T>
  bool Aerosol_SIREAM_SORGAM<T>::IsRequired(string field)
  {
    if (field == "wet_diameter_aer") return true;
    if (field == "rain") return true;
    return false;
  }


  //! Checks whether a given field is computed by this module.
  /*! Checks whether a given field is computed by this module and updated in
    the underlying model.
    \param field the field name.
    \return True is the field is computed, false otherwise.
  */
  template<class T>
  bool Aerosol_SIREAM_SORGAM<T>::IsComputed(string field)
  {
    if (field == "pH") return true;
    return false;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SIREAM_SORGAM_CXX
#endif
