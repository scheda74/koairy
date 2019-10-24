// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYRADM_CXX


#include "ChemistryRADM.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  ChemistryRADM<T>::ChemistryRADM()
  {
    _dimensions_radm(&Ns, &Nr, &Nr_photolysis);

    Ncycle = 1;
    photolysis_reaction_index.resize(Nr_photolysis);
    molecular_weight.resize(Ns);
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
  void ChemistryRADM<T>::Init(ClassModel& Model)
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

    const string& config_file = Model.GetConfigurationFile();
    const string& species_file = Model.GetSpeciesFile();
    ConfigStreams config(config_file, species_file);

    string molecular_weight_file = species_file;
    config.SetSection("[molecular_weight]");

    // Checks if there is an included file.
    if (config.GetElement() == "file")
      molecular_weight_file = config.GetElement();
    ConfigStream molecular_weight_config(molecular_weight_file);
    molecular_weight_config.SetSection("[molecular_weight]");
    for (int i = 0; i < Ns; i++)
      molecular_weight_config.PeekValue(species_list[i], molecular_weight(i));

    config.SetSection("[photolysis_reaction_index]");
    string species;
    for (int i = 0; i < Nr_photolysis; i++)
      {
        species = config.GetElement();
        photolysis_reaction_name[species] = convert<int>(config.GetElement());
      }

    for (int i = 0; i < Nr_photolysis; i++)
      photolysis_reaction_index(i) =
        photolysis_reaction_name[Model.GetPhotolysisReactionList()[i]];

    Ns_source = Model.GetNs_source();
    Nz_source = Model.GetNz_source();
    source_index.resize(Ns_source);
    for (int i = 0; i < Ns_source; i++)
      source_index(i) = Model.SourceGlobalIndex(i);

    // To be called even if there is no parallelization.
    BaseModuleParallel::Init(Model);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Partition along axis x.
    BuildPartition_x();
#endif

    // Photolysis options.
    ConfigStream config_main(Model.GetConfigurationFile());
    config_main.SetSection("[options]");
    config_main.PeekValue("Photolysis_tabulation_option", option_photolysis_tabulation);

  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i
    <li> D3("SpecificHumidity_i"),
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
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void ChemistryRADM<T>::Forward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (Model.source_splitting)
      {
        ScatterSlice_x_MPI(Model.GetSource_i());
        ScatterSlice_x_MPI(Model.GetSource_f());
      }
#endif

    Forward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
            Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
            Model.D3("Pressure_i"), Model.GetSource_i(),
            Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
            Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
            Model.D3("Temperature_f"), Model.D3("Pressure_f"),
            Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
            Model.GetGridXArray1D(), Model.GetGridYArray1D(),
            Model.GetConcentration());
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
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration concentrations.
  */
  template<class T>
  void ChemistryRADM<T>::Forward(T current_time,
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
                                 Data<T, 4>& Concentration)
  {
    int Nz = Concentration.GetLength(1);
    int Ny = Concentration.GetLength(2);
    int Nx = Concentration.GetLength(3);

    int Nthreads(1);
    int first_index_slice_x, last_index_slice_x;
    Array<int, 1> first_index_subslice_x(1), last_index_subslice_x(1);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GetEdgePartition_x(first_index_slice_x, last_index_slice_x);
    ScatterSlice_x_MPI(Concentration);
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
      {
        _chem_radm(first_index_subslice_x.data() + c,
                   last_index_subslice_x.data() + c,
                   &Nx, &Ny, &Nz, &Ns, &Nr,
                   &Nr_photolysis, photolysis_reaction_index.data(),
                   &Ns_source, source_index.data(), &Nz_source,
                   molecular_weight.data(),
                   &current_time, Attenuation_i.GetData(),
                   SpecificHumidity_i.GetData(), Temperature_i.GetData(),
                   Pressure_i.GetData(), Source_i.GetData(),
                   PhotolysisRate_i.GetData(), &delta_t,
                   Attenuation_f.GetData(),
                   SpecificHumidity_f.GetData(), Temperature_f.GetData(),
                   Pressure_f.GetData(), Source_f.GetData(),
                   PhotolysisRate_f.GetData(), &Ncycle, Longitude.data(),
                   Latitude.data(), Concentration.GetData(),
                   &option_photolysis_tabulation);
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(Concentration);
#endif
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYRADM_CXX
#endif
