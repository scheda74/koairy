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

// This file is part of the Eulerian model Castor.

// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).


#ifndef POLYPHEMUS_FILE_MODELS_CASTORCHEMISTRY_CXX


#include "CastorChemistry.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  CastorChemistry<T, ClassTransport, ClassChemistry>
  ::CastorChemistry(string config_file):
    CastorTransport<T, ClassTransport>(config_file)
  {

    /*** Managed data ***/

    this->option_manage["specific_humidity"] = true;
    this->option_manage["photolysis_rate"] = true;
    this->option_manage["attenuation"] = true;
    this->option_manage["liquid_water_content"] = true;

    /*** Pointers to 3D data ***/

    this->D3_map["SpecificHumidity"] = &SpecificHumidity;
    this->D3_map["SpecificHumidity_i"] = &SpecificHumidity;

    this->D2_map["Attenuation"] = &Attenuation;
    this->D2_map["Attenuation_i"] = &Attenuation;

    this->D3_map["LiquidWaterContent"] = &LiquidWaterContent;
    this->D3_map["LiquidWaterContent_i"] = &LiquidWaterContent;
  }


  //! Destructor.
  template<class T, class ClassTransport, class ClassChemistry>
  CastorChemistry<T, ClassTransport, ClassChemistry>
  ::~CastorChemistry()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::ReadConfiguration()
  {
    CastorTransport<T, ClassTransport>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[options]");
    this->config.PeekValue("With_chemistry",
                           this->option_process["with_chemistry"]);

    /*** Output units ***/

    this->config.SetSection("[output]");
    this->config.PeekValue("Unit", "mass | number | ratio", unit);

    if (unit == "mass")
      {
        // Molecular weight.
        ConfigStream config(this->GetSpeciesFile());
        config.SetSection("[molecular_weight]");
        for (int i = 0; i < this->Ns; i++)
          config.PeekValue(this->species_list[i],
                           molecular_weight[this->species_list[i]]);
      }
  }


  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::CheckConfiguration()
  {
    CastorTransport<T, ClassTransport>::CheckConfiguration();

    if (this->option_manage["specific_humidity"]
        && this->input_files["meteo"]("SpecificHumidity").empty())
      throw string("Specific humidity is needed but no input data file was")
        + " provided.";
    if (this->option_manage["attenuation"]
        && this->input_files["meteo"]("Attenuation").empty())
      throw "Attenuation is needed but no input data file was provided.";
    if (this->option_manage["liquid_water_content"]
        && this->input_files["meteo"]("LiquidWaterContent").empty())
      throw string("Liquid water content is needed but no")
        + " input data file was provided.";

    if (unit != "mass" && unit != "number" && unit != "ratio")
      throw string("Output unit \"") + unit + "\" is not supported.";
  }


  //! Return the index of a species among extended species.
  /*!
    \param species species name.
    \return The species index in extended species.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  int CastorChemistry<T, ClassTransport, ClassChemistry>
  ::GetSpeciesIndex_ext(string species) const
  {
    if (species == "M")
      return this->Ns;
    else if (species == "O2")
      return this->Ns + 1;
    else if (species == "N2")
      return this->Ns + 2;
    else if (species == "H2O")
      return this->Ns + 3;
    else
      return this->GetSpeciesIndex(species);
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::Allocate()
  {
    CastorTransport<T, ClassTransport>::Allocate();

    /*** Specific humidity ***/

    SpecificHumidity.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileSpecificHumidity_i.Resize(this->GridZ3D, this->GridY3D,
                                  this->GridX3D);
    FileSpecificHumidity_f.Resize(this->GridZ3D, this->GridY3D,
                                  this->GridX3D);

    /*** Attenuation ***/

    Attenuation.Resize(this->GridY2D, this->GridX2D);
    FileAttenuation_i.Resize(this->GridY2D, this->GridX2D);
    FileAttenuation_f.Resize(this->GridY2D, this->GridX2D);

    /*** Liquid water content ***/

    LiquidWaterContent.Resize(this->GridZ3D, this->GridY3D, this->GridX3D);
    FileLiquidWaterContent_i.Resize(this->GridZ3D, this->GridY3D,
                                    this->GridX3D);
    FileLiquidWaterContent_f.Resize(this->GridZ3D, this->GridY3D,
                                    this->GridX3D);
  }


  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::Init()
  {
    CastorTransport<T, ClassTransport>::Init();

    /*** Input data ***/

    CastorChemistry<T, ClassTransport, ClassChemistry>::InitAllData();

    /*** Chemical mechanism ***/

    if (this->option_process["with_chemistry"])
      Chemistry_.Init(*this);

    // Unit conversion.
    ConvertExternalUnit(this->Concentration);
  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::InitStep()
  {
    // Unit conversion.
    ConvertInternalUnit(this->Concentration);

    /*** Transport ***/

    CastorTransport<T, ClassTransport>::InitStep();

    /*** Specific humidity ***/

    if (this->option_manage["specific_humidity"])
      this->InitData("meteo", "SpecificHumidity", FileSpecificHumidity_i,
                     FileSpecificHumidity_f, this->intermediate_date,
                     SpecificHumidity);

    /*** Attenuation ***/

    if (this->option_manage["attenuation"])
      this->InitData("meteo", "Attenuation", FileAttenuation_i,
                     FileAttenuation_f, this->intermediate_date,
                     Attenuation);

    /*** Liquid water content ***/

    if (this->option_manage["liquid_water_content"])
      this->InitData("meteo", "LiquidWaterContent", FileLiquidWaterContent_i,
                     FileLiquidWaterContent_f, this->intermediate_date,
                     LiquidWaterContent);

    /*** Chemical mechanism ***/

    if (this->option_process["with_chemistry"])
      Chemistry_.InitStep(*this);

    // Unit conversion.
    ConvertExternalUnit(this->Concentration);
  }


  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::SetDate(Date date)
  {
    CastorTransport<T, ClassTransport>::SetDate(date);
    CastorChemistry<T, ClassTransport, ClassChemistry>::InitAllData();
  }


  /////////////////
  // INTEGRATION //
  /////////////////


  //! Performs one step forward.
  /*! It performs one advection step, then one diffusion step and finally
    integrates chemistry. The first two steps are split (operator
    splitting). The last step (chemistry) may be split or partially coupled
    through source splitting.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::Forward()
  {
    int s, h, k, j, i;

    // Unit conversion.
    ConvertInternalUnit(this->Concentration);

    Data<T, 4> Concentration_tmp(this->GridS4D, this->GridZ4D,
                                 this->GridY4D, this->GridX4D);

    T conc;
    for (s = 0; s < this->Ns; s++)
      for (k = 0; k < this->Nz; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            {
              Concentration_tmp(s, k, j, i)
                = (4. *  this->Concentration(s, k, j, i)
                   -  this->PreviousConcentration(s, k, j, i)) / 3.;
              conc = this->Concentration(s, k, j, i);
              this->Concentration(s, k, j, i) = 2. * conc
                - this->PreviousConcentration(s, k, j, i);
              this->PreviousConcentration(s, k, j, i) = conc;
            }
    this->Concentration.ThresholdMin(1.);
    Concentration_tmp.ThresholdMin(1.);

    int Niter = 1;
    // During the first hour (spin-up).
    if (T(this->step) * this->Delta_t < 3600.)
      Niter = 5;

    // To speed up tests.
    bool with_chemistry(this->option_process["with_chemistry"]);
    bool deposition_velocity(this->option_manage["deposition_velocity"]);
    bool with_transport(this->option_process["with_transport"]);

    T loss, production;
    T factor = 2. / 3. * this->Delta_t;
    for (h = 0; h < Niter; h++)
      for (k = 0; k < this->Nz; k++)
        for (j = 0; j < this->Ny; j++)
          for (i = 0; i < this->Nx; i++)
            for (s = 0; s < this->Ns; s++)
              {
                loss = 0.;
                production = 0.;

                // Chemistry.
                if (with_chemistry)
                  this->Chemistry_.LossProduction(*this, s, k, j, i,
                                                  loss, production);

                // Emissions.
                if (k < this->Nz_vol_emis && this->HasVolumeEmission(s))
                  {
                    int e = this->VolumeEmissionIndex(s);
                    production += this->VolumeEmission(e, k, j, i) / 100.
                      / (this->Altitude(k + 1, j, i)
                         - this->Altitude(k, j, i));
                  }

                // Deposition.
                if (deposition_velocity && k == 0
                    && this->HasDepositionVelocity(s))
                  {
                    int e = this->DepositionVelocityIndex(s);
                    loss += this->DepositionVelocity(e, j, i)
                      * this->Concentration(s, k, j, i);
                  }

                if (with_transport)
                  this->Transport_.LossProduction(*this, s, k, j, i,
                                                  loss, production);

                this->Concentration(s, k, j, i) =
                  (Concentration_tmp(s, k, j, i) + factor * production)
                  / (1. + factor * loss
                     / this->Concentration(s, k, j, i));
              }

    this->AddTime(this->Delta_t);
    this->step++;

    // Unit conversion.
    ConvertExternalUnit(this->Concentration);
  }


  //! Converts the concentrations Data into the output unit.
  /*!
    \param Concentration_ Concentrations.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::ConvertExternalUnit(Data<T, 4>& Concentration_)
  {
    const T avogadro = 6.02213e23;

    // From molecules / cm^3 to ppbv.
    if (unit == "ratio")
      for (int s = 0; s < this->Ns; s++)
        for (int k = 0; k < this->Nz; k++)
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              Concentration_(s, k, j, i) *= 1.e9 / this->AirDensity(k, j, i);

    // From molecules / cm^3 to microgramms / m^3.
    else if (unit == "mass")
      for (int s = 0; s < this->Ns; s++)
        for (int k = 0; k < this->Nz; k++)
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              Concentration_(s, k, j, i) *=
                1.e12 / avogadro * molecular_weight[this->species_list[s]];
  }


  //! Converts the concentrations Data into the internal unit.
  /*!
    \param Concentration_ Concentrations.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::ConvertInternalUnit(Data<T, 4>& Concentration_)
  {
    const T avogadro = 6.02213e23;

    // From ppbv to molecules / cm^3.
    if (unit == "ratio")
      for (int s = 0; s < this->Ns; s++)
        for (int k = 0; k < this->Nz; k++)
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              Concentration_(s, k, j, i) *= this->AirDensity(k, j, i) / 1.e9;

    // From microgramms / m^3 to molecules / cm^3.
    else if (unit == "mass")
      for (int s = 0; s < this->Ns; s++)
        for (int k = 0; k < this->Nz; k++)
          for (int j = 0; j < this->Ny; j++)
            for (int i = 0; i < this->Nx; i++)
              Concentration_(s, k, j, i) *=
                avogadro / 1.e12 / molecular_weight[this->species_list[s]];
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Moves model input-data to the current date.
  /*! This method prepares the model for a time integration from the current
    date. It reads input data to related to chemistry be read before InitStep
    and Forward.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  void CastorChemistry<T, ClassTransport, ClassChemistry>
  ::InitAllData()
  {

    /*** Specific humidity ***/

    if (this->option_manage["specific_humidity"])
      this->InitData("meteo", "SpecificHumidity",
                     FileSpecificHumidity_i, FileSpecificHumidity_f,
                     this->intermediate_date, SpecificHumidity);

    /*** Attenuation ***/

    if (this->option_manage["attenuation"])
      this->InitData("meteo", "Attenuation", FileAttenuation_i,
                     FileAttenuation_f, this->intermediate_date,
                     Attenuation);

    /*** Liquid water content ***/

    if (this->option_manage["liquid_water_content"])
      this->InitData("meteo", "LiquidWaterContent", FileLiquidWaterContent_i,
                     FileLiquidWaterContent_f, this->intermediate_date,
                     LiquidWaterContent);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CASTORCHEMISTRY_CXX
#endif
