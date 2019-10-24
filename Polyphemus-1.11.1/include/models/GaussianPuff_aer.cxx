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

// This file is part of a Gaussian puff model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_AER_CXX


#include "GaussianPuff_aer.hxx"


namespace Polyphemus
{


  /////////////////
  // CONSTRUCTOR //
  /////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T>
  GaussianPuff_aer<T>::GaussianPuff_aer(string config_file):
    GaussianPuff<T>(config_file), Npuff_aer(0)
  {
  }


  //! Destructor.
  template<class T>
  GaussianPuff_aer<T>::~GaussianPuff_aer()
  {
    ClearPuffList_aer();
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T>
  void GaussianPuff_aer<T>::ReadConfiguration()
  {
    GaussianPuff<T>::ReadConfiguration();

    // File containing the aerosol diameters.
    this->config.SetSection("[domain]");
    this->config.PeekValue("File_diameter", file_diameter);

    // File containing the source data for aerosol species.
    this->config.SetSection("[gaussian]");
    this->config.PeekValue("File_puff_aer", file_puff_aer);

    // Gets aerosol species names.
    ConfigStream species_stream(this->file_species);
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());

    // Gets the aerosol diameters.
    ConfigStream diameter_stream(file_diameter);
    diameter_stream.SetSection("[diameter]");
    while (!diameter_stream.IsEmpty())
      diameter_list.push_back(to_num<T>(diameter_stream.GetElement()));
    this->Nbin_aer = int(diameter_list.size());

    /*** Aerosol species data ***/

    // Reads the aerosol half-lives in the species file.
    if (this->option_radioactive_decay)
      {
        half_life_time_aer.resize(this->Ns_aer);
        for (int i = 0; i < this->Ns_aer; i++)
          {
            species_stream.SetSection("[half_life]");
            species_stream.PeekValue(this->species_list_aer[i],
                                     half_life_time_aer(i));
          }

        // Conversion from days to seconds.
        T factor = 24. * 3600.;
        for (int i = 0; i < this->Ns_aer; i++)
          half_life_time_aer(i) *= factor;
      }

    // Reads the aerosol biological half-lives in the species file.
    if (this->option_biological_decay)
      {
        T day_value, night_value;
        biological_half_life_time_aer.resize(this->Ns_aer);
        for (int i = 0; i < this->Ns_aer; i++)
          {
            species_stream.SetSection("[half_life_time]");
            species_stream.Find(this->species_list_aer[i]);
            species_stream.GetNumber(day_value);
            species_stream.GetNumber(night_value);
            if (this->option_day)
              biological_half_life_time_aer(i) = day_value;
            else
              biological_half_life_time_aer(i) = night_value;
          }
      }
  }


  //! Allocates memory.
  /*! Allocates the grids and the concentration Data for aerosol species.
   */
  template<class T>
  void GaussianPuff_aer<T>::Allocate()
  {
    GaussianPuff<T>::Allocate();

    // Grids.
    GridS5D_aer = RegularGrid<T>(this->Ns_aer);
    GridB5D_aer = RegularGrid<T>(this->Nbin_aer);
    GridX5D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);
    GridY5D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridZ5D = RegularGrid<T>(this->Nz);

    // Concentrations of aerosol species.
    Concentration_aer.Resize(GridS5D_aer, GridB5D_aer, GridZ5D,
                             GridY5D, GridX5D);
  }


  //! Clears the puff list.
  /*!
   */
  template<class T>
  void GaussianPuff_aer<T>::ClearPuffList_aer()
  {
    int puff, diam;
    for (puff = 0; puff < Npuff_aer; puff++)
      for (diam = 0; diam < this->Nbin_aer; diam++)
        delete PuffList_aer[puff][diam];

    PuffList_aer.clear();
    Npuff_aer = 0;
  }


  //! Model initialization.
  /*! It reads the configuration.
   */
  template<class T>
  void GaussianPuff_aer<T>::Init()
  {
    /*** Model initialization for gaseous species ***/

    GaussianPuff<T>::Init();

    // Vertical layers heights.
    for (int k = 0; k < this->Nz; k++)
      GridZ5D(k) = this->GridZ4D(k);

    /*** Puffs initialization for aerosol species ***/

    InitPuffSource_aer(file_puff_aer);
  }


  //! Puff sources initialization.
  /*! It sets all aerosol sources from a text file. Each source is described
    in a dedicated section "[source]" in which one finds the following
    entries:
    <ul>
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species_name: name of aerosol species.
    <li> Release_time: the release time (s),
    <li> Quantity: the released quantity (mass unit),
    </ul>
    \param file_puff_aer file that describes the sources.
  */
  template<class T>
  void GaussianPuff_aer<T>::InitPuffSource_aer(string file_puff_aer)
  {
    // Puffs definitions.
    T time_puff, quantity, abscissa, ordinate, height;
    T velocity, temperature, diameter;
    string line, species;
    // auto-generated source ids for now, this class will be superseded
    // with the GaussianPuffAerosol implementation.
    int source_id = 1;

    const T pi = 3.14159265358979323846264;

    ConfigStream puff_stream(file_puff_aer);

    while (!puff_stream.IsEmpty())
      {
        line = puff_stream.GetLine();
        if (split(line)[0] == "[aerosol_source]")
          {
            puff_stream.PeekValue("Abscissa", abscissa);
            puff_stream.PeekValue("Ordinate", ordinate);
            puff_stream.PeekValue("Altitude", "positive", height);
            puff_stream.PeekValue("Release_time", time_puff);
            puff_stream.PeekValue("Quantity", "positive", quantity);
            puff_stream.PeekValue("Species_name", species);
            puff_stream.PeekValue("Velocity", velocity);
            puff_stream.PeekValue("Temperature", temperature);
            puff_stream.PeekValue("Diameter", "positive", diameter);
            temperature += 273.15;

            vector< Puff<T>* > temp;
            for (int j = 0; j < this->Nbin_aer; j++)
              {
                Puff<T>* puff = new Puff<T>(time_puff, velocity,
                                            temperature, diameter,
                                            quantity, abscissa,
                                            ordinate, height,
                                            this
                                            ->GetSpeciesIndex_aer(species),
                                            to_str(source_id++));
                temp.push_back(puff);
              }
            PuffList_aer.push_back(temp);
          }
      }
    Npuff_aer = int(PuffList_aer.size());

  }


  //! Initializes meteorological conditions.
  /*! It sets the meteorological data from a configuration file.  The
    situation is described in a dedicated section "[situation]" in which one
    finds the following entries:
    <ul>
    <li> Temperature: the temperature of ambient air (Celsius degrees).
    <li> Wind_angle: the angle between positive x-axis and wind,
    counterclockwise (degrees).
    <li> Wind: the wind velocity (m/s).
    <li> Stability: stability class in [A, F].
    </ul>
    \param meteo ConfigStream instance through which all entries (Temperature,
    Wind_angle, ...) may be read to set the meteorological situation.
    \param show_meteo indicates whether the meteorological data is to be
    displayed on screen.
  */
  template<class T>
  void GaussianPuff_aer<T>::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {

    /*** Main meteorological variables and gaseous species data ***/

    GaussianPuff<T>::InitMeteo(meteo, show_meteo);

    /*** Initialization of puffs positions (aerosol species) ***/

    for (int i = 0; i < Npuff_aer; i++)
      for (int j = 0; j < this->Nbin_aer; j++)
        PuffList_aer[i][j]->InitPuff();

    /*** Scavenging coefficients for aerosol species ***/

    if (this->option_scavenging)
      {
        scavenging_coefficient_aer.resize(this->Ns_aer, this->Nbin_aer);
        meteo.Find("Scavenging_coefficient_aer");
        vector<string> scav_coefficient = split(meteo.GetLine());
        int Nscav = int(scav_coefficient.size());

        for (int  i = 0; i < this->Ns_aer; i++)
          for (int k = 0; k < this->Nbin_aer; k++)
            {
              T lambda = 0;
              for (int j = 0; j < Nscav - 1; j++)
                if (scav_coefficient[j] == this->species_list_aer[i]
                    + string("_") + to_str(k))
                  lambda = to_num<T>(scav_coefficient[j + 1]);
              scavenging_coefficient_aer(i, k) = lambda;
            }
      }

    /*** Deposition velocities for aerosol species ***/

    if (this->option_dry_deposition)
      {
        deposition_velocity_aer.resize(this->Ns_aer, this->Nbin_aer);
        meteo.Find("Deposition_velocity_aer");
        vector<string> dep_velocity = split(meteo.GetLine());
        int Ndep = int(dep_velocity.size());

        for (int  i = 0; i < this->Ns_aer; i++)
          for (int k = 0; k < this->Nbin_aer; k++)
            {
              T vd = 0;
              for (int j = 0; j < Ndep - 1; j++)
                if (dep_velocity[j] == this->species_list_aer[i]
                    + string("_") + to_str(k))
                  vd = to_num<T>(dep_velocity[j + 1]);
              deposition_velocity_aer(i, k) = vd;
            }
      }
  }


  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */ template<class T>
  void GaussianPuff_aer<T>::InitStep()
  {
    GaussianPuff<T>::InitStep();
  }


  //! Computes the new position of the puff center after advection.
  template<class T>
  void GaussianPuff_aer<T>::Advection_aer()
  {
    for (int i = 0; i < Npuff_aer; i++)
      for (int j = 0; j < this->Nbin_aer; j++)
        if (this->current_time >= PuffList_aer[i][j]->GetReleaseTime())
          {
            this->SetCurrentMeteo(PuffList_aer[i][j]);
            this->Advection(PuffList_aer[i][j]);
          }
  }


  //!  Computes the horizontal and vertical standard deviations.
  template<class T>
  void GaussianPuff_aer<T>::Diffusion_aer()
  {
    for (int i = 0; i < Npuff_aer; i++)
      for (int j = 0; j < this->Nbin_aer; j++)
        if (this->current_time >= PuffList_aer[i][j]->GetReleaseTime())
          {
            this->SetCurrentMeteo(PuffList_aer[i][j]);
            this->Diffusion(PuffList_aer[i][j]);
          }
  }


  //! Performs one step forward.
  template<class T>
  void GaussianPuff_aer<T>::Forward()
  {
    Advection_aer();
    Diffusion_aer();

    int puff, diam, i, j, k;
    Concentration_aer.SetZero();
    for (puff = 0; puff < Npuff_aer; puff++)
      for (diam = 0; diam < this->Nbin_aer; diam++)
        if (this->current_time >= PuffList_aer[puff][diam]->GetReleaseTime())
          {
            this->SetCurrentMeteo(PuffList_aer[puff][diam]);
            // Computing loss factors.
            ComputeLossFactor_aer(puff, diam);
            int s = PuffList_aer[puff][diam]->GetSpeciesIndex();

            // Loop on all points to compute concentration.
            for (k = 0; k < this->Nz; k++)
              for (j = 0; j < this->Ny; j++)
                for (i = 0; i < this->Nx; i++)
                  Concentration_aer(s, diam, k, j, i) +=
                    this->ComputePuffConcentration(PuffList_aer[puff][diam],
                                                   s, GridX5D(i), GridY5D(j),
                                                   GridZ5D(k));
          }
    GaussianPuff<T>::Forward();
  }

  //! Computes concentration at a given point for aerosols.
  /*!
    \return The concentration at the point.
  */
  template<class T>
  T GaussianPuff_aer<T>::GetConcentration_aer(int species, int diameter,
                                              T z, T y, T x)
  {
    T concentration = 0.;
    this->SubtractTime(this->Delta_t);
    for (int puff = 0; puff < Npuff_aer; puff++)
      for (int diam = 0; diam < this->Nbin_aer; diam++)
        if (this->current_time >= PuffList_aer[puff][diam]->GetReleaseTime())
          {
            this->SetCurrentMeteo(PuffList_aer[puff][diam]);
            concentration += this->
              ComputePuffConcentration(PuffList_aer[puff][diam], species,
                                       x, y, z);
          }
    this->AddTime(this->Delta_t);
    return concentration;
  }



  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  Data<T, 5>& GaussianPuff_aer<T>::GetConcentration_aer()
  {
    return Concentration_aer;
  }


  //!  Computes the loss factor.
  template<class T>
  void GaussianPuff_aer<T>::ComputeLossFactor_aer(int index, int diam)
  {
    Puff<T>* puff = PuffList_aer[index][diam];
    T loss_factor = 1.;
    T overcamp_factor = 1.;
    int s = puff->GetSpeciesIndex();
    T distance = puff->GetDistance();
    T transfer_time = puff->GetPuffTime();
    T quantity = puff->GetQuantity(s);
    T z = puff->GetZ();

    T rad, bio, scav, dep;

    if (this->option_radioactive_decay)
      rad = half_life_time_aer(s);
    else
      rad = 0.;
    if (this->option_biological_decay)
      bio = biological_half_life_time_aer(s);
    else
      bio = 0.;

    if (puff->HasMeteo())
      {
        if (this->option_scavenging)
          scav = puff->GetScavengingCoefficient(s);
        else
          scav = 0.;
        if (this->option_dry_deposition)
          dep = puff->GetDepositionVelocity(s);
        else
          dep = 0.;
      }
    else
      {
        if (this->option_scavenging)
          scav = scavenging_coefficient_aer(s, diam);
        else
          scav = 0.;
        if (this->option_dry_deposition)
          dep = deposition_velocity_aer(s, diam);
        else
          dep = 0.;
      }

    this->ComputeLossFactor(distance, transfer_time, z, rad, bio, scav, dep,
                            loss_factor, overcamp_factor);
    puff->SetQuantity(loss_factor * quantity, s);
    puff->SetReflectionFactor(overcamp_factor, s);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFF_AER_CXX
#endif
