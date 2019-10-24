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

// This file is part of a Gaussian plume model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_AER_CXX


#include "GaussianPlume_aer.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T>
  GaussianPlume_aer<T>::GaussianPlume_aer(string config_file):
    GaussianPlume<T>(config_file)
  {
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
  void GaussianPlume_aer<T>::ReadConfiguration()
  {
    GaussianPlume<T>::ReadConfiguration();

    /*** Configuration ***/

    // File containing the aerosol diameters.
    this->config.SetSection("[domain]");
    this->config.PeekValue("File_diameter", file_diameter);

    // File containing the source data for aerosol species.
    this->config.SetSection("[gaussian]");
    this->config.PeekValue("File_source_aer", file_source_aer);

    // Gets aerosol species names.
    ConfigStream species_stream(this->file_species);
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());

    // Gets the diameters of aerosol species.
    ConfigStream diameter_stream(file_diameter);
    diameter_stream.SetSection("[diameter]");
    while (!diameter_stream.IsEmpty())
      diameter_list.push_back(to_num<T>(diameter_stream.GetElement()));
    this->Nbin_aer = int(diameter_list.size());

    /*** Aerosol species data ***/

    // Reads the aerosol species half-lives in the species file.
    if (this->option_radioactive_decay)
      {
        half_life_time_aer.resize(this->Ns_aer);
        species_stream.SetSection("[half_life]");
        for (int i = 0; i < this->Ns_aer; i++)
          species_stream.PeekValue(this->species_list_aer[i],
                                   half_life_time_aer(i));

        // Conversion from days to seconds.
        T factor = 24. * 3600.;
        for (int i = 0; i < this->Ns_aer; i++)
          half_life_time_aer(i) *= factor;
      }

    // Reads the biological half-lives of aerosol species in the species file.
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
  void GaussianPlume_aer<T>::Allocate()
  {
    GaussianPlume<T>::Allocate();

    // Grids.
    GridS5D_aer = RegularGrid<T>(this->Ns_aer);
    GridB5D_aer = RegularGrid<T>(this->Nbin_aer);
    GridX5D = RegularGrid<T>(this->x_min, this->Delta_x, this->Nx);
    GridY5D = RegularGrid<T>(this->y_min, this->Delta_y, this->Ny);
    GridZ5D = RegularGrid<T>(this->Nz);

    // Aerosol concentrations.
    Concentration_aer.Resize(GridS5D_aer, GridB5D_aer, GridZ5D,
                             GridY5D, GridX5D);
  }


  //! Model initialization.
  /*! It reads the configuration.
    \param read_all_input Should all input data and options be read in the
    configuration files? If so, even if unnecessary data will be read,
    e.g. the Monin-Obukhov length while the Doury parameterization is
    used. This enables to switch et a later time to other parameterizations
    (e.g., to "similarity_theory").
  */
  template<class T>
  void GaussianPlume_aer<T>::Init(bool read_all_input)
  {
    /*** Model initialization for gaseous species ***/

    GaussianPlume<T>::Init(read_all_input);

    /*** Sources initialization for aerosol species ***/

    // Vertical layers heights.
    for (int k = 0; k < this->Nz; k++)
      GridZ5D(k) = this->GridZ4D(k);

    InitSource_aer(file_source_aer);
  }


  //! Sources initialization.
  /*! It sets all aerosol sources from a text file. Each source is described
    in a dedicated section "[source]" in which one finds the following
    entries:
    <ul>
    <li> Rate: the rate (mass unit per second),
    <li> Velocity: the efflux speed (m/s),
    <li> Temperature: the temperature of emissions (Celsius degrees),
    <li> Abscissa: abscissa (m),
    <li> Ordinate: ordinate (m),
    <li> Altitude: height (m),
    <li> Species_name: name of aerosol species.
    </ul>
    \param source_file_aer file that describes the sources for aerosol
    species.
  */
  template<class T>
  void GaussianPlume_aer<T>::InitSource_aer(string source_file_aer)
  {
    T rate, velocity, temperature, diameter, x, y, z;
    string line, species_name;

    ConfigStream source_stream(source_file_aer);

    SourceList_aer.clear();
    Nsource_aer = 0;
    while (!source_stream.IsEmpty())
      {
        line = source_stream.GetLine();
        if (split(line)[0] == "[aerosol_source]")
          {
            source_stream.PeekValue("Abscissa", x);
            source_stream.PeekValue("Ordinate", y);
            source_stream.PeekValue("Altitude", "positive", z);
            source_stream.PeekValue("Rate", "positive", rate);
            source_stream.PeekValue("Velocity", velocity);
            source_stream.PeekValue("Temperature", temperature);
            source_stream.PeekValue("Diameter", diameter);
            source_stream.PeekValue("Species_name", species_name);
            temperature += 273.15;

            SourceList_aer.push_back
              (PlumeSource<T>(rate, velocity, temperature, diameter, x, y,
                              z, this->GetSpeciesIndex_aer(species_name)));
            Nsource_aer++;

          }
      }
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
    <li> Inversion_height: the inversion height (m).
    <li> Stability: stability class in [A, F].
    </ul>
    \param meteo ConfigStream instance through which all entries (Temperature,
    Wind_angle, ...) may be read to set the meteorological situation.
    \param show_meteo indicates whether the meteorological data is to be
    displayed on screen.
  */
  template<class T>
  void GaussianPlume_aer<T>::InitMeteo(ConfigStream& meteo, bool show_meteo)
  {
    /*** Main meteorological variables and gaseous species data ***/

    GaussianPlume<T>::InitMeteo(meteo, show_meteo);

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


  //////////////////
  // COMPUTATIONS //
  //////////////////


  //! Computations performed before any concentration can be computed.
  /*!  After this method is called, one may call 'GetConcentration' or
    'Compute'. But never before!
  */
  template<class T>
  void GaussianPlume_aer<T>::InitCompute()
  {
    GaussianPlume<T>::InitCompute();
    if (this->option_plume_rise)
      ComputePlumeRise_aer();
  }


  //! Computes plume rise for aerosol species.
  template<class T>
  void GaussianPlume_aer<T>::ComputePlumeRise_aer()
  {
    // Loop over all sources.
    for (typename list<PlumeSource<T> >::iterator iter
           = this->SourceList_aer.begin();
         iter != this->SourceList_aer.end(); iter++)
      this->ComputeSourcePlumeRise(*iter);
  }


  //! Computes concentrations in the whole domain.
  template<class T>
  void GaussianPlume_aer<T>::Compute()
  {
    GaussianPlume<T>::Compute();

    int i, j, k, species_aer;

    for (species_aer = 0; species_aer < this->Ns_aer; species_aer++)
      for (int diam = 0; diam < this->Nbin_aer; diam++)
        for (k = 0; k < this->Nz; k++)
          for (j = 0; j < this->Ny; j++)
            for (i = 0; i < this->Nx; i++)
              Concentration_aer(species_aer, diam, k, j, i)
                = GetConcentration_aer(species_aer, diam, GridZ5D(k),
                                       GridY5D(j), GridX5D(i));
  }


  /*! \brief Computes concentration of a given aerosol species and a given
    aerosol diameter at a given point.
  */
  /*!
    \param species index of aerosol species.
    \param diam index of aerosol diameter.
    \param x abscissa (m).
    \param y ordinate (m).
    \param z height (m).
    \return The value of concentration at the given point (mass/m^3).
  */
  template<class T>
  T GaussianPlume_aer<T>::GetConcentration_aer(int species, int diameter,
                                               T z, T y, T x)
  {
    const T pi = 3.14159265358979323846264;

    // Distances downwind and crosswind from the source (m).
    T distance_x, distance_y, time;
    // Horizontal and vertical standard deviations (m).
    T sigma_y, sigma_z;
    // Effective height of release (m).
    T effective_height;

    // Output concentration.
    T concentration = 0.;

    // Loop over all sources for aerosol species.
    for (typename list<PlumeSource<T> >::const_iterator iter
           = this->SourceList_aer.begin();
         iter != this->SourceList_aer.end(); iter++)
      if (iter->GetSpeciesIndex() == species)
        {
          // Downwind distance from source.
          distance_x = (x - iter->GetX()) * this->cos_angle_
            + (y - iter->GetY()) * this->sin_angle_;

          // Distance from downwind axis.
          distance_y = (iter->GetX() - x) * this->sin_angle_
            + (y - iter->GetY()) * this->cos_angle_;
          distance_y = abs(distance_y);

          // If there is anything to compute.
          if (distance_x > 0)
            {
              effective_height = iter->GetHeight();

              // Computing standard deviations.
              time = distance_x / this->wind_;
              this->ComputeSigma(distance_x, time, effective_height,
                                 sigma_y, sigma_z);


              // Computing loss factors.
              T loss_factor = 1.;
              T overcamp_factor = 1.;
              ComputeLossFactor_aer(distance_x, effective_height,
                                    species, diameter,
                                    loss_factor, overcamp_factor);

              // Minimum volume.
              T minimum_volume = pi * iter->GetDiameter()
                * iter->GetDiameter() * iter->GetVelocity();

              concentration += loss_factor
                * ComputePlumeConcentration(this->wind_,
                                            this->inversion_height_,
                                            effective_height, iter->GetRate(),
                                            distance_y, z, sigma_y, sigma_z,
                                            overcamp_factor, minimum_volume);
            }
        }

    return concentration;
  }


  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  Data<T, 5>& GaussianPlume_aer<T>::GetConcentration_aer()
  {
    return Concentration_aer;
  }


  //!  Computes the loss factor.
  template<class T>
  void GaussianPlume_aer<T>::ComputeLossFactor_aer(T distance, T z,
                                                   int species, int diam,
                                                   T& loss_factor,
                                                   T& overcamp_factor)
  {
    T rad, bio, scav, dep;
    T transfer_time = distance / this->wind_;
    if (this->option_radioactive_decay)
      rad = half_life_time_aer(species);
    else
      rad = 0.;
    if (this->option_biological_decay)
      bio = biological_half_life_time_aer(species);
    else
      bio = 0.;
    if (this->option_scavenging)
      scav = scavenging_coefficient_aer(species, diam);
    else
      scav = 0.;
    if (this->option_dry_deposition)
      dep = deposition_velocity_aer(species, diam);
    else
      dep = 0.;

    this->ComputeLossFactor(distance, transfer_time, z, rad, bio, scav, dep,
                            loss_factor, overcamp_factor);
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPLUME_AER_CXX
#endif
