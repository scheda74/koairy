// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi, Vivien Mallet, Denis Quélo
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_DECAY_CXX


#include "Decay.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  Decay<T>::Decay()
  {
  }


  //! Initialization of the scheme.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConfigurationFile()
    <li> GetSpeciesFile()
    <li> GetNs()
    <li> GetSpeciesList()
    <li> GetNs_aer()
    <li> GetSpeciesList_aer()
    <li> GetNbin_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Decay<T>::Init(ClassModel& Model)
  {

    /*** Main options ***/

    ConfigStream general_config(Model.GetConfigurationFile());

    general_config.SetSection("[options]");

    general_config.PeekValue("With_filiation_matrix", with_filiation_matrix);
    general_config.PeekValue("With_time_dependence", with_time_dependence);

    // Throw an exception if a filiation matrix and time dependence are used
    // at the same time (incompatible models).
    if (with_filiation_matrix && with_time_dependence)
      throw string("Incompatible options: use of filiation matrix together ")
        + "with time-dependent half-lives.";


    /*** Gas ***/

    Ns = Model.GetNs();
    species_list = Model.GetSpeciesList();

    ConfigStream species_stream(Model.GetSpeciesFile());

    if (Ns != 0)
      if (!with_filiation_matrix && !with_time_dependence)
        {
          species_stream.SetSection("[half_life]");
          half_life.resize(Ns);
          for (int i = 0; i < Ns; i++)
            species_stream.PeekValue(species_list[i], half_life(i));
        }
      else if (with_time_dependence)
        {
          half_life_day.resize(Ns);
          half_life_night.resize(Ns);
          for (int i = 0; i < Ns; i++)
            {
              species_stream.SetSection("[half_life_time]");
              species_stream.Find(species_list[i]);
              species_stream.GetNumber(half_life_day(i));
              species_stream.GetNumber(half_life_night(i));
            }
        }
      else if (with_filiation_matrix)
        {
          species_stream.SetSection("[filiation_matrix]");
          ConfigStream matrix_stream(species_stream.PeekValue("File"));
          filiation_matrix.resize(Ns, Ns);
          for (int i = 0; i < Ns; i++)
            for (int j = 0; j < Ns; j++)
              filiation_matrix(i, j) = T(matrix_stream.GetNumber());
        }


    /*** Aerosol ***/

    Ns_aer = Model.GetNs_aer();
    species_list_aer = Model.GetSpeciesList_aer();
    Nbin_aer = Model.GetNbin_aer();

    if (Ns_aer != 0)
      if (!with_filiation_matrix && !with_time_dependence)
        {
          species_stream.SetSection("[half_life_aerosol]");
          half_life_aer.resize(Ns_aer);
          for (int i = 0; i < Ns_aer; i++)
            species_stream.PeekValue(species_list_aer[i], half_life_aer(i));
        }
      else if (with_time_dependence)
        {
          half_life_day_aer.resize(Ns_aer);
          half_life_night_aer.resize(Ns_aer);
          for (int i = 0; i < Ns_aer; i++)
            {
              species_stream.SetSection("[half_life_time_aerosol]");
              species_stream.Find(species_list_aer[i]);
              species_stream.GetNumber(half_life_day_aer(i));
              species_stream.GetNumber(half_life_night_aer(i));
            }
        }
      else if (with_filiation_matrix)
        {
          species_stream.SetSection("[filiation_matrix_aerosol]");
          ConfigStream matrix_stream(species_stream.PeekValue("File"));
          filiation_matrix_aer.resize(Ns_aer, Ns_aer);
          for (int i = 0; i < Ns_aer; i++)
            for (int j = 0; j < Ns_aer; j++)
              filiation_matrix_aer(i, j) = T(matrix_stream.GetNumber());
        }
  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetDelta_t()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> GetCurrentDate()
    <li> GetSource_i()
    <li> GetNs_source()
    <li> GetNz_source()
    <li> SourceGlobalIndex
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Decay<T>::Forward(ClassModel& Model)
  {
    int k, j, i, s, r;
    T decay;

    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    /*** Sources ***/

    int Ns_source = Model.GetNs_source();
    int Nz_source = Model.GetNz_source();

    for (s = 0; s < Ns_source; s++)
      {
        int gs = Model.SourceGlobalIndex(s);
        for (k = 0; k < Nz_source; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              Model.GetConcentration()(gs, k, j, i) += Model.GetDelta_t()
                * Model.GetSource_i()(s, k, j, i);
      }

    /*** Decay ***/

    if (!with_filiation_matrix && !with_time_dependence)
      for (s = 0; s < Ns; s++)
        {
          // If the species has no decay, its half-life is set to 0.
          if (half_life(s) == 0.)
            decay = 1.;
          else
            // Delta_t is in seconds while half_life is in days.
            decay = exp(- Model.GetDelta_t() * log(2.) /
                        (half_life(s) * 24 * 3600));
          for (k = 0; k < Nz; k++)
            for (j = 0; j < Ny; j++)
              for (i = 0; i < Nx; i++)
                Model.GetConcentration()(s, k, j, i) *= decay;
        }
    else if (with_time_dependence)
      for (s = 0; s < Ns; s++)
        {
          if (8 < Model.GetCurrentDate().GetHour() &&
              Model.GetCurrentDate().GetHour() < 19)
            if (half_life_day(s) == 0.)
              decay = 1.;
            else
              decay = exp(-Model.GetDelta_t() * log(2.) /
                          (half_life_day(s) * 24 * 3600));
          else if (half_life_night(s) == 0.)
            decay = 1.;
          else
            decay = exp(-Model.GetDelta_t() * log(2.) /
                        (half_life_night(s) * 24 * 3600));
          for (k = 0; k < Nz; k++)
            for (j = 0; j < Ny; j++)
              for (i = 0; i < Nx; i++)
                Model.GetConcentration()(s, k, j, i) *= decay;
        }
    else if (with_filiation_matrix)
      {
        Array<T, 1> Concentration_tmp(Ns);
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                Concentration_tmp = 0.;
                for (s = 0; s < Ns; s++)
                  for (r = 0; r < Ns; r++)
                    Concentration_tmp(s) += filiation_matrix(s, r)
                      * Model.GetConcentration()(r, k, j, i);
                for (s = 0; s < Ns; s++)
                  Model.GetConcentration()(s, k, j, i) = Concentration_tmp(s);
              }
      }
  }


  //! Performs an integration over one time step for aerosol species.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetDelta_t()
    <li> GetNz()
    <li> GetNy()
    <li> GetNx()
    <li> GetCurrentDate()
    <li> GetSource_aer_i()
    <li> SourceGlobalIndex_aer()
    <li> HasSource_aer()
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Decay<T>::Forward_aer(ClassModel& Model)
  {
    int i, j, k, s, b, r;
    T decay;

    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    /*** Decay for aerosol species ***/

    if (!with_filiation_matrix && !with_time_dependence)
      for (s = 0; s < Ns_aer; s++)
        {
          if (half_life_aer(s) == 0.)
            decay = 1.;
          else
            decay = exp(-Model.GetDelta_t() * log(2.) /
                        (half_life_aer(s) * 24 * 3600));
          for (b = 0; b < Nbin_aer; b++)
            for (k = 0; k < Nz; k++)
              for (j = 0; j < Ny; j++)
                for (i = 0; i < Nx; i++)
                  Model.GetConcentration_aer()(s, b, k, j, i) *= decay;
        }
    else if (with_time_dependence)
      for (s = 0; s < Ns_aer; s++)
        {
          if (8 < Model.GetCurrentDate().GetHour()  &&
              Model.GetCurrentDate().GetHour() < 19)
            if (half_life_day_aer(s) == 0.)
              decay = 1.;
            else
              decay =  exp(-Model.GetDelta_t() * log(2.) /
                           (half_life_day_aer(s) * 24 * 3600));
          else if (half_life_night_aer(s) == 0.)
            decay = 1.;
          else
            decay =  exp(-Model.GetDelta_t() * log(2.) /
                         (half_life_night_aer(s)) * 24 * 3600);
          for (b = 0; b < Nbin_aer; b++)
            for (k = 0; k < Nz; k++)
              for (j = 0; j < Ny; j++)
                for (i = 0; i < Nx; i++)
                  Model.GetConcentration_aer()(s, b, k, j, i) *= decay;
        }
    else if (with_filiation_matrix)
      {
        Array<T, 2> Concentration_tmp(Ns_aer, Nbin_aer);
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                Concentration_tmp = 0.;
                for (s = 0; s < Ns_aer; s++)
                  for (r = 0; r < Ns_aer; r++)
                    for (b = 0; b < Nbin_aer; b++)
                      Concentration_tmp(s, b) += filiation_matrix_aer(s, r)
                        * Model.GetConcentration_aer()(r, b, k, j, i);
                for (s = 0; s < Ns_aer; s++)
                  for (b = 0; b < Nbin_aer; b++)
                    Model.GetConcentration_aer()(s, b, k, j, i) =
                      Concentration_tmp(s, b);
              }
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_DECAY_CXX
#endif
