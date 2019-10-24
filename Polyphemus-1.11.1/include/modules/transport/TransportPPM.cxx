// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).


#ifndef POLYPHEMUS_FILE_MODULES_TRANSPORT_TRANSPORTPPM_CXX


#include "TransportPPM.hxx"


namespace Polyphemus
{


  //! Initialization of the scheme.
  /*!
    \param Model (input/output) model. Its interface must contain:
    <ul>
    <li> GetSpeciesFile()
    <li> IsSpecies()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::Init(ClassModel& Model)
  {

    /*** PPM species ***/

    ConfigStream species_stream(Model.GetSpeciesFile());

    species_stream.SetSection("[ppm_species]");

    string species;
    while (!species_stream.IsEmpty())
      {
        species = species_stream.GetElement();
        if (!Model.IsSpecies(species))
          throw string("Unknown species in section \"[ppm_species]\": ")
            + species + ".";
        species_list_ppm.push_back(species);
      }

    Ns_ppm = int(species_list_ppm.size());
    species_list_ppm_index.resize(Ns_ppm);
    for (int i = 0; i < int(species_list_ppm.size()); i++)
      {
        int index = Model.GetSpeciesIndex(species_list_ppm[i]);
        species_list_ppm_index(i) = index;
      }

    /*** Dimensions ***/

    Nz = Model.GetNz();
    Ny = Model.GetNy();
    Nx = Model.GetNx();

    /*** Allocations ***/

    WestFlux.Resize(Nz, Ny, Nx);
    EastFlux.Resize(Nz, Ny, Nx);
    SouthFlux.Resize(Nz, Ny, Nx);
    NorthFlux.Resize(Nz, Ny, Nx);

    VerticalFlux_out.Resize(Nz, Ny, Nx);
    VerticalFlux_in.Resize(Nz + 1, Nz, Ny, Nx);

    // Winds on interfaces.
    WestWind.Resize(Nz, Ny, Nx);
    EastWind.Resize(Nz, Ny, Nx);
    SouthWind.Resize(Nz, Ny, Nx);
    NorthWind.Resize(Nz, Ny, Nx);
    VerticalWind.Resize(Nz, Ny, Nx);

    AirDensity_ext.Resize(Nz, Ny + 2, Nx + 2);
  }


  //! Initialization for the current step.
  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::InitStep(ClassModel& Model)
  {
    int k, j, i;

    Data<T, 3>& ZonalWind = Model.ZonalWind;
    Data<T, 3>& MeridionalWind = Model.MeridionalWind;
    Data<T, 3>& VerticalDiffusionCoefficient
      = Model.VerticalDiffusionCoefficient;

    Data<T, 3>& AirDensity = Model.AirDensity;
    Data<T, 3>& FileAirDensity_i = Model.FileAirDensity_i;
    Data<T, 3>& FileAirDensity_f = Model.FileAirDensity_f;

    Array<T, 1>& CellWidth_x = Model.CellWidth_x;
    Array<T, 1>& CellWidth_y = Model.CellWidth_y;
    Data<T, 3>& Altitude = Model.Altitude;
    Data<T, 3>& Thickness = Model.Thickness;

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        {
          WestWind(k, j, 0) =
            1.5 * ZonalWind(k, j, 0) - 0.5 * ZonalWind(k, j, 1);
          WestFlux(k, j, 0) = WestWind(k, j, 0)
            * (1.5 * Thickness(k, j, 0) - 0.5 * Thickness(k, j, 1))
            / (Thickness(k, j, 0) * CellWidth_x(j));
          for (i = 1; i < Nx; i++)
            {
              WestWind(k, j, i) =
                0.5 * (ZonalWind(k, j, i) + ZonalWind(k, j, i - 1));
              WestFlux(k, j, i) = WestWind(k, j, i) * 0.5
                * (Thickness(k, j, i) + Thickness(k, j, i - 1))
                / (Thickness(k, j, i) * CellWidth_x(j));
            }
        }

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        {
          for (i = 0; i < Nx - 1; i++)
            {
              EastWind(k, j, i) =
                0.5 * (ZonalWind(k, j, i + 1) + ZonalWind(k, j, i));
              EastFlux(k, j, i) = EastWind(k, j, i) * 0.5
                * (Thickness(k, j, i + 1) + Thickness(k, j, i))
                / (Thickness(k, j, i) * CellWidth_x(j));
            }
          EastWind(k, j, Nx - 1) =
            1.5 * ZonalWind(k, j, Nx - 1)
            - 0.5 * ZonalWind(k, j, Nx - 2);
          EastFlux(k, j, Nx - 1) = EastWind(k, j, Nx - 1)
            * (1.5 * Thickness(k, j, Nx - 1)
               - 0.5 * Thickness(k, j, Nx - 2))
            / (Thickness(k, j, Nx - 1) * CellWidth_x(j));
        }

    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        {
          SouthWind(k, 0, i) =
            1.5 * MeridionalWind(k, 0, i) - 0.5 * MeridionalWind(k, 1, i);
          SouthFlux(k, 0, i) = SouthWind(k, 0, i)
            * (1.5 * Thickness(k, 0, i) - 0.5 * Thickness(k, 1, i))
            * (1.5 * CellWidth_x(0) - 0.5 * CellWidth_x(1))
            / (Thickness(k, 0, i) * CellWidth_x(0) * CellWidth_y(0));
          for (j = 1; j < Ny; j++)
            {
              SouthWind(k, j, i) =
                0.5 * (MeridionalWind(k, j, i)
                       + MeridionalWind(k, j - 1, i));
              SouthFlux(k, j, i) = SouthWind(k, j, i) * 0.5
                * (Thickness(k, j, i) + Thickness(k, j - 1, i))
                * 0.5 * (CellWidth_x(j) + CellWidth_x(j - 1))
                / (Thickness(k, j, i) * CellWidth_x(j) * CellWidth_y(j));
            }
        }

    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        {
          for (j = 0; j < Ny - 1; j++)
            {
              NorthWind(k, j, i) =
                0.5 * (MeridionalWind(k, j + 1, i) + MeridionalWind(k, j, i));
              NorthFlux(k, j, i) = NorthWind(k, j, i) * 0.5
                * (Thickness(k, j + 1, i) + Thickness(k, j, i))
                * 0.5 * (CellWidth_x(j + 1) + CellWidth_x(j))
                / (Thickness(k, j, i) * CellWidth_x(j) * CellWidth_y(j));
            }
          NorthWind(k, Ny - 1, i) =
            1.5 * MeridionalWind(k, Ny - 1, i)
            - 0.5 * MeridionalWind(k, Ny - 2, i);
          NorthFlux(k, Ny - 1, i) = NorthWind(k, Ny - 1, i)
            * (1.5 * Thickness(k, Ny - 1, i) - 0.5 * Thickness(k, Ny - 2, i))
            * (1.5 * CellWidth_x(Ny - 1) - 0.5 * CellWidth_x(Ny - 2))
            / (Thickness(k, Ny - 1, i) * CellWidth_x(Ny - 1)
               * CellWidth_y(Ny - 1));
        }

    // Extended densities.
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          AirDensity_ext(k, j + 1, i + 1) = AirDensity(k, j, i);
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        {
          AirDensity_ext(k, j + 1, 0) = AirDensity_ext(k, j + 1, 1);
          AirDensity_ext(k, j + 1, Nx + 1) = AirDensity_ext(k, j + 1, Nx);
        }
    for (k = 0; k < Nz; k++)
      for (i = 0; i < Nx; i++)
        {
          AirDensity_ext(k, 0, i + 1) = AirDensity_ext(k, 1, i + 1);
          AirDensity_ext(k, Ny + 1, i + 1) = AirDensity_ext(k, Ny, i + 1);
        }

    T budget(0.);
    T below_flux;
    T thickness, flux;
    VerticalFlux_in.SetZero();
    VerticalFlux_out.SetZero();
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        for (k = 0; k < Nz; k++)
          {
            below_flux = k == 0 ? 0. : VerticalFlux_in(k - 1, k, j, i)
              * AirDensity(k - 1, j, i);
            budget = below_flux
              + max(WestFlux(k, j, i), 0.) * AirDensity_ext(k, j + 1, i)
              + max(-EastFlux(k, j, i), 0.) * AirDensity_ext(k, j + 1, i + 2)
              + max(SouthFlux(k, j, i), 0.) * AirDensity_ext(k, j, i + 1)
              + max(-NorthFlux(k, j, i), 0.) * AirDensity_ext(k, j + 2, i + 1)
              - (max(-WestFlux(k, j, i), 0.) + max(EastFlux(k, j, i), 0.)
                 + max(-SouthFlux(k, j, i), 0.) + max(NorthFlux(k, j, i), 0.)
                 + VerticalFlux_out(k, j, i)) * AirDensity(k, j, i)
              - (FileAirDensity_f(k, j, i) - FileAirDensity_i(k, j, i))
              / Model.GetMeteoDelta_t();

            if (budget > 0.)
              {
                VerticalWind(k, j, i)
                  = budget * (Altitude(k + 1, j, i) - Altitude(k, j, i))
                  / AirDensity(k, j, i);
                VerticalFlux_out(k, j, i) += budget / AirDensity(k, j, i);
                if (k < Nz - 1)
                  VerticalFlux_in(k, k + 1, j, i) += VerticalWind(k, j, i)
                    / (Altitude(k + 2, j, i) - Altitude(k + 1, j, i));
              }
            else
              {
                VerticalWind(k, j, i)
                  = budget * (Altitude(k + 1, j, i) - Altitude(k, j, i))
                  / AirDensity(min(k + 1, Nz - 1), j, i);
                VerticalFlux_in(k + 1, k, j, i)
                  -= budget / AirDensity(min(k + 1, Nz - 1), j, i);
                if (k < Nz - 1)
                  VerticalFlux_out(k + 1, j, i) -= VerticalWind(k, j, i)
                    / (Altitude(k + 2, j, i) - Altitude(k + 1, j, i));
              }

            // Vertical diffusion.
            if (k != Nz - 1)
              thickness = (Altitude(k + 2, j, i) - Altitude(k, j, i)) / 2.;
            else
              thickness = Altitude(k + 1, j, i) - Altitude(k, j, i);
            flux = VerticalDiffusionCoefficient(k + 1, j, i)
              / (Altitude(k + 1, j, i) - Altitude(k, j, i)) / thickness;
            VerticalFlux_in(k + 1, k, j, i) += flux * AirDensity(k, j, i)
              / AirDensity(min(k + 1, Nz - 1), j, i);
            VerticalFlux_out(k, j, i) += flux;

            if (k != Nz - 1)
              {
                flux = VerticalDiffusionCoefficient(k + 1, j, i)
                  / (Altitude(k + 2, j, i) - Altitude(k + 1, j, i)) / thickness;
                VerticalFlux_in(k, k + 1, j, i) += flux;
                VerticalFlux_out(k + 1, j, i) += flux
                  * AirDensity(k, j, i) / AirDensity(min(k + 1, Nz - 1), j, i);
              }
          }
  }


  //! Performs an integration over one time step.
  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::LossProduction(ClassModel& Model,
                                       int s, int k, int j, int i,
                                       T& loss, T& production)
  {
    int h = 0;
    while (h < Ns_ppm && species_list_ppm_index(h) != s)
      h++;
    if (h == Ns_ppm)
      LossProductionUpwind(Model, s, k, j, i, loss, production);
    else
      LossProductionPPM(Model, s, k, j, i, loss, production);
  }


  //! Performs an integration over one time step.
  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::LossProductionUpwind(ClassModel& Model,
                                             int s, int k, int j, int i,
                                             T& loss, T& production)
  {
    Data<T, 4>& Concentration = Model.GetConcentration();

    Data<T, 3>& BoundaryCondition_z = Model.BoundaryCondition_z;
    Data<T, 4>& BoundaryCondition_y = Model.BoundaryCondition_y;
    Data<T, 4>& BoundaryCondition_x = Model.BoundaryCondition_x;

    // Western fluxes.
    if (WestWind(k, j, i) > 0. && i > 0)
      production += WestFlux(k, j, i) * Concentration(s, k, j, i - 1);
    else if (WestWind(k, j, i) > 0.)
      production += WestFlux(k, j, i) * BoundaryCondition_x(s, k, j, 0);
    else
      loss -= WestFlux(k, j, i) * Concentration(s, k, j, i);

    // Eastern fluxes.
    if (EastWind(k, j, i) < 0. && i < Nx - 1)
      production -= EastFlux(k, j, i) * Concentration(s, k, j, i + 1);
    else if (EastWind(k, j, i) < 0.)
      production -= EastFlux(k, j, i) * BoundaryCondition_x(s, k, j, 1);
    else
      loss += EastFlux(k, j, i) * Concentration(s, k, j, i);

    // Southern fluxes.
    if (SouthWind(k, j, i) > 0. && j > 0)
      production += SouthFlux(k, j, i) * Concentration(s, k, j - 1, i);
    else if (SouthWind(k, j, i) > 0.)
      production += SouthFlux(k, j, i) * BoundaryCondition_y(s, k, 0, i);
    else
      loss -= SouthFlux(k, j, i) * Concentration(s, k, j, i);

    // Northern side.
    if (NorthWind(k, j, i) < 0. && j < Ny - 1)
      production -= NorthFlux(k, j, i) * Concentration(s, k, j + 1, i);
    else if (NorthWind(k, j, i) < 0.)
      production -= NorthFlux(k, j, i) * BoundaryCondition_y(s, k, 1, i);
    else
      loss += NorthFlux(k, j, i) * Concentration(s, k, j, i);

    // Vertical.
    loss += VerticalFlux_out(k, j, i) * Concentration(s, k, j, i);
    for (int k_prod = 0; k_prod < Nz; k_prod++)
      production += VerticalFlux_in(k_prod, k, j, i)
        * Concentration(s, k_prod, j, i);
    production += VerticalFlux_in(Nz, k, j, i)
      * BoundaryCondition_z(s, j, i);
  }


  //! Performs an integration over one time step.
  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::LossProductionPPM(ClassModel& Model,
                                          int s, int k, int j, int i,
                                          T& loss, T& production)
  {
    Data<T, 4>& Concentration = Model.GetConcentration();

    Data<T, 3>& AirDensity = Model.AirDensity;

    Data<T, 3>& BoundaryCondition_z = Model.BoundaryCondition_z;

    T Delta_t = Model.GetDelta_t();

    // Vertical.
    loss += VerticalFlux_out(k, j, i) * Concentration(s, k, j, i);

    for (int k_prod = 0; k_prod < Nz; k_prod++)
      production += VerticalFlux_in(k_prod, k, j, i)
        * Concentration(s, k_prod, j, i);
    production += VerticalFlux_in(Nz, k, j, i)
      * BoundaryCondition_z(s, j, i);

    // Western fluxes.
    if (WestWind(k, j, i) > 0.)
      production += WestFlux(k, j, i)
        * MeanConcRight(Model, s, k, j, i - 1, WestWind(k, j, i) * Delta_t, true)
        * AirDensity(k, j, max(i - 1, 0));
    else
      loss -= WestFlux(k, j, i)
        * MeanConcLeft(Model, s, k, j, i, -WestWind(k, j, i) * Delta_t, true)
        * AirDensity(k, j, i);

    // Eastern fluxes.
    if (EastWind(k, j, i) < 0.)
      production -= EastFlux(k, j, i)
        * MeanConcLeft(Model, s, k, j, i + 1, -EastWind(k, j, i) * Delta_t, true)
        * AirDensity(k, j, min(i + 1, Nx - 1));
    else
      loss += EastFlux(k, j, i)
        * MeanConcRight(Model, s, k, j, i, EastWind(k, j, i) * Delta_t, true)
        * AirDensity(k, j, i);

    // Southern fluxes.
    if (SouthWind(k, j, i) > 0.)
      production += SouthFlux(k, j, i)
        * MeanConcRight(Model, s, k, j - 1, i,
                        SouthWind(k, j, i) * Delta_t, false)
        * AirDensity(k, max(j - 1, 0), i);
    else
      loss -= SouthFlux(k, j, i)
        * MeanConcLeft(Model, s, k, j, i, -SouthWind(k, j, i) * Delta_t, false)
        * AirDensity(k, j, i);

    // Northern side.
    if (NorthWind(k, j, i) < 0.)
      production -= NorthFlux(k, j, i)
        * MeanConcLeft(Model, s, k, j + 1, i,
                       -NorthWind(k, j, i) * Delta_t, false)
        * AirDensity(k, min(j + 1, Ny - 1), i);
    else
      loss += NorthFlux(k, j, i)
        * MeanConcRight(Model, s, k, j, i, NorthWind(k, j, i) * Delta_t, false)
        * AirDensity(k, j, i);
  }


  template<class T>
  template<class ClassModel>
  void TransportPPM<T>::Reconstruct(ClassModel& Model, int s, int k, int j,
                                    int i, T y, T& conc_left, T& conc_right,
                                    T& delta_conc, T& parab, T& x, bool zonal)
  {
    Data<T, 4>& Concentration = Model.GetConcentration();

    Data<T, 4>& BoundaryCondition_y = Model.BoundaryCondition_y;
    Data<T, 4>& BoundaryCondition_x = Model.BoundaryCondition_x;

    Data<T, 3>& AirDensity = Model.AirDensity;

    Array<T, 1>& CellWidth_x = Model.CellWidth_x;
    Array<T, 1>& CellWidth_y = Model.CellWidth_y;

    T conc;

    if (zonal && i == -1)
      conc = BoundaryCondition_x(s, k, j, 0) / AirDensity(k, j, 0);
    else if (zonal && i == Nx)
      conc = BoundaryCondition_x(s, k, j, 1) / AirDensity(k, j, i - 1);
    else if (!zonal && j == -1)
      conc = BoundaryCondition_y(s, k, 0, i) / AirDensity(k, 0, i);
    else if (!zonal && j == Ny)
      conc = BoundaryCondition_y(s, k, 1, i) / AirDensity(k, j - 1, i);
    else
      conc = Concentration(s, k, j, i) / AirDensity(k, j, i);

    T conc_ll, conc_l, conc_r, conc_rr;

    int jleft = zonal ? j : j - 1;
    int jright = zonal ? j : j + 1;
    int jlleft = zonal ? j : j - 2;
    int jrright = zonal ? j : j + 2;

    int ileft = zonal ? i - 1 : i;
    int iright = zonal ? i + 1 : i;
    int illeft = zonal ? i - 2 : i;
    int irright = zonal ? i + 2 : i;

    if ((zonal && (i <= 0 || i >= Nx - 1))
        || (!zonal && (j <= 0 || j >= Ny - 1)))
      {
        conc_left = conc;
        conc_right = conc;
        delta_conc = 0.;
        parab = 0.;
        x = 0.;
        return;
      }

    // Courant.
    if (zonal)
      x = y / CellWidth_x(j);
    else
      x = y / CellWidth_y(j);

    if ((zonal && i == 1) || (!zonal && j == 1))
      {
        conc_l = Concentration(s, k, jleft, ileft)
          / AirDensity(k, jleft, ileft);
        conc_r = Concentration(s, k, jright, iright)
          / AirDensity(k, jright, iright);
        conc_rr = Concentration(s, k, jrright, irright)
          / AirDensity(k, jrright, irright);
        conc_left = .5 * (conc_l + conc);
        conc_right = PPMInterpolation(conc_l, conc, conc_r, conc_rr);
      }
    else if ((zonal && i == Nx - 2) || (!zonal && j == Ny - 2))
      {
        conc_ll = Concentration(s, k, jlleft, illeft)
          / AirDensity(k, jlleft, illeft);
        conc_l = Concentration(s, k, jleft, ileft)
          / AirDensity(k, jleft, ileft);
        conc_r = Concentration(s, k, jright, iright)
          / AirDensity(k, jright, iright);
        conc_left = PPMInterpolation(conc_ll, conc_l, conc, conc_r);
        conc_right = .5 * (conc + conc_r);
      }
    else
      {
        conc_ll = Concentration(s, k, jlleft, illeft)
          / AirDensity(k, jlleft, illeft);
        conc_l = Concentration(s, k, jleft, ileft)
          / AirDensity(k, jleft, ileft);
        conc_r = Concentration(s, k, jright, iright)
          / AirDensity(k, jright, iright);
        conc_rr = Concentration(s, k, jrright, irright)
          / AirDensity(k, jrright, irright);
        conc_left = PPMInterpolation(conc_ll, conc_l, conc, conc_r);
        conc_right = PPMInterpolation(conc_l, conc, conc_r, conc_rr);
      }

    delta_conc = conc_right - conc_left;
    parab = 6. * (conc - .5 * (conc_right + conc_left));

    T test = (conc_right - conc) * (conc - conc_left);

    if (test <= 0.)
      {
        conc_left = conc;
        conc_right = conc;
        delta_conc = 0.;
        parab = 0.;
      }
    else
      {
        T dcpar = delta_conc * parab;
        T dc2 = delta_conc * delta_conc;
        if (dcpar > dc2)
          {
            conc_left = 3. * conc - 2. * conc_right;
            delta_conc = conc_right - conc_left;
            parab = 6. * (conc - .5 * (conc_left + conc_right));
          }
        if (dcpar < -dc2)
          {
            conc_right = 3. * conc - 2. * conc_left;
            delta_conc = conc_right - conc_left;
            parab = 6. * (conc - .5 * (conc_left + conc_right));
          }
      }
  }


  template<class T>
  T TransportPPM<T>::PPMInterpolation(T c1, T c2, T c3, T c4)
  {
    T delta_l = PPMDelta(c1, c2, c3);
    T delta_r = PPMDelta(c2, c3, c4);

    return .5 * (c3 + c2) + (delta_l - delta_r) / 6.;
  }


  template<class T>
  T TransportPPM<T>::PPMDelta(T c1, T c2, T c3)
  {
    T delta_front = c3 - c2;
    T delta_behind = c2 - c1;

    if (delta_front > 0. && delta_behind > 0.)
      return 2. * min(min(.25 * (c3 - c1), delta_front), delta_behind);
    if (delta_front < 0. && delta_behind < 0.)
      return 2. * max(max(.25 * (c3 - c1), delta_front), delta_behind);

    return 0.;
  }


  template<class T>
  template<class ClassModel>
  T TransportPPM<T>::MeanConcLeft(ClassModel& Model,
                                  int s, int k, int j, int i, T y, bool zonal)
  {
    T conc_left;
    T conc_right;
    T delta_conc;
    T parab;
    T x;

    Reconstruct(Model, s, k, j, i, y, conc_left, conc_right, delta_conc,
                parab, x, zonal);

    return conc_left + .5 * x * (delta_conc + (1. - 2. / 3. * x) * parab);
  }


  template<class T>
  template<class ClassModel>
  T TransportPPM<T>::MeanConcRight(ClassModel& Model, int s, int k, int j,
                                   int i, T y, bool zonal)
  {
    T conc_left;
    T conc_right;
    T delta_conc;
    T parab;
    T x;

    Reconstruct(Model, s, k, j, i, y, conc_left, conc_right, delta_conc,
                parab, x, zonal);

    return conc_right - .5 * x * (delta_conc - (1. - 2. / 3. * x) * parab);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_TRANSPORT_TRANSPORTPPM_CXX
#endif
