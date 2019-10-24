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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYCASTOR_CXX


#include "ChemistryCastor.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  ChemistryCastor<T>::ChemistryCastor()
  {
  }


  //! Initialization of the mechanism.
  template<class T>
  template<class ClassModel>
  void ChemistryCastor<T>::Init(ClassModel& Model)
  {

    /*** Configuration ***/

    Nz = Model.GetNz();
    Ny = Model.GetNy();
    Nx = Model.GetNx();

    Ns = Model.GetNs();
    species_list = Model.GetSpeciesList();

    ConfigStream config(Model.GetConfigurationFile());

    config.SetSection("[chemistry]");

    config.PeekValue("Reaction_file", reaction_file);
    config.PeekValue("Stoichiometry_file", stoichiometry_file);
    config.PeekValue("Photolysis_file", photolysis_file);
    config.PeekValue("Rates_file", rates_file);

    /*** Reads constants and allocates fields ***/

    int r, s, i, j, k;

    // Fixed values to be put in configuration files.
    Nr = 118;
    Ns_ext = this->Ns + 4;
    Ntemps = 4;

    kreacp.Resize(Ns_ext);
    kreacl.Resize(Ns_ext);

    ireacp.Resize(Ns_ext, Nr);
    ireacl.Resize(Ns_ext, Nr);

    nreactants.Resize(Nr);
    nprods.Resize(Nr);
    irctt.Resize(Nr, 10);

    stoi.Resize(Ns, Nr, Ntemps);

    kreacp.SetZero();
    kreacl.SetZero();

    /*** Reactions ***/

    ExtStream file(reaction_file);
    string species;
    for (r = 0; r < Nr; r++)
      {
        file >> nreactants(r);
        for (s = 0; s < nreactants(r); s++)
          {
            file >> species;
            i = Model.GetSpeciesIndex_ext(species);
            ireacl(i, ++kreacl(i) - 1) = r;
            irctt(r, s) = i;
          }

        file >> nprods(r);
        for (s = 0; s < nprods(r); s++)
          {
            file >> species;
            if (Model.IsSpecies(species))
              {
                i = Model.GetSpeciesIndex(species);
                ireacp(i, ++kreacp(i) - 1) = r;
              }
          }
      }

    stoi.Fill(1.);

    file.Close();

    /*** Stoichiometry ***/

    Array<T, 1> fr(Ntemps);
    file.Open(stoichiometry_file);
    while (!file.IsEmpty())
      {
        file >> species;
        for (j = 0; j < Ntemps; j++)
          file >> fr(j);
        file >> r;
        if (Model.IsSpecies(species))
          {
            s = Model.GetSpeciesIndex(species);
            for (j = 0; j < Ntemps; j++)
              stoi(s, r - 1, j) = fr(j);
          }
      }

    file.Close();

    /*** Photolysis ***/

    file.Open(photolysis_file);

    file >> Ntab_phot;
    file >> Nphot;
    file >> Nwave;
    file >> Nz_phot;

    altiphot.Resize(Nz_phot);
    zenang.Resize(Ntab_phot);
    zetaref.Resize(Ntab_phot);
    photoj.Resize(Ntab_phot, Nz_phot, Nphot);

    for (k = 0; k < Nz_phot; k++)
      file >> altiphot(k);

    for (i = 0; i < Ntab_phot; i++)
      {
        file.SkipLines(3);
        for (s = 0; s < Nphot; s++)
          {
            file >> species;
            file >> zenang(i);
            file.SkipElements(3);
            for (k = 0; k < Nz_phot; k++)
              file >> photoj(i, k, s);
          }
      }

    // ltabrate.
    ltabrate.Resize(24);

    ltabrate(1 - 1) = 1;
    ltabrate(2 - 1) = 2;
    ltabrate(3 - 1) = 3;
    ltabrate(4 - 1) = 7;
    ltabrate(5 - 1) = 1;
    ltabrate(6 - 1) = 4;
    ltabrate(7 - 1) = 4;
    ltabrate(8 - 1) = 8;
    ltabrate(9 - 1) = 8;
    ltabrate(10 - 1) = 4;
    ltabrate(11 - 1) = 2;
    ltabrate(12 - 1) = 8;
    ltabrate(13 - 1) = 1;
    ltabrate(14 - 1) = 7;
    ltabrate(15 - 1) = 1;
    ltabrate(16 - 1) = 6;
    ltabrate(17 - 1) = 3;
    ltabrate(18 - 1) = 4;
    ltabrate(19 - 1) = 2;
    ltabrate(20 - 1) = 5;
    ltabrate(21 - 1) = 4;
    ltabrate(22 - 1) = 5;
    ltabrate(23 - 1) = 5;
    ltabrate(24 - 1) = 5;

    // Cosine of zenith angle.
    const T pi = 3.14159265358979323846264;
    const T ratio = pi / 180.;
    for (i = 0; i < Ntab_phot; i++)
      zetaref(i) = cos(ratio * zenang(i));

    file.Close();

    /*** Reaction rates ***/

    file.Open(rates_file);

    tabrate.Resize(22, Nr);
    tabrate.SetZero();

    ityperate.Resize(Nr);

    iphoto.Resize(Nr);
    for (r = 0; r < Nr; r++)
      {
        iphoto(r) = 0;
        file >> i;
        file >> j;
        for (k = 0; k < ltabrate(j - 1); k++)
          file >> tabrate(k, i - 1);
        ityperate(i - 1) = j - 1;
        if (j == 5 || j == 13)
          iphoto(r) = int(tabrate(0, r) - 1);
      }

    /*** Other allocations ***/

    istoit.Resize(Nz, Ny, Nx);
    wgstl.Resize(Nz, Ny, Nx);
    wgsth.Resize(Nz, Ny, Nx);

    rate.Resize(Nr, Nz, Ny, Nx);
    photorate.Resize(Nphot, Nz, Ny, Nx);
  }


  //! Initialization of the mechanism for the current step.
  template<class T>
  template<class ClassModel>
  void ChemistryCastor<T>::InitStep(ClassModel& Model)
  {
    ComputePhotorate(Model);
    ComputeRate(Model);
  }


  template<class T>
  template<class ClassModel>
  void ChemistryCastor<T>::ComputePhotorate(ClassModel& Model)
  {
    Array<T, 3>& Altitude = Model.A3("Altitude");
    Array<T, 2>& Attenuation = Model.A2("Attenuation");

    int i, j, k, r, kk;

    const T pi = 3.14159265358979323846264;

    T wz1(0.), wz2(0.), z2, z1, ph1, ph2;
    int iz(0), k1(0), k2(0);

    Array<T, 1> w(Nz_phot);

    Date sol_date(Model.GetCurrentDate());
    sol_date.SetMonth(12);
    sol_date.SetDay(21);
    sol_date.SetHour(12);

    T sol_time, hour, loc_time, dec, zenith_angle;
    sol_time = Model.GetCurrentDate().GetSecondsFrom(sol_date) / 3600.;
    hour = T(Model.GetCurrentDate().GetHour())
      + T(Model.GetCurrentDate().GetMinutes()) / 60.
      + T(Model.GetCurrentDate().GetSeconds()) / 3600.
      + Model.GetDelta_t() / 7200.;
    dec = -0.409105176 * cos(7.168657935e-04 * sol_time);

    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
        {
          // Zenith angle.
          loc_time = hour - 12. + 0.066666666666667
            * (Model.GetGridXArray1D()(i));
          loc_time *= 0.261799387;
          zenith_angle =
            sin(Model.GetGridYArray1D()(j) / 180. * pi) * sin(dec)
            + cos(Model.GetGridYArray1D()(j) / 180. * pi)
            * cos(dec) * cos(loc_time);

          if (zenith_angle <= 0.)
            {
              iz = 0;
              wz2 = 0.;
              wz1 = 1.;
            }
          else
            for (int nt = 0; nt < Ntab_phot; nt++)
              if (zenith_angle > zetaref(nt)
                  && zenith_angle <= zetaref(nt + 1))
                {
                  iz = nt;
                  wz2 = (zenith_angle - zetaref(iz))
                    / (zetaref(iz + 1) - zetaref(iz));
                  wz1 = 1. - wz2;
                }

          for (k = 0; k < Nz; k++)
            {
              z2 = Altitude(k + 1, j, i);
              z1 = Altitude(k, j, i);

              for (kk = 0; kk < Nz_phot - 1; kk++)
                if (altiphot(kk) <= z1 && altiphot(kk + 1) >= z1)
                  k1 = kk;
              for (kk = 0; kk < Nz_phot - 1; kk++)
                if (altiphot(kk) <= z2 && altiphot(kk + 1) >= z2)
                  k2 = kk + 1;

              for (kk = 0; kk < Nz_phot - 1; kk++)
                w(kk) = 0.;

              if (k2 - k1 == 1)
                {
                  w(k1) = 0.5 * (2. * altiphot(k2) - z1 - z2) * (z2 - z1)
                    / (altiphot(k2) - altiphot(k1));
                  w(k2) = 0.5 * (z1 + z2 - 2. * altiphot(k1)) * (z2 - z1)
                    / (altiphot(k2) - altiphot(k1));
                }
              else
                {
                  w(k1) = 0.5 * (altiphot(k1 + 1) - z1) * (altiphot(k1 + 1) - z1)
                    / (altiphot(k1 + 1) - altiphot(k1));
                  w(k2) = 0.5 * (z2 - altiphot(k2 - 1)) * (z2 - altiphot(k2 - 1))
                    / (altiphot(k2) - altiphot(k2 - 1));
                  w(k1 + 1) = w(k1 + 1) + 0.5 * (altiphot(k1 + 1) - z1)
                    * (altiphot(k1 + 1) + z1 - 2. * altiphot(k1))
                    / (altiphot(k1 + 1) - altiphot(k1));
                  w(k2 - 1) = w(k2 - 1) + 0.5 * (z2 - altiphot(k2 - 1))
                    * (2. * altiphot(k2) - altiphot(k2 - 1) - z2)
                    / (altiphot(k2) - altiphot(k2 - 1));
                  for (kk = k1 + 1; kk < k2 - 1; kk++)
                    {
                      w(kk) = w(kk) + 0.5 * (altiphot(kk + 1) - altiphot(kk));
                      w(kk + 1) = w(kk + 1)
                        + 0.5 * (altiphot(kk + 1) - altiphot(kk));
                    }
                }

              for (kk = k1; kk < k2 + 1; kk++)
                w(kk) = w(kk) / (z2 - z1);

              for (r = 0; r < Nphot; r++)
                {
                  ph1 = 0.;
                  ph2 = 0.;
                  for (kk = k1; kk < k2 + 1; kk++)
                    {
                      ph1 += w(kk) * photoj(iz, kk, r);
                      ph2 += w(kk) * photoj(iz + 1, kk, r);
                    }
                  photorate(r, k, j, i) = Attenuation(j, i)
                    * (wz1 * ph1 + wz2 * ph2);
                }
            }
        }
  }


  template<class T>
  template<class ClassModel>
  void ChemistryCastor<T>::ComputeRate(ClassModel& Model)
  {
    Data<T, 3>& Temperature = Model.D3("Temperature");
    Data<T, 3>& Pressure = Model.D3("Pressure");
    Data<T, 3>& AirDensity = Model.D3("AirDensity");
    Data<T, 3>& SpecificHumidity = Model.D3("SpecificHumidity");
    Data<T, 2>& Attenuation = Model.D2("Attenuation");
    Data<T, 3>& LiquidWaterContent = Model.D3("LiquidWaterContent");

    Data<T, 3>& DepositionVelocity = Model.D3("DepositionVelocity");

    int i, j, k, nr, kk;

    int ity;
    T te, pe, ai, hu, at, cw, ph, c1, c2, c3, c4,
      f1, f2, f3, f4, ex, factor, t;
    T wl, kc, ve, gama, rhoa, diff, He, ft, tho, za, diam, ph1;

    tho  = 1. / Model.GetDelta_t();
    diff = 0.25e-4;

    const T Rg = 8.205e-2;
    const T conv = 1000.;
    const T an = 6.022045e23;
    const T pi = 3.14159265358979323846264;

    int ino2 =  Model.DepositionVelocityIndex("NO2");
    bool no2_deposition = Model.HasDepositionVelocity("NO2");

    for (nr = 0; nr < Nr; nr++)
      {
        ity = ityperate(nr) + 1;

        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                rate(nr, k, j, i) = 0.;

                te = Temperature(k, j, i);
                pe = Pressure(k, j, i);
                ai = AirDensity(k, j, i);
                hu = SpecificHumidity(k, j, i) * ai * 1.6;
                at = Attenuation(j, i);
                cw = LiquidWaterContent(k, j, i);
                ph = 5.2;

                if (ity == 1)
                  rate(nr, k, j, i) = tabrate(0, nr);

                else if (ity == 2)
                  rate(nr, k, j, i) = tabrate(0, nr)
                    * exp(-tabrate(1, nr) / te);

                else if (ity == 3)
                  rate(nr, k, j, i) =
                    tabrate(0, nr) * exp(-tabrate(1, nr) / te)
                    * pow(300. / te, tabrate(2, nr));

                else if (ity == 4)
                  {
                    c1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te)
                      * pow(300. / te, tabrate(2, nr));
                    c2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                      * pow(300. / te, tabrate(5, nr));
                    c3 = ai * c1;
                    c4 = c3 / c2;
                    ex = 1. / (1. + log10(c4) * log10(c4));
                    rate(nr, k, j, i) = c1 * pow(tabrate(6, nr), ex)
                      / (1. + c4);
                  }

                else if (ity == 15)
                  if (k == 0 && no2_deposition)
                    rate(nr, k, j, i) = tabrate(0, nr)
                      * DepositionVelocity(ino2, j, i);
                  else
                    rate(nr, k, j, i) = 0.;

                else if (ity == 14)
                  {
                    c1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te)
                      * pow(300. / te, tabrate(2, nr));
                    c2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                      * pow(300. / te, tabrate(5, nr));
                    c3 = ai * c1;
                    c4 = c3 / c2;
                    ex = (log10(c4) - 0.12) / 1.2;
                    ex = 1. / (1. + ex * ex);
                    rate(nr, k, j, i) = c1 * pow(tabrate(6, nr), ex)
                      / (1. + c4);
                  }

                else if (ity == 5)
                  rate(nr, k, j, i) = photorate(iphoto(nr), k, j, i);

                else if (ity == 13)
                  {
                    rate(nr, k, j, i) = rate(nr, k, j, i) * at;
                    factor = hu / (hu + ai * (0.02909 * exp(70. / te)
                                              + 0.06545 * exp(110. / te)));
                    rate(nr, k, j, i) = photorate(iphoto(nr), k, j, i)
                      * factor;
                  }

                else if (ity == 6)
                  {
                    f1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te);
                    f2 = tabrate(2, nr) * exp(-tabrate(3, nr) / te);
                    rate(nr, k, j, i) = f1 * f2 / (1. + f2);
                  }
                else if (ity == 7)
                  {
                    f1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te);
                    f2 = tabrate(2, nr) * exp(-tabrate(3, nr) / te);
                    rate(nr, k, j, i) = f1 / (1. + f2);
                  }
                else if (ity == 8)
                  {
                    f1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te);
                    f2 = tabrate(2, nr) * exp(-tabrate(3, nr) / te);
                    f3 = tabrate(4, nr) * exp(-tabrate(5, nr) / te);
                    f4 = tabrate(6, nr) * exp(-tabrate(7, nr) / te);
                    rate(nr, k, j, i) = 2.
                      * sqrt(f1 * f2 * f3 * f4 / ((1. + f3) * (1. + f4)));
                  }
                else if (ity == 9)
                  {
                    f1 = tabrate(0, nr) * exp(-tabrate(1, nr) / te);
                    f2 = tabrate(2, nr) * exp(-tabrate(3, nr) / te);
                    f3 = tabrate(4, nr) * exp(-tabrate(5, nr) / te);
                    f4 = tabrate(6, nr) * exp(-tabrate(7, nr) / te);
                    f3 = f3 / (1. + f3);
                    f4 = f4 / (1. + f4);
                    rate(nr, k, j, i) = 2. * sqrt(f1 * f2)
                      * (1. - sqrt(f3 * f4)) * (1. - f4) / (2. - f3 - f4);
                  }
                else if (ity == 16)
                  {
                    ft = exp(tabrate(2, nr) * (1. / te - 1. / 298.));
                    if (tabrate(3, nr) == 0)
                      {
                        c1 = pow(10., -ph);
                        c1 = pow(c1, tabrate(0, nr));
                      }
                    else
                      c1 = 1. + 13. * pow(10., -ph);

                    rate(nr, k, j, i) = 1.e-20;
                    if (cw < 1.e-11)
                      rate(nr, k, j, i) = 1.e3 * cw * (Rg * Rg * te * te / an)
                        * ft * tabrate(1, nr) * tabrate(4, nr)
                        * tabrate(5, nr) / c1;
                  }
                else if (ity == 23)
                  {
                    ft = exp(tabrate(3, nr) * (1. / te - 1. / 298.));
                    ph1 = min(ph, 5.0);
                    c1 = pow(10., -ph1);
                    rate(nr, k, j, i) = ft * cw * tabrate(0, nr)
                      * tabrate(1, nr) * tabrate(2, nr) * 8.314 * te
                      * 1.e6 / (1.013e5 * conv) / pow(c1, tabrate(4, nr));
                  }
                else if (ity == 24)
                  {
                    ft = exp(tabrate(3, nr) * (1. / te - 1. / 298.));
                    ph1 = min(ph, 6.0);
                    c1 = pow(10., -ph1);
                    rate(nr, k, j, i) = ft * cw * tabrate(0, nr)
                      * tabrate(1, nr) * tabrate(2, nr) * 8.314
                      * te * 1.e6 / (1.013e5 * conv)
                      / pow(c1, tabrate(4, nr));
                  }
                else if (ity == 20)
                  {
                    rate(nr, k, j, i) = 1.e-10;
                    ve = sqrt(8. + 3. * 8.314 * te / (pi * tabrate(1, nr)));
                    if (tabrate(3, nr) == 0.)
                      gama = tabrate(0, nr);
                    else
                      {
                        za = (log(1.) - log(tabrate(0, nr)))
                          / (log(273.) - log(298.));
                        gama = za * (log(te) - log(298.))
                          + log(tabrate(0, nr));
                        gama = exp(gama);
                        gama = max(min(1., gama), 0.001);
                      }
                    if (tabrate(4, nr) == 0.)
                      {
                        throw "Particles!";
                      }
                    else
                      {
                        // Cloud droplet diameter.
                        diam = 10.e-6;
                        c1 = diam / diff / 2.;
                        c2 = 4. / (gama * ve);
                        rate(nr, k, j, i) =
                          1.e-6 * (1. + 8. * cw / (diam * 1.e2)) / (c1 + c2);
                      }
                    rate(nr, k, j, i) = min(rate(nr, k, j, i), tho);
                  }
                else if (ity == 21)
                  {
                    rhoa = ai * 29. * 1000. / an;
                    if (tabrate(0, nr) == 1)
                      {
                        if (cw > 1.e-11)
                          wl = an * cw / (ai * 29.);
                        else
                          wl = 0.0;
                        ve = sqrt(8. + 3. * 8.314 * te
                                  / (pi * tabrate(2, nr)));
                        kc = tabrate(1, nr) / 2. / diff
                          + 4. / ve / tabrate(3, nr);
                        rate(nr, k, j, i) = 6. * wl * rhoa
                          / (1. + 3. * tabrate(1, nr) * kc);
                      }
                    else
                      throw "*** RATES: ERROR, ITY=21 and C1 != 1";
                    rate(nr, k, j, i) = min(rate(nr, k, j, i), tho);
                  }
                else if (ity == 22)
                  {
                    He = tabrate(1, nr)
                      * exp(-tabrate(2, nr) * (1. / te - 1. / 298.));
                    ve = sqrt(8. + 3. * 8.314 * te / (pi * tabrate(3, nr)));
                    kc = tabrate(0, nr) / 2. / diff
                      + 4. / ve / tabrate(4, nr);
                    rate(nr, k, j, i) = 6.e2
                      / (8.314 * He * te * tabrate(0, nr) * kc);
                    rate(nr, k, j, i) = min(rate(nr, k, j, i), tho);
                  }
                else
                  throw string("*** ERROR: Undefined rate type: ")
                    + to_str(ity);
              }
      }

    Data<T, 1> dtemp(Ntemps - 1);
    Data<T, 1> tabtemp(Ntemps);
    tabtemp(0) = 260.;
    tabtemp(1) = 280.;
    tabtemp(2) = 300.;
    tabtemp(3) = 320.;
    for (k = 0; k < Ntemps - 1; k++)
      dtemp(k) = tabtemp(k + 1) - tabtemp(k);

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            t = Temperature(k, j, i);
            if (t < tabtemp(0))
              {
                istoit(k, j, i) = 0;
                wgstl(k, j, i) = 1.;
                wgsth(k, j, i) = 0.;
              }
            else if (t < tabtemp(Ntemps - 1))
              {
                for (kk = 0; kk < Ntemps - 1; kk++)
                  if (t >= tabtemp(kk) && t < tabtemp(kk + 1))
                    {
                      istoit(k, j, i) = kk;
                      wgstl(k, j, i) = (tabtemp(kk + 1) - t) / dtemp(kk);
                      wgsth(k, j, i) = 1. - wgstl(k, j, i);
                    }
              }
            else
              {
                istoit(k, j, i) = Ntemps - 2;
                wgstl(k, j, i) = 0.;
                wgsth(k, j, i) = 1.;
              }
          }
  }


  template<class T>
  template<class ClassModel>
  void ChemistryCastor<T>
  ::LossProduction(ClassModel& Model, int s, int k, int j, int i,
                   T& loss, T& production)
  {
    Data<T, 3>& AirDensity = Model.AirDensity;
    Data<T, 3>& SpecificHumidity = Model.SpecificHumidity;

    Data<T, 4>& Concentration = Model.GetConcentration();

    int ir;
    T trat;

    int kr, it, nr, sp;
    for (kr = 0; kr < kreacl(s); kr++)
      {
        ir = ireacl(s, kr);
        trat = rate(ir, k, j, i);
        for (it = 0; it < nreactants(ir); it++)
          {
            sp = irctt(ir, it);
            if (sp < Ns)
              trat *= Concentration(sp, k, j, i);
            else if (sp == Ns)
              trat *= AirDensity(k, j, i);
            else if (sp == Ns + 1)
              trat *= 0.2095 * AirDensity(k, j, i);
            else if (sp == Ns + 2)
              trat *= 0.7905 * AirDensity(k, j, i);
            else if (sp == Ns + 3)
              trat *= SpecificHumidity(k, j, i) * AirDensity(k, j, i) * 1.6;
            else
              throw "Wrong irctt index.";
          }
        loss += trat;
      }

    int isto = istoit(k, j, i);
    T wgl = wgstl(k, j, i);
    T wgh = wgsth(k, j, i);

    for (kr = 0; kr < kreacp(s); kr++)
      {
        ir = ireacp(s, kr);
        trat = rate(ir, k, j, i);
        for (nr = 0; nr < nreactants(ir); nr++)
          {
            sp = irctt(ir, nr);
            if (sp < Ns)
              trat *= Concentration(sp, k, j, i);
            else if (sp == Ns)
              trat *= AirDensity(k, j, i);
            else if (sp == Ns + 1)
              trat *= 0.2095 * AirDensity(k, j, i);
            else if (sp == Ns + 2)
              trat *= 0.7905 * AirDensity(k, j, i);
            else if (sp == Ns + 3)
              trat *= SpecificHumidity(k, j, i) * AirDensity(k, j, i) * 1.6;
            else
              throw "Wrong irctt index.";
          }
        production += trat
          * (wgl * stoi(s, ir, isto) + wgh * stoi(s, ir, isto + 1));
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_CHEMISTRYCASTOR_CXX
#endif
