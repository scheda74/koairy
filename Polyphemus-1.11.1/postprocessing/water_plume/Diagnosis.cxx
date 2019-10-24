// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Meryem Ahmed de Biasi
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


//! Computes the saturation vapor pressure.
/*!
  \param temperature temperature (K).
  \return The saturation vapor pressure (Pa).
*/
template<class T>
inline T Pvsat(T temperature)
{
  return 610.78 * exp(17.2694 * (temperature - 273.15)
                      / (temperature - 35.86));
}


/*! \brief Computes the derivative of saturation vapor pressure with respect
  to temperature.
*/
/*!
  \param temperature temperature (K).
  \return The derivative of saturation vapor pressure with respect to
  temperature.
*/
template<class T>
inline T dPvsat_dT(T temperature)
{
  return Pvsat(temperature) * 17.2694 * (273.15 - 35.86)
    / ((temperature - 35.86) * (temperature - 35.86));
}


//! Computes the saturation humidity.
/*!
  \param pressure pressure (Pa).
  \param temperature temperature (K).
  \return The saturation humidity (kg/kg).
*/
template<class T>
inline T qsat(T pressure, T temperature)
{
  return 0.62197 * Pvsat(temperature)
    / (pressure - (1 - 0.62197) * Pvsat(temperature));
}


/*! \brief Computes the derivative of saturation humidity with respect to
  temperature.
*/
/*!
  \param pressure pressure (Pa).
  \param temperature temperature (K).
  \return The derivative of saturation humidity with respect to temperature.
*/
template<class T>
inline T dqsat_dT(T pressure, T temperature)
{
  return 0.62197 * dPvsat_dT(temperature) * pressure
    / ((pressure + (0.62197 - 1) * Pvsat(temperature))
       * (pressure + (0.62197 - 1) * Pvsat(temperature)));
}


//! Computes a function.
/*!
  \param pressure pressure (Pa).
  \param temperature temperature (K).
  \param thetal potential liquid temperature (K).
  \param qw specific content in water (kg/kg).
  \return The value of the function.
*/
template<class T>
inline T Function(T pressure, T temperature, T thetal, T qw)
{
  const T cp = 1005.;
  const T L = 2.501e6;
  const T P0 = 100000.;
  const T Rair = 287.;

  return log(temperature) - L / (cp * temperature)
    * (qw - qsat(pressure, temperature)) + Rair / cp * log(P0 / pressure)
    - log(thetal);
}


/*! \brief Computes the derivative of a function (with respect to the
  temperature).
*/
/*!
  \param pressure pressure (Pa).
  \param temperature temperature (K).
  \param thetal potential liquid temperature (K).
  \param qw specific content in water (kg/kg).
  \return The derivative of the function with respects to temperature.
*/
template<class T>
inline T dFunction_dT(T pressure, T temperature, T thetal, T qw)
{
  const T cp = 1005.;
  const T L = 2.501e6;

  return 1 / temperature + L / (cp * temperature)
    * (dqsat_dT(pressure, temperature)
       + (qw - qsat(pressure, temperature)) / temperature);
}


//! Computes the temperature by solving Function(T) = 0.
/*!
  \param pressure pressure (Pa).
  \param temperature_amb ambient temperature (first guess) (K).
  \param thetal potential liquid temperature (K).
  \param qw specific content in water (kg/kg).
  \return The temperature after evaporation and condensation of
  part of the water.
*/
template<class T>
T Newton(T pressure, T temperature_amb, T thetal, T qw)
{
  int i, Nmax(1000);
  T epsilon(1e-2), temperature(temperature_amb), temperature_prev;
  T diff;
  bool convergence =  false;

  for (i = 0; i < Nmax; i++)
    if (!convergence)
      {
        temperature_prev = temperature;
        temperature = temperature
          - Function(pressure, temperature, thetal, qw)
          / dFunction_dT(pressure, temperature, thetal, qw);
        diff = abs(temperature - temperature_prev);
        convergence = (i != 0) && (diff < epsilon);
      }

  if (!convergence)
    throw string("The computation of the temperature did not converge.\n");

  return temperature;
}


//! Computes air density on the basis of temperature and pressure.
/*! Formula: AirDensity = Pressure / (287. * Temperature).
  \param Temperature temperature.
  \param Pressure pressure.
  \param AirDensity (output) air density.
*/
template<class T>
void ComputeAirDensity(const Data<T, 4>& Temperature,
                       const Data<T, 4>& Pressure,
                       Data<T, 4>& AirDensity)
{
  const T r = 287.;

  int h, k, j, i;

  int Nt = AirDensity.GetLength(0);
  int Nz = AirDensity.GetLength(1);
  int Ny = AirDensity.GetLength(2);
  int Nx = AirDensity.GetLength(3);

  for (h = 0; h < Nt; h++)
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          AirDensity(h, k, j, i) = Pressure(h, k, j, i)
            / (r * Temperature(h, k, j, i));
}


//! Computes the liquid water potential temperature.
/*!
  \param Temperature temperature (K).
  \param Pressure pressure (Pa).
  \param LiquidWaterContent the liquid water content (kg.kg^{-1}).
  \param PotentialTemperatureLW (output) the liquid water potential
  temperature (K).
  \param P0 (optional) standard pressure. Default: 100000 Pa.
  \param L (optional) latent heat of vaporization for water (J.kg^{-1}).
  Default: 2.501e6 J.kg^{-1}.
  \param cp (optional) specific heat of dry air at constant pressure.
  Default: 1005 J.kg^{-1}.K^{-1}.
  \param r (optional) molar gas constant for air.
  Default: 287.0 J.kg^{-1}.K^{-1}.
*/
template<class T, class TG>
void ComputePotentialTemperatureLW(const Data<T, 4, TG>& Temperature,
                                   const Data<T, 4, TG>& Pressure,
                                   const Data<T, 4, TG>& LiquidWaterContent,
                                   Data<T, 4, TG>& PotentialTemperatureLW,
                                   T P0 = 100000., T L = 2.501e6,
                                   T cp = 1005., T r = 287.)
{
  int h, k, j, i;
  int Nt = PotentialTemperatureLW.GetLength(0);
  int Nz = PotentialTemperatureLW.GetLength(1);
  int Ny = PotentialTemperatureLW.GetLength(2);
  int Nx = PotentialTemperatureLW.GetLength(3);

  for (h = 0; h < Nt; h++)
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            T theta = Temperature(h, k, j, i)
              * pow(P0 / Pressure(h, k, j, i), r / cp);

            PotentialTemperatureLW(h, k, j, i) = theta
              * exp(- L / (cp * Temperature(h, k, j, i))
                    * LiquidWaterContent(h, k, j, i));
          }
}


/*! \brief Computes the liquid water content on the basis of ambient air
  characteristics and additional water.
*/
/*! It is assumed that water is added to ambient water. On the basis of the
  new water content and ambient air characteristics, the new liquid water
  content and the new specific humidity are computed.
  \param Pressure pressure (Pa).
  \param Temperature ambient temperature (K).
  \param AmbientLiquidWaterContent ambient liquid water content (kg.kg^{-1}).
  \param PlumePotentialTemperatureLW liquid water potential temperature in
  the plume (K).
  \param PlumeWaterContent total water content (liquid and vapor water)
  added to ambient air (kg.kg^{-1}).
  \param LiquidWaterContent (output) updated liquid water content (with
  additional water content if option is "total") (kg.kg^{-1}).
  \param option gives whether the plume liquid water content should be added
  to the ambient air liquid water content.
  \param P0 (optional) standard pressure. Default: 100000 Pa.
  \param L (optional) latent heat of vaporization for water (J.kg^{-1}).
  Default: 2.501e6 J.kg^{-1}.
  \param cp (optional) specific heat of dry air at constant pressure.
  Default: 1005 J.kg^{-1}.K^{-1}.
  \param r (optional) molar gas constant for air.
  Default: 287.0 J.kg^{-1}.K^{-1}.
*/
template<class T, class TG>
void WaterDiagnosis(const Data<T, 4, TG>& Pressure,
                    const Data<T, 4, TG>& Temperature,
                    const Data<T, 4, TG>& AmbientLiquidWaterContent,
                    const Data<T, 4, TG>& PlumePotentialTemperatureLW,
                    const Data<T, 4, TG>& PlumeWaterContent,
                    Data<T, 4, TG>& LiquidWaterContent, string option,
                    T P0 = 100000., T L = 2.501e6,
                    T cp = 1005., T r = 287.)
{
  int h, k, j, i;
  int Nt = LiquidWaterContent.GetLength(0);
  int Nz = LiquidWaterContent.GetLength(1);
  int Ny = LiquidWaterContent.GetLength(2);
  int Nx = LiquidWaterContent.GetLength(3);

  T pressure, temperature, temperature_amb, thetal, qw, ql;

  for (h = 0; h < Nt; h++)
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            // Computes the temperature with Newton method.
            // temperature_amb, the ambient temperature, is the first guess.
            temperature_amb = Temperature(h, k, j, i);
            pressure = Pressure(h, k, j, i);
            thetal = PlumePotentialTemperatureLW(h, k, j, i);
            qw = PlumeWaterContent(h, k, j, i);
            temperature = Newton(pressure, temperature_amb,
                                 thetal, qw);

            // Computes the liquid water content in the plume and the total
            // liquid water content.
            ql = max(0., qw - qsat(pressure, temperature));
            LiquidWaterContent(h, k, j, i) = ql;
            if (option == "total")
              LiquidWaterContent(h, k, j, i) +=
                max(0., AmbientLiquidWaterContent(h, k, j, i));
          }
}
