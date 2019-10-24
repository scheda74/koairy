// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Irène Korsakissok
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

// This file is part of the Gaussian models for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_BRIGGSFORMULA_HXX


namespace Polyphemus
{


  //! Estimates the potential temperature gradient.
  /*! Estimates the potential temperature gradient.
    \param stability the stability in [0, 5].
    \param dTdz the potential temperature gradient.
  */
  template<class T>
  void ComputePotentialTemperatureGradient(int stability, T& dTdz)
  {
    // Vertical gradient of the potential temperature.
    if (stability == 0)
      dTdz = -0.015;
    else if (stability == 1)
      dTdz = -0.008;
    else if (stability == 2)
      dTdz = -0.006;
    else if (stability == 3)
      dTdz = 0.;
    else if (stability == 4)
      dTdz = 0.02;
    else if (stability == 5)
      dTdz = 0.035;
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Estimates Briggs static stability parameter.
  /*! Estimates Briggs static stability parameter.
    \param stability stability in [0, 5].
    \param temperature temperature of ambient air (K).
    \param sp_briggs Briggs static stability parameter.
  */
  template<class T>
  void ComputeBriggsParameter(T temperature, int stability, T& sp_briggs)
  {
    const T g = 9.81;
    T dTdz;
    ComputePotentialTemperatureGradient(stability, dTdz);
    sp_briggs = (g * dTdz) / temperature ;
  }


  //! Estimates an effective source height of a hot source.
  /*! It uses the Briggs' plume rise formulation modified in HPDM (Hanna 87).
    \param temperature temperature of ambient air (K).
    \param wind wind velocity (m/s).
    \param inversion_height inversion height (m).
    \param stability stability in [0, 5].
    \param source_velocity efflux speed of gases (m/s).
    \param source_temperature temperature of emitted gases (K).
    \param diameter source diameter (m).
    \param source_height source height (m).
    \param option_breakup true if formulae for plume breakup are used.
    \return The plume rise (m).
  */
  template<class T>
  T ComputeHPDMPlumeRise(T temperature, T wind, T inversion_height,
                         T convective_velocity,
                         T friction_velocity,
                         int stability, T source_velocity,
                         T source_temperature, T diameter,
                         T source_height,
                         bool option_breakup = 0)
  {
    const T g = 9.81;

    T sp_briggs, Fb;
    T Dhz = 0.;

    // Briggs static stability parameter.
    ComputeBriggsParameter(temperature, stability, sp_briggs);

    // Initial buoyancy flux parameter.
    Fb = g * source_velocity * 0.25 * diameter * diameter
      * abs(source_temperature - temperature) / source_temperature;

    if (stability == 4 || stability == 5)
      {
        // Stable.
        T dhz1, dhz2;
        dhz1 = 2.6 * pow(Fb / (wind * sp_briggs), 1. / 3.);
        dhz2 = 4. * pow(Fb, 1. / 4.) * pow(sp_briggs, - 3. / 8.);
        Dhz = min(dhz1, dhz2);
      }
    else
      {
        T dhz3, w_star_2, u_star_2;
        // Unstable or neutral.
        if (Fb < 55)
          Dhz = 21.4 / wind * pow(Fb, 3. / 4.);
        else
          Dhz = 38.71 / wind * pow(Fb, 3. / 5.);
        if (option_breakup)
          if (stability == 0 || stability == 1)
            {
              w_star_2 = convective_velocity * convective_velocity;
              dhz3 = 4.3 * pow((Fb / (wind * w_star_2)), 3. / 5.)
                * pow(inversion_height, 2. / 5.);
              Dhz = min(Dhz, dhz3);
            }
          else if (stability == 3)
            {
              u_star_2 = friction_velocity * friction_velocity;
              dhz3 = 1.54 * pow(Fb / (wind * u_star_2), 2. / 3.)
                * pow(source_height, 1. / 3.);
              Dhz = min(Dhz, dhz3);
            }
      }
    return Dhz;
  }


  //! Estimates plume rise.
  /*! It uses Holland (53) formulation modified by Stumke (63).
    \param temperature temperature of ambient air (K).
    \param wind wind velocity (m/s).
    \param source_velocity efflux speed of gases (m/s).
    \param source_temperature temperature of emitted gases (K).
    \param diameter source diameter (m).
    \return The plume rise (m).
  */
  template<class T>
  T ComputeHollandPlumeRise(T temperature, T wind,
                            T source_velocity,
                            T source_temperature,
                            T diameter)
  {
    if (source_temperature < temperature)
      {
	cout << "Warning, Stack exit temperature lower than ambiant temperature, plume rise set to 0" << endl;
	return 0;
      }
    else
      {
    return 1.5 * diameter * source_velocity / wind
      + 65 * pow(diameter, 1.5) * pow((source_temperature - temperature)
                                      / (source_temperature), 0.25) / wind;
      }
  }


  //! Estimates plume rise.
  /*! It uses Concawe (66) formulation.
    \param temperature temperature of ambient air (K).
    \param wind wind velocity (m/s).
    \param source_velocity efflux speed of gases (m/s).
    \param source_temperature temperature of emitted gases (K).
    \param diameter source diameter (m).
    \return The plume rise (m).
  */
  template<class T>
  T ComputeConcawePlumeRise(T temperature, T wind, T inversion_height,
                            T convective_velocity,
                            T friction_velocity,
                            int stability, T source_velocity,
                            T source_temperature, T diameter,
                            T source_height,
                            bool option_breakup = 0)
  {
    if (wind > 1.)
      {
        T D;
	if (source_temperature < temperature)
	  {
	    cout << "Warning, Stack exit temperature lower than ambiant temperature, plume rise set to 0" << endl;
	    D = 0;
	  }
	else
	  {
        D = 228.19 * source_velocity * diameter * diameter
          * (source_temperature - temperature);
      }
        return 0.071 * pow(D, 0.55) / pow(wind, 0.67);
      }
    else
      return ComputeHPDMPlumeRise(temperature, wind, inversion_height,
                                  convective_velocity,
                                  friction_velocity,
                                  stability, source_velocity,
                                  source_temperature, diameter,
                                  source_height,
                                  option_breakup);
  }


  //! Computes horizontal diffusion parameter over rural area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeRuralPlumeHorizontalSigma(T distance, int stability)
  {
    if (stability == 0)
      return 0.22 * distance / sqrt(1. + 1.e-4 * distance);
    else if (stability == 1)
      return 0.16 * distance / sqrt(1. + 1.e-4 * distance);
    else if (stability == 2)
      return 0.11 * distance / sqrt(1. + 1.e-4 * distance);
    else if (stability == 3)
      return 0.08 * distance / sqrt(1. + 1.e-4 * distance);
    else if (stability == 4)
      return 0.06 * distance / sqrt(1. + 1.e-4 * distance);
    else if (stability == 5)
      return 0.04 * distance / sqrt(1. + 1.e-4 * distance);
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  /*! \brief Computes distance corresponding to a given horizontal diffusion
    parameter over rural area.*/
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param sigma horizontal plume-dispersion parameter (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The distance distance downwind of the source (m).
  */
  template<class T>
  T ComputeRuralPlumeHorizontalDistance(T sigma, int stability)
  {
    if (stability == 0)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.22 * 0.22 / (sigma * sigma)))
        / (2 * 0.22 * 0.22 / (sigma * sigma));
    else if (stability == 1)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.16 * 0.16 / (sigma * sigma)))
        / (2 * 0.16 * 0.16 / (sigma * sigma));
    else if (stability == 2)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.11 * 0.11 / (sigma * sigma)))
        / (2 * 0.11 * 0.11 / (sigma * sigma));
    else if (stability == 3)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.08 * 0.08 / (sigma * sigma)))
        / (2 * 0.08 * 0.08 / (sigma * sigma));
    else if (stability == 4)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.06 * 0.06 / (sigma * sigma)))
        / (2 * 0.06 * 0.06 / (sigma * sigma));
    else if (stability == 5)
      return (1.e-4 + sqrt(1.e-8 + 4 * 0.04 * 0.04 / (sigma * sigma)))
        / (2 * 0.04 * 0.04 / (sigma * sigma));
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Computes vertical diffusion parameter over rural area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The vertical plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeRuralPlumeVerticalSigma(T distance, int stability)
  {
    if (stability == 0)
      return 0.2 * distance;
    else if (stability == 1)
      return 0.12 * distance;
    else if (stability == 2)
      return 0.08 * distance / sqrt(1. + 2.e-4 * distance);
    else if (stability == 3)
      return 0.06 * distance / sqrt(1. + 1.5e-3 * distance);
    else if (stability == 4)
      return 0.03 * distance / (1. + 3.e-4 * distance);
    else if (stability == 5)
      return 0.016 * distance / (1. + 3.e-4 * distance);
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Computes derivative of vertical diffusion parameter over rural area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The derivative of vertical plume-dispersion parameter
    with respect to the distance downwind of the source.
  */
  template<class T>
  T DifferentiateRuralPlumeVerticalSigma(T distance, int stability)
  {
    if (stability == 0)
      return 0.2;
    else if (stability == 1)
      return 0.12;
    else if (stability == 2)
      {
        T diff_factor = sqrt(1. + 2.e-4 * distance) * (1. + 2.e-4 * distance);
        return 0.08 * (1. + 1.e-4 * distance) / diff_factor;
      }
    else if (stability == 3)
      {
        T diff_factor = sqrt(1. + 1.5e-3 * distance) *
          (1. + 1.5e-3 * distance);
        return 0.03 * (2. + 1.5e-3 * distance) / diff_factor;
      }
    else if (stability == 4)
      {
        T diff_factor = (1. + 3.e-4 * distance) * (1. + 3.e-4 * distance);
        return 0.03 * (1. + 6.e-3 *  distance) / diff_factor;
      }
    else if (stability == 5)
      {
        T diff_factor = (1. + 3.e-4 * distance) * (1. + 3.e-4 * distance);
        return 0.016 * (1. + 6.e-3 *  distance) / diff_factor;
      }
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Computes horizontal diffusion parameter over urban area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeUrbanPlumeHorizontalSigma(T distance, int stability)
  {
    if (stability == 0)
      return 0.32 * distance / sqrt(1. + 4.e-4 * distance);
    else if (stability == 1)
      return 0.32 * distance / sqrt(1. + 4.e-4 * distance);
    else if (stability == 2)
      return 0.22 * distance / sqrt(1. + 4.e-4 * distance);
    else if (stability == 3)
      return 0.16 * distance / sqrt(1. + 4.e-4 * distance);
    else if (stability == 4)
      return 0.11 * distance / sqrt(1. + 4.e-4 * distance);
    else if (stability == 5)
      return 0.11 * distance / sqrt(1. + 4.e-4 * distance);
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  /*! \brief Computes time corresponding to a given horizontal diffusion
    parameter over urban area.*/
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param sigma horizontal plume-dispersion parameter (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The  distance downwind of the source (m).
  */
  template<class T>
  T ComputeUrbanPlumeHorizontalDistance(T sigma, int stability)
  {
    if (stability == 0 || stability == 1)
      return (4.e-4 + sqrt(16.e-8 + 4 * 0.32 * 0.32 / (sigma * sigma)))
        / (2 * 0.32 * 0.32 / (sigma * sigma));
    else if (stability == 2)
      return (4.e-4 + sqrt(16.e-8 + 4 * 0.22 * 0.22 / (sigma * sigma)))
        / (2 * 0.22 * 0.22 / (sigma * sigma));
    else if (stability == 3)
      return (4.e-4 + sqrt(16.e-8 + 4 * 0.16 * 0.16 / (sigma * sigma)))
        / (2 * 0.16 * 0.16 / (sigma * sigma));
    else if (stability == 4 || stability == 5)
      return (4.e-4 + sqrt(16.e-8 + 4 * 0.11 * 0.11 / (sigma * sigma)))
        / (2 * 0.11 * 0.11 / (sigma * sigma));
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Computes vertical diffusion parameter over urban area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The vertical plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeUrbanPlumeVerticalSigma(T distance, int stability)
  {
    if (stability == 0)
      return 0.24 * distance * sqrt(1. + 1.e-3 * distance);
    else if (stability == 1)
      return 0.24 * distance * sqrt(1. + 1.e-3 * distance);
    else if (stability == 2)
      return 0.2 * distance;
    else if (stability == 3)
      return 0.14 * distance / sqrt(1. + 3.e-4 * distance);
    else if (stability == 4)
      return 0.08 * distance / sqrt(1. + 1.5e-3 * distance);
    else if (stability == 5)
      return 0.08 * distance / sqrt(1. + 1.5e-3 * distance);
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  //! Computes derivative of vertical diffusion parameter over urban area.
  /*! It uses Briggs' dispersion parameterizations for open country with
    Pasquill stability classes.
    \param distance distance downwind of the source (m).
    \param stability Pasquill stability classes in [0, 5].
    \return The derivative of vertical plume-dispersion parameter
    with respect to the distance downwind of the source.
  */
  template<class T>
  T DifferentiateUrbanPlumeVerticalSigma(T distance, int stability)
  {
    if (stability == 0 || stability == 1)
      return 0.12 * (2. + 3.e-3 * distance) / sqrt(1. + 1.e-3 * distance);
    else if (stability == 2)
      return 0.20;
    else if (stability == 3)
      {
        T diff_factor = sqrt(1. + 3e-4 * distance) * (1. + 3.e-4 * distance);
        return 0.07 * (2. + 3.e-4 * distance) / diff_factor;
      }
    else if (stability == 4 || stability == 5)
      {
        T diff_factor = sqrt(1. + 1.5e-3 * distance) *
          (1. + 1.5e-3 * distance);
        return 0.04 * (2. + 1.5e-3 *  distance) / diff_factor;
      }
    else
      throw string("Stability index should be in [0, 5], but ")
        + to_str(stability) + " was provided.";
  }


  /*! Computes horizontal diffusion parameter for situations of normal
    diffusion.
  */
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeDouryNormalPlumeHorizontalSigma(T t)
  {
    if (t >= 0. && t < 240.)
      return pow(0.405 * t, 0.859);
    else if (t >= 240. && t < 97000.)
      return pow(0.135 * t, 1.130);
    else if (t >= 97000. && t < 508000.)
      return pow(0.463 * t, 1.000);
    else if (t >= 508000. && t < 1300000.)
      return pow(6.5 * t, 0.824);
    else
      return pow(2.e+05 * t, 0.5);
  }


  /*! Computes transfer time corresponding to a given horizontal diffusion
    parameter for situations of normal diffusion.
  */
  /*! It uses Doury's dispersion parameterizations.
    \param sigma horizontal plume-dispersion parameter (m).
    \return The transfert time from source (s).
  */
  template<class T>
  T ComputeDouryNormalPlumeHorizontalTime(T sigma)
  {
    if (sigma >= 0. && sigma < 50.9)
      return pow(sigma, 1. / 0.859) / 0.405;
    else if (sigma >= 50.9 && sigma < 44908.)
      return pow(sigma, 1. / 1.130) / 0.135;
    else if (sigma >= 44908. && sigma < 235204.)
      return sigma / 0.463;
    else if (sigma >= 235204. && sigma < 510187)
      return pow(sigma, 1. / 0.824) / 6.5;
    else
      return sigma * sigma / 2.e+05;
  }


  /*! Computes vertical diffusion parameter for situations of normal
    diffusion.
  */
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The vertical plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeDouryNormalPlumeVerticalSigma(T t)
  {
    if (t >= 0. && t < 240.)
      return pow(0.42 * t, 0.814);
    else if (t >= 240. && t < 3280.)
      return pow(t, 0.685);
    else
      return sqrt(20. * t);
  }


  /*! Computes derivative pf vertical diffusion parameter for situations of
    normal diffusion.
  */
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The derivative of vertical plume-dispersion parameter (m)
    with respect to the distance from the source.
  */
  template<class T>
  T DifferentiateDouryNormalPlumeVerticalSigma(T t)
  {
    if (t >= 0. && t < 240.)
      return (0.42 * 0.814) * pow(0.42 * t, (0.814 - 1.));
    if (t >= 240. && t < 3280.)
      return 0.685 * pow(t, (0.685 - 1.));
    else
      return 10. / sqrt(20. * t);
  }


  //! Computes horizontal diffusion parameter for situations of low diffusion.
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeDouryLowPlumeHorizontalSigma(T t)
  {
    if (t >= 0. && t < 240.)
      return pow(0.405 * t, 0.859);
    else if (t >= 240. && t < 97000.)
      return pow(0.135 * t, 1.130);
    else if (t >= 97000. && t < 508000.)
      return pow(0.463 * t, 1.000);
    else if (t >= 508000. && t < 1300000.)
      return pow(6.5 * t, 0.824);
    else
      return pow(2.e+05 * t, 0.5);
  }


  //! Computes vertical diffusion parameter for situations of low diffusion.
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The vertical plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeDouryLowPlumeVerticalSigma(T t)
  {
    return sqrt(0.20 * t);
  }


  /*! Computes derivative of vertical diffusion parameter for situations of
    low diffusion.
  */
  /*! It uses Doury's dispersion parameterizations.
    \param t transfert time from source (s).
    \return The derivative of vertical plume-dispersion parameter (m)
    with respect to the distance from the source.
  */
  template<class T>
  T DifferentiateDouryLowPlumeVerticalSigma(T t)
  {
    return 0.10 / sqrt(0.20 * t);
  }


  //! Computes the stability class corresponding to a Monin-Obukhov length.
  template<class T>
  void ComputeStabilityClass(T L, int& stability, string& stability_class)
  {
    if (-100. < L && L < 0.)
      {
        stability_class = "B";
        stability = 1;
      }
    else if (-1e+05 <= L && L <= -100.)
      {
        stability_class = "C";
        stability = 2;
      }
    else if (abs(L) > 1.e+05)
      {
        stability_class = "D";
        stability = 3;
      }
    else if (10. <= L && L <= 1.e+05)
      {
        stability_class = "E";
        stability = 4;
      }
    else
      {
        stability_class = "F";
        stability = 5;
      }
  }

  //! Computes horizontal downwind wind standard deviation.
  /*! It uses Hanna (82) formula.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The horizontal downwind wind standard deviation (m/s).
  */
  template<class T>
  T ComputeSigma_u(T z, T u_star, T w_star, T h, T L, T f, int stability)
  {
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          return u_star * pow((12. - 0.5 * h / L), 1. / 3.);
        // Stable conditions.
        else if (stability == 4 || stability == 5)
          return 2. * u_star * (1. - z / h);
        // Neutral conditions.
        else
          return 2. * u_star * exp(-3. * f * z / u_star);
      }
    else
      return 0.1 * sqrt(3.6 * u_star * u_star + 0.35 * w_star * w_star);
  }


  //! Computes horizontal crosswind wind standard deviation.
  /*! It uses Hanna (82) formula.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The horizontal wind standard deviation (m/s).
  */
  template<class T>
  T ComputeSigma_v(T z, T u_star, T w_star, T h, T L, T f, int stability)
  {
    T sigma_v;
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          sigma_v = u_star * pow((12. - 0.5 * h / L), 1. / 3.);
        // Stable conditions.
        else if (stability == 4 || stability == 5)
          sigma_v =  max(1.3 * u_star * (1. - z / h), 0.2);
        // Neutral conditions.
        else
          sigma_v =  1.3 * u_star * exp(-2. * f * z / u_star);
      }
    else
      sigma_v = 0.1 * sqrt(3.6 * u_star * u_star + 0.35 * w_star * w_star);
    return sigma_v;
  }


  //! Computes crosswind wind standard deviation.
  /*! It uses Weil (88) formula for convective cases
    and HPDM formulae (Hanna 88) for stable/neutral cases.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The vertical wind standard deviation (m/s).
  */
  template<class T>
  T ComputeHPDMSigma_v(T z, T u_star, T w_star, T h, T L, int stability)
  {
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          return u_star * pow((12. - 0.5 * h / L), 1. / 3.);
        // Neutral/slightly stable conditions.
        else if (L >= 100)
          return 0.7 * sqrt(3.6 * u_star * u_star + 0.35 * w_star * w_star);
        // Stable conditions.
        else
          return max(1.5 * sqrt(3.6 * u_star * u_star
                                + 0.35 * w_star * w_star), 0.5);
      }
    else
      return 0.1 * sqrt(3.6 * u_star * u_star + 0.35 * w_star * w_star);
  }

  //! Computes Lagrangian time scale for downwind dispersion.
  /*! It uses Hanna (84) formula.
    \param z height (m).
    \param h boundary layer height (m).
    \param L Monin-Obukhov length (m).
    \param f Coriolis parameter.
    \param sigma_u downwind wind standard deviation (m/s).
    \param sigma_v crosswind wind standard deviation (m/s).
    \return The Lagrangian time scale for downwind dispersion (s).
  */
  template<class T>
  T ComputeDownwindTimeScale(T z, T u_star, T h, T L, T f,
			     T sigma_u, T sigma_v, int stability)
  {
    // Unstable conditions.
    if (stability == 0 || stability == 1 || stability == 2)
      return 0.15 * h / sigma_u;
    // Stable conditions.
    else if (stability == 4 || stability == 5)
      return 0.15 * h / sigma_u * sqrt(z / h);
    // Neutral conditions.
    else
      return 0.5 * z / (sigma_v * (1 + 15 * f * z / u_star) );
  }

  //! Computes downwind diffusion parameter with similarity theory.
  /*! It uses Weil (88) formula for the general form.
    \param t transfert time from source (s).
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \param stability the stability class.
    \return The downwind dispersion parameter (m).
  */
  template<class T>
  T ComputeDownwindSigma(T t, T z, T u_star, T w_star, T h,
                         T L, T f, bool hpdm, int stability)
  {
    T sigma_u = ComputeSigma_u(z, u_star, w_star, h, L, f, stability);
    T sigma_v;
    if (!hpdm)
      sigma_v = ComputeSigma_v(z, u_star, w_star, h, L, f, stability);
    else
      sigma_v = ComputeHPDMSigma_v(z, u_star, w_star, h, L, stability);
    T tlx = ComputeDownwindTimeScale(z, u_star, h, L, f,
                                     sigma_u, sigma_v, stability);
    T sigma_x = t *  sigma_u / sqrt(1 + 0.5 * t  / tlx);
    return sigma_x;
  }


  //! Computes Lagrangian time scale for crosswind dispersion.
  /*! It uses Hanna (84) formula.
    \param z height (m).
    \param h boundary layer height (m).
    \param L Monin-Obukhov length (m).
    \param f Coriolis parameter.
    \param sigma_v crosswind wind standard deviation (m/s).
    \return The Lagrangian time scale for crosswind dispersion (s).
  */
  template<class T>
  T ComputeCrosswindTimeScale(T z, T u_star, T h, T L, T f, T sigma_v,
                              int stability)
  {
    // Unstable conditions.
    if (stability == 0 || stability == 1 || stability == 2)
      return 0.15 * h / sigma_v;
    // Stable conditions.
    else if (stability == 4 || stability == 5)
      return 0.07 * h / sigma_v * sqrt(z / h);
    // Neutral conditions.
    else
      return 0.5 * z / (sigma_v * (1 + 15 * f * z / u_star));
  }


  //! Computes crosswind diffusion parameter with similarity theory.
  /*! It uses Weil (88) formula for the general form.
    \param t transfert time from source (s).
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The crosswind dispersion parameter (m).
  */
  template<class T>
  T ComputeCrosswindSigma(T t, T z, T u_star, T w_star,
                          T h, T L, T f, bool hpdm, int stability)
  {
    T sigma_v, tly, sigma_y;
    if (!hpdm)
      sigma_v = ComputeSigma_v(z, u_star, w_star, h, L, f, stability);
    else
      sigma_v = ComputeHPDMSigma_v(z, u_star, w_star, h, L, stability);
    tly = ComputeCrosswindTimeScale(z, u_star, h, L, f, sigma_v, stability);
    sigma_y = t *  sigma_v / sqrt(1 + 0.5 * t  / tly);
    return sigma_y;
  }


  /*! Computes time corresponding to a given horizontal crosswind
    diffusion parameter with similarity theory.*/
  /*! It uses Weil (88) formula for the general form.
    \param sigma horizontal plume-dispersion parameter (m).
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The transfert time from source (s).
  */
  template<class T>
  T ComputeHorizontalTime(T sigma, T z, T u_star, T w_star,
                          T h, T L, T f, bool hpdm, int stability)
  {
    T sigma_v;
    if (!hpdm)
      sigma_v = ComputeSigma_v(z, u_star, w_star, h, L, f, stability);
    else
      sigma_v = ComputeHPDMSigma_v(z, u_star, w_star, h, L, stability);
    T tly = ComputeCrosswindTimeScale(z, u_star, h, L, f, sigma_v, stability);
    T delta = 0.25 / (tly * tly) + 4 * sigma_v * sigma_v / (sigma * sigma);
    return (0.5 / tly + sqrt(delta))
      / (2 * sigma_v * sigma_v / (sigma * sigma));
  }


  //! Computes vertical wind standard deviation.
  /*! It uses Weil (88) formula for convective cases
    and Venkatram (84) for stable/neutral cases.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The vertical wind standard deviation (m/s).
  */
  template<class T>
  T ComputeWindVerticalSigma(T z, T u_star, T w_star, T h, T L, int stability)
  {
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          return 0.6 * w_star;
        // Neutral/Stable conditions.
        else
          return 1.3 * u_star * pow((1. - z / h), 3. / 4.);
      }
    else
      return 0.1 * sqrt(1.2 * u_star * u_star + 0.35 * w_star * w_star);
  }


  //! Computes vertical wind standard deviation.
  /*! It uses Weil (88) formula for convective cases
    and HPDM formulae (Hanna 88) for stable/neutral cases.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The vertical wind standard deviation (m/s).
  */
  template<class T>
  T ComputeHPDMWindVerticalSigma(T z, T u_star, T w_star, T h, T L,
                                 int stability)
  {
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          return 0.6 * w_star;
        // Neutral/slightly stable conditions.
        else if (L >= 100)
          return 0.5 * sqrt(1.2 * u_star * u_star + 0.35 * w_star * w_star);
        // Stable conditions.
        else
          return 1.3 * u_star;
      }
    else
      return 0.1 * sqrt(1.2 * u_star * u_star + 0.35 * w_star * w_star);
  }


  /*! Computes Lagrangian time scale for vertical dispersion in unstable
    cases.*/
  /*! It uses Hanna (82) formula.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param h boundary layer height (m).
    \param L Monin-Obukhov length (m).
    \param sigma_w vertical wind standard deviation (m/s).
    \param f Coriolis parameter.
    \return The Lagrangian time scale for vertical dispersion (s).
  */
  template<class T>
  T ComputeVerticalTimeScale(T z, T u_star, T h, T L, T sigma_w, T f,
                             int stability)
  {
    // Unstable conditions.
    if (stability == 0 || stability == 1 || stability == 2)
      return 0.15 * h * (1. - exp(-5. * z / h)) / sigma_w;
    // Stable conditions.
    else if (stability == 4 || stability == 5)
      return 0.10 * h / sigma_w * pow(z / h, 0.8);
    // Neutral conditions.
    else
      return 0.5 * z / (sigma_w * (1 + 15 * f * z / u_star));
  }


  /*! Computes Lagrangian time scale for vertical dispersion in unstable
    cases.*/
  /*! It uses Hanna (82) formula for unstable cases
    and HPDM formulae (Hanna 88) for stable/neutral cases.
    \param z height (m).
    \param h boundary layer height (m).
    \param sigma_w vertical wind standard deviation (m/s).
    \param u_star friction velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param f Coriolis parameter.
    \return The Lagrangian time scale for vertical dispersion (s).
  */
  template<class T>
  T ComputeHPDMVerticalTimeScale(T z, T h, T sigma_w, T u_star,
                                 T L, T f, T temperature, int stability)
  {
    // Stable conditions.
    if (L > 0.)
      {
        T sp_briggs;
        ComputeBriggsParameter(temperature, stability, sp_briggs);
        if (z <= L)
          return z / sigma_w;
        else if (L <= 10.)
          return 0.27 / sqrt(sp_briggs);
        else
          return (z / sigma_w) * (L - 10.) / (z - 10.)
            + (0.27 / sqrt(sp_briggs))  * (z - L) / (z - 10.);
      }
    else if (abs(L) < 100)
      {
        if (z <= abs(L))
          return 0.27 * (z / sigma_w) * (0.55 - 0.38 * z / abs(L));
        else
          return (0.3 * h / sigma_w) * (1. - exp(-5 * z / h))
            - 0.0003 * exp(8 * z / h);
      }
    else
      return 0.15 * h * (1. - exp(-5. * z / h)) / sigma_w;
  }


  //! Computes vertical diffusion parameter with similarity theory.
  /*! It uses Weil (88) formula for the general form in unstable
    (convective) cases, and Irwin (79) for stable/neutral cases.
    \param t transfert time from source (s).
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T ComputeVerticalSigma(T t, T z, T u_star, T w_star,
                         T h, T L, T f, bool hpdm,
                         T temperature, int stability)
  {
    T Fz, sigma_z, tlz, sigma_w;
    if (!hpdm)
      sigma_w = ComputeWindVerticalSigma(z, u_star, w_star, h, L, stability);
    else
      sigma_w = ComputeHPDMWindVerticalSigma(z, u_star, w_star, h,
                                             L, stability);
    // Unstable conditions.
    if (stability == 0 || stability == 1 || stability == 2 || t >= 240.)
      {
        if (!hpdm)
          tlz = ComputeVerticalTimeScale(z, u_star, h, L, sigma_w,
                                         f, stability);
        else
          tlz = ComputeHPDMVerticalTimeScale(z, h, sigma_w, u_star,
                                             L, f, temperature, stability);
        sigma_z = t *  sigma_w / sqrt(1 + 0.5 * t / tlz);
      }
    // Neutral/Stable conditions.
    else
      {
        if (z < 50)
          Fz = 1. / (1. + 0.9 * sqrt(t / 50.));
        else
          Fz = 1. / (1. + 0.945 * pow(0.1 * t, 0.806));
        sigma_z = t * sigma_w * Fz;
      }
    return sigma_z;
  }


  /*! Computes derivative of the vertical diffusion parameter with similarity
    theory.*/
  /*! It uses Weil (88) formula for the general form in unstable
    (convective) cases, and Irwin (79) for stable/neutral cases.
    \param t transfert time from source (s).
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param w_star convective velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \return The horizontal plume-dispersion parameter (m).
  */
  template<class T>
  T DifferentiateVerticalSigma(T t, T z, T u, T u_star,
                               T w_star, T h, T L, T f, bool hpdm,
                               T temperature, int stability)
  {
    T Fz, dFz, tlz, sigma_w;
    if (!hpdm)
      sigma_w = ComputeWindVerticalSigma(z, u_star, w_star, h, L, stability);
    else
      sigma_w = ComputeHPDMWindVerticalSigma(z, u_star, w_star,
                                             h, L, stability);
    // Unstable conditions.
    if (h / L < - 0.3)
      {
        if (!hpdm)
          tlz = ComputeVerticalTimeScale(z, u_star, h, L,
                                         sigma_w, f, stability);
        else
          tlz = ComputeHPDMVerticalTimeScale(z, h, sigma_w, u_star,
                                             L, f, temperature, stability);
        T dfz = - 0.25 / (u * tlz * pow((1. + 0.5 * t / tlz), 3. / 2.));
        T fz = 1. / sqrt(1 + 0.5 * t / tlz);
        return sigma_w * t * dfz + sigma_w * fz / u;
      }
    // Neutral/Stable conditions.
    else
      {
        if (z < 50)
          {
            Fz = 1. / (1. + 0.9 * sqrt(t / 50.));
            dFz = - Fz * Fz * 0.09 / (u * sqrt(t / 50.));
          }
        else
          {
            Fz = 1. / (1. + 0.945 * pow(0.1 * t, 0.806));
            dFz = - Fz * Fz * 0.762 / (100. * u * pow((t / 100), 0.194));
          }
        return t * sigma_w * dFz + sigma_w * Fz / u;
      }
  }


  //! Computes vertical wind standard deviation.
  /*! It uses Hanna (82) formulae.
    \param z height (m).
    \param u_star friction velocity (m/s).
    \param L Monin-Obukhov length (m).
    \param h boundary layer height (m).
    \param f Coriolis parameter.
    \return The vertical wind standard deviation (m/s).
  */
  template<class T>
  T ComputeHannaWindVerticalSigma(T z, T u_star, T w_star, T h, T L, T f,
                                  int stability)
  {
    if (z / h < 1.)
      {
        // Unstable conditions.
        if (stability == 0 || stability == 1 || stability == 2)
          {
            if (z / h  < 0.03)
              return 0.96 * w_star * pow((3. * z / h - L / h), 1. / 3.);
            else if (z / h < 0.4)
              return min(0.96 * w_star * pow((3. * z / h - L / h), 1. / 3.),
                         0.763 * w_star * pow(z / h, 0.175));
            else if (z / h < 0.96)
              return 0.722 * w_star * pow((1. - z / h), 0.207);
            else
              return 0.37 * w_star;
          }
        // Stable conditions.
        else if (stability == 4 || stability == 5)
          return 1.3 * u_star * (1. - z / h);
        // Neutral conditions.
        else
          return 1.3 * u_star * exp(-2. * f * z / u_star);
      }
    else
      return 0.1 * sqrt(1.2 * u_star * u_star + 0.35 * w_star * w_star);
  }


  //! Computes a concentration at a given point.
  /*! Gaussian plume model for a continuous point source in uniform flow with
    homogeneous turbulence.
    \param wind mean transport velocity across the plume (m/s).
    \param inversion_height inversion height (m).
    \param effective_height effective source height (m).
    \param rate source strength or emission rate (mass unit/s).
    \param y ordinate (crosswind) (m).
    \param z height (m).
    \param sigma_y horizontal diffusion parameter (m).
    \param sigma_z vertical diffusion parameter (m).
    \param alpha coefficient for Overcamp model (in [-1, 1])..
    \param minimum_volume (optional) minimum volume in which the tracer is
    mixed within one second (m^3/s). It is used to set an upper bound on the
    concentration. Put zero to avoid this threshold. Default: 1.e-6.
    \return The concentration at (\a y, \a z).
    \warning The wind should not be zero.
  */
  template<class T>
  T ComputePlumeConcentration(T wind, T inversion_height, T effective_height,
                              T rate, T y, T z, T sigma_y, T sigma_z, T alpha,
                              T minimum_volume = 1.e-6)
  {
    const T pi = 3.14159265358979323846264;

    // Crosswind contribution.
    T Fy = 0.;
    Fy = exp(-y * y  / (2. * sigma_y * sigma_y));

    // Vertical contribution.
    T Fz = 0.;
    T dzr = z - effective_height;

    T sigma_z_2 = sigma_z * sigma_z;

    Fz = exp(-dzr * dzr / (2. * sigma_z_2));

    // Ground reflection: plume is below the inversion layer.
    T FzG = 0.;
    if (inversion_height == 0.
        || effective_height <= inversion_height)
      if (effective_height - sigma_z <= 0.)
        {
          T dzg = z + effective_height;
          FzG = alpha * exp(-dzg * dzg / (2. * sigma_z_2));
        }

    T concentration = 0.;

    // If there is an inversion layer.
    if (inversion_height > 0.)
      {
        T FzI = 0.;

        // When plume is below inversion layer.
        if (effective_height <= inversion_height)
          {
            // Reflection occurs only when plume touches the inversion layer.
            if (effective_height + sigma_z >= inversion_height)
              {
                T dzi = z + effective_height - 2. * inversion_height;
                T dzig = z - effective_height + 2. * inversion_height;
                T dzgi = z - effective_height - 2. * inversion_height;

                FzI += exp(-dzi * dzi / (2. * sigma_z_2))
                  + exp(-dzgi * dzgi / (2. * sigma_z_2))
                  + exp(-dzig * dzig / (2. * sigma_z_2));
              }
            // Plume only impacts below the inversion layer.
            if (z <= inversion_height)
              {
                // Far field.
                if (sigma_z > 1.5 * inversion_height)
                  concentration = rate * Fy
                    / (sqrt(2. * pi) * wind * sigma_y * inversion_height);
                else
                  concentration = rate * Fy * (Fz + FzG + FzI)
                    / (2. * pi * wind * sigma_y * sigma_z);
              }
            else
              concentration = 0.;
          }

        // When plume is above inversion layer.
        if (effective_height > inversion_height)
          {
            // Reflection on the inversion layer.
            T dzi = z + effective_height - 2. * inversion_height;
            FzI += exp(-dzi * dzi / (2. * sigma_z_2));
            // Plume only impacts above the inversion layer.
            if (z >= inversion_height)
              concentration = rate * Fy * (Fz + FzG + FzI)
                / (2. * pi * wind * sigma_y * sigma_z);
            else
              concentration = 0.;
          }
      }
    else
      {
        concentration = rate * Fy * (Fz + FzG)
          / (2. * pi * wind * sigma_y * sigma_z);
      }

    return minimum_volume != 0. ? min(concentration, rate / minimum_volume)
      : concentration;
  }


  //! Computes a concentration at a given point from a line source.
  /*! Gaussian plume model for a continuous line source in uniform flow with
    homogeneous turbulence.
    \param[in] wind mean transport velocity across the plume (m/s).
    \param[in] inversion_height inversion height (m).
    \param[in] effective_height effective source height (m).
    \param[in] rate source strength or emission rate (mass unit/s).
    \param[in] x abscisa (m).
    \param[in] y ordinate (crosswind) (m).
    \param[in] z height (m).
    \param[in] sigma_y1 horizontal diffusion parameter (m) from the first
    source extremity.
    \param[in] sigma_y1 horizontal diffusion parameter (m) from the second
    source extremity.
    \param[in] sigma_z vertical diffusion parameter (m).
    \param[in] alpha coefficient for Overcamp model (in [-1, 1])..
    \param[in] y1 ordinate of the source first extremity (m).
    \param[in] y2 ordinate of the source second extremity (m).
    \param[in] cos_theta cosine of the wind angle.
    \param[in] sin_theta sine of the wind angle.
    \param[in] minimum_volume (optional) minimum volume in which the tracer is
    mixed within one second (m^3/s). It is used to set an upper bound on the
    concentration. Put zero to avoid this threshold. Default: 1.e-6.
    \return The concentration at (\a y, \a z).
    \warning The wind and the wind angle cosine should not be zero.
  */
  template<class T>
  T ComputePlumeLineConcentration(T wind, T inversion_height,
                                  T effective_height, T rate, T x, T y, T z,
                                  T sigma_y1, T sigma_y2, T sigma_z, T alpha,
                                  T y1, T y2, T cos_theta, T sin_theta,
                                  T minimum_volume = 1.e-6)
  {
    const T pi = 3.14159265358979323846264;

    // Constant factor
    T C = rate / (2 * sqrt(2. * pi) * wind * cos_theta * sigma_z);

    // Crosswind contribution.
    T Fy = 0.;

    //Sign of cos_theta * sin_theta
    T sign_cos_sin_theta = 1;

    if (cos_theta * sin_theta < 0)
      sign_cos_sin_theta = -1;

    /* Modif VR */
//    T Fy1 = - sign_cos_sin_theta;
//    T Fy2 = - sign_cos_sin_theta;
//    if (sigma_y2 != 0)
//      Fy2 = erf(((y - y2) * cos_theta - x * sin_theta) / (sqrt(2) * sigma_y2));
//    if (sigma_y1 != 0)
//      Fy1 = erf(((y - y1) * cos_theta - x * sin_theta) / (sqrt(2) * sigma_y1));
//    Fy = Fy1 - Fy2;

    if(sigma_y1 != 0 && sigma_y2 != 0)
      Fy = erf(((y - y1) * cos_theta - x * sin_theta) / (sqrt(2) * sigma_y1))
        - erf(((y - y2) * cos_theta - x * sin_theta) / (sqrt(2) * sigma_y2));
    else if (sigma_y1 == 0)
      Fy = - sign_cos_sin_theta - erf(((y - y2) * cos_theta - x * sin_theta)
                                      / (sqrt(2) * sigma_y2));
    else if (sigma_y2 == 0)
      Fy = erf(((y - y1) * cos_theta - x * sin_theta) / (sqrt(2) * sigma_y1))
        + sign_cos_sin_theta;

    /* Modif VR */

    // Vertical contribution.
    T Fz = 0.;
    T dzr = z - effective_height;

    T sigma_z_2 = sigma_z * sigma_z;

    Fz = exp(-dzr * dzr / (2. * sigma_z_2));

    // Ground reflection: plume is below the inversion layer.
    T FzG = 0.;
    if (inversion_height == 0.
        || effective_height <= inversion_height)
      if (effective_height - sigma_z <= 0.)
        {
          T dzg = z + effective_height;
          FzG = alpha * exp(-dzg * dzg / (2. * sigma_z_2));
        }

    T concentration = 0.;

    // If there is an inversion layer.
    if (inversion_height > 0.)
      {
        T FzI = 0.;

        // When plume is below inversion layer.
        if (effective_height <= inversion_height)
          {
            // Reflection occurs only when plume touches the inversion layer.
            if (effective_height + sigma_z >= inversion_height)
              {
                T dzi = z + effective_height - 2. * inversion_height;
                T dzig = z - effective_height + 2. * inversion_height;
                T dzgi = z - effective_height - 2. * inversion_height;

                FzI += exp(-dzi * dzi / (2. * sigma_z_2))
                  + exp(-dzgi * dzgi / (2. * sigma_z_2))
                  + exp(-dzig * dzig / (2. * sigma_z_2));
              }
            // Plume only impacts below the inversion layer.
            if (z <= inversion_height)
              {
                // Far field.
                if (sigma_z > 1.5 * inversion_height)
                  concentration = Fy * sigma_z * sqrt(2 * pi) /
                    inversion_height;
                else
                  concentration = Fy * (Fz + FzG + FzI);
              }
            else
              concentration = 0.;
          }

        // When plume is above inversion layer.
        if (effective_height > inversion_height)
          {
            // Reflection on the inversion layer.
            T dzi = z + effective_height - 2. * inversion_height;
            FzI += exp(-dzi * dzi / (2. * sigma_z_2));
            // Plume only impacts above the inversion layer.
            if (z >= inversion_height)
              concentration = Fy * (Fz + FzG + FzI);
            else
              concentration = 0.;
          }
      }
    else
      concentration = Fy * (Fz + FzG);

    // Multilication by the constant.
    concentration = C * concentration;

    return minimum_volume != 0. ? min(concentration, rate / minimum_volume)
      : concentration;
  }
}

#define POLYPHEMUS_FILE_MODELS_BRIGGSFORMULA_HXX
#endif
