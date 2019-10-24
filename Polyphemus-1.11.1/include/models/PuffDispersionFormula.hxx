#ifndef POLYPHEMUS_FILE_MODELS_PUFFDISPERSIONFORMULA_HXX


namespace Polyphemus
{
  /* Compute buoyancy frequency 
  temperature in K */
  template<class T>
  void ComputeBuoyancyFrequency(T temperature,
				int stability,
				T& N2)
  {
    const T g = 9.81;
    T dTdz;
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
      throw string("In puffdispersion Stability index should be in [0, 5], but ")
	+ to_str(stability) + " was provided.";

    N2 = (g * dTdz) / temperature;
  }

  /* Compute Plume rise in stable/neutral conditions.
  x: downwind distance to the source
  N2 buoyancy frequency */
  template<class T>
  void ComputeStablePlumeRise(T source_velocity,
			      T source_diameter,
			      T source_temperature,
			      T source_height,
			      T temperature,
			      T x,
			      T lmo,
			      T wind_speed,
			      T friction_velocity,
			      T N2,
			      int stability,
			      T& plume_rise)
  {
    T Fb, Fm;
    // T Dh_max;
    // T Dh_max1, Dh_max2;
    T xf;
    
    if (wind_speed < 1.)
      wind_speed = 1.;

    const T g = 9.81;
    Fb = 0.25 * g * source_velocity * source_diameter
      * source_diameter * (source_temperature - temperature) / source_temperature;
    Fm = (temperature / source_temperature) * wind_speed * wind_speed * 
      source_diameter * source_diameter / 4.;

    T Bi = 0.6;
    T Bj = 1./3. + wind_speed / source_velocity;

    if (Fb > 55.)
      xf = 34. * pow(Fb, 0.4);
    else
      xf = 14. * pow(Fb, 0.625);
    
    if (Fb > 0.)
      xf = 3.5 * xf;
    else
      xf = 4. * source_diameter * pow((source_velocity + 3. * wind_speed), 2.) / (wind_speed * source_velocity);

    T zsf, zsf1, zsf2;
    zsf1 = 4. * pow(Fb, 1./4.) / pow(N2, 3./8.);
    zsf2 = 2.6 * pow(Fb / (wind_speed * N2), 1./3.);
    zsf = min(zsf1, zsf2);

    T dz_tmp;
    dz_tmp = pow( (3. * Fm * x / (Bj * Bj * wind_speed * wind_speed)) + 
		  (3. * Fb * x * x / (2. * Bi * Bi * pow(source_velocity, 3.))), 1./3.);
    plume_rise = min(dz_tmp, zsf);
  }

  
  /* Compute Plume rise in unstable conditions.
  x: downwind distance to the source
  N2 buoyancy frequency */
  template<class T>
  void ComputeUnstablePlumeRise(T source_velocity,
				T source_diameter,
				T source_temperature,
				T source_height,
				T temperature,
				T x,
				T lmo,
				T wind_speed,
				T friction_velocity,
				T N2,
				T convective_velocity,
				T inversion_height,
				int stability,
				T& plume_rise)
  {
    T Fb, Fm;
    T Dh_max;
    T dhf;
    const T g = 9.81;
    if (wind_speed < 1.)
      wind_speed = 1.;

    Fb = 0.25 * g * source_velocity * source_diameter
      * source_diameter * abs(source_temperature - temperature) / source_temperature;
    Fm = (temperature / source_temperature) * wind_speed * wind_speed 
      * source_diameter * source_diameter / 4.;

    T xf;
    if (Fb > 55.)
      xf = 34. * pow(Fb, 0.4);
    else
      xf = 14. * pow(Fb, 0.625);

    if (Fb  > 0.)
      xf = 3.5 * xf;
    else
      xf = 4. * source_diameter * pow((source_velocity + 3. * wind_speed), 2.)
	/ (wind_speed * source_velocity);

    if (x > xf)
      x = xf;

    T Bj;
    Bj = 1. / 3. + wind_speed / source_velocity;

    T dz_tmp, H, dz_breakup;
    dz_tmp = 1.6 * pow(Fb * x * x, 1./3.) / wind_speed + 
      pow((3. * Fm * x / (Bj * Bj * pow(wind_speed, 3.))), 1./3.);
 
    H = inversion_height;

    if (stability == 0 || stability == 1)
      dz_breakup = 4.3 * pow(Fb / (wind_speed * convective_velocity * convective_velocity), 3. / 5.)
	* pow(H, 2./5.);
    else if (stability == 3)
      dz_breakup = 1.54 * pow(Fb / (wind_speed * friction_velocity * friction_velocity), 2. / 3.)
	* pow(source_height, 1. / 3.);
    else
      dz_breakup = dz_tmp;
    plume_rise = min(dz_tmp, dz_breakup);

  }

  /* Compute Distance corresponding to a given 
     Plume rise in unstable conditions.*/
  template<class T>
  T ComputeUnstablePlumeRiseDistance(T plume_rise,
				     T source_velocity,
				     T source_diameter,
				     T source_temperature,
				     T temperature,
				     T wind_speed)
  {
    T Fb, Fm;
    T Dh_max;
    T dhf;
    const T g = 9.81;
    if (wind_speed < 1.)
      wind_speed = 1.;

    Fb = 0.25 * g * source_velocity * source_diameter
      * source_diameter * abs(source_temperature - temperature) / source_temperature;
    Fm = (temperature / source_temperature) * wind_speed * wind_speed 
      * source_diameter * source_diameter / 4.;

    T Bj;
    Bj = 1. / 3. + wind_speed / source_velocity;

    T x_rise, delta, x3;
    delta = pow(3. * Fm / (Bj * Bj * pow(wind_speed, 3.)), 2. / 3.)
      + 4. * (1.6 / wind_speed * pow(Fb, 1./3.) * plume_rise);
    x3 = (-1. * pow(3. * Fm / (Bj * Bj * pow(wind_speed, 3.)), 1. / 3.)
	  + sqrt(delta)) / (2. * 1.6 / wind_speed * pow(Fb, 1./3.));
    x_rise = pow(x3, 3.);    		    

    return x_rise;

  }

  /* Compute distance corresponding to a given
     Plume rise in stable/neutral conditions.*/
  template<class T>
  T ComputeStablePlumeRiseDistance(T plume_rise,
				   T source_velocity,
				   T source_diameter,
				   T source_temperature,
				   T temperature,
				   T wind_speed)
  {
    T Fb, Fm;
    T xf;
    
    if (wind_speed < 1.)
      wind_speed = 1.;

    const T g = 9.81;
    Fb = 0.25 * g * source_velocity * source_diameter
      * source_diameter * (source_temperature - temperature) / source_temperature;
    Fm = (temperature / source_temperature) * wind_speed * wind_speed * 
      source_diameter * source_diameter / 4.;

    T Bi = 0.6;
    T Bj = 1./3. + wind_speed / source_velocity;

    T x_rise, delta;
    delta = pow(3. * Fm / (Bj * Bj * wind_speed * wind_speed), 2.)
      + 4. * 3. * Fb * pow(plume_rise, 3.) / (2. * Bi * Bi * pow(source_velocity, 3.));

    x_rise = (-3. * Fm / (Bj * Bj * wind_speed * wind_speed) + sqrt(delta))
      * Bi * Bi * pow(source_velocity, 3.) / (3. * Fb);

    return x_rise;
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
  T ComputeLagrangianVerticalTimeScale(T z, T h, T sigma_w, T u_star,
				       T L, T f, T temperature, int stability)
  {
    // Stable conditions.
    if (L > 0.)
      {
	T sp_briggs;
	ComputeBuoyancyFrequency(temperature, stability, sp_briggs);
	if (z <= L)
	  return z / sigma_w;
	else if (L <= 10.)
	  return 0.27 / sqrt(sp_briggs);
	else
	  return (z / sigma_w) * (L - 10.) / (z - 10.)
	    + (0.27 / sqrt(sp_briggs))  * (z - L) / (z - 10.);
      }
    else if (abs(L) < 100. && z <= abs(L))
      return 0.27 * (z / sigma_w) * (0.55 - 0.38 * z / abs(L));
    else if (abs(L) < 100. && z > abs(L) && z < h)
      return (0.3 * h / sigma_w) * (1. - exp(-5. * z / h))
	- 0.0003 * exp(8. * z / h);
    else
      return 0.15 * h * (1. - exp(-5. * z / h)) / sigma_w;
  }

  template<class T>
  T ComputeSigmaMin(T lmo)
  {
    T sigma_min;

    if (lmo <= 10.)
      sigma_min = 0.1;
    else if (lmo > 10. && lmo <= 30.)
      sigma_min = 0.01 * (lmo - 10.);
    else
      sigma_min = 0.2;
    return sigma_min;    
  }

  /* Compute zi, zim and zic.
  Zi: inversion height 
  Zim: mixing depth
  Zic: Convection depth*/
  template<class T>
  void ComputeZi(T w_star, T u_star, T lmo,
		 T& zi, T& zim, T& zic)
  {
    T w_star3, u_star3;
    w_star3 = pow(w_star, 3.);
    u_star3 = pow(u_star, 3.);

    zim = 2300. * u_star3;

    T inv_lmo;
    if (lmo != 0.)
      inv_lmo = 1. / lmo;
    else
      inv_lmo = 999.;

    if (inv_lmo != 0.)
      zic = -1. * w_star3 * 0.4 / (inv_lmo * u_star3);
    else
      zic = 0.;
    zi = max(zim, zic);
  }

  template<class T>
  T ComputeTlv(T zi, T sigmav)
  {
    T Tlv;
    Tlv = (1/78.) * zi / sigmav;
    return Tlv;
  }

  template<class T>
  T ComputeTlw(T zi, T sigmaw, T z)
  {
    T Tlw;
    
    if (zi == 0.)
      Tlw = 0.;
    else
      Tlw = 0.2 * (zi / sigmaw) * (1 - exp(-8. * z / zi) -
				   0.0003 * exp(8.5 * z / zi));
    if (Tlw < 0.)
      Tlw = 0.;
    
    return Tlw;
  }

  template<class T>
  T ComputeN2Adms(T z, T lmo, T u_star, T boundary_height)
  {
    T u_star2 = pow(u_star, 2.);
    T zsu, N2, dz;

    zsu = min(100., boundary_height);
    if (lmo == 0.)
      N2 = 0.;
    else
      {
	if (z <= zsu)
	  {
	    dz = z + 0.3;
	    N2 = u_star2 / (0.4 * 0.4 * lmo)
	      * (1 / dz + 0.7 / lmo + (0.75 / lmo - 0.75 * (dz / lmo - 0.5 / 0.35)
				       * 0.35 / lmo) * exp(-0.35 * dz / lmo));
	  }
	else
	  {
	    if (z <= boundary_height)
	      {
		dz = zsu + 0.3;
		N2 = u_star2 / (0.4 * 0.4 * lmo)
		  * (1 / dz + 0.7 / lmo + (0.75 / lmo - 0.75 * (dz / lmo - 0.5 / 0.35)
					   * 0.35 / lmo) * exp(-0.35 * dz / lmo));
	      }
	    else
	      {
		dz = boundary_height + 0.3;
		N2 = u_star2 / (0.4 * 0.4 * lmo)
		  * (1 / dz + 0.7 / lmo + (0.75 / lmo - 0.75 * (dz / lmo - 0.5 / 0.35)
					   * 0.35 / lmo) * exp(-0.35 * dz / lmo));
	      }
	  }
      }
    return N2;
    
  }
    
  /* Computes Puff lateral dispersion in stable/neutral conditions.
  Based on AMDS gaussian dispersion formula.
  time: time since release (s)
  sigma_v: horizontal crosswind standard deviation
  x: downwind distance to the source
  u10: 10m wind speed*/
  template<class T>
  T ComputeADMSSigmaY(T t,
		      T sigma_v,
		      T boundary_height,
		      T friction_velocity,
		      T source_diameter,
		      T x,
		      T u10,
		      T plume_rise,
		      T Tlv)
  {
    //T sigma_yt, sigma_yw, sigma_ypr;
    T sigma_yw, sigma_ypr;
    T sigma_y;
    T sigma_yi;

    sigma_yw = sigma_v * t / pow((1 + t / (2 * Tlv)), 0.3);
    sigma_ypr = plume_rise * 0.4 / pow(2., 1./2.);
    sigma_yi = source_diameter * source_diameter / 4.;

    sigma_y = pow((sigma_yw * sigma_yw +
		   sigma_ypr * sigma_ypr +
		   sigma_yi * sigma_yi), 1./2.);

    return sigma_y;
  }

  /* Computes Puff lateral dispersion in unstable conditions.
  Based on AMDS gaussian dispersion formula.
  time: time since release (s)
  sigma_v: horizontal crosswind standard deviation
  x: downwind distance to the source
  z: puff center height
  u10: 10m wind speed*/
  template<class T>
  T ComputeADMSConvectiveSigmaY(T t,
				T boundary_height,
				T friction_velocity,
				T convective_velocity,
				T source_diameter,
				T x,
				T u10,
				T wind_,
				T z,
				T plume_rise,
				T Tlv)
  {
    T z0, Twn;
    T sigma_vc2, sigma_vn2;
    T sigma_vc, sigma_vn, sigma_v;
    T sigma_yi;
    T sigma_yc, sigma_yn, sigma_yw, sigma_ypr, sigma_y;

    z0 = 0.5;
    Twn = 1. - 0.8 * (z + z0) / boundary_height;

    sigma_vc2 = 0.3 * convective_velocity * convective_velocity;
    sigma_vn2 = 4. * Twn * Twn * friction_velocity * friction_velocity;

    sigma_vc = sqrt(sigma_vc2);
    sigma_vn = sqrt(sigma_vn2);
    sigma_v = pow(sigma_vc2 + sigma_vn2, 1./2.);

    sigma_yc = sigma_vc * t *
      pow((1 + 0.91 * convective_velocity * t / boundary_height), -1./2.);
    sigma_yn = sigma_vn * t *
      pow((1 + 2.5 * friction_velocity * t / boundary_height), -1./2.);

    T ti;
    ti = t * wind_ / u10;
    sigma_yw = sigma_v * ti / pow((1. + ti / (2. * Tlv)), 0.3);
    sigma_ypr = plume_rise * 0.4 / sqrt(2.);

    sigma_yi = source_diameter * source_diameter / 4.;

    sigma_y = sqrt(sigma_yc * sigma_yc +
		   sigma_yn * sigma_yn +
		   sigma_yw * sigma_yw +
		   sigma_ypr * sigma_ypr +
		   sigma_yi * sigma_yi);
    return sigma_y;
  }

   /* Computes Puff vertical dispersion in unstable conditions.
  Based on AMDS gaussian dispersion formula.
  time: time since release (s)
  tl: Lagrangian time step for vertical diffusion
  sigma_w: vertical crosswind standard deviation*/
  template<class T>
  T ComputeADMSConvectiveSigmaZ(T t,
				T sigma_w,
				T boundary_height,
				T convective_velocity,
				T source_diameter,
				T plume_rise,
				T Tlw)
  {
    T sigma_zpr, sigma_z1, sigma_z2, sigma_z;

    sigma_zpr = plume_rise / 2.;
    sigma_z1 = sigma_w * t * pow((1. + t / (2. * Tlw)), -1. / 2.);
    sigma_z2 = source_diameter * source_diameter / 4.;

    sigma_z = sqrt(sigma_z1 * sigma_z1 +
		   sigma_z2 * sigma_z2 +
		   sigma_zpr * sigma_zpr);

    return sigma_z;    
  }

  /* Computes Puff vertical dispersion in stable/neutral conditions.
  Based on AMDS gaussian dispersion formula.
  time: time since release (s)
  sigma_w: vertical crosswind standard deviation
  N2: Buoyancy frequency*/
  template<class T>
  T ComputeADMSSigmaZ(T t,
		      T sigma_w,
		      T boundary_height,
		      T friction_velocity,
		      T N2,
		      T source_diameter,
		      T source_height,
		      T plume_rise)
  {
    T b, b2;

    if ((source_height / boundary_height) <= 0.05)
      b = (1. + 0.4 * friction_velocity * t / source_height)
	/ (1. + friction_velocity * t / source_height);
    else if ((source_height / boundary_height) > 0.05
      && (source_height / boundary_height) < 0.15)
      b = (1. - (source_height / boundary_height - 0.05) / 0.1)
	* ((1. + 0.4 * friction_velocity * t / source_height) /
	   (1. + friction_velocity *t / source_height))
	+ (source_height / boundary_height - 0.05) / 0.1;
    else
      b = 1.;

    T N;
    N = sqrt(N2);
    b2 = pow((1. / (b * b) + N2 * t * t / (1. + 2. * N * t)), -1./2.);

    T sigma_zpr, sigma_z1, sigma_z2;
    T sigma_z;
    sigma_zpr = plume_rise * 0.4 / sqrt(2.);
    sigma_z1 = sigma_w * t * b2;
    sigma_z2 = source_diameter * source_diameter / 4.;

    sigma_z = sqrt(sigma_z1 * sigma_z1 +
		   sigma_z2 * sigma_z2 +
		   sigma_zpr * sigma_zpr);

    return sigma_z;
  }

  template<class T>
  T ComputeSigmaWADMS(T z,
		      T boundary_height,
		      T friction_velocity,
		      T lmo,
		      T convective_velocity,
		      T zic,
		      T zi)
  {
    T z0;
    T sigma_wc2, sigma_wm2;
    T sigma_w2, sigma_w;
    // T Twc, Twn;

    if (zic > 0.)
      {
	if (z <= 0.1 * zic)
	  sigma_wc2 = 1.6 * pow((z/zic), 2./3.) * pow(convective_velocity, 2.);
	else if (z > (0.1 * zic) && z <= zic)
	  sigma_wc2 = 0.35 * pow(convective_velocity, 2.);
	else
	  sigma_wc2 = 0.35 * pow(convective_velocity, 2.) *
	    exp(-6. * (z - zic) / zic);	
      }
    else
      sigma_wc2 = 0.;

    if (z < zi)
      sigma_wm2 = 1.3 * friction_velocity * pow((1 - z / zi), 1./2.);
    else
      sigma_wm2 = 0.;

    sigma_w2 = sigma_wc2 + sigma_wm2;
    sigma_w = pow(sigma_w2, 1./2.);

    return sigma_w;
  }

  template<class T>
  T ComputeSigmaVADMS(T z,
		      T boundary_height,
		      T friction_velocity,
		      T lmo,
		      T convective_velocity,
		      T zim)
  {
    T z0;
    T sigma_vc2, sigma_vm2;
    T sigma_v02, sigma_vm2zim;
    T sigma_v2, sigma_v;

    
    if (zim > 0.)
      {
	if (z <= zim)
	  sigma_vc2 = 0.35 * pow(convective_velocity, 2.);
	else
	  sigma_vc2 = 0.5;
      }
    else
      sigma_vc2 = 0.;

    sigma_v02 = 3.6 * pow(friction_velocity, 2.);
    sigma_vm2zim = min(sigma_v02, 0.25);

    if (z <= zim)
      sigma_vm2 = ((sigma_vm2zim - sigma_v02) / zim) * z + sigma_v02;
    else
      sigma_vm2 = sigma_vm2zim;

    sigma_v2 = sigma_vc2 + sigma_vm2;
    sigma_v = pow(sigma_v2, 1./2.);
    return sigma_v;

  }

 /* Computes time corresponding to a griven sigmay 
 Formula used in stable conditions */
  template<class T>
  T ComputeADMSStableTime(T sigma_y,
			  T friction_velocity,
			  T convective_velocity,
			  T plume_rise,
			  T Tlv_stable,
			  T source_diameter,
			  T sigma_v)
  {
    T sigma_yw, sigma_ypr;
    T sigma_yi;

    sigma_ypr = plume_rise * 0.4 / pow(2., 1./2.);
    sigma_yi = source_diameter * source_diameter / 4.;
    
    T time_tmp, sigma_tmp;
    time_tmp = 0.;
    sigma_tmp = 0.;

    while (sigma_tmp < (0.9995 * sigma_y))
      {
	sigma_yw = sigma_v * time_tmp 
	  / pow((1 + time_tmp / (2 * Tlv_stable)), 0.3);
		
	sigma_tmp = pow((sigma_yw * sigma_yw +
		       sigma_ypr * sigma_ypr +
		       sigma_yi * sigma_yi), 1./2.);
	time_tmp += 1.;
      }
    time_tmp -= 1.;
    return time_tmp;    
  }

 /* Computes time corresponding to a griven sigmay 
 Formula used in convective conditions */
  template<class T>
  T ComputeADMSConvectiveTime(T sigma_y,
			      T boundary_height,
			      T friction_velocity,
			      T convective_velocity,
			      T u10,
			      T u,
			      T z,
			      T plume_rise,
			      T Tlv_convective,
			      T source_diameter)
  {
    T z0, Twn;
    T sigma_vc2, sigma_vn2;
    T sigma_vc, sigma_vn, sigma_v;
    T sigma_yi;
    T sigma_yc, sigma_yn, sigma_yw, sigma_ypr;

    z0 = 0.5;
    Twn = 1. - 0.8 * (z + z0) / boundary_height;

    sigma_vc2 = 0.3 * convective_velocity * convective_velocity;
    sigma_vn2 = 4. * Twn * Twn * friction_velocity * friction_velocity;

    sigma_vc = sqrt(sigma_vc2);
    sigma_vn = sqrt(sigma_vn2);
    sigma_v = pow(sigma_vc2 + sigma_vn2, 1./2.);

    sigma_ypr = plume_rise * 0.4 / sqrt(2.);

    sigma_yi = source_diameter * source_diameter / 4.;

    T time_tmp, sigma_tmp, ti;
    time_tmp = 0.;
    sigma_tmp = 0.;
    ti = 0.;

    while (sigma_tmp < (0.9995 * sigma_y))
      {
	ti = time_tmp * (u / u10);
	sigma_yc = sigma_vc * time_tmp *
	  pow((1 + 0.91 * convective_velocity * time_tmp / boundary_height), -1./2.);
	sigma_yn = sigma_vn * time_tmp *
	  pow((1 + 2.5 * friction_velocity * time_tmp / boundary_height), -1./2.);
	
	sigma_yw = sigma_v * ti / pow((1. + ti / (2. * Tlv_convective)), 0.3);
	
	
	sigma_tmp = sqrt(sigma_yc * sigma_yc +
			 sigma_yn * sigma_yn +
			 sigma_yw * sigma_yw +
			 sigma_ypr * sigma_ypr +
			 sigma_yi * sigma_yi);
	time_tmp += 1.;
      }

    time_tmp -= 1.;

    return time_tmp;
  }


  /* Computes time corresponding to a griven sigmay */
  template<class T>
  T ComputeADMSTime(T sigma_y,
		    T boundary_height,
		    T friction_velocity,
		    T convective_velocity,
		    T u,
		    T z,
		    T lmo,
		    T coriolis,
		    T plume_rise,
		    T source_diameter,
		    int stability)
  {
    T time_adms;

    T u10;
    u10 = friction_velocity / 0.4 * log(10. / 0.5);

    T zi, zim, zic;
    zi = 0.;
    zim = 0.;
    zic = 0.;
    ComputeZi(convective_velocity, friction_velocity, lmo,
	      zi, zim, zic);
    
    // Computes Crosswind horizontal dispersion parameters
    T sigma_v;
    sigma_v = ComputeSigmaVADMS(z, boundary_height, friction_velocity,
				lmo, convective_velocity, zim);

    T sigma_min;
    sigma_min = ComputeSigmaMin(lmo);

    if (sigma_v < sigma_min)
      sigma_v = sigma_min;

    // Compute Lagrangian timestep for horizontal dispersion
    T Tlv_convective, Tlv_stable;
    Tlv_convective = ComputeTlv(zi, sigma_v);
    Tlv_stable = ComputeCrosswindTimeScale(z, friction_velocity,
					   boundary_height, lmo,
					   coriolis, sigma_v, stability);

    // cout << "Compute time adms" << endl;
    // cout << "Stability = " << stability << " sigmav: " << sigma_v << " tlv_conv: " << Tlv_convective << endl;

    if (stability <= 3)
      time_adms = ComputeADMSConvectiveTime(sigma_y, boundary_height, friction_velocity,
					    convective_velocity, u10, u, z, plume_rise, 
					    Tlv_convective, source_diameter);
    else
      time_adms = ComputeADMSStableTime(sigma_y, friction_velocity,
					convective_velocity, plume_rise, 
					Tlv_stable, source_diameter, sigma_v);
    return time_adms;
      
  }

   /* Compute time corresponding to a given sigmaz in convective conditions*/
  template<class T>
  T ComputeADMSConvectiveVerticalTime(T sigma_z,
				      T sigma_w,
				      T boundary_height,
				      T convective_velocity,
				      T source_diameter,
				      T plume_rise,
				      T Tlw)
  {
    T sigma_zpr, sigma_z1, sigma_z2;

    sigma_zpr = plume_rise / 2.;
    sigma_z2 = source_diameter * source_diameter / 4.;

    T time_tmp, sigma_tmp;
    time_tmp = 0.;
    sigma_tmp = 0.;

    while (sigma_tmp < (0.9995 * sigma_z))
      {
	sigma_z1 = sigma_w * time_tmp * 
	  pow((1. + time_tmp / (2. * Tlw)), -1. / 2.);
	
	sigma_tmp = sqrt(sigma_z1 * sigma_z1 +
		       sigma_z2 * sigma_z2 +
		       sigma_zpr * sigma_zpr);
	time_tmp += 1.;	
      }

    time_tmp -= 1.;

    return time_tmp;    
  }

   /* Compute time corresponding to a given sigmaz in stable conditions*/
  template<class T>
  T ComputeADMSStableVerticalTime(T sigma_z,
				  T sigma_w,
				  T boundary_height,
				  T friction_velocity,
				  T N2,
				  T source_diameter,
				  T source_height,
				  T plume_rise)
  {
    T sigma_zpr, sigma_z1, sigma_z2;
    sigma_zpr = plume_rise * 0.4 / sqrt(2.);
    sigma_z2 = source_diameter * source_diameter / 4.;

    T N;
    N = sqrt(N2);

    T b, b2;
    T time_tmp, sigma_tmp;
    time_tmp = 0.;
    sigma_tmp = 0.;

    while (sigma_tmp < (0.9995 * sigma_z))
      {	
	if ((source_height / boundary_height) <= 0.05)
	  b = (1. + 0.4 * friction_velocity * time_tmp / source_height)
	    / (1. + friction_velocity * time_tmp / source_height);
	else if ((source_height / boundary_height) > 0.05
		 && (source_height / boundary_height) < 0.15)
	  b = (1. - (source_height / boundary_height - 0.05) / 0.1)
	    * ((1. + 0.4 * friction_velocity * time_tmp / source_height) /
	       (1. + friction_velocity *time_tmp / source_height))
	    + (source_height / boundary_height - 0.05) / 0.1;
	else
	  b = 1.;

	b2 = pow((1. / (b * b) + N2 * time_tmp * time_tmp
		  / (1. + 2. * N * time_tmp)), -1./2.);

	sigma_z1 = sigma_w * time_tmp * b2;


	sigma_tmp = sqrt(sigma_z1 * sigma_z1 +
		       sigma_z2 * sigma_z2 +
		       sigma_zpr * sigma_zpr);
	time_tmp += 1.;
      }
    time_tmp -= 1.;
    return time_tmp;
  }

  /* Computes time corresponding to a griven sigmaz */
  template<class T>
  T ComputeADMSVerticalTime(T sigma_z,
			    T boundary_height,
			    T friction_velocity,
			    T convective_velocity,
			    T lmo,
			    T source_diameter,
			    T source_height,
			    T plume_rise,
			    T z,
			    int stability)
  {
    // Computes Buoyancy frequency
    T N2;
    N2 = ComputeN2Adms(z, lmo, friction_velocity, boundary_height);

    T zi, zim, zic;
    zi = 0.;
    zim = 0.;
    zic = 0.;
    ComputeZi(convective_velocity, friction_velocity, lmo,
	      zi, zim, zic);

    // Computes Crosswind vertical dispersion parameters
    T sigma_w;
    sigma_w = ComputeSigmaWADMS(z, boundary_height, friction_velocity,
				lmo, convective_velocity, zic, zi);

    T sigma_min;
    sigma_min = ComputeSigmaMin(lmo);
    if (sigma_w < sigma_min)
      sigma_w = sigma_min;

    // Compute Lagrangian timestep for vertical dispersion
    T Tlw;
    Tlw = ComputeTlw(zi, sigma_w, z);
    T time_adms;

    if (stability <= 3)
      time_adms = ComputeADMSConvectiveVerticalTime(sigma_z, sigma_w, boundary_height,
						    convective_velocity, source_diameter,
						    plume_rise, Tlw);
    else
      time_adms = ComputeADMSStableVerticalTime(sigma_z, sigma_w, boundary_height, 
						friction_velocity, N2, source_diameter,
						source_height, plume_rise);
    return time_adms;
      
  }

}

#define POLYPHEMUS_FILE_MODELS_PUFFDISPERSIONFORMULA_HXX
#endif
