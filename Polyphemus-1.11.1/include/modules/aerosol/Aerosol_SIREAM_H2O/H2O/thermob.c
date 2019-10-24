/* Copyright (C) 2013, ENPC
 *    Author(s): Florian Couvidat
 *
 * This file is part of the air quality modeling system Polyphemus.
 *
 * Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
 * the ENPC - EDF R&D joint laboratory CEREA.
 *
 * Polyphemus is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * Polyphemus is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * For more information, visit the Polyphemus web site:
 *      http://cerea.enpc.fr/polyphemus/
 */


#include <stdio.h>
#include <stdlib.h>

#include "glodef.h"


extern int thermoflag;
extern void bunidriver(double xpass[], double gamma[], int size);


/**
 * \brief Calls the thermodynamic model.
 *
 * The thermodynamic model is chosen with the global variable thermoflag.
 */
void ThermoB(double xpass[], double gamma[], int size)
{
  static const double gamma_inf[NBSP + NBSPAOM] =
    {
      1.0, 1.0, 1.0, 1.0, 1.0, 5.66, 2.29, 1.16, 1.16, 5.99, 3.81
    };
  static const double gamma_f[NBSP + NBSPAOM] =
    {
      1.0, 1.0, 1.0, 1.0, 1.0, 5.66, 2.29, 1.16, 1.16, 5.99, 3.81
    };

  const double tiny = 1e-6;
  int i;

  /* Gamma is constant. */
  if (thermoflag == 0)
    {
      for (i = 0; i < size; i++)
        if (xpass[i] <= tiny) /* Infinite dilution. */
          gamma[i] = gamma_inf[i];
        else
          gamma[i] = gamma_f[i];
    }
  /* Computes Gamma with Unifac model. */
  else if (thermoflag == 1)
    {
      bunidriver(xpass, gamma, size);
    }
  /* The system is ideal. */
  else if (thermoflag == 2)
    {
      for (i = 0; i < size; i++)
        gamma[i] = 1.0;
    }
}
