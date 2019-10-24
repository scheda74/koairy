// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_MODELS_NORMS_HXX


namespace Polyphemus
{

  //! Checks convergence with 1-norm.
  template<class T>
  bool CheckConvergenceOne(Data<T, 4> Dp, Data<T, 4> Dc, T epsilon,
                           string method)
  {
    int h, i, j, k;

    int Nt =  Dp.GetLength(0);
    int Nz =  Dp.GetLength(1);
    int Ny =  Dp.GetLength(2);
    int Nx =  Dp.GetLength(3);

    T norm1(0.);
    double mean;
    Dp.Mean(mean);
    T maximum = Dp.GetMax();

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            norm1 += abs(Dc(h, k, j, i) - Dp(h, k, j, i));
    norm1 = norm1 / T(Nt * Nz * Ny * Nx);

    if (method == "mean")
      return (norm1 <= (T(mean) * epsilon));
    else if (method == "max")
      return (norm1 <= (maximum * epsilon));
    else
      throw string("Wrong option.");
  }


  //! Checks convergence with 2-norm.
  template<class T>
  bool CheckConvergenceTwo(Data<T, 4> Dp, Data<T, 4> Dc, T epsilon,
                           string method)
  {
    int h, i, j, k;

    int Nt =  Dp.GetLength(0);
    int Nz =  Dp.GetLength(1);
    int Ny =  Dp.GetLength(2);
    int Nx =  Dp.GetLength(3);

    T norm2(0.);
    double mean;
    Dp.Mean(mean);
    T maximum = Dp.GetMax();

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            norm2 += pow(Dc(h, k, j, i) - Dp(h, k, j, i), 2);
    norm2 = sqrt(norm2) / T(Nt * Nz * Ny * Nx);

    if (method == "mean")
      return (norm2 <= (T(mean) * epsilon));
    else if (method == "max")
      return (norm2 <= (maximum * epsilon));
    else
      throw string("Wrong option.");
  }


  //! Checks convergence with infinity-norm.
  template<class T>
  bool CheckConvergenceInfinity(Data<T, 4> Dp, Data<T, 4> Dc, T epsilon,
                                string method)
  {
    int h, i, j, k;

    int Nt =  Dp.GetLength(0);
    int Nz =  Dp.GetLength(1);
    int Ny =  Dp.GetLength(2);
    int Nx =  Dp.GetLength(3);

    T norminf(0.);
    double mean;
    Dp.Mean(mean);
    T maximum = Dp.GetMax();

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            norminf = max(norminf, abs(Dc(h, k, j, i) - Dp(h, k, j, i)));

    if (method == "mean")
      return (norminf <= (T(mean) * epsilon));
    else if (method == "max")
      return (norminf <= (maximum * epsilon));
    else
      throw string("Wrong option.");
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_NORMS_HXX
#endif
