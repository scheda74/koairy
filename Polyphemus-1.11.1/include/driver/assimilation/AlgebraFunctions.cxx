// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_DRIVER_ALGEBRAFUNCTIONS_CXX


/////////////////////////////
// BLAS & LAPACK INTERFACE //
/////////////////////////////


extern "C"
{
#include "cblas.h"
#include "clapack.h"
}


namespace Polyphemus
{


  //////////////////////
  // ALGEBRAFUNCTIONS //
  //////////////////////


  //! Truncates columns of a matrix.
  /*!
    \param Matrix the matrix to be truncated.
    \param Nreduce the expected number of columns for the truncation.
    \param TruncateMatrix (output) the truncated matrix.
    \return The number of columns in \a TruncateMatrix.
  */

  template<class T>
  int TruncateMatrixByColumn(Array<T, 2>& Matrix, int Nreduce,
                             Array<T, 2>& TruncateMatrix)
  {
    // Reports an error if Matrix or TruncateMatrix is not column-majored.
    if (Matrix.ordering(0) != 0 || TruncateMatrix.ordering(0) != 0)
      throw Error("TruncateMatrixByColumn",
                  "A matrix is not column-majored.");

    int Nstate, Nmode;
    Nstate = Matrix.rows();
    Nmode = Matrix.columns();

    Array<T, 2> MTM(Nmode, Nmode, ColumnMajorArray<2>());
    double one = 1.;
    double zero = 0;
    char Transpose = 'T';
    char NoTranspose = 'N';

    dgemm_(&Transpose, &NoTranspose, &Nmode, &Nmode, &Nstate, &one,
           Matrix.data(), &Nstate, Matrix.data(), &Nstate, &zero,
           MTM.data(), &Nmode);

    Array<T, 1> SingularValue(Nmode);
    Array<T, 2> V(Nmode, Nmode);
    Array<T, 2> VT(Nmode, Nmode);
    Array<T, 1> work(1);

    int lwork = -1;
    char JobV = 'A';
    char JobVT = 'A';
    int info;

    dgesvd_(&JobV, &JobVT, &Nmode, &Nmode, MTM.data(), &Nmode,
            SingularValue.data(), V.data(), &Nmode, VT.data(), &Nmode,
            work.data(),  &lwork, &info);

    if (info != 0)
      throw Error("TruncateMatrixByColumn",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESVD.");

    lwork = int(ceil(work(0)));
    work.resize(lwork);
    dgesvd_(&JobV, &JobVT, &Nmode, &Nmode, MTM.data(), &Nmode,
            SingularValue.data(), V.data(), &Nmode, VT.data(), &Nmode,
            work.data(), &lwork, &info);

    if (info != 0)
      throw Error("TruncateMatrixByColumn",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESVD.");

    int Ntruncate = 0;
    for (int i = 0; i < Nmode; i++)
      if (SingularValue(i) > 1.e-10)
        Ntruncate++;
    if (Ntruncate >= Nreduce)
      Ntruncate = Nreduce;

    TruncateMatrix.resize(Nstate, Ntruncate);

    dgemm_(&NoTranspose, &NoTranspose, &Nstate, &Ntruncate, &Nmode, &one,
           Matrix.data(), &Nstate, V.data(), &Nmode, &zero,
           TruncateMatrix.data(), &Nstate);

    return Ntruncate;
  }


  //! Computes the square root of a square matrix.
  /*! Only the first \a num_mode leading columns of the square root matrix are
    selected. The selection is performed by SVD.
    \param Matrix the square matrix.
    \param num_mode the number of modes (that is, columns) expected in the
    resulting square root.
    \param SqrtMatrix (output) the first \a num_mode leading columns of the
    square root of \a Matrix.
  */
  template<class T>
  void ComputeSqrtMatrix(Array<T, 2>& Matrix, int num_mode,
                         Array<T, 2>& SqrtMatrix)
  {
    // Reports an error if Matrix or SqrtMatrix is not column-majored.
    if (Matrix.ordering(0) != 0 || SqrtMatrix.ordering(0) != 0)
      throw Error("ComputeSqrtMatrix",
                  "A matrix is not column-majored.");

    if (Matrix.rows() != Matrix.columns())
      throw Error("ComputeSqrtMatrix",
                  "Input should be square matrix.");

    int num_row = Matrix.rows();

    SqrtMatrix.resize(num_row, num_mode);

    Array<T, 2> working_matrix(Matrix.rows(), Matrix.columns(),
                               ColumnMajorArray<2>());
    working_matrix = Matrix;

    Array<T, 1> SingularValue(num_row);
    Array<T, 2> U(num_row, num_row, ColumnMajorArray<2>());
    Array<T, 2> VT(num_row, num_row, ColumnMajorArray<2>());
    Array<T, 1> work(1);

    int lwork = -1;
    char JobV = 'A';
    char JobVT = 'A';
    int info;

    dgesvd_(&JobV, &JobVT, &num_row, &num_row, working_matrix.data(),
            &num_row, SingularValue.data(), U.data(), &num_row, VT.data(),
            &num_row, work.data(),  &lwork, &info);
    if (info != 0)
      throw Error("ComputeSqrtMatrix",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESVD.");

    lwork = int(ceil(work(0)));
    work.resize(lwork);
    dgesvd_(&JobV, &JobVT, &num_row, &num_row, working_matrix.data(),
            &num_row, SingularValue.data(), U.data(), &num_row, VT.data(),
            &num_row, work.data(),  &lwork, &info);
    if (info != 0)
      throw Error("ComputeSqrtMatrix",
                  string("Lapack has returned code error #") + to_str(info)
                  + " from DGESVD.");

    for (int j = 0; j < num_mode; j++)
      {
        if (SingularValue(j) < 0.)
          throw Error("ComputeSqrtMatrix",
                      "Matrix singular value less than 0.");
        for (int i = 0; i < num_row; i++)
          SqrtMatrix(i, j) = U(i, j) * sqrt(SingularValue(j));
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ALGEBRAFUNCTIONS_CXX
#endif
