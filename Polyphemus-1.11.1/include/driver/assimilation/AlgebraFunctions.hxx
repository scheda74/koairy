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


#ifndef POLYPHEMUS_FILE_DRIVER_ALGEBRAFUNCTIONS_HXX


#include "AtmoData.hxx"


namespace Polyphemus
{


  using namespace AtmoData;


  //////////////////////
  // ALGEBRAFUNCTIONS //
  //////////////////////

  template<class T>
  int TruncateMatrixByColumn(Array<T, 2>& Matrix, int num_truncate,
                             Array<T, 2>& TruncateMatrix);

  template<class T>
  void ComputeSqrtMatrix(Array<T, 2>& Matrix, int num_mode,
                         Array<T, 2>& SqrtMatrix >);


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ALGEBRAFUNCTIONS_HXX
#endif
