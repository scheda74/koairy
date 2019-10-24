// Copyright (C) 2012, INRIA
// Author(s): Vivien Mallet
//
// This file is part of AtmoPy library, a tool for data processing and
// visualization in atmospheric sciences.
//
// AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in the
// ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoPy is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// AtmoPy is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoPy home page:
//     http://cerea.enpc.fr/polyphemus/atmopy.html


#ifndef ATMOPY_FILE_MANAGER_CPP

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "Verdandi.hxx"
#include "method/NewranPerturbationManager.cxx"

#include "AtmoDataHeader.hxx"
#include "SeldonDataHeader.hxx"

#include "../../models/Polair3DVerdandi.cxx"

#include "../../models/BaseModel.cxx"
#include "../../models/Polair3DChemistry.cxx"

#include "../../modules/transport/SplitAdvectionDST3.cxx"
#include "../../modules/transport/DiffusionROS2.cxx"
#include "../../modules/chemistry/Photochemistry.cxx"

#include "../../driver/common/output_saver/BaseOutputSaver.cxx"


namespace Seldon
{
  template class MallocAlloc<int>;
  template class Vector_Base<int, MallocAlloc<int> >;
  template class Vector<int, VectFull, MallocAlloc<int> >;
  template class Vector < Vector<int, VectFull, MallocAlloc<int> >, Collection,
                          MallocAlloc<Vector<int, VectFull, MallocAlloc<int> > > >;
}


namespace Polyphemus
{
  template class DiffusionROS2<double>;
  template class BaseModel<double>;
  template class Polair3DTransport<double, SplitAdvectionDST3<double>, DiffusionROS2<double> >;
  template class Polair3DChemistry<double, SplitAdvectionDST3<double>, DiffusionROS2<double>, Photochemistry<double> >;
  template class Polair3DVerdandi<double, Polair3DTransport<double, SplitAdvectionDST3<double>, DiffusionROS2<double> >, BaseOutputSaver<double, Polair3DTransport<double, SplitAdvectionDST3<double>, DiffusionROS2<double> > > >;
  template class Polair3DVerdandi<double, Polair3DChemistry<double, SplitAdvectionDST3<double>, DiffusionROS2<double>, Photochemistry<double> >, BaseOutputSaver<double, Polair3DChemistry<double, SplitAdvectionDST3<double>, DiffusionROS2<double>, Photochemistry<double> > > >;
}


#define ATMOPY_FILE_MANAGER_CPP
#endif
