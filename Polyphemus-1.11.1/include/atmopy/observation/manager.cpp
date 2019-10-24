// Copyright (C) 2011, INRIA
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

#include "../../observation/GroundNetworkObservationManager.cxx"
#include "../../observation/EnsembleObservationManager.cxx"
#include "share/Functions_Vector2.cxx"

namespace Seldon
{
  template class MallocAlloc<int>;
  template class Vector_Base<int, MallocAlloc<int> >;
  template class Vector<int, VectFull, MallocAlloc<int> >;
  template class MallocAlloc<double>;
  template class Vector_Base<double, MallocAlloc<double> >;
  template class Vector<double, VectFull, MallocAlloc<double> >;

  template class MallocObject<Vector<int> >;
  template class Vector2<int>;
  template bool Vector2<int>::HasSameShape<Vector2<int> >(const Vector2<int>& V) const;
  template bool Vector2<int>::HasSameShape<Vector2<double> >(const Vector2<double>& V) const;
  template void Vector2<int>::Flatten<int, MallocAlloc<int> >
  (Vector<int, VectFull, MallocAlloc<int> >& data) const;
  template void Vector2<int>::Flatten<int, MallocAlloc<int> >
  (int beg, int end, Vector<int, VectFull, MallocAlloc<int> >& data) const;

  template class MallocObject<Vector<double> >;
  template class Vector2<double>;
  template bool Vector2<double>::HasSameShape<Vector2<double> >(const Vector2<double>& V) const;
  template bool Vector2<double>::HasSameShape<Vector2<int> >(const Vector2<int>& V) const;
  template void Vector2<double>::Flatten<double, MallocAlloc<double> >
  (Vector<double, VectFull, MallocAlloc<double> >& data) const;
  template void Vector2<double>::Flatten<double, MallocAlloc<double> >
  (int beg, int end, Vector<double, VectFull, MallocAlloc<double> >& data) const;

  template class MallocObject < Vector < Vector<double, Vect_Full, MallocAlloc<double> >,
                                         Vect_Full, MallocObject<Vector<double> > > >;
  template class Vector3<double>;
  template void Vector3<double>::Flatten<double, MallocAlloc<double> >
  (Vector<double, VectFull, MallocAlloc<double> >& data) const;
  template void Vector3<double>::Flatten<double, MallocAlloc<double> >
  (int beg, int end, Vector<double, VectFull, MallocAlloc<double> >& data) const;
  template void Vector3<double>::Flatten<double, MallocAlloc<double> >
  (int beg0, int end0, int beg1, int end1, Vector<double, VectFull, MallocAlloc<double> >& data) const;
}

namespace Verdandi
{
  template void SelectLocation<MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > >, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > > >
  (const Vector2<int, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > > >&,
   Vector2<int, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > > >&);
  template void SelectLocation<MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > >, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > >, double, MallocAlloc<double>, MallocObject<Vector<double, VectFull, MallocAlloc<double> > > >
  (const Vector2<int, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > > >&,
   Vector2<int, MallocAlloc<int>, MallocObject<Vector<int, VectFull, MallocAlloc<int> > > >&,
   Vector2<double, MallocAlloc<double>, MallocObject<Vector<double, VectFull, MallocAlloc<double> > > >&);
}

namespace Polyphemus
{
  template class GroundNetworkObservationManager<double>;
  template void GroundNetworkObservationManager<double>::ApplyOperator(const Seldon::Vector<double>& x, Seldon::Vector<double>& y) const;

  template class EnsembleObservationManager<double>;
}


#define ATMOPY_FILE_MANAGER_CPP
#endif
