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


%module polair3d
%{
  // TODO: getting rid of the following inclusions.
#include "seldon/SeldonHeader.hxx"
#include "verdandi/VerdandiHeader.hxx"
#include "model/QuadraticModel.hxx"
#include "model/ClampedBar.hxx"
#include "model/PythonModel.hxx"
#include "observation_manager/GridToNetworkObservationManager.hxx"
#include "observation_manager/LinearObservationManager.hxx"
#include "observation_manager/PythonObservationManager.hxx"
#include "method/OptimalInterpolation.hxx"
#include "method/ForwardDriver.hxx"
#include "method/ReducedOrderExtendedKalmanFilter.hxx"
#include "share/Functions_Vector2.hxx"

#include "method/NewranPerturbationManager.hxx"

#include "SeldonDataHeader.hxx"
#include "AtmoDataHeader.hxx"

#include "../../models/BaseModel.hxx"
#include "../../models/Polair3DTransport.hxx"
#include "../../models/Polair3DChemistry.hxx"
#include "../../models/Polair3DVerdandi.hxx"

#include "../../modules/transport/SplitAdvectionDST3.hxx"
#include "../../modules/transport/DiffusionROS2.hxx"
#include "../../modules/chemistry/Photochemistry.hxx"

#include "../../driver/common/output_saver/BaseOutputSaver.hxx"
  %}

%include "std_string.i"
%include "std_vector.i"
using namespace std;

namespace std
{
  %template(vector_string) vector<string>;
}

%import "../../verdandi/python/verdandi.i"

%include "../../verdandi/VerdandiHeader.hxx"
%include "../../verdandi/share/VerdandiBase.hxx"
%include "../../verdandi/share/UsefulFunction.hxx"
%include "../../verdandi/share/Functions_Vector2.hxx"

%include "../../Talos/TalosHeader.hxx"

%include "../../verdandi/include/seldon/vector/Vector2.hxx"
%include "../../verdandi/include/seldon/vector/Vector3.hxx"

%include "/usr/include/blitz/array.h"

%include "../../SeldonData/SeldonDataHeader.hxx"
%include "../../AtmoData/AtmoDataHeader.hxx"
%include "../../models/BasePointEmission.hxx"
%include "../../models/BaseModel.hxx"
%include "../../models/Polair3DTransport.hxx"
%include "../../models/Polair3DChemistry.hxx"
%include "../../models/Polair3DVerdandi.hxx"

namespace Polyphemus
{
  %template(BasePointEmissionDouble) BasePointEmission<double>;
  %template(BaseModelDouble) BaseModel<double>;

  %template(Polair3DTransportDouble) Polair3DTransport<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double> >;
  %template(Polair3DChemistryDouble) Polair3DChemistry<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double>, Polyphemus::Photochemistry<double> >;

  %template(Polair3DVerdandiTransport) Polair3DVerdandi<double, Polyphemus::Polair3DTransport<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double> >, Polyphemus::BaseOutputSaver<double, Polyphemus::Polair3DTransport<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double> > > >;
  %template(Polair3DVerdandiChemistry) Polair3DVerdandi<double, Polyphemus::Polair3DChemistry<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double>, Polyphemus::Photochemistry<double> >, Polyphemus::BaseOutputSaver<double, Polyphemus::Polair3DChemistry<double, Polyphemus::SplitAdvectionDST3<double>, Polyphemus::DiffusionROS2<double>, Polyphemus::Photochemistry<double> > > >;
}
