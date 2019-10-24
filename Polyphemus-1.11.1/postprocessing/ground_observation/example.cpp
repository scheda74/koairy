// Copyright (C) 2011, INRIA
// Author(s): Vivien Mallet
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


//////////////
// INCLUDES //

#define VERDANDI_WITH_ABORT
#define SELDON_DEBUG_LEVEL_2

#include "observation/GroundNetworkObservationManager.cxx"

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  Polyphemus::GroundNetworkObservationManager<double> om("observation.lua");

  Seldon::Vector<double> observation;
  om.SetDate("2001-06-15 09:00");
  om.GetObservation(observation);
  observation.Print();

  END;

  return 0;

}
