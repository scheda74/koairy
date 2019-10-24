// Copyright (C) 2007, ENPC - INRIA - EDF R&D
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

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"

#include "Gradient4DVarDriver.cxx"
#include "Polair3DChemistryAssimConc.cxx"
#include "AdvectionDST3.cxx"
#include "DiffusionROS2.cxx"
#include "Photochemistry.cxx"
#include "BaseOutputSaver.cxx"
#include "GroundObservationManager.cxx"

using namespace Polyphemus;

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  typedef double real;
  typedef Polair3DChemistryAssimConc < real, AdvectionDST3<real>,
                                       DiffusionROS2<real>,
                                       Photochemistry<real> > ClassModel;

  Gradient4DVarDriver < real, ClassModel,
                        BaseOutputSaver<real, ClassModel>,
                        GroundObservationManager<real> >
    Driver(argv[1]);

  Driver.Run();

  END;

  return 0;

}
