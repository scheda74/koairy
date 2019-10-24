// Copyright (C) 2016, CEREA - ENPC - EDF R&D
// Author(s): Youngseob Kim
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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

#include "StreetInGridChemistry.cxx"
#include "BaseDriver.cxx"
#include "StreetNetworkChemistry.cxx"
#include "Polair3DChemistrySplit.cxx"
#include "Photochemistry.cxx"
#include "SplitAdvectionDST3.cxx"
#include "DiffusionROS2Split.cxx"
#include "BaseOutputSaver.cxx"

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
  typedef Polair3DChemistry<real, SplitAdvectionDST3<real>,
                            DiffusionROS2<real>,
                            Photochemistry<real> > ClassEulerianModel;
  typedef StreetNetworkChemistry<real, Photochemistry<real> > ClassLocalModel;
  typedef StreetInGridChemistry<real, ClassEulerianModel, ClassLocalModel> ClassModel;

  BaseDriver<real, ClassModel, BaseOutputSaver<double, ClassModel> >
    Driver(argv[1]);

  Driver.Run();

  END;

  return 0;

}
