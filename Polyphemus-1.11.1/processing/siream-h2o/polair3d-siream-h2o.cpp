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

#include "BaseDriver.cxx"
#include "Polair3DAerosol.cxx"
#include "SplitAdvectionDST3.cxx"
#include "DiffusionROS2.cxx"
#include "Aerosol_SIREAM_H2O.cxx"
#include "BaseOutputSaver.cxx"

using namespace Polyphemus;

// INCLUDES //
//////////////


int main(int argc, char** argv)
{
  char s[256];
  gethostname(s, sizeof(s));
  ofstream fout("machine.pid");
  fout << s << "." << getpid() << endl;
  fout.close();

  TRY;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  typedef double real;
  typedef Polair3DAerosol < real, SplitAdvectionDST3<real>,
                            DiffusionROS2<real>,
                            Aerosol_SIREAM_AEC<real> > ClassModel;

  BaseDriver<real, ClassModel, BaseOutputSaver<real, ClassModel> >
    Driver(argv[1]);

  Driver.Run();

  END;

  return 0;

}
