// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Hervé Njomgang
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


#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;


int main(int argc, char** argv)
{
  TRY;

  cout << endl;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [file]\n\n";
      mesg += "Argument:\n";
      mesg += "  [file]: path of the file (absolute).\n";
      cout << mesg << endl;
      return 1;
    }

  string filename = argv[1];
  if (!exists(filename))
    throw string("Unable to find the file ") + filename;

  ifstream file(filename.c_str());
  FormatMM5 MM5;
  MM5SubHeader SH;


  // Big-header (general metadata).
  MM5.ReadFlag(file);

  Array<int, 2> BHI;
  Array<float, 2> BHR;
  Array<string, 2> BHIC;
  Array<string, 2> BHRC;
  MM5.ReadBigHeader(file, BHI, BHR, BHIC, BHRC);

  cout << "Metadata (-999 means unknown):\n" << endl;

  // Domain.
  for (int j = 0; j < 24; j++)
    cout << BHIC(0, j) << ": " << BHI(0, j) << endl;
  for (int j = 0; j < 15; j++)
    cout << BHRC(0, j) << ": " << BHR(0, j) << endl;

  cout << endl;

  // Simulation.
  for (int j = 1; j < 12; j++)
    cout << BHIC(10, j) << ": " << BHI(10, j) << endl;

  for (int j = 0; j < 4; j++)
    if (j != 2)
      cout << BHRC(11, j) << ": " << BHR(11, j) << endl;


  // Sub-headers (outputs metadata).
  cout << "\n\nOutputs:\n" << endl;
  cout << "Name      Dim. 1   2   3   4   Stag. Ord. Units"
       << "                    Description" << endl;

  int steps = 0;

  while (!is_empty(file))
    {
      while (MM5.ReadFlag(file) == 1)
        {
          MM5.ReadSubHeader(file, SH);

          if (steps == 0)
            {
              cout << fill(trim(SH.name), 10, ' ', ios_base::left);
              cout << to_str_fill(SH.ndim, 5, ' ', ios_base::left);

              for (int i = 0 ; i < SH.ndim ; i++)
                cout << to_str_fill(SH.end_index(i) - SH.start_index(i) + 1,
                                    4, ' ', ios_base::left);

              for (int i = SH.ndim ; i < 4 ; i++)
                cout << "    ";

              cout << fill(trim(SH.staggering), 6, ' ', ios_base::left);
              cout << fill(trim(SH.ordering), 5, ' ', ios_base::left);
              cout << fill(trim(SH.unit), 25, ' ', ios_base::left);
              cout << trim(SH.description) << endl;
            }

          MM5.ReadField(file);
        }

      steps++;
    }

  cout << endl;
  cout << "Total number of time steps read in the file: " << steps << endl;

  file.close();


  END;

  return 0;
}
