// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Hervé Njomgang, Vivien Mallet
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


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////

  if (argc != 3)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [file] [field]\n\n";
      mesg += "Arguments:\n";
      mesg += "  [file]: path to the MM5 file.\n";
      mesg += "  [field]: field to be extracted.\n";
      cout << mesg << endl;
      return 1;
    }

  string filename = argv[1];
  if (!exists(filename))
    throw string("Unable to open file \"") + filename + "\".";

  // Field to be extracted from the MM5 file.
  string fieldname = argv[2];


  //////////////////////
  // PARSING MM5 FILE //
  //////////////////////

  ifstream file(filename.c_str());

  FormatMM5 MM5;

  Array<int, 2> BHI;
  Array<float, 2> BHR;
  Array<string, 2> BHIC;
  Array<string, 2> BHRC;

  MM5SubHeader SH;

  MM5.ReadFlag(file);  // Skips the first flag.
  MM5.ReadBigHeader(file, BHI, BHR, BHIC, BHRC);

  // Searches for the variable 'fieldname'.
  string name;
  while (MM5.ReadFlag(file) == 1 && name != fieldname)
    {
      MM5.ReadSubHeader(file, SH);
      MM5.ReadField(file);  // Skips the field.
      name = trim(SH.name);
    }

  if (name != fieldname)
    throw string("Unable to find the field \"") + fieldname
      + string("\" in file \"") + filename + "\".";

  // Number of time steps.
  int Nt;
  if (BHR(11, 4) < BHR(11, 3))
    Nt = int(BHR(11, 0) / BHR(11, 3)) + 1;
  else
    Nt = int(BHR(11, 4) / BHR(11, 3));

  // Other dimensions are retrieved from SH.
  // SH is currently the sub-header associated with 'fieldname'.
  int offset = 0;
  if (trim(SH.staggering) == "D")
    offset = 1;
  int Ny = max(SH.end_index(0) - SH.start_index(0) + offset, 1);
  int Nx = max(SH.end_index(1) - SH.start_index(1) + offset, 1);
  int Nz = max(SH.end_index(2) - SH.start_index(2) + 1, 1);


  //////////////////
  // READING DATA //
  //////////////////

  Data<float, 4> Field(Nt, Nz, Nx, Ny);
  MM5.ReadWholeField(filename, name, Field);

  double mean(0), std(0);
  Field.Mean(mean);
  Field.StandardDeviation(std);

  cout << "Min: " << Field.GetMin() << endl;
  cout << "Max: " << Field.GetMax() << endl;
  cout << "Mean: " << mean << endl;
  cout << "Std. dev.: " << std << endl;
  cout << endl;
  END;

  return 0;
}
