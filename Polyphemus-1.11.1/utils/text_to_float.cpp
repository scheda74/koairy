// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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


#include<iostream>
#include<fstream>
#include<sstream>

using namespace std;

//! Converts string to most types, specially numbers.
/*!
  \param s string to be converted.
  \param num 'input' converted to 'T'.
*/
template <class T>
void to_num(std::string s, T& num)
{
  std::istringstream str(s);
  str >> num;
}

int main(int argc,  char* argv[])
{

  cout << endl;

  typedef float value_type;

  if (argc < 3)
    {
      cout << string("Usage: \"") + argv[0]
           << " [input text file] [output binary file]\".\n"
           << "It converts a text file to a single-precision binary file.\n"
           << endl;
      return 1;
    }

  string InputFile = argv[1];
  string OutputFile = argv[2];

  ifstream InputStream;
  InputStream.open(InputFile.c_str());

  // Checks if the input file was opened.
  if (!InputStream.is_open())
    {
      cout << "Unable to open input file \"" + InputFile + "\".\n" << endl;
      return 1;
    }

  ofstream OutputStream;
  OutputStream.open(OutputFile.c_str(), ofstream::binary);

  // Checks if the output file was opened.
  if (!OutputStream.is_open())
    {
      cout << "Unable to open output file \"" + OutputFile + "\"." << endl;
      return 1;
    }

  // Data.
  value_type* data;

  value_type value;

  long int nb_values(0), i(0);

  while (InputStream >> value)
    nb_values++;
  InputStream.close();
  InputStream.clear();

  data = new value_type[nb_values];

  InputStream.open(InputFile.c_str());
  while (InputStream >> data[i++]);
  InputStream.close();

  OutputStream.write(reinterpret_cast<char*>(data),
                     nb_values * sizeof(value_type));
  OutputStream.close();

  cout << endl;

  return 0;
}
