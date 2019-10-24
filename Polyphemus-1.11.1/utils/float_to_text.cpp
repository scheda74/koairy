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

  typedef float value_type_in;

  if (argc < 3)
    {
      cout << string("Usage: \"") + argv[0]
           << " [input binary file] [output text file]\".\n"
           << "It converts a binary file to a text file.\n" << endl;
      return 1;
    }

  string InputFile = argv[1];
  string OutputFile = argv[2];

  ifstream InputStream;
  InputStream.open(InputFile.c_str(), ifstream::binary);

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
  value_type_in* data_in;

  // Gets file length.
  streampos position;
  position = InputStream.tellg();
  InputStream.seekg(0, ios::end);
  unsigned long file_size = InputStream.tellg() - position;

  cout << "\"" << InputFile << "\" size: " << file_size << " bytes." << endl;

  // Allocates memory.
  try
    {
      data_in = new value_type_in[file_size / sizeof(value_type_in)];
    }
  catch (bad_alloc)
    {
      cout << "Out of memory.\n" << endl;
      return 1;
    }
  catch (...)
    {
      cout << "Error in memory allocation.\n" << endl;
      return 1;
    }

  // Go home.
  InputStream.seekg(position);

  // Reads data.
  InputStream.read(reinterpret_cast<char*>(data_in), file_size);
  InputStream.close();

  // Writes data.
  OutputStream.precision(8);
  for (unsigned long i = 0; i < file_size / sizeof(value_type_in) - 1; i++)
    OutputStream << scientific << data_in[i] << "\t";
  OutputStream << scientific
               << data_in[file_size / sizeof(value_type_in) - 1];

  if (!OutputStream.good())
    {
      cout << "Unable to write in file \"" + OutputFile + "\"." << endl;
      return 1;
    }

  OutputStream.close();

  cout << "Number of values in \"" << OutputFile << "\": "
       << file_size / sizeof(value_type_in) << " values." << endl;

  cout << endl;

  return 0;
}
