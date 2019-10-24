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

  if (argc != 4)
    {
      cout << string("Program \"") + argv[0]
           << "\" requires three arguments:\n"
           << "[input file] [number] [output file]\n"
           << "It multiplies the single-precision binary file [input file] by"
           << " a given number, and it saves it in [output file].\n" << endl;
      return 1;
    }

  string InputFile = argv[1];
  value_type number;
  to_num<value_type>(argv[2], number);
  string OutputFile = argv[3];

  ifstream InputStream;
  InputStream.open(InputFile.c_str(), ifstream::binary);

  // Checks if the input file was opened.
  if (!InputStream.is_open())
    {
      cout << "Unable to open input file \"" + InputFile + "\"." << endl;
      return 1;
    }

  ofstream OutputStream;
  OutputStream.open(OutputFile.c_str(), ofstream::binary);

  // Checks if the input file was opened.
  if (!OutputStream.is_open())
    {
      cout << "Unable to open output file \"" + OutputFile + "\"." << endl;
      return 1;
    }

  // Data.
  char* char_data(NULL);
  value_type* data;

  // Checks files length.
  streampos position;

  position = InputStream.tellg();
  InputStream.seekg(0, ios::end);
  unsigned long file_size = InputStream.tellg() - position;

  cout << "Files size: " << file_size << " bytes." << endl;

  // Allocates memory.
  try
    {
      char_data = new char[file_size];
    }
  catch (bad_alloc)
    {
      cout << "Out of memory." << endl;
    }
  catch (...)
    {
      cout << "Error in memory allocation." << endl;
    }

  data = reinterpret_cast<value_type*>(char_data);

  // Go home.
  InputStream.seekg(position);

  cout << "Reading file \"" << InputFile << "\"...";
  cout.flush();
  // Reads the input file.
  InputStream.read(reinterpret_cast<char*>(data), file_size);

  if (InputStream.bad())
    {
      InputStream.close();
      OutputStream.close();
      cout << "Error while reading input file \"" + InputFile + "\"!" << endl
           << "It may have been removed or it may be too large." << endl;
      return 1;
    }

  InputStream.close();
  cout << " done." << endl;

  // Transformation.
  for (unsigned long i = 0; i < file_size * sizeof(char) / sizeof(value_type); i++)
    data[i] *= number;

  // Writes data.
  cout << "Writing in file \"" << OutputFile << "\"...";
  cout.flush();
  OutputStream.write(reinterpret_cast<char*>(data), file_size);

  if (OutputStream.bad())
    {
      OutputStream.close();
      cout << "Error while writing in output file \"" + OutputFile + "\"!\n"
           << "Maybe there is no space left on device." << endl;
      return 1;
    }

  OutputStream.close();
  cout << " done." << endl;

  cout << endl;

  return 0;
}
