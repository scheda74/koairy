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

  if (argc != 3)
    {
      cout << string("Program \"") + argv[0]
           << "\" requires two arguments:\n"
           << "[input file] [output file]\n"
           << "It switches a file from "
           << "big-endian to little-endian, or vice versa.\n" << endl;
      return 1;
    }

  string InputFile = argv[1];
  string OutputFile = argv[2];

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

  // Checks if the output file was opened.
  if (!OutputStream.is_open())
    {
      cout << "Unable to open output file \"" + OutputFile + "\"." << endl;
      return 1;
    }

  // Data.
  char* data;

  // Checks file length.
  streampos position;
  position = InputStream.tellg();
  InputStream.seekg(0, ios::end);
  unsigned long file_size = InputStream.tellg() - position;

  cout << "File size: " << file_size << " bytes." << endl;

  if (file_size % 4 != 0)
    {
      InputStream.close();
      OutputStream.close();
      cout << "The length of \"" << InputFile
           << "\" is not a multiple of 4! Cannot proceed." << endl;
      return 1;
    }

  // Allocates memory.
  try
    {
      data = new char[file_size];
    }
  catch (bad_alloc)
    {
      cout << "Out of memory." << endl;
      return 1;
    }
  catch (...)
    {
      cout << "Error in memory allocation." << endl;
      return 1;
    }

  // Go home.
  InputStream.seekg(position);

  cout << "Reading file \"" << InputFile << "\"...";
  cout.flush();
  // Reads the input file.
  InputStream.read(data, file_size);

  if (InputStream.bad())
    {
      InputStream.close();
      OutputStream.close();
      cout << "Error while reading input file \"" + InputFile + "\"!" << endl
           << "It may have been removed or it may be too big." << endl;
      return 1;
    }

  InputStream.close();
  cout << " done." << endl;

  // Transformation.
  char a, b;
  for (unsigned long i = 0; i < file_size / 4; i++)
    {
      a = data[4 * i];
      b = data[4 * i + 1];
      data[4 * i] = data[4 * i + 3];
      data[4 * i + 1] = data[4 * i + 2];
      data[4 * i + 2] = b;
      data[4 * i + 3] = a;
    }

  // Writes data.
  cout << "Writing in file \"" << OutputFile << "\"...";
  cout.flush();
  OutputStream.write(data, file_size);

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
