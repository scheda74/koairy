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
           << "[input file #1] [input file #2] [output file]\n"
           << "It computes [input file #1] - [input file #2]"
           << " and saves it into [output file].\n" << endl;
      return 1;
    }

  string InputFile0 = argv[1];
  string InputFile1 = argv[2];
  string OutputFile = argv[3];

  ifstream InputStream0;
  InputStream0.open(InputFile0.c_str(), ifstream::binary);

  // Checks if the input file was opened.
  if (!InputStream0.is_open())
    {
      cout << "Unable to open input file \"" + InputFile0 + "\"." << endl;
      return 1;
    }

  ifstream InputStream1;
  InputStream1.open(InputFile1.c_str(), ifstream::binary);

  // Checks if the input file was opened.
  if (!InputStream1.is_open())
    {
      cout << "Unable to open input file \"" + InputFile1 + "\"." << endl;
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
  char* char_data0(NULL);
  char* char_data1(NULL);
  value_type* data0;
  value_type* data1;

  // Checks files length.
  streampos position0, position1;

  position0 = InputStream0.tellg();
  InputStream0.seekg(0, ios::end);
  unsigned long file_size0 = InputStream0.tellg() - position0;

  position1 = InputStream1.tellg();
  InputStream1.seekg(0, ios::end);
  unsigned long file_size1 = InputStream1.tellg() - position1;

  if (file_size0 != file_size1)
    {
      cout << "\"" + InputFile0 + "\" size: " << file_size0
           << " bytes." << endl;
      cout << "\"" + InputFile1 + "\" size: " << file_size1
           << " bytes." << endl;
      cout << "Sizes differ. Cannot proceed." << endl;
      return 1;
    }

  cout << "Files size: " << file_size0 << " bytes." << endl;

  // Allocates memory.
  try
    {
      char_data0 = new char[file_size0];
      char_data1 = new char[file_size0];
    }
  catch (bad_alloc)
    {
      cout << "Out of memory." << endl;
    }
  catch (...)
    {
      cout << "Error in memory allocation." << endl;
      return 1;
    }

  data0 = reinterpret_cast<value_type*>(char_data0);
  data1 = reinterpret_cast<value_type*>(char_data1);

  // Go home.
  InputStream0.seekg(position0);

  cout << "Reading file \"" << InputFile0 << "\"...";
  cout.flush();
  // Reads the input file.
  InputStream0.read(reinterpret_cast<char*>(data0), file_size0);

  if (InputStream0.bad())
    {
      InputStream0.close();
      InputStream1.close();
      OutputStream.close();
      cout << "Error while reading input file \"" + InputFile0 + "\"!" << endl
           << "It may have been removed or it may be too large." << endl;
      return 1;
    }

  InputStream0.close();
  cout << " done." << endl;

  // Go home.
  InputStream1.seekg(position0);

  cout << "Reading file \"" << InputFile1 << "\"...";
  cout.flush();
  // Reads the input file.
  InputStream1.read(reinterpret_cast<char*>(data1), file_size1);

  if (InputStream1.bad())
    {
      InputStream0.close();
      InputStream1.close();
      OutputStream.close();
      cout << "Error while reading input file \"" + InputFile1 + "\"!" << endl
           << "It may have been removed or it may be too large." << endl;
      return 1;
    }

  InputStream1.close();
  cout << " done." << endl;

  // Transformation.
  for (unsigned long i = 0; i < file_size0 * sizeof(char) / sizeof(value_type); i++)
    data0[i] -= data1[i];

  // Writes data.
  cout << "Writing in file \"" << OutputFile << "\"...";
  cout.flush();
  OutputStream.write(reinterpret_cast<char*>(data0), file_size0);

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
