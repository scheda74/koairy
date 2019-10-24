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

  typedef double value_type;

  if (argc < 2)
    {
      cout << string("Program \"") + argv[0]
           << "\" takes at least one argument.\n" << endl;
      return 1;
    }

  for (int file_index = 1; file_index < argc; file_index++)
    {

      string InputFile = argv[file_index];

      ifstream InputStream;
      InputStream.open(InputFile.c_str(), ifstream::binary);

      // Checks if the input file was opened.
      if (!InputStream.is_open())
        {
          cout << "Unable to open input file \"" + InputFile + "\".\n"
               << endl;
          return 1;
        }

      // Gets file length.
      streampos position;
      position = InputStream.tellg();
      InputStream.seekg(0, ios::end);
      unsigned long file_size = InputStream.tellg() - position;

      // Data.
      char* char_data;
      value_type* data;

      // Allocates memory.
      try
        {
          char_data = new char[file_size];
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

      data = reinterpret_cast<value_type*>(char_data);

      // Go home.
      InputStream.seekg(position);

      // Reads the input file.
      InputStream.read(reinterpret_cast<char*>(data), file_size);

      if (!InputStream.good())
        {
          InputStream.close();
          cout << "Error while reading input file \""
               << InputFile + "\"!" << endl
               << "It may have been removed.\n" << endl;
          return 1;
        }

      InputStream.close();

      // Computations.
      value_type maximum(data[0]);
      value_type minimum(data[0]);
      double mean(data[0]);

      for (unsigned long i = 1;
           i < file_size * sizeof(char) / sizeof(value_type); i++)
        {
          maximum = max(maximum, data[i]);
          minimum = min(minimum, data[i]);
          mean += data[i];
        }

      mean /= value_type(file_size * sizeof(char) / sizeof(value_type));

      if (argc > 2)
        cout << "-- File \"" + InputFile + "\"" << endl;
      cout << "Minimum: " << minimum << endl;
      cout << "Maximum: " << maximum << endl;
      cout << "Mean: " << mean << endl;

      cout << endl;

    }

  return 0;
}
