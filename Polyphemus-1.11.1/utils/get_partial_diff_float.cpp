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
#include<algorithm>
#include<sstream>
#include<cmath>

using namespace std;

//! Converts most types to string.
/*!
  \param input variable to be converted.
  \return A string containing 'input'.
*/
template<typename T>
std::string to_str(const T& input)
{
  std::ostringstream output;
  output << input;
  return output.str();
}

//! Converts string to most types, specially numbers.
/*!
  \param input string to be converted.
  \return 'input' converted to 'T'.
*/
template <class T>
void to_num(std::string s, T& num)
{
  std::istringstream str(s);
  str >> num;
}

//! Converts most types to string of a given size.
/*!
  \param input variable to be converted.
  \param size length of the string.
  \return A string of length 'size' and containing 'input' padded left
  (and filled with spaces).
*/
template<typename T>
std::string to_str(const T& input, int size)
{
  std::ostringstream output;
  output.width(size);
  output.fill(' ');
  output.flags(ios_base::left);
  output << input;
  return output.str();
}


int main(int argc,  char* argv[])
{

  cout << endl;

  typedef float value_type;

  if (argc != 3)
    {
      cout << string("Program \"") + argv[0]
           << "\" requires two arguments:\n"
           << "    [input file #1] [input file #2]\n"
           << "It computes the subtraction: [input file #1] "
           << "- [input file #2].\n" << endl;
      return 1;
    }

  string InputFile0 = argv[1];
  string InputFile1 = argv[2];

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

  // Data.
  char* char_data0(NULL);
  char* char_data1(NULL);
  value_type* data0(NULL);
  value_type* data1(NULL);

  // Checks files length.
  streampos position0, position1;

  position0 = InputStream0.tellg();
  InputStream0.seekg(0, ios::end);
  unsigned long file_size0 = InputStream0.tellg() - position0;

  position1 = InputStream1.tellg();
  InputStream1.seekg(0, ios::end);
  unsigned long file_size1 = InputStream1.tellg() - position1;

  file_size0 = min(file_size0, file_size1);

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

  // Reads the input file.
  InputStream0.read(reinterpret_cast<char*>(data0), file_size0);

  if (InputStream0.bad())
    {
      InputStream0.close();
      InputStream1.close();
      cout << "Error while reading input file \"" + InputFile0 + "\"!" << endl
           << "It may have been removed or it may be too large." << endl;
      return 1;
    }

  InputStream0.close();

  // Go home.
  InputStream1.seekg(position0);

  // Reads the input file.
  InputStream1.read(reinterpret_cast<char*>(data1), file_size0);

  if (InputStream1.bad())
    {
      InputStream0.close();
      InputStream1.close();
      cout << "Error while reading input file \"" + InputFile1 + "\"!" << endl
           << "It may have been removed or it may be too large." << endl;
      return 1;
    }

  InputStream1.close();

  unsigned long length = file_size0 * sizeof(char) / sizeof(value_type);

  // Statistics for file 0.
  value_type maximum0(data0[0]);
  value_type minimum0(data0[0]);
  double mean0(data0[0]);

  for (unsigned long i = 1; i < length; i++)
    {
      maximum0 = max(maximum0, data0[i]);
      minimum0 = min(minimum0, data0[i]);
      mean0 += data0[i];
    }

  mean0 /= double(length);

  double var0(0);

  for (unsigned long i = 0; i < length; i++)
    var0 += (data0[i] - mean0) * (data0[i] - mean0);
  var0 /= double(length - 1);

  // Statistics for file 1.
  value_type maximum1(data1[0]);
  value_type minimum1(data1[0]);
  double mean1(data1[0]);

  for (unsigned long i = 1; i < length; i++)
    {
      maximum1 = max(maximum1, data1[i]);
      minimum1 = min(minimum1, data1[i]);
      mean1 += data1[i];
    }

  mean1 /= double(length);

  double var1(0);

  for (unsigned long i = 0; i < length; i++)
    var1 += (data1[i] - mean1) * (data1[i] - mean1);
  var1 /= double(length - 1);

  // Computes the correlation.
  double corr(0);
  for (unsigned long i = 0; i < length; i++)
    corr += (data0[i] - mean0) * (data1[i] - mean1);
  corr /= double(length - 1);
  corr /= sqrt(var0 * var1);

  // Computes the difference.
  for (unsigned long i = 0; i < length; i++)
    data0[i] -= data1[i];

  // Statistics for the difference.
  value_type maximum(data0[0]);
  value_type minimum(data0[0]);
  double mean(data0[0]);

  for (unsigned long i = 1; i < length; i++)
    {
      maximum = max(maximum, data0[i]);
      minimum = min(minimum, data0[i]);
      mean += data0[i];
    }

  mean /= double(length);

  double var(0);

  for (unsigned long i = 0; i < length; i++)
    var += (data0[i] - mean) * (data0[i] - mean);
  var /= double(length - 1);

  var = sqrt(var);
  var0 = sqrt(var0);
  var1 = sqrt(var1);

  size_t max_length(0);
  max_length = max(max_length, to_str(minimum).size());
  max_length = max(max_length, to_str(minimum0).size());
  max_length = max(max_length, to_str(minimum1).size());
  max_length = max(max_length, to_str(maximum).size());
  max_length = max(max_length, to_str(maximum0).size());
  max_length = max(max_length, to_str(maximum1).size());
  max_length = max(max_length, to_str(mean).size());
  max_length = max(max_length, to_str(mean0).size());
  max_length = max(max_length, to_str(mean1).size());
  max_length = max(max_length, to_str(var).size());
  max_length = max(max_length, to_str(var0).size());
  max_length = max(max_length, to_str(var1).size());
  max_length = max(max_length, to_str("File #0").size());
  max_length = max(max_length, to_str("Difference").size());

  cout << "               \t" << to_str("File #0", max_length)
       << '\t' << to_str("File #1", max_length) << endl;
  cout << "Minima:        \t" << to_str(minimum0, max_length)
       << '\t' << to_str(minimum1, max_length) << endl;
  cout << "Maxima:        \t" << to_str(maximum0, max_length)
       << '\t' << to_str(maximum1, max_length) << endl;
  cout << "Means:         \t" << to_str(mean0, max_length)
       << '\t' << to_str(mean1, max_length) << endl;
  cout << "Standard dev.: \t" << to_str(var0, max_length)
       << '\t' << to_str(var1, max_length) << endl;

  cout << endl;

  cout << "               \t" << to_str("Difference", max_length) << endl;
  cout << "Minimum:       \t" << to_str(minimum, max_length) << endl;
  cout << "Maximum:       \t" << to_str(maximum, max_length) << endl;
  cout << "Mean:          \t" << to_str(mean, max_length) << endl;
  cout << "Standard dev.: \t" << to_str(var, max_length) << endl;

  cout << endl;

  cout << "Correlation between files #0 and #1: " << corr << endl;

  cout << endl;

  return 0;
}
