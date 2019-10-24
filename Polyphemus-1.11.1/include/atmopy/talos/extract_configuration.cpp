// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
//     Author(s): Vivien Mallet
//
// This file is part of AtmoPy library, a tool for data processing and
// visualization in atmospheric sciences.
//
// AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in the
// ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoPy is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// AtmoPy is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoPy home page:
//     http://cerea.enpc.fr/polyphemus/atmopy.html


// This program extracts data from a configuration file.
// It is based on Talos.


//////////////
// INCLUDES //

#include <iostream>
#include <list>
#include <vector>
using namespace std;

#undef TALOS_DEBUG
#include "Talos.hxx"
using namespace Talos;

// INCLUDES //
//////////////


int main(int argc, char** argv)
{

  TRY;

  if (argc == 1)
    {
      string mesg  = "\nUsage:\n";
      mesg += string("  ") + argv[0] + " [configuration file] {{options} or {elements}}\n\n";
      mesg += "Arguments:\n";
      mesg += "  [configuration file]: configuration file.\n";
      mesg += "  {options}:\n";
      mesg += "    -ls (default): lists all sections.\n";
      mesg += "    -ll: lists all lines (discarding empty lines and comments).\n";
      mesg += "    -t: test the configuration file for existence.\n";
      mesg += "    -s [section] {elements}: if {elements} is empty, it lists all lines\n";
      mesg += "      (discarding empty lines and comments) in the section; if {elements}\n";
      mesg += "      is not empty, it lists the first values of the elements\n";
      mesg += "      (in the section or after this section).\n";
      mesg += "  {elements}: elements whose values are to be returned.\n";
      mesg += "              No option should be provided.\n\n";
      mesg += "Examples:\n";
      mesg += "  \"extract_configuration input.cfg -s [Section]\"\n";
      mesg += "      Returns all lines in section [Section].\n";
      mesg += "  \"extract_configuration input.cfg -s [Section] length date\"\n";
      mesg += "      Returns the values of \"length\" and \"date\" as set in section [Section].\n";
      mesg += "  \"extract_configuration input.cfg -s length date\"\n";
      mesg += "      Returns the values of \"length\" and \"date\" as first set in the file.\n";
      cout << mesg << endl;
      return 1;
    }

  string configuration_file = argv[1];
  if (!exists(configuration_file))
    throw string("Unable to open file \"") + configuration_file + "\".";

  ConfigStream configuration(configuration_file, "#");

  string option(""), value, line;
  vector<list<string> > keys;
  vector<string> sections;
  vector<string> extract_section;

  // With no option: lists all sections.
  if (argc == 2)
    {
      configuration.Rewind();
      while (configuration.GetLine(line))
        if (line[0] == '[')
          {
            extract_section = split(line, "[]");
            if (!extract_section.empty())
              cout << "[" << split(line, "[]")[0] << "]" << endl;
          }
    }

  int i = 2;
  while (i < argc)
    {
      // Is there an option?
      if (option == "" && argv[i][0] == '-')
        option = &argv[i][1];

      // Stores the option tag and skips it.
      if (option != "")
        {
          ++i;
          if (i == argc && option != "ls" && option != "ll")
            throw string("Option -") + option + " should be followed by a value.";
          else if (option != "ls" && option != "ll")
            value = argv[i];
        }

      // Lists all lines.
      if (option == "ll")
        {
          configuration.Rewind();
          while (configuration.GetLine(line))
            cout << line << endl;
          option = "";
          --i;
        }
      // Lists all sections.
      else if (option == "ls")
        {
          configuration.Rewind();
          while (configuration.GetLine(line))
            if (line[0] == '[')
              {
                extract_section = split(line, "[]");
                if (!extract_section.empty())
                  cout << "[" << split(line, "[]")[0] << "]" << endl;
              }
          option = "";
          --i;
        }
      // Section tag.
      else if (option == "s")
        {
          sections.push_back(value);
          keys.push_back(list<string>());
          option = "";
        }
      // Just a test for configuration existence.
      else if (option == "t")
        {
          return 0;
        }
      // No option: just a key.
      else if (option == "")
        {
          value = argv[i];
          // Not in a given section yet.
          if (sections.empty())
            {
              sections.push_back("");
              keys.push_back(list<string>(1, value));
            }
          else // In a given section.
            keys[keys.size() - 1].push_back(value);
        }
      else
        throw string("Option -") + option + " unrecognized.";
      ++i;
    }

  configuration.Rewind();
  for (i = 0; i < int(sections.size()); i++)
    {
      if (sections[i] != "")
        configuration.SetSection(sections[i]);
      if (!keys[i].empty())
        {
          list<string>::iterator j;
          for (j = keys[i].begin(); j != keys[i].end(); j++)
            cout << configuration.PeekValue(*j) << endl;
        }
      else // Prints lines of the current section.
        while (configuration.PeekLine(line) && line[0] != '[')
          {
            cout << line << endl;
            configuration.GetLine(line);
          }
    }

  END;

  return 0;
}
