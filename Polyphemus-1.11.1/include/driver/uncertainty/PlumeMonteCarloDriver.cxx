// Copyright (C) 2008, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEMONTECARLODRIVER_CXX


#include "PlumeMonteCarloDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  PlumeMonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::PlumeMonteCarloDriver(string config_file):
    BaseDriver<T, ClassModel, ClassOutputSaver>(config_file),
    config(config_file)
  {
    int i;

    /*** Display options ***/

    config.SetSection("[display]");
    // Should meteorological data be displayed on screen?
    config.PeekValue("Show_meteorological_data",
                     option_display["show_meteo"]);
    // Should iterations be displayed on screen?
    config.PeekValue("Show_iterations", option_display["show_iterations"]);

    /*** Reads meteorological conditions ***/

    config.SetSection("[gaussian]");
    config.GetValue("File_meteo", file_meteo);

    /*** Reads the Monte Carlo parameters ***/

    config.SetSection("[uncertainty]");
    config.GetValue("File_perturbation", file_perturbation);
    config.GetValue("Number_samples", "positive", Nsample);
    config.GetValue("Random_seed", seed);

    ConfigStream perturbation(file_perturbation);
    string line;
    int line_number;
    vector<string> line_element;
    bool parse;

    perturbation.SetSection("[boolean_option]");
    while (!perturbation.IsEmpty())
      {
        line = perturbation.GetLine();
        split(line, line_element, " \n\t:=,;");
        if (line_element.size() != 2 || !is_num(line_element[1]))
          throw "Unable to parse this line:\n" + line
            + "\nin section [boolean_option] of file \""
            + file_perturbation
            + "\".\nThere should be one name followed by one number.";
        boolean_option.push_back(line_element[0]);
        boolean_option_probability.push_back(to_num<T>(line_element[1]));
      }

    perturbation.SetSection("[string_option]");
    line_number = 0;
    while (!perturbation.IsEmpty())
      {
        line = perturbation.GetLine();
        split(line, line_element, " \n\t:=,;()|[]");
        parse = line_element.size() >= 3 && (line_element.size() - 1) % 2 == 0;
        if (parse)
          for (int i = 0; i < int(line_element.size() - 1) / 2; i++)
            parse = parse && is_num(line_element[2 + 2 * i]);
        if (!parse)
          throw "Unable to parse this line:\n\"\"\"\n" + line
            + "\n\"\"\"\nin section [string_option] of file \""
            + file_perturbation
            + "\".\nRefer to Polyphemus user's guide for explanations.";
        string_option.push_back(line_element[0]);
        string_option_value.push_back(vector<string>());
        string_option_probability.push_back(vector<T>());
        for (i = 0; i < int(line_element.size() - 1) / 2; i++)
          {
            string_option_value[line_number]
              .push_back(line_element[1 + 2 * i]);
            if (i == 0)
              string_option_probability[line_number]
                .push_back(to_num<T>(line_element[2 + 2 * i]));
            else
              string_option_probability[line_number]
                .push_back(to_num<T>(line_element[2 + 2 * i])
                           + string_option_probability[line_number][i - 1]);
          }
        i = (line_element.size() - 1) / 2 - 1;
        if (string_option_probability[line_number][i] != 1.)
          throw "The sum of the probabilities for option \"" + line_element[0]
            + "\" is " + to_str(string_option_probability[line_number][i])
            + " instead of 1.";
        line_number++;
      }

    perturbation.SetSection("[numerical_value]");
    while (!perturbation.IsEmpty())
      {
        line = perturbation.GetLine();
        split(line, line_element, " \n\t:=,;");
        if (line_element.size() < 3 || !is_num(line_element[2]))
          throw "Unable to parse this line:\n\"\"\"\n" + line
            + "\n\"\"\"\nin section [numerical_value] of file \""
            + file_perturbation
            + "\".\nAwaited a value name and its PDF description.";

        numerical_value.push_back(line_element[0]);

        if (line_element[1] == "Uniform_relative")
          numerical_value_pdf.push_back(0);
        else if (line_element[1] == "Uniform")
          numerical_value_pdf.push_back(1);
        else if (line_element[1] == "Normal")
          numerical_value_pdf.push_back(2);
        else if (line_element[1] == "Log-normal")
          numerical_value_pdf.push_back(3);
        else
          throw "In section [numerical_value] of file \""
            + file_perturbation + "\", unknown PDF: " + line_element[1];

        vector<T> parameter;
        parameter.push_back(to_num<T>(line_element[2]));

        if (line_element[1] == "Uniform_relative"
            || line_element[1] == "Uniform")
          if (line_element.size() != 4)
            throw "Badly formatted line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [numerical_value] of file \""
              + file_perturbation
              + "\".\nAwaited two parameters for the PDF description.";
          else if (!is_num(line_element[3]))
            throw "Unable to parse this line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [numerical_value] of file \""
              + file_perturbation
              + "\".\nAwaited a value name and its PDF description.";
          else
            parameter.push_back(to_num<T>(line_element[3]));
        else // Normal or Log-normal.
          if (line_element.size() != 3)
            throw "Badly formatted line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [numerical_value] of file \""
              + file_perturbation
              + "\".\nAwaited one parameter for the PDF description.";

        numerical_value_parameter.push_back(parameter);
      }
    numerical_value_reference.resize(numerical_value.size());

    perturbation.SetSection("[source_data]");
    while (!perturbation.IsEmpty())
      {
        line = perturbation.GetLine();
        split(line, line_element, " \n\t:=,;");
        if (line_element.size() < 3 || !is_num(line_element[2]))
          throw "Unable to parse this line:\n\"\"\"\n" + line
            + "\n\"\"\"\nin section [source_data] of file \""
            + file_perturbation
            + "\".\nAwaited a data name and its PDF description.";

        source_data.push_back(line_element[0]);

        if (line_element[1] == "Uniform_relative")
          source_data_pdf.push_back(0);
        else if (line_element[1] == "Uniform")
          source_data_pdf.push_back(1);
        else if (line_element[1] == "Normal")
          source_data_pdf.push_back(2);
        else if (line_element[1] == "Log-normal")
          source_data_pdf.push_back(3);
        else
          throw "In section [source_data] of file \""
            + file_perturbation + "\", unknown PDF: " + line_element[1];

        vector<T> parameter;
        parameter.push_back(to_num<T>(line_element[2]));

        if (line_element[1] == "Uniform_relative"
            || line_element[1] == "Uniform")
          if (line_element.size() != 4)
            throw "Badly formatted line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [numerical_value] of file \""
              + file_perturbation
              + "\".\nAwaited two parameters for the PDF description.";
          else if (!is_num(line_element[3]))
            throw "Unable to parse this line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [source_data] of file \""
              + file_perturbation
              + "\".\nAwaited a data name and its PDF description.";
          else
            parameter.push_back(to_num<T>(line_element[3]));
        else // Normal or Log-normal.
          if (line_element.size() != 3)
            throw "Badly formatted line:\n\"\"\"\n" + line
              + "\n\"\"\"\nin section [numerical_value] of file \""
              + file_perturbation
              + "\".\nAwaited one parameter for the PDF description.";

        source_data_parameter.push_back(parameter);
      }
    source_data_reference.resize(source_data.size());

    /*** Random number generation ***/

    if (is_num(seed))
      {
        double seed_number;
        to_num(seed, seed_number);
        if (seed_number <= 0. || seed_number >= 1.)
          throw "Error: seed number must be in ]0, 1[.";
        urng = new NEWRAN::MotherOfAll(seed_number);
        NEWRAN::Random::Set(*urng);
      }
    else if (seed == "current_time")
      {
        srand(static_cast<unsigned long>(time(0)));
        double seed_number = rand() / double(RAND_MAX);
        urng = new NEWRAN::MotherOfAll(seed_number);
        NEWRAN::Random::Set(*urng);
      }
    else
      {
        NEWRAN::Random::SetDirectory(seed.c_str());
        urng = new NEWRAN::MotherOfAll;
        NEWRAN::Random::Set(*urng);
        NEWRAN::Random::CopySeedFromDisk(true);
      }
  }


  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  PlumeMonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::~PlumeMonteCarloDriver()
  {
    delete urng;
  }


  //! Performs Monte Carlo simulations.
  /*! Initializes the model, the output saver, and then performs the time
    loop.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void PlumeMonteCarloDriver<T, ClassModel, ClassOutputSaver>::Run()
  {
    int i, s, v;
    string line;
    double random_number;
    Array<double, 1> random_number_array;

    /*** Initializations ***/

    this->Model.Init(true);

    this->OutputSaver.Init(this->Model);
    this->OutputSaver.SetGroup("ensemble_forecast");

    /*** Meteorological conditions ***/

    if (option_display["show_meteo"])
      cout << "\tTemperature\tWind angle\tWind velocity\tStability" << endl;

    ConfigStream meteo(file_meteo);
    Nmeteo = 0;
    while (!meteo.IsEmpty())
      {
        line = meteo.GetLine();
        if (split(line)[0] == "[situation]")
          {
            if (option_display["show_iterations"])
              cout << "Case #" << Nmeteo << endl;

            this->Model.InitMeteo(meteo, option_display["show_meteo"]);
            this->OutputSaver.InitStep(this->Model);

            /*** Prepares the Monte Carlo simulations ***/

            // Saves the reference numerical values.
            for (v = 0; v < int(numerical_value.size()); v++)
              numerical_value_reference[v]
                = this->Model.GetModelParameter(numerical_value[v]);
            // Saves the reference sources data.
            for (v = 0; v < int(source_data.size()); v++)
              this->Model.GetSourceData(source_data[v],
                                        source_data_reference[v]);

            /*** Monte Carlo loop ***/

            for (s = 0; s < Nsample; s++)
              {
                for (v = 0; v < int(boolean_option.size()); v++)
                  this->Model.SetOption(boolean_option[v],
                                        uniform.Next() < boolean_option_probability[v]);
                for (v = 0; v < int(string_option.size()); v++)
                  {
                    random_number = uniform.Next();
                    i = 0;
                    while (random_number > string_option_probability[v].at(i))
                      i++;
                    this->Model.SetStringOption(string_option[v],
                                                string_option_value[v].at(i));
                  }
                for (v = 0; v < int(numerical_value.size()); v++)
                  {
                    random_number
                      = RandomNumber(numerical_value_pdf[v],
                                     numerical_value_parameter[v],
                                     numerical_value_reference[v]);
                    this->Model.SetModelParameter(numerical_value[v],
                                                  random_number);
                  }
                for (v = 0; v < int(source_data.size()); v++)
                  {
                    random_number_array
                      .resize(source_data_reference[v].size());
                    for (i = 0; i < (int) source_data_reference[v].size(); i++)
                      random_number_array(i)
                        = RandomNumber(source_data_pdf[v],
                                       source_data_parameter[v],
                                       source_data_reference[v](i));
                    this->Model.SetSourceData(source_data[v],
                                              random_number_array);
                  }

                this->Model.Compute();
                this->OutputSaver.Save(this->Model);
              }

            Nmeteo++;
          }
      }
  }


  //! Returns a random number.
  /*!
    \param pdf the probability density function from which to generate the
    random number; it is encoded in an integer: (1) Uniform_relative, (2)
    Uniform, (3) Normal and (4) Log-normal.
    \param parameter the vector of parameters that define the probability
    density function.
    \param reference the reference value to perturb.
    \return The perturbed value.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  T PlumeMonteCarloDriver<T, ClassModel, ClassOutputSaver>
  ::RandomNumber(int pdf, const vector<T>& parameter, T reference)
  {
    if (pdf == 0) // Uniform_relative.
      {
        double random_number = uniform.Next();
        return reference * (parameter[0] * (1. - random_number)
                            + parameter[1] * random_number);
      }
    else if (pdf == 1) // Uniform.
      {
        double random_number = uniform.Next();
        return reference + (parameter[0] * (1. - random_number)
                            + parameter[1] * random_number);
      }
    else if (pdf == 2) // Normal.
      return reference + parameter[0] / 2. * normal.Next();
    else if (pdf == 3) // Log-normal.
      return reference * pow(sqrt(parameter[0]), normal.Next());
    else
      throw "Unknown PDF index: " + to_str(pdf) + ".";
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEMONTECARLODRIVER_CXX
#endif
