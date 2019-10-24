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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEMONTECARLODRIVER_HXX


#include "BaseDriver.cxx"
#include "newran.h"


namespace Polyphemus
{


  //////////////////////
  // MONTECARLODRIVER //
  //////////////////////


  /*! \brief This driver performs Monte Carlo simulations.
   */
  template<class T, class ClassModel, class ClassOutputSaver>
  class PlumeMonteCarloDriver:
    public BaseDriver<T, ClassModel, ClassOutputSaver>
  {

  protected:

    /*** Configuration ***/

    //! Configuration stream.
    ConfigStream config;

    //! Display options.
    map<string, bool> option_display;

    //! Meteorological conditions.
    string file_meteo;
    //! Number of meteorological conditions.
    int Nmeteo;

    /*** Perturbations ***/

    //! Perturbation to apply.
    string file_perturbation;
    //! Number of Monte Carlo samples.
    int Nsample;
    //! Random seed (directory, number or "current_time").
    string seed;

    //! The base uniform random number generator.
    NEWRAN::MotherOfAll* urng;
    //! Uniform random generator.
    NEWRAN::Uniform uniform;
    //! Normal random generator.
    NEWRAN::Normal normal;

    //! List of Boolean options.
    vector<string> boolean_option;
    //! Probability to set a Boolean option.
    vector<T> boolean_option_probability;

    //! List of string options.
    vector<string> string_option;
    //! List of possible values for the string options.
    vector<vector<string> > string_option_value;
    //! Probability to set a string option to a given value.
    vector<vector<T> > string_option_probability;

    //! List of uncertain numerical values.
    vector<string> numerical_value;
    //! PDFs of the uncertain numerical values.
    vector<int> numerical_value_pdf;
    //! PDFs parameter for the uncertain numerical values.
    vector<vector<T> > numerical_value_parameter;
    //! Reference numerical values (without perturbation).
    vector<T> numerical_value_reference;

    //! List of uncertain sources data.
    vector<string> source_data;
    //! PDFs of the uncertain sources data.
    vector<int> source_data_pdf;
    //! PDFs parameter for the uncertain sources data.
    vector<vector<T> > source_data_parameter;
    //! Reference sources data (without perturbation).
    vector<Array<T, 1> > source_data_reference;

  public:

    /*** Constructor and destructor ***/

    PlumeMonteCarloDriver(string config_file);
    ~PlumeMonteCarloDriver();

    /*** Methods ***/

    void Run();

  protected:

    T RandomNumber(int pdf, const vector<T>& parameter, T reference);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEMONTECARLODRIVER_HXX
#endif
