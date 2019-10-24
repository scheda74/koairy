// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Mohamed Aissaoui, Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_PERTURBATION_RANDGENERATOR_CXX


#include "RandGenerator.hxx"


namespace Polyphemus
{


  //! Main constructor.
  template <class T, class Y>
  RandGenerator<T, Y>::RandGenerator(string seed_directory_):
    seed_directory(seed_directory_), seed_number(0.),
    with_seed_number(is_num(seed_directory_)), urng_(NULL)
  {
    if (with_seed_number)
      {
        to_num(seed_directory_, seed_number);
        if (seed_number <= 0. || seed_number >= 1.)
          throw "Error: seed number must be in ]0, 1[.";
      }

    if (!with_seed_number)
      {
        if (seed_directory_ == "current_time")
          {
            srand((unsigned long)time(0));
            seed_number = rand() / double(RAND_MAX);
            if (seed_number <= 0. || seed_number >= 1.)
              throw "Error: seed number must be in ]0, 1[.";
            urng_ = new T(seed_number);
            NEWRAN::Random::Set(*urng_);
          }
        else
          {
            NEWRAN::Random::SetDirectory(seed_directory.c_str());
            urng_ = new T;
            NEWRAN::Random::Set(*urng_);
            NEWRAN::Random::CopySeedFromDisk();
          }
      }
    else
      {
        urng_ = new T(seed_number);
        NEWRAN::Random::Set(*urng_);
      }

    maximum_spread = 2.;
  }


  //! Destructor.
  template <class T, class Y>
  RandGenerator<T, Y>::~RandGenerator()
  {
    if (urng_ != NULL)
      delete urng_;
  }


  ////////////
  // ACCESS //
  ////////////


  //! Returns maximum spread of the random numbers to be generated.
  /*!
    \return The maximum spread of the random numbers to be generated.
  */
  template <class T, class Y>
  Y RandGenerator<T, Y>::GetMaximumSpread()
  {
    return maximum_spread;
  }


  //! Sets maximum spread.
  /*!
    \param spread the maximum spread to be set.
  */
  template <class T, class Y>
  void RandGenerator<T, Y>::SetMaximumSpread(Y spread)
  {
    maximum_spread = spread;
  }


  //////////////////////////////
  // RANDOM NUMBER GENERATION //
  //////////////////////////////


  //! Initialization.
  template <class T, class Y>
  void RandGenerator<T, Y>::Init(string seed_directory_)
  {
    with_seed_number = is_num(seed_directory_);
    if (urng_ != NULL)
      delete urng_;
    urng_ = NULL;

    if (with_seed_number)
      {
        to_num(seed_directory_, seed_number);
        if (seed_number <= 0. || seed_number >= 1.)
          throw "Error: seed number must be in ]0, 1[.";
      }

    if (!with_seed_number)
      {
        NEWRAN::Random::SetDirectory(seed_directory.c_str());
        urng_ = new T;
        NEWRAN::Random::Set(*urng_);
        NEWRAN::Random::CopySeedFromDisk();
      }
    else
      {
        urng_ = new T(seed_number);
        NEWRAN::Random::Set(*urng_);
      }

    maximum_spread = 2.;
  }


  //! Generates the random number of a given parameter.
  /*!
    Generates the random number of a given parameter, taking into account
    its correlation with other fields.
  */
  template <class T, class Y>
  Y RandGenerator<T, Y>::GenRandomNumber(PolairParam<Y>& param)
  {

    vector<PolairParam<Y>* >& corr_param = param.GetCorrelatedParam();

    NEWRAN::Normal N1;
    NEWRAN::Normal N2;
    Y alpha(0.);

    if (param.GetPdf() == "LN")
      {
        if (corr_param.empty())
          {
            while (abs(double(alpha = N1.Next())) > maximum_spread) ;
            alpha = pow(double(sqrt(param.GetStd())), double(alpha));
            param.SetRandomValue(alpha);
            if (!with_seed_number && seed_directory != "current_time")
              NEWRAN::Random::CopySeedToDisk();
            return param.GetRandomValue();
          }

        if (!corr_param[0]->IsSet())
          {
            while (abs(double(alpha = N1.Next())) > maximum_spread) ;
            alpha = pow(double(sqrt(corr_param[0]->GetStd())), double(alpha));
            corr_param[0]->SetRandomValue(alpha);
          }
        else
          alpha = corr_param[0]->GetRandomValue();

        param.SetRandomValue(alpha);
      }
    else if (param.GetPdf() == "N")
      {
        if (corr_param.empty())
          {
            while (abs(double(alpha = N1.Next())) > maximum_spread) ;
            param.SetRandomValue(param.GetStd() * alpha);
            if (!with_seed_number && seed_directory != "current_time")
              NEWRAN::Random::CopySeedToDisk();
            return param.GetRandomValue();
          }

        if (!corr_param[0]->IsSet())
          {
            while (abs(double(alpha = N1.Next())) > maximum_spread) ;
            corr_param[0]->SetRandomValue(corr_param[0]->GetStd() * alpha);
          }
        else
          alpha = corr_param[0]->GetRandomValue();

        Y alpha_bis(0.), beta(0.);
        while (abs(double(alpha_bis = N2.Next())) > maximum_spread) ;
        beta = param.GetStd() * param.GetCorrelation()[0] * alpha
          / corr_param[0]->GetStd() + param.GetStd() *
          sqrt(1 - param.GetCorrelation()[0] * param.GetCorrelation()[0])
          * alpha_bis;

        param.SetRandomValue(beta);
      }

    if (!with_seed_number && seed_directory != "current_time")
      NEWRAN::Random::CopySeedToDisk();

    return param.GetRandomValue();
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_RANDGENERATOR_CXX
#endif
