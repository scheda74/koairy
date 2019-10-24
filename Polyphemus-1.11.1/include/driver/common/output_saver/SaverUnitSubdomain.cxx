// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Nikki Vercauteren.
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_CXX


#include "SaverUnitSubdomain.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  ////////////////////////
  // SAVERUNITSUBDOMAIN //
  ////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitSubdomain<T, ClassModel>
  ::SaverUnitSubdomain(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitSubdomain<T, ClassModel>::~SaverUnitSubdomain()
  {
  }


  //! Type of saver.
  /*!
    \return The string "subdomain".
  */
  template<class T, class ClassModel>
  string SaverUnitSubdomain<T, ClassModel>::GetType()  const
  {
    return "subdomain";
  }


  //! First initialization.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex(string)
    <li> GetConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain<T, ClassModel>::Init(ConfigStream& config_stream,
                                               ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Horizontal domain to be saved.
    config_stream.PeekValue("i_min", i_min);
    config_stream.PeekValue("i_max", i_max);
    config_stream.PeekValue("j_min", j_min);
    config_stream.PeekValue("j_max", j_max);

    // Vertical levels to be saved.
    config_stream.Find("Levels");
    split(config_stream.GetLine(), levels);
    Nlevels = int(levels.size());

    // Species to be saved.
    this->Ns = int(this->species_list.size());
    for (int s = 0; s < this->Ns; s++)
      this->species_index.push_back(Model.
                                    GetSpeciesIndex(this->species_list[s]));

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(this->Ns);
    for (unsigned int i = 0; i < this->species_list.size(); i++)
      output_file[i] = find_replace(filename, "&f", this->species_list[i]);

    // Empties output files.
    for (unsigned int s = 0; s < this->species_list.size(); s++)
      ofstream tmp_stream(output_file[s].c_str());

    if (this->averaged)
      {
        int s, k, j, i;
        int sizei = i_max - i_min + 1;
        int sizej = j_max - j_min + 1;
        Concentration_.Resize(this->Ns, Nlevels, sizej, sizei);
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < Nlevels; k++)
            for (j = j_min; j < j_max + 1; j++)
              for (i = i_min; i < i_max + 1; i++)
                Concentration_(s, k, j - j_min, i - i_min) = 0.5
                  * Model.GetConcentration()(this->species_index[s],
                                             levels[k], j, i);
      }

    if (this->initial_concentration && !this->averaged)
      this->Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetCurrentDate()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain<T, ClassModel>::Save(ClassModel& Model)
  {
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
        {
          int s, k, j, i;
          for (s = 0; s < this->Ns; s++)
            for (k = 0; k < Nlevels; k++)
              for (j = j_min; j < j_max + 1; j++)
                for (i = i_min; i < i_max + 1; i++)
                  Concentration_(s, k, j - j_min, i - i_min) += 0.5
                    * Model.GetConcentration()(this->species_index[s],
                                               levels[k], j, i);

          Concentration_.GetArray() /= T(this->interval_length);

          if (Model.GetCurrentDate() >= this->date_beg
              && Model.GetCurrentDate() <= this->date_end)
            for (s = 0; s < this->Ns; s++)
              {
                Data<T, 3>
                  Concentration_tmp(&Concentration_(s, 0, 0, 0),
                                    shape(Nlevels, j_max - j_min + 1,
                                          i_max - i_min + 1));
                FormatBinary<float>().Append(Concentration_tmp,
                                             output_file[s]);
              }

          for (s = 0; s < this->Ns; s++)
            for (k = 0; k < Nlevels; k++)
              for (j = j_min; j < j_max + 1; j++)
                for (i = i_min; i < i_max + 1; i++)
                  Concentration_(s, k, j - j_min, i - i_min) = 0.5
                    * Model.GetConcentration()(this->species_index[s],
                                               levels[k], j, i);

          this->counter = 0;
        }
      else
        {
          int s, k, j, i;
          for (s = 0; s < this->Ns; s++)
            for (k = 0; k < Nlevels; k++)
              for (j = j_min; j < j_max + 1; j++)
                for (i = i_min; i < i_max + 1; i++)
                  Concentration_(s, k, j - j_min, i - i_min) +=
                    Model.GetConcentration()(this->species_index[s],
                                             levels[k], j, i);
        }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (int s = 0; s < this->Ns; s++)
        for (int k = 0; k < Nlevels; k++)
          {
            int sizei = i_max - i_min + 1;
            int sizej = j_max - j_min + 1;

            Data<T, 2>
              Concentration_tmp2(&Model.GetConcentration()
                                 (this->species_index[s], levels[k], 0, 0),
                                 shape(this->base_Ny, this->base_Nx));
            Data<T, 2>
              Concentration_tmp(shape(sizej, sizei));

            Concentration_tmp.SubData(Concentration_tmp2, Range(j_min, j_max),
                                      Range(i_min, i_max));


            FormatBinary<float>().Append(Concentration_tmp, output_file[s]);
          }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_CXX
#endif
