// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Nikki Vercauteren, Meryem Ahmed de Biasi.
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_AER_CXX


#include "SaverUnitSubdomain_aer.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  /////////////////////////////
  // SAVERUNITSUBDOMAIN_AER //
  ///////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitSubdomain_aer<T, ClassModel>
  ::SaverUnitSubdomain_aer(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitSubdomain_aer<T, ClassModel>::~SaverUnitSubdomain_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "subdomain_aer".
  */
  template<class T, class ClassModel>
  string SaverUnitSubdomain_aer<T, ClassModel>::GetType()  const
  {
    return "subdomain_aer";
  }


  //! First initialization.
  /*!
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex_aer(string)
    <li> GetNs_aer()
    <li> GetNbin_aer()
    <li> GetSpeciesList_aer()
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain_aer<T, ClassModel>::Init(ConfigStream&
                                                   config_stream,
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

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species and bins.
    string field, species, file;
    vector<string> vsplit, bounds;
    int first, last, i, j, k;
    if (this->species_list[0] == "all")
      {
        int Ns_aer = Model.GetNs_aer();
        vector<string> list_aer = Model.GetSpeciesList_aer();
        vector<int> bins;
        pair<string, vector<int> > tmp;
        for (i = 0; i < Model.GetNbin_aer(); i++)
          bins.push_back(i);
        tmp.second = bins;
        for (i = 0; i < Ns_aer; i++)
          {
            tmp.first = list_aer[i];
            species_list_aer.push_back(tmp);
            // Adds a buffer to compute averaged concentrations.
            Concentration_.push_back(vector<Data<T, 3> >());
            output_file.push_back(vector<string>());
            for (k = 0; k < int(bins.size()); k++)
              {
                file = find_replace(filename, "&f", tmp.first);
                file = find_replace(file, "&n", to_str(bins[k]));
                output_file[i].push_back(file);
                Concentration_[i].push_back(Data<T, 3>());
              }
          }
      }
    else
      for (unsigned int i = 0; i < this->species_list.size(); i++)
        {
          field = this->species_list[i];
          vsplit = split(field, "{}");
          if (field[0] == '{')
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: it cannot be parsed by the output saver.";
          if (vsplit.size() == 1)
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: bins are needed by the output saver.";
          if (vsplit.size() > 2)
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: it cannot be parsed by the output saver.";

          // Searches for species index in 'species_list_aer' and
          // 'output_file'.  If the species is not found, it is added in
          // 'species_list_aer' and 'output_file'.
          species = split(field, "_")[0];
          j = 0;
          while (j < int(species_list_aer.size())
                 && species_list_aer[j].first != species)
            j++;
          if (j == int(species_list_aer.size()))
            {
              species_list_aer.push_back(pair<string, vector<int> >
                                         (species, vector<int>()));
              output_file.push_back(vector<string>());
              Concentration_.push_back(vector<Data<T, 3> >());
            }
          bounds = split(vsplit[1], "-");
          // First bound.
          first = convert<int>(bounds[0]);
          // Last bound.
          if (bounds.size() != 1)
            last = convert<int>(bounds[1]);
          else
            last = first;
          for (k = first; k < last + 1; k++)
            {
              // Adds the bin (associated to species #j).
              species_list_aer[j].second.push_back(k);
              // Adds a buffer to compute averaged concentrations.
              Concentration_[j].push_back(Data<T, 3>());

              // Field base name is replaced in the generic file name.
              file = find_replace(filename, "&f",
                                  species);

              // Field number is also replaced.
              file = find_replace(file, "&n", to_str(k));
              // Adds the corresponding output file.
              output_file[j].push_back(file);
            }
        }

    // Indices in the model.
    int base_s, base_b;

    int s, b;
    // Empties output files.
    for (s = 0; s < int(this->species_list_aer.size()); s++)
      for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
        ofstream tmp_stream(output_file[s][b].c_str());

    if (this->averaged)
      {
        int sizei = i_max - i_min + 1;
        int sizej = j_max - j_min + 1;
        for (s = 0; s < int(this->species_list_aer.size()); s++)
          for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
            {
              base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
              base_b = species_list_aer[s].second[b];
              Concentration_[s][b].Resize(Nlevels, sizej, sizei);
              for (k = 0; k < Nlevels; k++)
                for (j = j_min; j < j_max + 1; j++)
                  for (i = i_min; i < i_max + 1; i++)
                    Concentration_[s][b](k, j - j_min, i - i_min) = 0.5
                      * Model.GetConcentration_aer()(base_s, base_b,
                                                     levels[k], j, i);
            }
      }

    if (this->initial_concentration && !this->averaged)
      this->Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetSpeciesIndex_aer()
    <li> GetCurrentDate()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitSubdomain_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int s, b, k, j, i;
    int base_s, base_b;
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
        for (s = 0; s < int(this->species_list_aer.size()); s++)
          for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
            {
              base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
              base_b = species_list_aer[s].second[b];
              for (k = 0; k < Nlevels; k++)
                for (j = j_min; j < j_max + 1; j++)
                  for (i = i_min; i < i_max + 1; i++)
                    Concentration_[s][b](k, j - j_min, i - i_min) += 0.5
                      * Model.GetConcentration_aer()(base_s, base_b,
                                                     levels[k], j, i);

              Concentration_[s][b].GetArray() /= T(this->interval_length);

              if (Model.GetCurrentDate() >= this->date_beg
                  && Model.GetCurrentDate() <= this->date_end)
                FormatBinary<float>().Append(Concentration_[s][b],
                                             output_file[s][b]);

              for (k = 0; k < Nlevels; k++)
                for (j = j_min; j < j_max + 1; j++)
                  for (i = i_min; i < i_max + 1; i++)
                    Concentration_[s][b](k, j - j_min, i - i_min) = 0.5
                      * Model.GetConcentration_aer()(base_s, base_b,
                                                     levels[k], j, i);

              this->counter = 0;
            }
      else
        for (s = 0; s < int(this->species_list_aer.size()); s++)
          for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
            {
              base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
              base_b = species_list_aer[s].second[b];
              for (k = 0; k < Nlevels; k++)
                for (j = j_min; j < j_max + 1; j++)
                  for (i = i_min; i < i_max + 1; i++)
                    Concentration_[s][b](k, j - j_min, i - i_min) +=
                      Model.GetConcentration_aer()(base_s, base_b,
                                                   levels[k], j, i);
            }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (s = 0; s < int(this->species_list_aer.size()); s++)
        for (b = 0; b < int(this->species_list_aer[s].second.size()); b++)
          {
            base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
            base_b = species_list_aer[s].second[b];
            int sizei = i_max - i_min + 1;
            int sizej = j_max - j_min + 1;
            for (int k = 0; k < Nlevels; k++)
              {

                Data<T, 2>
                  Concentration_tmp2(&Model.GetConcentration_aer()
                                     (base_s, base_b, levels[k], 0, 0),
                                     shape(this->base_Ny, this->base_Nx));
                Data<T, 2>
                  Concentration_tmp(shape(sizej, sizei));

                Concentration_tmp.SubData(Concentration_tmp2,
                                          Range(j_min, j_max),
                                          Range(i_min, i_max));

                FormatBinary<float>().Append(Concentration_tmp,
                                             output_file[s][b]);
              }
          }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_AER_CXX
#endif
