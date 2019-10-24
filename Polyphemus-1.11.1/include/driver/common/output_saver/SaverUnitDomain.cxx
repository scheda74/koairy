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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_CXX


#include "SaverUnitDomain.hxx"

#include "AtmoData.hxx"


namespace Polyphemus
{


  /////////////////////
  // SAVERUNITDOMAIN //
  /////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitDomain<T, ClassModel>
  ::SaverUnitDomain(): BaseSaverUnit<T, ClassModel>(),
    Ncall(0)
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitDomain<T, ClassModel>::~SaverUnitDomain()
  {
  }


  //! Type of saver.
  /*!
    \return The string "domain".
  */
  template<class T, class ClassModel>
  string SaverUnitDomain<T, ClassModel>::GetType()  const
  {
    return type;
  }


  //! Group of the saver unit.
  /*!
    \return The group of the saver unit, that is, "forecast" or "ensemble".
  */
  template<class T, class ClassModel>
  string SaverUnitDomain<T, ClassModel>::GetGroup()  const
  {
    if (with_ensemble)
      if (type == "domain_ensemble_forecast")
        return "ensemble_forecast";
      else if (type == "domain_ensemble_analysis")
        return "ensemble_analysis";
      else
        throw string("Error: type of an ensemble saver is not supported");
    else
      return "forecast";
  }


  //! First initialization.
  /*! Reads the configuration.
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex(string)
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNz()
    <li> GetConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitDomain<T, ClassModel>::Init(ConfigStream& config_stream,
                                            ClassModel& Model)
  {
    config_stream.PeekValue("Type", type);
    with_ensemble = type == "domain_ensemble_forecast"
      || type == "domain_ensemble_analysis";

    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Ensemble restriction.
    if (with_ensemble && this->averaged)
      throw string("Error: saving averaged concentrations ")
        + "from an ensemble is not supported.";
    if (with_ensemble && this->initial_concentration)
      throw string("Error: saving initial concentrations ")
        + "from an ensemble is not supported.";

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
    if (!with_ensemble)
      for (unsigned int s = 0; s < this->species_list.size(); s++)
        ofstream tmp_stream(output_file[s].c_str());

    if (this->averaged)
      {
        int s, k, j, i;
        Concentration_.Resize(this->Ns, Nlevels,
                              this->base_Ny, this->base_Nx);
        for (s = 0; s < this->Ns; s++)
          for (k = 0; k < Nlevels; k++)
            for (j = 0; j < this->base_Ny; j++)
              for (i = 0; i < this->base_Nx; i++)
                Concentration_(s, k, j, i) = 0.5
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
  void SaverUnitDomain<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetCurrentDate()
    <li> ComputeConcentration(const vector<int>&, const vector<int>&)
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitDomain<T, ClassModel>::Save(ClassModel& Model)
  {
    if (with_ensemble)
      {
        // If the saver has already been called at the current date, the model
        // number in the ensemble is increased. Otherwise a new time step is
        // saved, starting with the first model.
        if (previous_date == Model.GetCurrentDate())
          Nmodel++;
        else
          {
            Nmodel = 0;
            previous_date = Model.GetCurrentDate();
          }

        // Empties output files.
        if (Nmodel == Ncall)
          {
            for (unsigned int s = 0; s < this->species_list.size(); s++)
              {
                string filename =
                  find_replace(output_file[s], "&m", to_str(Nmodel));
                ofstream tmp_stream(filename.c_str());
              }
            Ncall++;
          }
      }

    if (this->averaged)
      {
        Model.ComputeConcentration(this->species_index, levels);

        if (this->counter % this->interval_length == 0)
          {
            int s, k, j, i;
            for (s = 0; s < this->Ns; s++)
              for (k = 0; k < Nlevels; k++)
                for (j = 0; j < this->base_Ny; j++)
                  for (i = 0; i < this->base_Nx; i++)
                    Concentration_(s, k, j, i) += 0.5
                      * Model.GetConcentration()(this->species_index[s],
                                                 levels[k], j, i);

            Concentration_.GetArray() /= T(this->interval_length);

            if (Model.GetCurrentDate() >= this->date_beg
                && Model.GetCurrentDate() <= this->date_end)
              for (s = 0; s < this->Ns; s++)
                {
                  Data<T, 3>
                    Concentration_tmp(&Concentration_(s, 0, 0, 0),
                                      shape(Nlevels, this->base_Ny,
                                            this->base_Nx));
                  FormatBinary<float>().Append(Concentration_tmp,
                                               output_file[s]);
                }

            for (s = 0; s < this->Ns; s++)
              for (k = 0; k < Nlevels; k++)
                for (j = 0; j < this->base_Ny; j++)
                  for (i = 0; i < this->base_Nx; i++)
                    Concentration_(s, k, j, i) = 0.5
                      * Model.GetConcentration()(this->species_index[s],
                                                 levels[k], j, i);

            this->counter = 0;
          }
        else
          {
            int s, k, j, i;
            for (s = 0; s < this->Ns; s++)
              for (k = 0; k < Nlevels; k++)
                for (j = 0; j < this->base_Ny; j++)
                  for (i = 0; i < this->base_Nx; i++)
                    Concentration_(s, k, j, i) +=
                      Model.GetConcentration()(this->species_index[s],
                                               levels[k], j, i);
          }
      }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      {
        Model.ComputeConcentration(this->species_index, levels);

        // Instantaneous concentrations.
        for (int s = 0; s < this->Ns; s++)
          for (int k = 0; k < Nlevels; k++)
            {
              Data<T, 2>
                Concentration_tmp(&Model.GetConcentration()
                                  (this->species_index[s], levels[k], 0, 0),
                                  shape(this->base_Ny, this->base_Nx));
              if (!with_ensemble)
                FormatBinary<float>().Append(Concentration_tmp,
                                             output_file[s]);
              else
                {
                  string filename =
                    find_replace(output_file[s], "&m", to_str(Nmodel));
                  FormatBinary<float>().Append(Concentration_tmp, filename);
                }
            }
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_CXX
#endif
