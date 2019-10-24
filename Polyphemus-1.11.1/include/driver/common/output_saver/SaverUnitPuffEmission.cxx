// Copyright (C) 2012, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPUFFEMISSION_CXX


#include "SaverUnitPuffEmission.hxx"


namespace Polyphemus
{


  ////////////////////////////
  // SAVERUNITDRYDEPOSITION //
  ////////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitPuffEmission<T, ClassModel>
  ::SaverUnitPuffEmission(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitPuffEmission<T, ClassModel>::~SaverUnitPuffEmission()
  {
  }


  //! Type of saver.
  /*!
    \return The string "puff_emission".
  */
  template<class T, class ClassModel>
  string SaverUnitPuffEmission<T, ClassModel>::GetType()  const
  {
    return "puff_emission";
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
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitPuffEmission<T, ClassModel>
  ::Init(ConfigStream& config_stream, ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Check that memory has been allocated for puff emission.
    if (Model.HasField("PuffEmission"))
      cout << "The model has the field PuffEmission!!" << endl;

    cout << Model.D4("PuffEmission").GetLength(1) << " " << Model.GetNz()
         << " "  <<  Model.D4("PuffEmission").GetLength(2) << " "
         << Model.GetNy() << " " << Model.D4("PuffEmission").GetLength(3)
         << " " << Model.GetNx() << endl;
    if (!Model.HasField("PuffEmission") ||
        //       Model.D4("PuffEmission").GetLength(1) != Model.GetNz_source() ||
        Model.D4("PuffEmission").GetLength(2) != Model.GetNy() ||
        Model.D4("PuffEmission").GetLength(3) != Model.GetNx())
      throw string("The model does not collect the puff emissions.");

    // Species to be saved among dry deposited ones.
    this->Ns = int(this->species_list.size());
    vector<string> model_species_list = Model.GetSpeciesList();
    vector<string> emission_species_list =
      Model.GetSpeciesList("VolumeEmission");
    Ns_emis = emission_species_list.size();
    Nz_emis = Model.D4("PuffEmission").GetLength(1);
    int emission_index;
    for (int s = 0; s < Ns_emis; s++)
      if (find(model_species_list.begin(), model_species_list.end(),
               emission_species_list[s]) == model_species_list.end())
        throw string("Species \"") + emission_species_list[s] + "\" unknown.";
      else
        {
          emission_index = 0;
          while (emission_species_list[s] != this->species_list[emission_index])
            emission_index++;
          this->species_index.push_back(emission_index);
        }

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(Ns_emis);
    for (unsigned int i = 0; i < emission_species_list.size(); i++)
      output_file[i] = find_replace(filename, "&f",
                                    emission_species_list[i]);

    // Empties output files.
    for (unsigned int s = 0; s < emission_species_list.size(); s++)
      ofstream tmp_stream(output_file[s].c_str());

    if (this->averaged)
      {
        int s, k, j, i;
        PuffEmission_.Resize(Ns_emis, Nz_emis, this->base_Ny, this->base_Nx);
        for (s = 0; s < Ns_emis; s++)
          for (k = 0; k < Nz_emis; k++)
            for (j = 0; j < this->base_Ny; j++)
              for (i = 0; i < this->base_Nx; i++)
                PuffEmission_(s, k, j, i) = 0.5
                  * Model.D4("PuffEmission")(this->species_index[s], k, j, i);
      }
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitPuffEmission<T, ClassModel>::InitStep(ClassModel& Model)
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
  void SaverUnitPuffEmission<T, ClassModel>::Save(ClassModel& Model)
  {
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
        {
          int s, k, j, i;
          for (s = 0; s < Ns_emis; s++)
            for (k = 0; k < Nz_emis; k++)
              for (j = 0; j < this->base_Ny; j++)
                for (i = 0; i < this->base_Nx; i++)
                  PuffEmission_(s, k, j, i) += 0.5
                    * Model.D4("PuffEmission")(this->species_index[s], k, j, i);

          PuffEmission_.GetArray() /= T(this->interval_length);

          if (Model.GetCurrentDate() >= this->date_beg
              && Model.GetCurrentDate() <= this->date_end)
            for (s = 0; s < Ns_emis; s++)
              {
                Data<T, 3>
                  PuffEmission_tmp(&PuffEmission_(s, 0, 0, 0),
                                   shape(Nz_emis, this->base_Ny, this->base_Nx));
                FormatBinary<float>().Append(PuffEmission_tmp, output_file[s]);
              }

          for (s = 0; s < Ns_emis; s++)
            for (k = 0; k < Nz_emis; k++)
              for (j = 0; j < this->base_Ny; j++)
                for (i = 0; i < this->base_Nx; i++)
                  PuffEmission_(s, k, j, i) = 0.5
                    * Model.D4("PuffEmission")(this->species_index[s], k, j, i);

          this->counter = 0;
        }
      else
        {
          int s, k, j, i;
          for (s = 0; s < Ns_emis; s++)
            for (k = 0; k < Nz_emis; k++)
              for (j = 0; j < this->base_Ny; j++)
                for (i = 0; i < this->base_Nx; i++)
                  PuffEmission_(s, k, j, i) +=
                    Model.D4("PuffEmission")(this->species_index[s], k, j, i);
        }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (int s = 0; s < Ns_emis; s++)
        for (int k = 0; k < Nz_emis; k++)
          {
            Data<T, 2>
              PuffEmission_tmp(&Model.D4("PuffEmission")
                               (this->species_index[s], k, 0, 0),
                               shape(this->base_Ny, this->base_Nx));
            FormatBinary<float>().Append(PuffEmission_tmp, output_file[s]);
          }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPUFFEMISSION_CXX
#endif
