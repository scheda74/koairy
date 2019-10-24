// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPUFFEMISSION_AER_CXX


#include "SaverUnitPuffEmission_aer.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // SAVERUNITDRYDEPOSITION_AER //
  ////////////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitPuffEmission_aer<T, ClassModel>
  ::SaverUnitPuffEmission_aer(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitPuffEmission_aer<T, ClassModel>::~SaverUnitPuffEmission_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "puff_emission_aer".
  */
  template<class T, class ClassModel>
  string SaverUnitPuffEmission_aer<T, ClassModel>::GetType()  const
  {
    return "puff_emission_aer";
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
    if (Model.HasField("PuffEmission_aer"))
      cout << "The model has the field PuffEmission!!" << endl;

    cout << Model.D5("PuffEmission_aer").GetLength(2) << " " << Model.GetNz()
         << " "  <<  Model.D5("PuffEmission_aer").GetLength(3) << " "
         << Model.GetNy() << " " << Model.D5("PuffEmission_aer").GetLength(4)
         << " " << Model.GetNx() << endl;
    if (!Model.HasField("PuffEmission_aer") ||
        //       Model.D4("PuffEmission").GetLength(1) != Model.GetNz_source() ||
        Model.D5("PuffEmission_aer").GetLength(3) != Model.GetNy() ||
        Model.D5("PuffEmission_aer").GetLength(4) != Model.GetNx())
      throw string("The model does not collect the puff aerosol emissions.");

    // // Species to be saved among dry deposited ones.
    // int Ns_aer = int(this->species_list_aer.size());
    // vector<string> model_species_list_aer = Model.GetSpeciesList_aer();
    // vector<string> emission_species_list_aer =
    //   Model.GetSpeciesList("VolumeEmission_aer");
    // Ns_aer_emis = emission_species_list_aer.size();
    // Nz_emis = Model.D5("PuffEmission_aer").GetLength(2);
    // int emission_index_aer;
    // cout << "Ns_aer_emis: " << Ns_aer_emis << endl;
    // for (int s = 0; s < Ns_aer_emis; s++)
    //   if (find(model_species_list_aer.begin(), model_species_list_aer.end(),
    //            emission_species_list_aer[s]) == model_species_list_aer.end())
    //     throw string("Species \"") + emission_species_list_aer[s] + "\" unknown.";
    //   else
    //     {
    //       emission_index_aer = 0;
    //       while (emission_species_list_aer[s] != this->species_list_aer[emission_index_aer])
    //         emission_index_aer++;
    //       cout << "emission_index_aer: " << emission_index_aer << endl;
    //       this->species_index_aer.push_back(emission_index);
    //     }

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(Ns_emis);
    for (unsigned int i = 0; i < emission_species_list.size(); i++)
      {
        output_file[i] = find_replace(filename, "&f",
                                      emission_species_list[i]);
        cout << "emission_species_list[i]: " << emission_species_list[i] << endl;
      }
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
    cout << "End of Init of SaverUnit() " << endl;
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
