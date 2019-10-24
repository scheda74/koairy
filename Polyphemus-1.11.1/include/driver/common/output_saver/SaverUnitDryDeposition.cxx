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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDRYDEPOSITION_CXX


#include "SaverUnitDryDeposition.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  ////////////////////////////
  // SAVERUNITDRYDEPOSITION //
  ////////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitDryDeposition<T, ClassModel>
  ::SaverUnitDryDeposition(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitDryDeposition<T, ClassModel>::~SaverUnitDryDeposition()
  {
  }


  //! Type of saver.
  /*!
    \return The string "dry_deposition".
  */
  template<class T, class ClassModel>
  string SaverUnitDryDeposition<T, ClassModel>::GetType()  const
  {
    return "dry_deposition";
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
  void SaverUnitDryDeposition<T, ClassModel>
  ::Init(ConfigStream& config_stream, ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Check that memory has been allocated for dry deposition flux.
    if (!Model.HasField("DryDepositionFlux") ||
        Model.D3("DryDepositionFlux").GetLength(1) != Model.GetNy() ||
        Model.D3("DryDepositionFlux").GetLength(2) != Model.GetNx())
      throw string("The model does not collect the dry deposition fluxes.");

    // Species to be saved among dry deposited ones.
    this->Ns = int(this->species_list.size());
    vector<string> model_species_list = Model.GetSpeciesList();
    vector<string> deposition_species_list =
      Model.GetSpeciesList("DepositionVelocity");
    int deposition_index;
    for (int s = 0; s < this->Ns; s++)
      if (find(model_species_list.begin(), model_species_list.end(),
               this->species_list[s]) == model_species_list.end())
        throw string("Species \"") + this->species_list[s] + "\" unknown.";
      else if (find(deposition_species_list.begin(),
                    deposition_species_list.end(), this->species_list[s])
               == deposition_species_list.end())
        throw string("Species \"") + this->species_list[s]
          + "\" is not dry deposited.";
      else
        {
          deposition_index = 0;
          while (deposition_species_list[deposition_index] != this->species_list[s])
            deposition_index++;
          this->species_index.push_back(deposition_index);
        }

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
        int s, j, i;
        DryDepositionFlux_.Resize(this->Ns, this->base_Ny, this->base_Nx);
        for (s = 0; s < this->Ns; s++)
          for (j = 0; j < this->base_Ny; j++)
            for (i = 0; i < this->base_Nx; i++)
              DryDepositionFlux_(s, j, i) = 0.5
                * Model.D3("DryDepositionFlux")(this->species_index[s], j, i);
      }
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitDryDeposition<T, ClassModel>::InitStep(ClassModel& Model)
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
  void SaverUnitDryDeposition<T, ClassModel>::Save(ClassModel& Model)
  {
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
        {
          int s, j, i;
          for (s = 0; s < this->Ns; s++)
            for (j = 0; j < this->base_Ny; j++)
              for (i = 0; i < this->base_Nx; i++)
                DryDepositionFlux_(s, j, i) += 0.5
                  * Model.D3("DryDepositionFlux")(this->species_index[s], j, i);

          DryDepositionFlux_.GetArray() /= T(this->interval_length);

          if (Model.GetCurrentDate() >= this->date_beg
              && Model.GetCurrentDate() <= this->date_end)
            for (s = 0; s < this->Ns; s++)
              {
                Data<T, 2>
                  DryDepositionFlux_tmp(&DryDepositionFlux_(s, 0, 0),
                                        shape(this->base_Ny, this->base_Nx));
                FormatBinary<float>().Append(DryDepositionFlux_tmp, output_file[s]);
              }

          for (s = 0; s < this->Ns; s++)
            for (j = 0; j < this->base_Ny; j++)
              for (i = 0; i < this->base_Nx; i++)
                DryDepositionFlux_(s, j, i) = 0.5
                  * Model.D3("DryDepositionFlux")(this->species_index[s], j, i);

          this->counter = 0;
        }
      else
        {
          int s, j, i;
          for (s = 0; s < this->Ns; s++)
            for (j = 0; j < this->base_Ny; j++)
              for (i = 0; i < this->base_Nx; i++)
                DryDepositionFlux_(s, j, i) +=
                  Model.D3("DryDepositionFlux")(this->species_index[s], j, i);
        }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (int s = 0; s < this->Ns; s++)
        {
          Data<T, 2>
            DryDepositionFlux_tmp(&Model.D3("DryDepositionFlux")
                                  (this->species_index[s], 0, 0),
                                  shape(this->base_Ny, this->base_Nx));
          FormatBinary<float>().Append(DryDepositionFlux_tmp, output_file[s]);
        }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDRYDEPOSITION_CXX
#endif
