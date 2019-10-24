// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_AER_CXX


#include "SaverUnitDomain_aer.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  /////////////////////////
  // SAVERUNITDOMAIN_AER //
  /////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitDomain_aer<T, ClassModel>
  ::SaverUnitDomain_aer(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitDomain_aer<T, ClassModel>::~SaverUnitDomain_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "domain_aer".
  */
  template<class T, class ClassModel>
  string SaverUnitDomain_aer<T, ClassModel>::GetType()  const
  {
    return "domain_aer";
  }


  //! First initialization.
  /*!
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex_aer(string)
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNz()
    <li> GetNs_aer()
    <li> GetNbin_aer()
    <li> GetSpeciesList_aer()
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitDomain_aer<T, ClassModel>::Init(ConfigStream& config_stream,
                                                ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Checks whether the number is computed in the model.
    // Use the virtual function "HasNumberConcentration_aer" of BaseModel


    bool model_with_number = Model.HasNumberConcentration_aer();
	

    // Vertical levels to be saved.
    config_stream.Find("Levels");
    split(config_stream.GetLine(), levels);
    Nlevels = int(levels.size());

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species and bins.
    string field, species, file;
    vector<string> vsplit, bounds;
    int first, last, j, k;
    if (this->species_list[0] == "all")
      {
        int i;
        int Ns_aer = Model.GetNs_aer();
        vector<string> list_aer = Model.GetSpeciesList_aer();
        if (model_with_number)
          {
            list_aer.push_back("Number");
            Ns_aer++;
          }
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
      {
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

        if (model_with_number)
	  {
	    output_file.push_back(vector<string>());
	    Concentration_.push_back(vector<Data<T, 3> >());
		
	    vector<int> bins;
	    pair<string, vector<int> > tmp;
	    for (int i = 0; i < Model.GetNbin_aer(); i++)
	      bins.push_back(i);
	    tmp.second = bins;
	    tmp.first = "Number";
	    species_list_aer.push_back(tmp);
			
	    int rear = output_file.size()-1;
	    for (k = 0; k < int(bins.size()); k++)
	      {
		file = find_replace(filename, "&f", "Number");
		file = find_replace(file, "&n", to_str(bins[k]));
		output_file[rear].push_back(file);
		Concentration_[rear].push_back(Data<T, 3>());
	      }
	  }
      }

    // Indices in the model.
    int base_s, base_b;

    int s, b, i;
    // Empties output files.
    for (s = 0; s < int(species_list_aer.size()); s++)
      for (b = 0; b < int(species_list_aer[s].second.size()); b++)
        ofstream tmp_stream(output_file[s][b].c_str());

    if (this->averaged)
      for (s = 0; s < int(species_list_aer.size()); s++)
	for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	  {
	    if (species_list_aer[s].first != "Number")
	      {
		base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
		base_b = species_list_aer[s].second[b];
		Concentration_[s][b].Resize(Nlevels,
					    this->base_Ny, this->base_Nx);
		for (k = 0; k < Nlevels; k++)
		  for (j = 0; j < this->base_Ny; j++)
		    for (i = 0; i < this->base_Nx; i++)
		      Concentration_[s][b](k, j, i) = 0.5
			* Model.GetConcentration_aer()(base_s, base_b,
						       levels[k], j, i);
	      }
	    else if (model_with_number)
	      {
		base_b = species_list_aer[s].second[b];
		Concentration_[s][b].Resize(Nlevels,
					    this->base_Ny, this->base_Nx);
		for (k = 0; k < Nlevels; k++)
		  for (j = 0; j < this->base_Ny; j++)
		    for (i = 0; i < this->base_Nx; i++)
		      Concentration_[s][b](k, j, i) = 0.5
			* Model.GetNumberConcentration_aer()(base_b,
							     levels[k], j, i);
	      }
	  }
	
    if (this->initial_concentration && !this->averaged)
      Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitDomain_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration_aer()
    <li> GetSpeciesIndex_aer()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitDomain_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int s, b, k, j, i;
    // Indices in the model.
    int base_s, base_b;
    base_s = -999;
    bool model_with_number = Model.HasNumberConcentration_aer();
	
	
    if (this->averaged)
      {
        if (this->counter % this->interval_length == 0)
          for (s = 0; s < int(species_list_aer.size()); s++)
            for (b = 0; b < int(species_list_aer[s].second.size()); b++)
              {
                base_b = species_list_aer[s].second[b];
                
                if (species_list_aer[s].first != "Number")
                  {
                    base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) += 0.5
                            * Model.GetConcentration_aer()(base_s, base_b,
                                                           levels[k], j, i);
                  }
                else if (model_with_number)		
                  {
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) += 0.5
                            * Model.GetNumberConcentration_aer()(base_b,
                                                                 levels[k], j, i);
                  }

                Concentration_[s][b].GetArray() /= T(this->interval_length);

                if (Model.GetCurrentDate() >= this->date_beg
                    && Model.GetCurrentDate() <= this->date_end)
                  FormatBinary<float>().Append(Concentration_[s][b],
					     output_file[s][b]);

                if (species_list_aer[s].first != "Number")
                  {
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) = 0.5
                            * Model.GetConcentration_aer()(base_s, base_b,
                                                           levels[k], j, i);
                  }
                else if (model_with_number)
                  {
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) = 0.5
                            * Model.GetNumberConcentration_aer()(base_b,
                                                                 levels[k], j, i);
                  }

                this->counter = 0;
              }
        
        else
          for (s = 0; s < int(species_list_aer.size()); s++)
            for (b = 0; b < int(species_list_aer[s].second.size()); b++)
              {
                base_b = species_list_aer[s].second[b];
                if (species_list_aer[s].first != "Number")
                  {
                    base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) +=
                            Model.GetConcentration_aer()(base_s, base_b,
                                                         levels[k], j, i);
                  }
                else if (model_with_number)
                  {
                    for (k = 0; k < Nlevels; k++)
                      for (j = 0; j < this->base_Ny; j++)
                        for (i = 0; i < this->base_Nx; i++)
                          Concentration_[s][b](k, j, i) +=
                            Model.GetNumberConcentration_aer()(base_b,
                                                               levels[k], j, i);
                  }
              }
      }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (s = 0; s < int(species_list_aer.size()); s++)
        for (b = 0; b < int(species_list_aer[s].second.size()); b++)
          {
            base_b = species_list_aer[s].second[b];

            for (int k = 0; k < Nlevels; k++)
              {
		if (species_list_aer[s].first != "Number")
		  {
                    base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);

                    Data<T, 2> Concentration_tmp(&Model.GetConcentration_aer()
                                                 (base_s, base_b,
                                                  levels[k], 0, 0),
                                                 shape(this->base_Ny,
                                                       this->base_Nx));
                    FormatBinary<float>().Append(Concentration_tmp,
                                                 output_file[s][b]);
                  }
		else if (model_with_number)
		  {
		    Data<T, 2> NumberConcentration_tmp(&Model.GetNumberConcentration_aer()
						       (base_b,
							levels[k], 0, 0),
						       shape(this->base_Ny,
							     this->base_Nx));
		    FormatBinary<float>().Append(NumberConcentration_tmp,
						 output_file[s][b]);
		  }
              }
          }
  }
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_AER_CXX
#endif
