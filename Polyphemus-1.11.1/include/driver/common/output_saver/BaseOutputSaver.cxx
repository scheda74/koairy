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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_BASEOUTPUTSAVER_CXX


#include "BaseOutputSaver.hxx"

#include "SaverUnitDomain.cxx"
#include "SaverUnitDryDeposition.cxx"
#include "SaverUnitWetDeposition.cxx"
#include "SaverUnitDomain_aer.cxx"
#include "SaverUnitDryDeposition_aer.cxx"
#include "SaverUnitWetDeposition_aer.cxx"
#include "SaverUnitDomain_assimilation.cxx"
#include "SaverUnitDomain_prediction.cxx"
#include "SaverUnitNesting.cxx"
#include "SaverUnitNesting_aer.cxx"
#include "SaverUnitPoint.cxx"
#include "SaverUnitPoint_aer.cxx"
#include "SaverUnitPuffEmission.cxx"
#include "SaverUnitSubdomain.cxx"
#include "SaverUnitSubdomain_aer.cxx"
#include "SaverUnitBackup.cxx"
#include "SaverUnitBackup_aer.cxx"
#include "SaverUnitStreet.cxx"

#include "AtmoData.hxx"


namespace Polyphemus
{


  /////////////////////
  // BASEOUTPUTSAVER //
  /////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  BaseOutputSaver<T, ClassModel>::BaseOutputSaver()
    : group("all")
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  BaseOutputSaver<T, ClassModel>::~BaseOutputSaver()
  {
    for (unsigned int i = 0; i < Units.size(); i++)
      delete Units[i];
  }


  //! Sets the group of active saver units.
  /*! Only saver units with the same group as \a grp will be taken into
    account (in 'Init', 'InitStep' and 'Save'). These saver units are called
    active savers.
    \param grp saver group to be activated. Put "all" to activate all saver
    units.
  */
  template<class T, class ClassModel>
  void BaseOutputSaver<T, ClassModel>::SetGroup(string grp)
  {
    group = grp;
  }


  //! Gets the name of the group of active saver units.
  /*!
    \return The group of active saver units.
  */
  template<class T, class ClassModel>
  string BaseOutputSaver<T, ClassModel>::GetGroup() const
  {
    return group;
  }


  //! Check if there is only coordinates_list in the Unit list.
  /*!
    \return True The type of saver.
  */
  template<class T, class ClassModel>
  bool BaseOutputSaver<T, ClassModel>::OnlyCoordinatesList() const
  {
    bool OnlyCoordinatesList = 1;
    for (int i = 0; i < int(Units.size()); i++)
      if (Units[i]->GetType() != "coordinates_list")
        OnlyCoordinatesList = 0;
    return OnlyCoordinatesList;
  }


  //! First initializations.
  /*! Reads the configuration.
    \param Model model with the following interface:
    <ul>
    <li> GetConfigurationFile()
    </ul>
  */
  template<class T, class ClassModel>
  void BaseOutputSaver<T, ClassModel>::Init(ClassModel& Model)
  {
    // Clear saver units.
    for (unsigned int i = 0; i < Units.size(); i++)
      delete Units[i];
    Units.clear();

    string config_file = Model.GetConfigurationFile();
    // Retrieves the path to the configuration file for output savers.
    ConfigStream main_config_stream(config_file);
    string saver_config_file;
    main_config_stream.Find("[output]");
    main_config_stream.GetValue("Configuration_file", saver_config_file);
    ConfigStream config_stream(saver_config_file);

    // Browses all sections "[save]" that define saver units.
    string line, type;
    int i = 0;
    while (!config_stream.IsEmpty())
      {
        line = config_stream.GetLine();
        if (split(line)[0] == "[save]")
          {
            config_stream.PeekValue("Type", type);
            if (type == "domain" || type == "domain_ensemble_forecast"
                || type == "domain_ensemble_analysis")
              Units.push_back(new SaverUnitDomain<T, ClassModel>());
            else if (type == "dry_deposition")
              Units.push_back(new SaverUnitDryDeposition<T, ClassModel>());
            else if (type == "wet_deposition")
              Units.push_back(new SaverUnitWetDeposition<T, ClassModel>());
            else if (type == "domain_aer")
              Units.push_back(new SaverUnitDomain_aer<T, ClassModel>());
            else if (type == "dry_deposition_aer")
              Units.push_back(new
                              SaverUnitDryDeposition_aer<T, ClassModel>());
            else if (type == "wet_deposition_aer")
              Units.push_back(new
                              SaverUnitWetDeposition_aer<T, ClassModel>());
            else if (type == "domain_assimilation")
              Units.
                push_back(new SaverUnitDomain_assimilation<T, ClassModel>());
            else if (type == "domain_prediction")
              Units.
                push_back(new SaverUnitDomain_prediction<T, ClassModel>());
            else if (type == "nesting")
              Units.push_back(new SaverUnitNesting<T, ClassModel>());
            else if (type == "nesting_aer")
              Units.push_back(new SaverUnitNesting_aer<T, ClassModel>());
            else if (type == "subdomain")
              Units.push_back(new SaverUnitSubdomain<T, ClassModel>());
            else if (type == "subdomain_aer")
              Units.push_back(new SaverUnitSubdomain_aer<T, ClassModel>());
            else if (type == "indices_list" || type == "coordinates_list")
              Units.push_back(new SaverUnitPoint<T, ClassModel>());
            else if (type == "indices_list_aer"
                     || type == "coordinates_list_aer")
              Units.push_back(new SaverUnitPoint_aer<T, ClassModel>());
            else if (type == "backup")
              Units.push_back(new SaverUnitBackup<T, ClassModel>());
            else if (type == "backup_aer")
              Units.push_back(new SaverUnitBackup_aer<T, ClassModel>());
            else if (type == "puff_emission")
              Units.push_back(new SaverUnitPuffEmission<T, ClassModel>());
     	    else if (type == "street")
	          Units.push_back(new SaverUnitStreet<T, ClassModel>());
            else
              throw string("Unknown saver unit: ") + type + ".";
            Units[i]->Init(config_stream, Model);
            i++;
          }
      }
  }


  //! Initializes the savers at the beginning of each step.
  /*!
    \param Model model.
  */
  template<class T, class ClassModel>
  void BaseOutputSaver<T, ClassModel>::InitStep(ClassModel& Model)
  {
    for (unsigned int i = 0; i < Units.size(); i++)
      Units[i]->InitStep(Model);
  }


  //! Calls the methods 'Save' of all unit savers.
  /*!
    \param Model model.
  */
  template<class T, class ClassModel>
  void BaseOutputSaver<T, ClassModel>::Save(ClassModel& Model)
  {
    for (unsigned int i = 0; i < Units.size(); i++)
      if (group == "all" || Units[i]->GetGroup() == group)
        Units[i]->Save(Model);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_BASEOUTPUTSAVER_CXX
#endif
