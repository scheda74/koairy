// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_DRIVER_ASSIMILATIONDRIVER_CXX


#include "AssimilationDriver.hxx"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::AssimilationDriver(string config_file):
    BaseDriver<T, ClassModel, ClassOutputSaver>(config_file),
    ObsManager(config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~AssimilationDriver()
  {
  }


  //! Read configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    this->config.SetSection("[data_assimilation]");
    this->config.PeekValue("Nt_assimilation", Nt_assim);
    if (Nt_assim == -1)
      Nt_assim = this->Model.GetNt();
  }


  //! Driver initialization.
  /*! It reads configurations */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    this->ReadConfiguration();

    Nt_predict = this->Model.GetNt() - Nt_assim;

    if (Nt_predict < 0)
      throw string("Error: the assimilation window is longer")
        + " than the simulation period.";

    Nstate = this->Model.GetNstate();
  }


  //! Empty method.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    throw string("No assimilation method is implemented in class")
      + string(" AssimilationDriver.");
  }


  //////////////////////////////////////////
  // STORAGE MANAGEMENTS FOR PERTURBATION //
  //////////////////////////////////////////


  //! Allocates field arrays.
  /*! It allocates the arrays for the safeguard of the unperturbed field data.
   */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::AllocateFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      {
        if (this->Model.GetFieldDimension(field_list[i]) == 2)
          A2_perturb_map[field_list[i]] =
            new Array<T, 2>(this->Model.A2(field_list[i]).shape());
        if (this->Model.GetFieldDimension(field_list[i]) == 3)
          A3_perturb_map[field_list[i]] =
            new Array<T, 3>(this->Model.A3(field_list[i]).shape());
        if (this->Model.GetFieldDimension(field_list[i]) == 4)
          A4_perturb_map[field_list[i]] =
            new Array<T, 4>(this->Model.A4(field_list[i]).shape());
        if (this->Model.GetFieldDimension(field_list[i]) == 5)
          A5_perturb_map[field_list[i]] =
            new Array<T, 5>(this->Model.A5(field_list[i]).shape());
      }
  }


  //! Deallocates field arrays.
  /*! It frees the arrays allocated for the safeguard of the unperturbed field
    data.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::DeallocateFieldArray()
  {
    typename map<string, Array<T, 2>* >::iterator iter2;
    if (!A2_perturb_map.empty())
      for (iter2 = A2_perturb_map.begin();
           iter2 != A2_perturb_map.end(); iter2++)
        delete iter2->second;
    typename map<string, Array<T, 3>* >::iterator iter3;
    if (!A3_perturb_map.empty())
      for (iter3 = A3_perturb_map.begin();
           iter3 != A3_perturb_map.end(); iter3++)
        delete iter3->second;
    typename map<string, Array<T, 4>* >::iterator iter4;
    if (!A4_perturb_map.empty())
      for (iter4 = A4_perturb_map.begin();
           iter4 != A4_perturb_map.end(); iter4++)
        delete iter4->second;
    typename map<string, Array<T, 5>* >::iterator iter5;
    if (!A5_perturb_map.empty())
      for (iter5 = A5_perturb_map.begin();
           iter5 != A5_perturb_map.end(); iter5++)
        delete iter5->second;
  }


  //! Sets field arrays to model field data.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::SetFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      {
        if (this->Model.GetFieldDimension(field_list[i]) == 2)
          *A2_perturb_map[field_list[i]]
            = this->Model.A2(field_list[i]);
        if (this->Model.GetFieldDimension(field_list[i]) == 3)
          *A3_perturb_map[field_list[i]]
            = this->Model.A3(field_list[i]);
        if (this->Model.GetFieldDimension(field_list[i]) == 4)
          *A4_perturb_map[field_list[i]]
            = this->Model.A4(field_list[i]);
        if (this->Model.GetFieldDimension(field_list[i]) == 5)
          *A5_perturb_map[field_list[i]]
            = this->Model.A5(field_list[i]);
      }
  }


  //! Sets model field data from data stored in field arrays.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::SetModelFieldArray()
  {
    unsigned long i;
    vector<string> field_list_tmp = PerturbManager.GetPerturbDataList();

    vector<string> field_list;
    for (i = 0; i < field_list_tmp.size(); i++)
      if (this->Model.HasField(field_list_tmp[i] + "_i")
          && this->Model.HasField(field_list_tmp[i] + "_f"))
        {
          field_list.push_back(field_list_tmp[i] + "_i");
          field_list.push_back(field_list_tmp[i] + "_f");
        }
      else
        field_list.push_back(field_list_tmp[i]);

    for (i = 0; i < field_list.size(); i++)
      {
        if (this->Model.GetFieldDimension(field_list[i]) == 2)
          this->Model.A2(field_list[i])
            = *A2_perturb_map[field_list[i]];
        if (this->Model.GetFieldDimension(field_list[i]) == 3)
          this->Model.A3(field_list[i])
            = *A3_perturb_map[field_list[i]];
        if (this->Model.GetFieldDimension(field_list[i]) == 4)
          this->Model.A4(field_list[i])
            = *A4_perturb_map[field_list[i]];
        if (this->Model.GetFieldDimension(field_list[i]) == 5)
          this->Model.A5(field_list[i])
            = *A5_perturb_map[field_list[i]];
      }
  }


  //! Generates random numbers for the field data perturbation.
  /*! It generates random numbers for the field data perturbation. For each
    member in the ensemble, the generated randoms numbers are stored in a map
    that maps parameter names in their field to the corresponding
    random numbers.
    \param Nens ensemble number.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::SetPerturbNumEns(int Nens)
  {
    map<string, vector<T> > PerturbNum_map;
    vector<string> field_list = PerturbManager.GetPerturbFieldList();

    PerturbNum_ens.clear();

    for (int s = 0; s < Nens; s++)
      {
        PerturbNum_map = PerturbManager.GeneratePerturbNumMap();
        PerturbNum_ens.push_back(PerturbNum_map);
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ASSIMILATIONDRIVER_CXX
#endif
