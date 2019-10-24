// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Edouard Debry
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITBACKUP_HXX


#include "BaseSaverUnit.hxx"

#include <vector>
#include <fstream>

#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // SAVERUNITBACKUP //
  /////////////////////


  /*! \brief This class is used to perform a backup saving of concentrations
    over the whole domain of the underlying model in order to restart it. */
  template<class T, class ClassModel>
  class SaverUnitBackup: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! Type of saver.
    string type;

    //! List of output files.
    vector<string> output_file;

    //! Date file.
    string date_file;

  public:

    SaverUnitBackup();
    virtual ~SaverUnitBackup();

    virtual string GetType() const;
    virtual string GetGroup() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITBACKUP_HXX
#endif
