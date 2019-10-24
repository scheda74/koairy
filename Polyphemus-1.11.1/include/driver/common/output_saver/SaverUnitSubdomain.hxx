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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_HXX


#include "BaseSaverUnit.hxx"

#include <vector>
#include <fstream>

#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////////
  // SAVERUNITSUBDOMAIN //
  ///////////////////////


  /*! \brief This class is used to save concentrations over a subdomain
    of the underlying model. Chemical species and vertical levels may be
    selected. */
  template<class T, class ClassModel>
  class SaverUnitSubdomain: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    int i_min;
    int i_max;
    int j_min;
    int j_max;

    //! List of vertical levels to be saved.
    vector<int> levels;
    //! Number of levels to be saved.
    int Nlevels;

    //! List of output files.
    vector<string> output_file;

    //! Buffer used to compute averaged concentrations.
    Data<T, 4> Concentration_;

  public:

    SaverUnitSubdomain();
    virtual ~SaverUnitSubdomain();

    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSUBDOMAIN_HXX
#endif
