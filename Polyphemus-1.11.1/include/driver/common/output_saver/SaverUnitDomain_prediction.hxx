// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_PREDICTION_HXX


#include "SaverUnitDomain.hxx"

#include <vector>
#include <fstream>

#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // SAVERUNITDOMAIN //
  /////////////////////


  /*! \brief This class is used to save concentrations over the whole domain
    of the underlying model, in the context of data assimilation. Chemical
    species and vertical levels may be selected. */
  template<class T, class ClassModel>
  class SaverUnitDomain_prediction: public SaverUnitDomain<T, ClassModel>
  {

  protected:

    /*! Name of the file in which dates associated with the analyses are
      saved. */
    string date_file;

  public:

    SaverUnitDomain_prediction();
    virtual ~SaverUnitDomain_prediction();

    virtual string GetType() const;

    virtual string GetGroup() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);

    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDOMAIN_PREDICTION_HXX
#endif
