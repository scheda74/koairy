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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPUFFEMISSION_HXX


#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


#include <vector>
#include <fstream>
  using namespace std;

#include "AtmoData.hxx"
  using namespace AtmoData;


  ////////////////////////////
  // SAVERUNITDRYDEPOSITION //
  ////////////////////////////


  /*! \brief This class is used to save puff emissions in the Gaussian puff models. Chemical species may be selected. */
  template<class T, class ClassModel>
  class SaverUnitPuffEmission: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! List of output files.
    vector<string> output_file;

    //! Buffer used to compute averaged puff emissions.
    Data<T, 4> PuffEmission_;

    int Ns_emis, Nz_emis;

  public:

    SaverUnitPuffEmission();
    virtual ~SaverUnitPuffEmission();

    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPUFFEMISSION_HXX
#endif
