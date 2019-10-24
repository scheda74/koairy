// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDRYDEPOSITION_AER_HXX


#include "BaseSaverUnit.hxx"
#include <vector>
#include <fstream>
#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////////
  // SAVERUNITDRYDEPOSITION_AER //
  ////////////////////////////////


  /*! \brief This class is used to save aerosol dry deposition fluxes over the
    whole domain of the underlying model. Species and bins may be selected. */
  template<class T, class ClassModel>
  class SaverUnitDryDeposition_aer: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! List of aerosol species (with their bins) to be saved.
    vector<pair<string, vector<int> > > species_list_aer;

    //! List of output files.
    vector<vector<string> > output_file;

    /*! Buffer used to compute averaged concentrations. It is first indexed by
      the species name and the bin. */
    vector<vector<Data<T, 2> > > DryDepositionFlux_aer_;

  public:

    SaverUnitDryDeposition_aer();
    virtual ~SaverUnitDryDeposition_aer();

    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITDRYDEPOSITION_AER_HXX
#endif
