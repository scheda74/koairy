// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_AER_HXX


#include "BaseSaverUnit.hxx"
#include <vector>
#include <fstream>
#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////
  // SAVERUNITPOINT_AER //
  ////////////////////////


  /*! \brief This class is used to save concentrations over the whole domain
    of the underlying model. Chemical species and vertical levels may be
    selected. */
  template<class T, class ClassModel>
  class SaverUnitPoint_aer: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! Type of saver (indices or coordinates).
    string type;
    //! Boolean: true if the type of saver is "indices_list", false otherwise.
    bool with_indices_list;

    //! List of aerosol species (with their bins) to be saved.
    vector<pair<string, vector<int> > > species_list_aer;

    //! List of vertical levels to be saved (indices).
    vector<int> levels;
    //! List of vertical levels to be saved (coordinates).
    vector<T> levels_coord;
    //! Number of levels to be saved.
    int Nlevels;

    //! List of indices of points along which concentrations are saved.
    vector<vector<int> > point_list;
    //! List of coordinates of points along which concentrations are saved.
    vector<vector<T> > coord_list;
    //! Number of points to be saved.
    int Npoint_aer;

    //! List of output files.
    vector<vector<string> > output_file;

    /*! Buffer used to compute averaged concentrations. It is first indexed by
      the species name and the bin. */
    vector<vector<Data<T, 1> > > Concentration_;

  public:

    SaverUnitPoint_aer();
    virtual ~SaverUnitPoint_aer();

    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITPOINT_AER_HXX
#endif
