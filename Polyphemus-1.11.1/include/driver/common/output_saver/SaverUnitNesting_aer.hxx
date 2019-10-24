// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITNESTING_AER_HXX


#include "BaseSaverUnit.hxx"
#include <vector>
#include <fstream>
#include <utility>
#include "AtmoDataHeader.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////
  // SAVERUNITNESTING_AER  //
  //////////////////////////


  /*! \brief This class is used to save concentrations for aerosol species
    over a subdomain in order to perform nested simulations. What is saved is
    then used as a boundary condition for the smaller scale simulation. */
  template<class T, class ClassModel>
  class SaverUnitNesting_aer: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! List of aerosol species (with their bins) to be saved.
    vector<pair<string, vector<int> > > species_list_aer;

    //! Origin of the subdomain along x.
    T x_min_sub;
    //! Step along x for the subdomain.
    T Delta_x_sub;
    //! Number of points along x for the subdomain.
    int Nx_sub;
    //! Origin of the subdomain along y.
    T y_min_sub;
    //! Step along y for the subdomain.
    T Delta_y_sub;
    //! Number of points along y for the subdomain.
    int Ny_sub;
    //! Number of points along z for the subdomain.
    int Nz_sub;
    //! File storing the vertical levels.
    string file_vertical_levels_sub;
    //! Array  storing the vertical levels.
    Array<T, 1> LayerInterface;

    //! Grid for the subdomain along x.
    RegularGrid<T> GridX_sub;
    //! Grid for the subdomain along y.
    RegularGrid<T> GridY_sub;
    //! Grid for the subdomain along z.
    RegularGrid<T> GridZ_sub;
    //! Grid for the domain along x.
    RegularGrid<T> GridX;
    //! Grid for the domain along y.
    RegularGrid<T> GridY;
    //! Grid for the domain along z.
    RegularGrid<T> GridZ;
    //! Grid for boundary conditions along x.
    RegularGrid<T> GridX_bc_x;
    //! Grid for boundary conditions along y.
    RegularGrid<T> GridY_bc_y;
    //! Grid for boundary conditions along z.
    RegularGrid<T> GridZ_bc_z;

    //! Data used to store the concentrations interpolated along x.
    Data<T, 3> Concentration_bc_x_aer;
    //! Data used to store the concentrations interpolated along y.
    Data<T, 3> Concentration_bc_y_aer;
    //! Data used to store the concentrations interpolated along z.
    Data<T, 3> Concentration_bc_z_aer;

    //! List of the output files.
    vector<vector<vector<string> > > output_file;

  public:

    SaverUnitNesting_aer();
    virtual ~SaverUnitNesting_aer();

    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITNESTING_AER_HXX
#endif
