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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITNESTING_AER_CXX


#include "SaverUnitNesting_aer.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  ////////////////////////////
  // SAVERUNITNESTING_AER  //
  //////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitNesting_aer<T, ClassModel>
  ::SaverUnitNesting_aer(): BaseSaverUnit<T, ClassModel>()
  {
  }

  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitNesting_aer<T, ClassModel>::~SaverUnitNesting_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "nesting_aer".
  */
  template<class T, class ClassModel>
  string SaverUnitNesting_aer<T, ClassModel>::GetType()  const
  {
    return "nesting_aer";
  }


  //! First initialization.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetGridZArray1D()
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNs_aer()
    <li> GetNbin_aer()
    <li> GetSpeciesList_aer()
    <li> GetSpeciesIndex_aer(string)
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitNesting_aer<T, ClassModel>::Init(ConfigStream& config_stream,
                                                 ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    int s, i, b;

    // Grids for the domain.
    GridX = Model.GetGridXArray1D();
    GridY = Model.GetGridYArray1D();
    GridZ = Model.GetGridZArray1D();

    // Subdomain dimensions.
    config_stream.PeekValue("x_min", x_min_sub);
    config_stream.PeekValue("Delta_x", Delta_x_sub);
    config_stream.PeekValue("Nx", Nx_sub);

    config_stream.PeekValue("y_min", y_min_sub);
    config_stream.PeekValue("Delta_y", Delta_y_sub);
    config_stream.PeekValue("Ny", Ny_sub);

    config_stream.PeekValue("Vertical_levels", file_vertical_levels_sub);
    config_stream.PeekValue("Nz", Nz_sub);
    LayerInterface.resize(Nz_sub + 1);
    FormatText().Read(file_vertical_levels_sub, LayerInterface);

    // Grids for the subdomain.
    GridX_sub = RegularGrid<T>(x_min_sub, Delta_x_sub, Nx_sub);
    GridY_sub = RegularGrid<T>(y_min_sub, Delta_y_sub, Ny_sub);
    GridZ_sub = RegularGrid<T>(Nz_sub);
    for (i = 0; i < Nz_sub; i++)
      GridZ_sub(i) =  0.5 * (LayerInterface(i + 1) + LayerInterface(i));

    // Grids for boundary conditions.
    RegularGrid<T> GridX_bc_x(2);
    GridX_bc_x(0) = x_min_sub - Delta_x_sub;
    GridX_bc_x(1) = x_min_sub + T(Nx_sub) * Delta_x_sub;

    RegularGrid<T> GridY_bc_y(2);
    GridY_bc_y(0) = y_min_sub - Delta_y_sub;
    GridY_bc_y(1) = y_min_sub + T(Ny_sub) * Delta_y_sub;

    RegularGrid<T> GridZ_bc_z(1);
    GridZ_bc_z(0) = min(GridZ(GridZ.GetLength() - 1),
                        1.5 * LayerInterface(Nz_sub) -
                        0.5 * LayerInterface(Nz_sub - 1));

    // Checks that interpolation is feasible.
    T x_min = Model.GetX_min();
    T Delta_x = Model.GetDelta_x();
    int Nx = Model.GetNx();
    if (GridX_bc_x(0) < x_min || GridX_bc_x(1) > x_min + T(Nx) * Delta_x)
      throw string("Error in SaverUnitNesting_aer. Interpolation impossible")
        + "along x: the subdomain is not included in the domain.";

    T y_min = Model.GetY_min();
    T Delta_y = Model.GetDelta_y();
    int Ny = Model.GetNy();
    if (GridY_bc_y(0) < y_min || GridY_bc_y(1) > y_min + T(Ny) * Delta_y)
      throw string("Error in SaverUnitNesting_aer.Interpolation impossible")
        + "along y: the subdomain is not included in the domain.";

    Concentration_bc_x_aer.Resize(GridZ_sub, GridY_sub, GridX_bc_x);
    Concentration_bc_y_aer.Resize(GridZ_sub, GridY_bc_y, GridX_sub);
    Concentration_bc_z_aer.Resize(GridZ_bc_z, GridY_sub, GridX_sub);

    // Species to be saved.
    string field, species;
    vector<string> vsplit, bins;
    string filename_ref = config_stream.PeekValue("Output_file");
    string filename;
    if (this->species_list[0] == "all")
      {
        vector<string> list_aer = Model.GetSpeciesList_aer();
        vector<int> bins;
        pair<string, vector<int> > tmp;
        for (b = 0; b < Model.GetNbin_aer(); b++)
          bins.push_back(b);
        tmp.second = bins;
        for (s = 0; s < Model.GetNs_aer(); s++)
          {
            tmp.first = list_aer[s];
            species_list_aer.push_back(tmp);
          }
      }
    else
      for (s = 0; s < int(this->species_list.size()); s++)
        {
          field = this->species_list[s];
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

          species = split(field, "_")[0];
          int j = 0;
          while (j < int(species_list_aer.size())
                 && species_list_aer[j].first != species)
            j++;
          if (j == int(species_list_aer.size()))
            species_list_aer.push_back(pair<string, vector<int> >
                                       (species, vector<int>()));

          bins = split(vsplit[1], "-");
          if (bins.size() > 2)
            throw string("Species \"") + field + string("\" is badly ")
              + "formatted: it cannot be parsed by the output saver.";
          int bins_beg = convert<int>(bins[0]);
          int bins_end = convert<int>(bins[bins.size() - 1]);
          for (i = bins_beg; i < bins_end + 1; i++)
            species_list_aer[j].second.push_back(i);
        }
    this->Ns = int(species_list_aer.size());

    // Output filenames.
    output_file.resize(this->Ns);
    for (s = 0; s < this->Ns; s++)
      {
        int Nb = int(species_list_aer[s].second.size());
        output_file[s].resize(Nb);
        for (b = 0; b < Nb; b++)
          {
            filename = find_replace(filename_ref, "&f",
                                    species_list_aer[s].first);
            output_file[s][b].resize(3);
            string bin = to_str(species_list_aer[s].second[b]);
            filename = find_replace(filename, "&n", bin);
            output_file[s][b][0] = find_replace(filename, "&c", "x");
            output_file[s][b][1] = find_replace(filename, "&c", "y");
            output_file[s][b][2] = find_replace(filename, "&c", "z");
          }
      }

    // Empties output files.
    for (s = 0; s < this->Ns; s++)
      for (b = 0; b < int(species_list_aer[s].second.size()); b++)
        for (i = 0; i < 3; i++)
          ofstream tmp_stream(output_file[s][b][i].c_str());

    Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitNesting_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetNx()
    <li> GetNy()
    <li> GetNz()
    <li> GetCurrentDate()
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitNesting_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int Nz = Model.GetNz();
    int Ny = Model.GetNy();
    int Nx = Model.GetNx();

    if (this->counter % this->interval_length == 0 &&
        Model.GetCurrentDate() >= this->date_beg &&
        Model.GetCurrentDate() <= this->date_end)
      for (int s = 0; s < this->Ns;  s++)
        {
          int gs = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
          for (int b = 0; b < int(species_list_aer[s].second.size()); b++)
            {
              int lb = species_list_aer[s].second[b];
              Data<T, 3> Concentration_in(&Model.GetConcentration_aer()
                                          (gs, lb, 0, 0, 0),
                                          shape(Nz, Ny, Nx));
              Concentration_in.Resize(GridZ, GridY, GridX);


              LinearInterpolationRegular(Concentration_in,
                                         Concentration_bc_x_aer);
              LinearInterpolationRegular(Concentration_in,
                                         Concentration_bc_y_aer);
              LinearInterpolationRegular(Concentration_in,
                                         Concentration_bc_z_aer);

              FormatBinary<float>().Append(Concentration_bc_x_aer,
                                           output_file[s][b][0]);
              FormatBinary<float>().Append(Concentration_bc_y_aer,
                                           output_file[s][b][1]);
              FormatBinary<float>().Append(Concentration_bc_z_aer,
                                           output_file[s][b][2]);
            }
        }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITNESTING_AER_CXX
#endif
