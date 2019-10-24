// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Damien Garaud
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


#ifndef POLYPHEMUS_FILE_DRIVER_PERTURBATIONDRIVER_CXX


#include "PerturbationDriver.hxx"


namespace Polyphemus
{


  //! Main constructor.
  /*! Builds the driver and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  PerturbationDriver<T, ClassModel, ClassOutputSaver>
  ::PerturbationDriver(string config_file):
    Model(config_file), config(config_file)
  {

    /*** Display options ***/

    config.SetSection("[display]");
    // Should iterations be displayed on screen?
    config.PeekValue("Show_iterations", option_display["show_iterations"]);
    // Should current date be displayed on screen?
    config.PeekValue("Show_date", option_display["show_date"]);

    /*** Perturbations ***/

    config.SetSection("[data]");

    // Path to the file that describes the perturbations.
    string perturbation_file;
    config.PeekValue("Perturbation_file", perturbation_file);

    ConfigStream config_perturbation(perturbation_file);

    // Reads the list of fields that do not depend on the species.
    config_perturbation.SetSection("[AdditionalField]");
    string field_name;
    string perturbation_type;
    T perturbation_number;
    while (!config_perturbation.IsEmpty())
      {
        config_perturbation.GetElement(field_name);
        config_perturbation.GetElement(perturbation_type);
        if (perturbation_type != "multiply" &&  perturbation_type != "add")
          throw string("Error in \"PerturbationDriver::PertubationDriver\":")
            + "\n  The perturbation type must be \"add\" or \"multiply\", "
            + "but the type \"" + perturbation_type + "\" was given.";
        config_perturbation.GetElement(perturbation_number);
        additional_field.push_back(field_name);
        additional_field_type.push_back(perturbation_type);
        additional_field_perturbation.push_back(perturbation_number);
      }

    // Reads the list of fields that depend on the species.
    config_perturbation.SetSection("[general]");
    config_perturbation.Find("Field_list");
    string pollutant;
    string fields_string = config_perturbation.GetLine();
    // If there is any field.
    if (trim(fields_string) != "---" && trim(fields_string) != "--"
        && trim(fields_string) != "-")
      field_list = split(fields_string);
    for (int i = 0; i < int(field_list.size()); i++)
      {
        field[field_list[i]] = vector<string>();
        field_type[field_list[i]] = vector<string>();
        field_perturbation[field_list[i]] = vector<T>();

        config_perturbation.SetSection("[" + field_list[i] + "]");
        while (!config_perturbation.IsEmpty())
          {
            config_perturbation.GetElement(pollutant);
            config_perturbation.GetElement(perturbation_type);
            if (perturbation_type != "multiply"
                &&  perturbation_type != "add")
              throw string("Error in \"PerturbationDriver")
                + "::PertubationDriver\":\n"
                + "  The perturbation type must be \"add\" or \"multiply\", "
                + "but the type \"" + perturbation_type + "\" was given.";
            config_perturbation.GetElement(perturbation_number);
            field[field_list[i]].push_back(pollutant);
            field_type[field_list[i]].push_back(perturbation_type);
            field_perturbation[field_list[i]].push_back(perturbation_number);
          }
      }
  }


  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  PerturbationDriver<T, ClassModel, ClassOutputSaver>::~PerturbationDriver()
  {
  }


  //! Performs the simulation.
  /*! Initializes the model and the output saver, and then performs the time
    loop with calls to the model and to the output saver.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void PerturbationDriver<T, ClassModel, ClassOutputSaver>::Run()
  {
    int i, j, s, s_local;

    int rank;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Init();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    /*** Initializations ***/

    Model.Init();
    if (rank == 0)
      OutputSaver.Init(Model);

    /*** Time loop ***/

    for (i = 0; i < Model.GetNt(); i++)
      {
        if (rank == 0 && option_display["show_iterations"])
          cout << "Performing iteration #" << i << endl;

        if (rank == 0 && option_display["show_date"])
          cout << "Current date: "
               << Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        Model.InitStep();
        if (rank == 0)
          OutputSaver.InitStep(Model);

        /*** Perturbs the fields that do not depend on the species ***/

        for (j = 0; j < int(additional_field.size()); j++)
          if (additional_field[j] == "WindAngle")
            WindAngle(additional_field_type[j],
                      additional_field_perturbation[j]);
          else if (additional_field[j] == "WindModule")
            WindModule(additional_field_type[j],
                       additional_field_perturbation[j]);
          else if (Model.GetFieldDimension(additional_field[j]) == 2)
            Apply(Model.A2(additional_field[j]), additional_field_type[j],
                  additional_field_perturbation[j]);
          else if (Model.GetFieldDimension(additional_field[j]) == 3)
            Apply(Model.A3(additional_field[j]), additional_field_type[j],
                  additional_field_perturbation[j]);
          else if (Model.GetFieldDimension(additional_field[j]) == 4)
            Apply(Model.A4(additional_field[j]), additional_field_type[j],
                  additional_field_perturbation[j]);
          else if (Model.GetFieldDimension(additional_field[j]) == 5)
            Apply(Model.A5(additional_field[j]), additional_field_type[j],
                  additional_field_perturbation[j]);

        /*** Perturbs the fields that depend on the species ***/

        for (j = 0; j < int(field_list.size()); j++)
          for (s = 0; s < int(field[field_list[j]].size()); s++)
            {
              s_local = Model.GetSpeciesIndex(field_list[j],
                                              field[field_list[j]][s]);
              if (Model.GetFieldDimension(field_list[j]) == 2)
                {
                  int Nx = Model.A2(field_list[j]).extent(1);
                  Array<T, 1> array(&Model.A2(field_list[j])(s_local, 0),
                                    shape(Nx), neverDeleteData);
                  Apply(array, field_type[field_list[j]][s],
                        field_perturbation[field_list[j]][s]);
                }
              else if (Model.GetFieldDimension(field_list[j]) == 3)
                {
                  int Ny = Model.A3(field_list[j]).extent(1);
                  int Nx = Model.A3(field_list[j]).extent(2);
                  Array<T, 2> array(&Model.A3(field_list[j])(s_local, 0, 0),
                                    shape(Ny, Nx), neverDeleteData);
                  Apply(array, field_type[field_list[j]][s],
                        field_perturbation[field_list[j]][s]);
                }
              else if (Model.GetFieldDimension(field_list[j]) == 4)
                {
                  int Nz = Model.A4(field_list[j]).extent(1);
                  int Ny = Model.A4(field_list[j]).extent(2);
                  int Nx = Model.A4(field_list[j]).extent(3);
                  Array<T, 3> array(&Model.A4(field_list[j])(s_local,
                                                             0, 0, 0),
                                    shape(Nz, Ny, Nx), neverDeleteData);
                  Apply(array, field_type[field_list[j]][s],
                        field_perturbation[field_list[j]][s]);
                }
              else if (Model.GetFieldDimension(field_list[j]) == 5)
                {
                  int Nt = Model.A4(field_list[j]).extent(1);
                  int Nz = Model.A4(field_list[j]).extent(2);
                  int Ny = Model.A4(field_list[j]).extent(3);
                  int Nx = Model.A4(field_list[j]).extent(4);
                  Array<T, 4> array(&Model.A5(field_list[j])(s_local,
                                                             0, 0, 0, 0),
                                    shape(Nt, Nz, Ny, Nx), neverDeleteData);
                  Apply(array, field_type[field_list[j]][s],
                        field_perturbation[field_list[j]][s]);
                }
            }

        Model.Forward();
        if (rank == 0)
          OutputSaver.Save(Model);
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI::Finalize();
#endif
  }


  //! Applies a given operation on a multidimensional array.
  /*!
    \param A the array to be modified.
    \param type operation to apply.
    \param value number related to the operation.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  template<int N>
  void PerturbationDriver<T, ClassModel, ClassOutputSaver>
  ::Apply(Array<T, N>& A, string type, T value)
  {
    if (type == "multiply")
      A *= value;
    else if (type == "add")
      A += value;
    else
      throw "Error in PerturbationDriver::Apply: type \""
        + type + "\" is not supported.";
  }


  //! Applies a given operation on the wind angle.
  /*!
    \param type operation to apply.
    \param value number related to the operation.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void PerturbationDriver<T, ClassModel, ClassOutputSaver>
  ::WindAngle(string type, T value)
  {
    const T pi = 3.14159265358979323846264;
    int Nz = Model.A3("MeridionalWind").extent(0);
    Data<T, 3>& meridional_wind(Model.D3("MeridionalWind"));
    Data<T, 3>& zonal_wind(Model.D3("ZonalWind"));
    Grid<T>* GridYMeridional = meridional_wind[1].Duplicate();
    Grid<T>* GridXMeridional = meridional_wind[2].Duplicate();
    Grid<T>* GridYZonal = zonal_wind[1].Duplicate();
    Grid<T>* GridXZonal = zonal_wind[2].Duplicate();
    GridYMeridional->SetVariable(0);
    GridXMeridional->SetVariable(1);
    GridYZonal->SetVariable(0);
    GridXZonal->SetVariable(1);
    Data<T, 2> zonalw_interpolated(*GridYMeridional, *GridXMeridional);
    Data<T, 2> meridionalw_interpolated(*GridYZonal, *GridXZonal);
    Data<T, 2> U(*GridYZonal, *GridXZonal);
    Data<T, 2> V(*GridYMeridional, *GridXMeridional);
    int i, j;
    double angle;
    double module;

    for (int k = 0; k < Nz; k++)
      {
        U.SubData(zonal_wind, k, Range::all(), Range::all());
        V.SubData(meridional_wind, k, Range::all(), Range::all());
        LinearInterpolationRegular(U, zonalw_interpolated);
        LinearInterpolationRegular(V, meridionalw_interpolated);
        // Zonal wind.
        for (j = 0; j < zonal_wind.GetLength(1); j++)
          for (i = 0; i < zonal_wind.GetLength(2); i++)
            {
              module = sqrt(U(j, i) * U(j, i)
                            + meridionalw_interpolated(j, i)
                            * meridionalw_interpolated(j, i));
              if (type == "multiply")
                {
                  angle = acos(U(j, i) / module);
                  if (meridionalw_interpolated(j, i) < 0)
                    angle = - angle;
                  angle *= value;
                  U(j, i) = cos(angle) * module;
                }
              else if (type == "add")
                {
                  U(j, i) = U(j, i) * cos(value * pi / 180.) -
                    meridionalw_interpolated(j, i) * sin(value * pi / 180.);
                }
              else
                throw "Error in PerturbationDriver::WindAngle: type \""
                  + type + "\" is not supported.";
            }
        // Meridional wind.
        for (j = 0; j < meridional_wind.GetLength(1); j++)
          for (i = 0; i < meridional_wind.GetLength(2); i++)
            {
              module = sqrt(V(j, i) * V(j, i)
                            + zonalw_interpolated(j, i)
                            * zonalw_interpolated(j, i));
              if (type == "multiply")
                {
                  angle = asin(V(j, i) / module);
                  if (zonalw_interpolated(j, i) < 0)
                    angle = pi - angle;
                  angle *= value;
                  V(j, i) = sin(angle) * module;
                }
              else if (type == "add")
                {
                  V(j, i) = zonalw_interpolated(j, i) * sin(value * pi / 180.)
                    + V(j, i) * cos(value * pi / 180.);
                }
              else
                throw "Error in PerturbationDriver::WindAngle: type \""
                  + type + "\" is not supported.";
            }

        zonal_wind()(k, Range::all(), Range::all())
          = U()(Range::all(), Range::all());
        meridional_wind()(k, Range::all(), Range::all())
          = V()(Range::all(), Range::all());
      }
  }


  //! Applies a given operation on the wind module.
  /*!
    \param type operation to apply.
    \param value number related to the operation.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void PerturbationDriver<T, ClassModel, ClassOutputSaver>
  ::WindModule(string type, T value)
  {
    if (type == "multiply")
      {
        Apply(Model.A3("MeridionalWind"), type, value);
        Apply(Model.A3("ZonalWind"), type, value);
      }
    else if (type == "add")
      {
        int Nz = Model.A3("MeridionalWind").extent(0);
        Data<T, 3>& meridional_wind(Model.D3("MeridionalWind"));
        Data<T, 3>& zonal_wind(Model.D3("ZonalWind"));
        Grid<T>* GridYMeridional = meridional_wind[1].Duplicate();
        Grid<T>* GridXMeridional = meridional_wind[2].Duplicate();
        Grid<T>* GridYZonal = zonal_wind[1].Duplicate();
        Grid<T>* GridXZonal = zonal_wind[2].Duplicate();
        GridYMeridional->SetVariable(0);
        GridXMeridional->SetVariable(1);
        GridYZonal->SetVariable(0);
        GridXZonal->SetVariable(1);
        Data<T, 2> zonalw_interpolated(*GridYMeridional, *GridXMeridional);
        Data<T, 2> meridionalw_interpolated(*GridYZonal, *GridXZonal);
        Data<T, 2> U(*GridYZonal, *GridXZonal);
        Data<T, 2> V(*GridYMeridional, *GridXMeridional);
        int i, j;
        double module;
        double cos_angle, sin_angle;
        for (int k = 0; k < Nz; k++)
          {
            U.SubData(zonal_wind, k, Range::all(), Range::all());
            V.SubData(meridional_wind, k, Range::all(), Range::all());
            LinearInterpolationRegular(U, zonalw_interpolated);
            LinearInterpolationRegular(V, meridionalw_interpolated);
            // Zonal wind.
            for (j = 0; j < zonal_wind.GetLength(1); j++)
              for (i = 0; i < zonal_wind.GetLength(2); i++)
                {
                  module = sqrt(U(j, i) * U(j, i)
                                + meridionalw_interpolated(j, i)
                                * meridionalw_interpolated(j, i));
                  cos_angle = U(j, i) / module;
                  module += value;
                  U(j, i) = module * cos_angle;
                }
            // Meridional wind.
            for (j = 0; j < meridional_wind.GetLength(1); j++)
              for (i = 0; i < meridional_wind.GetLength(2); i++)
                {
                  module = sqrt(V(j, i) * V(j, i)
                                + zonalw_interpolated(j, i)
                                * zonalw_interpolated(j, i));
                  sin_angle = V(j, i) / module;
                  module += value;
                  V(j, i) = module * sin_angle;
                }
            zonal_wind()(k, Range::all(), Range::all())
              = U()(Range::all(), Range::all());
            meridional_wind()(k, Range::all(), Range::all())
              = V()(Range::all(), Range::all());
          }
      }
    else
      throw "Error in PerturbationDriver::WindModule: type \""
        + type + "\" is not supported.";
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PERTURBATIONDRIVER_CXX
#endif
