// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet, Mohamed Aissaoui
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


#ifndef POLYPHEMUS_FILE_PERTURBATION_PERTURBATIONMANAGER_CXX


#include "PerturbationManager.hxx"


namespace Polyphemus
{


  //! Main constructor.
  template <class T>
  PerturbationManager<T>::PerturbationManager():
    randomization(NULL)
  {
  }


  //! Destructor.
  template <class T>
  PerturbationManager<T>::~PerturbationManager()
  {
    if (randomization != NULL)
      delete randomization;
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads perturbation field configurations. A map named "Param_map" is
    established to associate the parameter vector with the name of its field.
    \param Model model with the following interface:
    <ul>
    <li> GetConfigurationFile()
    </ul>
  */
  template <class T>
  template <class ClassModel>
  void PerturbationManager<T>::ReadConfiguration(ClassModel& Model)
  {

    /*** Reads configuration file name ***/

    ConfigStream model_config(Model.GetConfigurationFile());
    string perturb_manager_file;
    model_config.SetSection("[perturbation_management]");
    model_config.PeekValue("Configuration_file", perturb_manager_file);

    ConfigStream perturb_config(perturb_manager_file);

    /*** Retrieves perturbation parameters ***/

    perturb_config.SetSection("[general]");
    perturb_config.Find("Fields");
    perturb_field_list = split(perturb_config.GetLine());
    if (perturb_field_list.size() == 1 && perturb_field_list[0] == "---")
      perturb_field_list.clear();
    // For additional fields.
    perturb_field_list.push_back("AdditionalField");

    perturb_config.PeekValue("Rand_seed", rand_seed_directory);
    perturb_config.PeekValue("Field_maximum_spread", field_maximum_spread);
    perturb_config.PeekValue("Observation_maximum_spread",
                             obs_maximum_spread);

    if (randomization != NULL)
      delete randomization;
    randomization =
      new RandGenerator<NEWRAN::LGM_mixed, float>(rand_seed_directory);
    randomization->SetMaximumSpread(field_maximum_spread);

    vector<string> param_name, param_pdf;
    vector<double> param_std;
    vector<string> param_dependence;
    vector<float> param_correlation;

    string line;
    vector<string> vtmp;

    for (unsigned long l = 0; l < perturb_field_list.size(); l++)
      {
        param_name.clear();
        param_std.clear();
        param_pdf.clear();
        param_dependence.clear();
        param_correlation.clear();

        perturb_config.SetSection(string("[") + perturb_field_list[l] + "]");

        while (!perturb_config.IsEmpty()
               && perturb_config.PeekElement()[0] != '[')
          {
            param_name.push_back(perturb_config.GetElement());
            param_std.push_back(perturb_config.GetNumber());
            param_pdf.push_back(perturb_config.GetElement());

            // At most one correlated species is possible.
            line = perturb_config.GetFullLine();
            if (!perturb_config.Discard(line))
              {
                split(line, vtmp);
                if (vtmp.size() != 2)
                  throw string("\"") + perturb_config.GetFileName() +
                    string("\" is badly formatted: correlation should be") +
                    + " specified in this way: [parameter] [correlation].";
                param_dependence.push_back(vtmp[0]);
                param_correlation.push_back(convert<float>(vtmp[1]));
              }
            else
              {
                param_dependence.push_back("");
                param_correlation.push_back(0.);
              }
          }


        // Sets parameters for each field.
        if (perturb_field_list[l] == "BoundaryCondition")
          {
            string fieldname_BC_x = perturb_field_list[l] + "_x";
            string fieldname_BC_y = perturb_field_list[l] + "_y";
            string fieldname_BC_z = perturb_field_list[l] + "_z";

            Param_map[fieldname_BC_x].resize(param_name.size());
            Param_map[fieldname_BC_y].resize(param_name.size());
            Param_map[fieldname_BC_z].resize(param_name.size());
            for (unsigned long p = 0; p < param_name.size(); p++)
              {

                Param_map[fieldname_BC_x][p]
                  = PolairParam<float>(fieldname_BC_x,
                                       param_name[p],
                                       param_std[p],
                                       param_pdf[p]);
                Param_map[fieldname_BC_y][p]
                  = PolairParam<float>(fieldname_BC_y,
                                       param_name[p],
                                       param_std[p],
                                       param_pdf[p]);
                Param_map[fieldname_BC_z][p]
                  = PolairParam<float>(fieldname_BC_z,
                                       param_name[p],
                                       param_std[p],
                                       param_pdf[p]);

                // Boundary conditions are correlated.
                Param_map[fieldname_BC_y][p].GetCorrelatedParam()
                  .push_back(&Param_map[fieldname_BC_x][p]);
                Param_map[fieldname_BC_y][p].GetCorrelation().push_back(1.);
                Param_map[fieldname_BC_z][p].GetCorrelatedParam()
                  .push_back(&Param_map[fieldname_BC_x][p]);
                Param_map[fieldname_BC_z][p].GetCorrelation().push_back(1.);
              }
          }
        else
          {
            Param_map[perturb_field_list[l]].resize(param_name.size());
            for (unsigned long p = 0; p < param_name.size(); p++)
              {
                Param_map[perturb_field_list[l]][p]
                  = PolairParam<float>(perturb_field_list[l],
                                       param_name[p],
                                       param_std[p],
                                       param_pdf[p]);

                // At most one correlated species is possible.
                if (!param_dependence[p].empty())
                  {
                    unsigned long j = 0;
                    while (j < p &&
                           param_dependence[p] !=
                           Param_map[perturb_field_list[l]][j].GetName())
                      j++;
                    if (j == p)
                      throw string("\"") + perturb_config.GetFileName()
                        + string("\" is badly formatted!")
                        + string(" Unable to find \"") + param_dependence[p]
                        + string("\" before \"") + param_name[p] + "\"";
                    Param_map[perturb_field_list[l]][p].GetCorrelatedParam()
                      .push_back(&Param_map[perturb_field_list[l]][j]);
                    Param_map[perturb_field_list[l]][p].GetCorrelation()
                      .push_back(param_correlation[p]);
                  }
              }
          }
      }
  }


  ////////////////////
  // INITIALIZATION //
  ////////////////////


  //! Initialization of observation manager.
  /*! It reads perturbed field configurations. It generates a name list of the
    fields to be perturbed, in which "BoundaryCondition" field is expanded
    with respect to three spatial coordinates. "AdditionalField" field is also
    replaced by its parameter names to generate a complete name list of all
    field data to be perturbed.
    \param Model model with the following interface:
    <ul>
    <li> GetConfigurationFile()
    <li> HasField(string)
    </ul>
  */
  template <class T>
  template <class ClassModel>
  void PerturbationManager<T>::Init(ClassModel& Model)
  {
    // Reads configurations.
    ReadConfiguration(Model);

    // Replaces "BoundaryCondition" in field list with "BoundaryCondition_x",
    // "BoundaryCondition_y", "BoundaryCondition_z".
    vector<string> tmp;
    for (unsigned long i = 0; i < perturb_field_list.size(); i++)
      {
        if (perturb_field_list[i] == "BoundaryCondition")
          {
            tmp.push_back("BoundaryCondition_x");
            tmp.push_back("BoundaryCondition_y");
            tmp.push_back("BoundaryCondition_z");
          }
        else
          tmp.push_back(perturb_field_list[i]);
      }
    perturb_field_list = tmp;

    // Replaces "AdditionalField" with its parameter names, for instance
    // "Attenuation" or "VerticalDiffusionCoefficient".
    for (unsigned long i = 0; i < perturb_field_list.size(); i++)
      {
        if (perturb_field_list[i] == "AdditionalField")
          for (unsigned long j = 0;
               j < Param_map["AdditionalField"].size();
               j++)
            {
              string str = Param_map["AdditionalField"].at(j).GetName();
              perturb_data_list.push_back(str);
            }
        else
          perturb_data_list.push_back(perturb_field_list[i]);
      }

    // Checks perturbed data names.
    for (unsigned long i = 0; i < perturb_data_list.size(); i++)
      if (perturb_data_list[i] == "WindAngle"
          && (!Model.HasField("MeridionalWind")
              || !Model.HasField("ZonalWind")))
        throw string("Field \"WindAngle\" cannot be perturbed because")
          + " \"MeridionalWind\" or \"ZonalWind\" or is not part of the"
          " model interface.";
      else if (!Model.HasField(perturb_data_list[i])
               && perturb_data_list[i] != "WindAngle")
        throw string("Field \"") + perturb_data_list[i]
          + string("\" to be perturbed is not part of the model interface.");
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Returns the map from field names to their parameter lists.
  /*!
    \return Map from field names to their parameter lists.
  */
  template <class T>
  map<string, vector<PolairParam<float> > >&
  PerturbationManager<T>::GetPerturbationMap()
  {
    return Param_map;
  }


  //! Return the name list of complete field data to be perturbed.
  /*!
    \return Name list of complete field data to be perturbed.
  */
  template <class T>
  vector<string>& PerturbationManager<T>::GetPerturbDataList()
  {
    return perturb_data_list;
  }


  //! Return the name list of the fields to be perturbed.
  /*!
    \return Name list of the fields to be perturbed.
  */
  template <class T>
  vector<string>& PerturbationManager<T>::GetPerturbFieldList()
  {
    return perturb_field_list;
  }


  //! Returns the field maximum spread.
  /*! Every perturbed value of the field cannot exceed the original value of
    the field plus or minus 'maximum_spread' times the standard deviation.
    \return The field maximum spread.
  */
  template <class T>
  T PerturbationManager<T>::GetFieldMaximumSpread()
  {
    return field_maximum_spread;
  }


  //! Return the observation maximum spread.
  /*! Every perturbed value cannot exceed the original value of the
    obsvervation plus or minus 'maximum_spread' times the standard
    deviation.
    \return The observation maximum spread.
  */
  template <class T>
  T PerturbationManager<T>::GetObsMaximumSpread()
  {
    return obs_maximum_spread;
  }


  //////////////////////////
  // PERTURBATION METHODS //
  //////////////////////////


  //! Generation of perturbation random numbers.
  /*! It generates perturbation random numbers for each parameters according
    to the parameter configurations, e.g. probability distribution function
    and standard derivation. Then a map is established to map the random
    number vectors for the parameters to their field name.
    \return Map from field name to its vector of random numbers.
  */
  template <class T>
  map<string, vector<T> > PerturbationManager<T>::GeneratePerturbNumMap()
  {
    unsigned long l, p;
    double pert;

    map<string, vector<T> > PerturbNum_map;
    for (l = 0; l < perturb_field_list.size(); l++)
      for (p = 0; p < Param_map[perturb_field_list[l]].size(); p++)
        {
          PolairParam<float>& param = Param_map[perturb_field_list[l]][p];
          pert = randomization->GenRandomNumber(param);
          PerturbNum_map[perturb_field_list[l]].push_back(T(pert));
        }

    return PerturbNum_map;
  }


  //! Generates random column vectors with given covariance.
  /*! It generates random vectors with given covariance and with 0-mean
    Gaussian distribution.
    \param Matrix (input/output) on entry, the matrix with its row number as
    the dimension of the random vectors to be generated, and with its column
    number as the total number of random vectors to be generated; on exit,
    each matrix column 'j' is filled with random numbers sampled from a
    Gaussian distribution with covariance \a Covariance.
    \param Covariance covariance matrix with dimensions as \a Matrix row
    number times \a Matrix row number.
    \param maximum_spread number which defines the maximum spread, that is, \a
    maximum_spread times the corresponding standard deviation. No random
    number can exceed the maximum spread.
    \warning Works only for diagonal covariance matrix \a Covariance
  */
  template <class T>
  void
  PerturbationManager<T>
  ::GenerateRandomColumnVector(Array<T, 2>& Covariance, Array<T, 2>& Matrix,
                               T maximum_spread)
  {
    if (randomization == NULL)
      throw "Error: perturbation manager is not initialized.";

    // Gaussian generator with 0 mean and unitary variance.
    NEWRAN::Normal N1;

    int m = Matrix.rows();
    int n = Matrix.cols();
    T random_number;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        {
          while (abs(double(random_number = N1.Next()))
                 > double(maximum_spread)) ;
          Matrix(i, j) = random_number * sqrt(Covariance(i, i));
        }
  }


  //! Generates random row vectors with given covariance.
  /*! It generates random vectors with given covariance and with 0-mean
    Gaussian distribution.
    \param Matrix (input/output) on entry, the matrix with its column number
    as the dimension of the random vectors to be generated, and with its row
    number as the total number of random vectors to be generated; on exit,
    each matrix row 'j' is filled with random numbers sampled from a
    Gaussian distribution with covariance \a Covariance.
    \param Covariance covariance matrix with dimensions as \a Matrix column
    number times \a Matrix column number.
    \param maximum_spread number which defines the maximum spread, that is, \a
    maximum_spread times the corresponding standard deviation. No random
    number can exceed the maximum spread.
    \warning Works only for diagonal covariance matrix \a Covariance.
  */
  template <class T>
  void
  PerturbationManager<T>
  ::GenerateRandomRowVector(Array<T, 2>& Covariance, Array<T, 2>& Matrix,
                            T maximum_spread)
  {
    if (randomization == NULL)
      throw "Error: perturbation manager is not initialized.";

    // Gaussian generator with 0 mean and unitary variance.
    NEWRAN::Normal N1;

    int m = Matrix.rows();
    int n = Matrix.cols();
    T random_number;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        {
          while (abs(double(random_number = N1.Next()))
                 > double(maximum_spread)) ;
          Matrix(j, i) = random_number * sqrt(Covariance(i, i));
        }
  }


  //! Array perturbations.
  /*! It perturbs array data with a given random number and a given
    probability density function.
    \param array the array to be perturbed.
    \param random_number the value of the given random number for the
    perturbation.
    \param pdf the name of the given probability density function.
  */
  template <class T>
  template <int N>
  void
  PerturbationManager<T>::PerturbArray(Array<T, N>& array,
                                       double random_number, string pdf)
  {
    int Nx, Ny, Nz, Ns, Np, i, j, k, s, p;
    int Ndim = array.dimensions();

    if (Ndim == 1)
      {
        Nx = array.extent(0);
        for (i = 0; i < Nx; i++)
          if (pdf == "LN")
            array(i) = array(i) * random_number;
          else if (pdf == "N")
            array(i) *= 1. + random_number;
      }
    else if (Ndim == 2)
      {
        Ny = array.extent(0);
        Nx = array.extent(1);
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            if (pdf == "LN")
              array(j, i) = array(j, i) * random_number;
            else if (pdf == "N")
              array(j, i) *= 1. + random_number;
      }
    else if (Ndim == 3)
      {
        Nz = array.extent(0);
        Ny = array.extent(1);
        Nx = array.extent(2);
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              if (pdf == "LN")
                array(k, j, i) = array(k, j, i) * random_number;
              else if (pdf == "N")
                array(k, j, i) *= 1. + random_number;
      }
    else if (Ndim == 4)
      {
        Ns = array.extent(0);
        Nz = array.extent(1);
        Ny = array.extent(2);
        Nx = array.extent(3);
        for (s = 0; s < Ns; s++)
          for (k = 0; k < Nz; k++)
            for (j = 0; j < Ny; j++)
              for (i = 0; i < Nx; i++)
                if (pdf == "LN")
                  array(s, k, j, i) = array(s, k, j, i) * random_number;
                else if (pdf == "N")
                  array(s, k, j, i) *= 1. + random_number;
      }
    else if (Ndim == 5)
      {
        Np = array.extent(0);
        Ns = array.extent(1);
        Nz = array.extent(2);
        Ny = array.extent(3);
        Nx = array.extent(4);
        for (p = 0; p < Np; p++)
          for (s = 0; s < Ns; s++)
            for (k = 0; k < Nz; k++)
              for (j = 0; j < Ny; j++)
                for (i = 0; i < Nx; i++)
                  if (pdf == "LN")
                    array(p, s, k, j, i)
                      = array(p, s, k, j, i) * random_number;
                  else if (pdf == "N")
                    array(p, s, k, j, i) *= 1. + random_number;
      }
  }


  //! Field perturbations.
  /*! It generates perturbed fields according to the configurations and
    random numbers associated with the parameters.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex(string, string)
    <li> GetFieldDimension(string)
    <li> A2(string)
    <li> A3(string)
    <li> A4(string)
    <li> A5(string)
    </ul>
    \param PerturbNum_map the map that maps field name to its parameter
    perturbation numbers.
  */
  template <class T>
  template <class ClassModel>
  void
  PerturbationManager<T>::GenerateField(ClassModel& Model,
                                        map<string, vector<T> >
                                        PerturbNum_map)
  {
    unsigned long l, p;
    // Are there two arrays per step ("_i" and "_f")?
    bool double_step;

    for (l = 0; l < perturb_field_list.size(); l++)
      {
        bool additional_field = perturb_field_list[l] == "AdditionalField";

        for (p = 0; p < Param_map[perturb_field_list[l]].size(); p++)
          {
            PolairParam<float>& param = Param_map[perturb_field_list[l]][p];

            string fieldname;
            if (additional_field)
              fieldname = param.GetName();
            else
              fieldname = param.GetFieldName();

            double_step = Model.HasField(fieldname + "_i")
              && Model.HasField(fieldname + "_f");

            T pert = PerturbNum_map[perturb_field_list[l]][p];

            if (fieldname == "WindAngle") // Special case because it requires
              // the perturbation of the wind
              // components.
              {
                Array<T, 3> U = Model.A3("ZonalWind");
                Array<T, 3> V = Model.A3("MeridionalWind");
                T zonal_wind;
                // 'pert' should be in degrees.
                T angle = pert * 3.14159265358979323846264 / 180.;
                // Applies the rotation.
                for (int i = 0; i < Model.A3("ZonalWind").extent(0); i++)
                  for (int j = 0; j < Model.A3("ZonalWind").extent(1); j++)
                    for (int k = 0; k < Model.A3("ZonalWind").extent(2); k++)
                      {
                        zonal_wind = U(i, j , k);
                        U(i, j, k) = U(i, j, k) * cos(angle)
                          - V(i, j, k) * sin(angle);
                        V(i, j, k) = zonal_wind * sin(angle)
                          + V(i, j, k) * cos(angle);
                      }
              }

            else if (Model.GetFieldDimension(fieldname) == 2)
              {
                int nx = Model.A2(fieldname).extent(1);

                // the parameter name is an additional field name.
                if (additional_field)
                  if (double_step)
                    {
                      PerturbArray<2>(Model.A2(fieldname + "_i"), pert,
                                      param.GetPdf());
                      PerturbArray<2>(Model.A2(fieldname + "_f"), pert,
                                      param.GetPdf());
                    }
                  else
                    PerturbArray<2>(Model.A2(fieldname), pert,
                                    param.GetPdf());
                // the parameter name is a species name.
                else if (double_step)
                  {
                    int index = Model.GetSpeciesIndex(fieldname + "_i",
                                                      param.GetName());
                    Array<T, 1> array1Di(&Model.A2(fieldname + "_i")
                                         (index, 0),
                                         shape(nx), neverDeleteData);
                    PerturbArray<1>(array1Di, pert, param.GetPdf());
                    index = Model.GetSpeciesIndex(fieldname + "_f",
                                                  param.GetName());
                    Array<T, 1> array1Df(&Model.A2(fieldname + "_f")
                                         (index, 0),
                                         shape(nx), neverDeleteData);
                    PerturbArray<1>(array1Df, pert, param.GetPdf());
                  }
                else
                  {
                    int index = Model.GetSpeciesIndex(fieldname,
                                                      param.GetName());
                    Array<T, 1> array1D(&Model.A2(fieldname)(index, 0),
                                        shape(nx), neverDeleteData);
                    PerturbArray<1>(array1D, pert, param.GetPdf());
                  }
              }

            else if (Model.GetFieldDimension(fieldname) == 3)
              {
                int ny = Model.A3(fieldname).extent(1);
                int nx = Model.A3(fieldname).extent(2);

                // the parameter name is an additional field name.
                if (additional_field)
                  if (double_step)
                    {
                      PerturbArray<3>(Model.A3(fieldname + "_i"), pert,
                                      param.GetPdf());
                      PerturbArray<3>(Model.A3(fieldname + "_f"), pert,
                                      param.GetPdf());
                    }
                  else
                    PerturbArray<3>(Model.A3(fieldname), pert,
                                    param.GetPdf());
                // the parameter name is a species name.
                else if (double_step)
                  {
                    int index = Model.GetSpeciesIndex(fieldname + "_i",
                                                      param.GetName());
                    Array<T, 2> array2Di(&Model.A3(fieldname + "_i")
                                         (index, 0, 0),
                                         shape(ny, nx), neverDeleteData);
                    PerturbArray<2>(array2Di, pert, param.GetPdf());
                    index = Model.GetSpeciesIndex(fieldname + "_f",
                                                  param.GetName());
                    Array<T, 2> array2Df(&Model.A3(fieldname + "_f")
                                         (index, 0, 0),
                                         shape(ny, nx), neverDeleteData);
                    PerturbArray<2>(array2Df, pert, param.GetPdf());
                  }
                else
                  {
                    int index = Model.GetSpeciesIndex(fieldname,
                                                      param.GetName());
                    Array<T, 2> array2D(&Model.A3(fieldname)(index, 0, 0),
                                        shape(ny, nx), neverDeleteData);
                    PerturbArray<2>(array2D, pert, param.GetPdf());
                  }
              }

            else if (Model.GetFieldDimension(fieldname) == 4)
              {
                int nz = Model.A4(fieldname).extent(1);
                int ny = Model.A4(fieldname).extent(2);
                int nx = Model.A4(fieldname).extent(3);

                // the parameter name is an additional field name.
                if (additional_field)
                  if (double_step)
                    {
                      PerturbArray<4>(Model.A4(fieldname + "_i"), pert,
                                      param.GetPdf());
                      PerturbArray<4>(Model.A4(fieldname + "_f"), pert,
                                      param.GetPdf());
                    }
                  else
                    PerturbArray<4>(Model.A4(fieldname), pert,
                                    param.GetPdf());
                // the parameter name is a species name.
                else if (double_step)
                  {
                    int index = Model.GetSpeciesIndex(fieldname + "_i",
                                                      param.GetName());
                    Array<T, 3> array3Di(&Model.A4(fieldname + "_i")
                                         (index, 0, 0, 0),
                                         shape(nz, ny, nx),
                                         neverDeleteData);
                    PerturbArray<3>(array3Di, pert, param.GetPdf());
                    index = Model.GetSpeciesIndex(fieldname + "_f",
                                                  param.GetName());
                    Array<T, 3> array3Df(&Model.A4(fieldname + "_f")
                                         (index, 0, 0, 0),
                                         shape(nz, ny, nx),
                                         neverDeleteData);
                    PerturbArray<3>(array3Df, pert, param.GetPdf());
                  }
                else
                  {
                    int index = Model.GetSpeciesIndex(fieldname,
                                                      param.GetName());
                    Array<T, 3> array3D(&Model.A4(fieldname)
                                        (index, 0, 0, 0),
                                        shape(nz, ny, nx), neverDeleteData);
                    PerturbArray<3>(array3D, pert, param.GetPdf());
                  }
              }

            else if (Model.GetFieldDimension(fieldname) == 5)
              {
                int nt = Model.A5(fieldname).extent(1);
                int nz = Model.A5(fieldname).extent(2);
                int ny = Model.A5(fieldname).extent(3);
                int nx = Model.A5(fieldname).extent(4);

                // the parameter name is an additional field name.
                if (additional_field)
                  if (double_step)
                    {
                      PerturbArray<5>(Model.A5(fieldname + "_i"), pert,
                                      param.GetPdf());
                      PerturbArray<5>(Model.A5(fieldname + "_f"), pert,
                                      param.GetPdf());
                    }
                  else
                    PerturbArray<5>(Model.A5(fieldname), pert,
                                    param.GetPdf());
                // the parameter name is a species name.
                else if (double_step)
                  {
                    int index = Model.GetSpeciesIndex(fieldname + "_i",
                                                      param.GetName());
                    Array<T, 4>
                      array4Di(&Model.A5(fieldname + "_i")
                               (index, 0, 0, 0, 0),
                               shape(nt, nz, ny, nx), neverDeleteData);
                    PerturbArray<4>(array4Di, pert, param.GetPdf());
                    index = Model.GetSpeciesIndex(fieldname + "_f",
                                                  param.GetName());
                    Array<T, 4>
                      array4Df(&Model.A5(fieldname + "_f")
                               (index, 0, 0, 0, 0),
                               shape(nt, nz, ny, nx), neverDeleteData);
                    PerturbArray<4>(array4Df, pert, param.GetPdf());
                  }
                else
                  {
                    int index = Model.GetSpeciesIndex(fieldname,
                                                      param.GetName());
                    Array<T, 4>
                      array4D(&Model.A5(fieldname)(index, 0, 0, 0, 0),
                              shape(nt, nz, ny, nx), neverDeleteData);
                    PerturbArray<4>(array4D, pert, param.GetPdf());
                  }
              }

          } // end for: parameter in field.
      } // end for: field.
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_PERTURBATIONMANAGER_CXX
#endif
