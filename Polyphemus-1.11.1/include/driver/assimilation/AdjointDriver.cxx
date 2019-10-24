// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Lin Wu, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_DRIVER_ADJOINTDRIVER_CXX


#include "AdjointDriver.hxx"
#include "newran.h"


namespace Polyphemus
{


  //! Constructor.
  /*!
    \param config_file configuration file.
  */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  AdjointDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::AdjointDriver(string config_file):
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
    (config_file)
  {
  }


  //! Destructor.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  AdjointDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::~AdjointDriver()
  {
  }


  //! Reads configurations.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AdjointDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::ReadConfiguration()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::ReadConfiguration();

    // Reads "adjoint" section.
    this->config.SetSection("[adjoint]");

    this->config.PeekValue("Display_precision", display_precision);
    this->config.PeekValue("With_trajectory_management",
                           with_trajectory_management);
    if (with_trajectory_management)
      {
        this->config.PeekValue("Trajectory_delta_t", Trajectory_delta_t);
        this->config.PeekValue("Trajectory_file", Trajectory_file_name);
      }

    this->config.PeekValue("Point_species_name", point_species_name);
    this->config.PeekValue("Point_nx", point_nx);
    this->config.PeekValue("Point_ny", point_ny);
    this->config.PeekValue("Point_nz", point_nz);

    this->config.PeekValue("Norm_perturbation_vector", norm_scale);
    this->config.PeekValue("With_random_perturbation",
                           with_random_perturbation);

    this->config.PeekValue("Decreasing_root", decreasing_root);
    this->config.PeekValue("Start_index", start_index);
    if (start_index < 0) start_index = 0;
    this->config.PeekValue("End_index", end_index);

    this->config.PeekValue("With_left_finite_difference_checking",
                           with_left_finite_difference_checking);
    this->config.PeekValue("Display_sensitivity", display_sensitivity);
  }


  //! Driver initialization.
  /*! It reads configurations */
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AdjointDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Init()
  {
    AssimilationDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
      ::Init();
  }


  //! Checks the sensitivity calculation using adjoint model.
  template < class T, class ClassModel, class ClassOutputSaver,
             class ClassObsManager >
  void AdjointDriver<T, ClassModel, ClassOutputSaver, ClassObsManager>
  ::Run()
  {
    int t, s, z, y, x;


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    this->Model.Init();
    Init();
    if (with_trajectory_management)
      TrajManager.Init(this->Model.GetConcentration().GetArray().shape(),
                       this->Model.GetDate_min(),
                       double(Trajectory_delta_t),
                       Trajectory_file_name,
                       true);
    else
      TrajData.Resize(this->Model.GetNt(), this->Model.GetNs(),
                      this->Model.GetNz(), this->Model.GetNy(),
                      this->Model.GetNx());

    // Sets output precision.
    cout.precision(streamsize(display_precision));


    //////////////////
    // Forward loop //
    //////////////////


    /*** Integrations ***/

    cout << "  Forward loop ..." << endl;
    for (t = 0; t < this->Model.GetNt(); t++)
      {
        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        if (with_trajectory_management)
          TrajManager.Append(this->Model.GetConcentration(),
                             this->Model.GetCurrentDate());
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                             shape(this->Model.GetNs(),
                                   this->Model.GetNz(),
                                   this->Model.GetNy(),
                                   this->Model.GetNx()),
                             neverDeleteData);

            Traj = this->Model.GetConcentration().GetArray();
          }


        this->Model.InitStep();
        this->Model.Forward();
      }

    if (this->option_display["show_date"])
      cout << "Current date: " <<
        this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

    // Stores reference sensitivity value.
    int point_ns = this->Model.GetSpeciesIndex(point_species_name);
    T sens = this->Model.GetConcentration()
      .GetArray()(point_ns, point_nz, point_ny, point_nx);


    ///////////////////
    // Backward loop //
    ///////////////////


    if (with_trajectory_management)
      TrajManager.InitLoad();

    this->Model.SetBackward(true);
    this->Model.GetConcentration_ccl().GetArray() = T(0.);
    this->Model.GetConcentration_ccl().GetArray()
      (point_ns, point_nz, point_ny, point_nx) = T(1.);

    Date cur_date = this->Model.GetDate_min();
    cur_date.AddSeconds((this->Model.GetNt() - 1) * this->Model.GetDelta_t());

    cout << "  Backward loop ..." << endl;
    for (t = this->Model.GetNt() - 1; t >= 0; t--)
      {
        if (with_trajectory_management)
          TrajManager.Load(cur_date, this->Model.GetConcentration(), true);
        else
          {
            Array<T, 4> Traj(&TrajData.GetArray()(t, 0, 0, 0, 0),
                             shape(this->Model.GetNs(),
                                   this->Model.GetNz(),
                                   this->Model.GetNy(),
                                   this->Model.GetNx()),
                             neverDeleteData);

            this->Model.GetConcentration().GetArray() = Traj;
          }

        this->Model.SetDate(cur_date);
        this->Model.InitStep();

        if (this->option_display["show_date"])
          cout << "Current date: " <<
            this->Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;

        this->Model.Backward();

        this->Model.SubtractTime(2. * this->Model.GetDelta_t());
        cur_date = this->Model.GetCurrentDate();
      }

    // Saves sensitivities of model output at given point with respect to
    // initial conditions.
    Data<T, 4> Sensitivity(this->Model.GetNs(), this->Model.GetNz(),
                           this->Model.GetNy(), this->Model.GetNx());
    Sensitivity.GetArray() = this->Model.GetConcentration_ccl().GetArray();


    //////////////////////////
    // Checks sensitivities //
    //////////////////////////


    Data<T, 4> conc_ic(this->Model.GetNs(), this->Model.GetNz(),
                       this->Model.GetNy(), this->Model.GetNx());

    if (with_trajectory_management)
      TrajManager.Load(this->Model.GetDate_min(), conc_ic, true);
    else
      {
        Array<T, 4> Traj(&TrajData.GetArray()(0, 0, 0, 0, 0),
                         shape(this->Model.GetNs(), this->Model.GetNz(),
                               this->Model.GetNy(), this->Model.GetNx()),
                         neverDeleteData);

        conc_ic.GetArray() = Traj;
      }


    /*** Generates initial perturbation vector ***/

    Data<T, 4> delta_h(this->Model.GetNs(), this->Model.GetNz(),
                       this->Model.GetNy(), this->Model.GetNx());

    srand((unsigned long)time(0));
    double seed = rand() / double(RAND_MAX);
    NEWRAN::LGM_mixed* urng = new NEWRAN::LGM_mixed(seed);
    NEWRAN::Random::Set(*urng);

    // Gaussian generator with 0 mean and unitary variance.
    NEWRAN::Normal N1;

    if (with_random_perturbation)
      for (s = 0; s < this->Model.GetNs(); s++)
        for (z = 0; z < this->Model.GetNz(); z++)
          for (y = 0; y < this->Model.GetNy(); y++)
            for (x = 0; x < this->Model.GetNx(); x++)
              delta_h(s, z, y, x) = N1.Next()
                * conc_ic.GetArray()(s, z, y, x);
    else
      delta_h.GetArray() = conc_ic.GetArray();

    delta_h.GetArray() = norm_scale / delta_h.Norm2() * delta_h.GetArray();

    // Temporary variables.
    Data<T, 4> pert_h(this->Model.GetNs(), this->Model.GetNz(),
                      this->Model.GetNy(), this->Model.GetNx());
    Data<T, 4> conc_pert(this->Model.GetNs(), this->Model.GetNz(),
                         this->Model.GetNy(), this->Model.GetNx());
    Array<T, 1> ratio_right(end_index - start_index);
    Array<T, 1> sens_right(end_index - start_index);
    Array<T, 1> ratio_left(end_index - start_index);
    Array<T, 1> sens_left(end_index - start_index);

    // Forward mode.
    this->Model.SetBackward(false);

    cout << "  Checks sensitivities ..." << endl;
    for (int i = start_index; i < end_index; i++)
      {

        // Decreases perturbations.
        T alpha = pow(double(decreasing_root), double(-i));
        pert_h.GetArray() = alpha * delta_h.GetArray();

        /*** right-side finite difference ***/

        Date cur_date = this->Model.GetDate_min();
        this->Model.SetDate(cur_date);

        // Calculates perturbed output at selected point.
        conc_pert.GetArray() = pert_h.GetArray() + conc_ic.GetArray();
        this->Model.GetConcentration().GetArray() = conc_pert.GetArray();

        for (t = 0; t < this->Model.GetNt(); t++)
          {
            this->Model.InitStep();
            this->Model.Forward();
          }

        sens_right(i - start_index)
          = this->Model.GetConcentration()(point_ns, point_nz,
                                           point_ny, point_nx);

        // Calculates inner product of sensitivity vector and initial
        // perturbation vector.
        T tmp = T(0.);
        for (s = 0; s < this->Model.GetNs(); s++)
          for (z = 0; z < this->Model.GetNz(); z++)
            for (y = 0; y < this->Model.GetNy(); y++)
              for (x = 0; x < this->Model.GetNx(); x++)
                tmp += Sensitivity(s, z, y, x) * delta_h(s, z, y, x);

        // Ratio of the results obtained by finite difference and by adjoint
        // model.
        ratio_right(i - start_index) = (sens_right(i - start_index) - sens)
          / (alpha * tmp);

        /*** left-side finite difference ***/

        if (with_left_finite_difference_checking)
          {
            cur_date = this->Model.GetDate_min();
            this->Model.SetDate(cur_date);
            conc_pert.GetArray() = conc_ic.GetArray() - pert_h.GetArray();
            this->Model.GetConcentration().GetArray() = conc_pert.GetArray();

            for (t = 0; t < this->Model.GetNt(); t++)
              {
                this->Model.InitStep();
                this->Model.Forward();
              }

            sens_left(i - start_index)
              = this->Model.GetConcentration()(point_ns, point_nz,
                                               point_ny, point_nx);

            ratio_left(i - start_index) = (sens - sens_left(i - start_index))
              / (alpha * tmp);
          }
      }

    /*** Display checking results ***/

    cout << "  Sensitivity checking results : " << endl;

    if (display_sensitivity)
      {
        cout << "  Right-side sensitivity values :" << endl;
        for (int i = 0; i < end_index - start_index; i++)
          cout << sens_right(i) << endl;
        if (with_left_finite_difference_checking)
          {
            cout << "  Left-side sensitivity values :" << endl;
            for (int i = 0; i < end_index - start_index; i++)
              cout << sens_left(i) << endl;
          }
      }

    cout << "  Right-side ratio values :" << endl;
    for (int i = 0; i < end_index - start_index; i++)
      cout << ratio_right(i) << endl;

    if (with_left_finite_difference_checking)
      {
        cout << "  Left-side ratio values :" << endl;
        for (int i = 0; i < end_index - start_index; i++)
          cout << ratio_left(i) << endl;
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_ADJOINTDRIVER_CXX
#endif
