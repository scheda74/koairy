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


#ifndef POLYPHEMUS_FILE_DRIVER_TRAJECTORY_HXX


#include "SeldonData.hxx"


namespace Polyphemus
{


  using namespace SeldonData;


  ////////////////
  // TRAJECTORY //
  ////////////////


  /*! \brief 'Trajectory' manages concentration trajectory. The trajectory
    tape stores concentrations of all species at regular time steps.
  */
  template<class T, int N>
  class Trajectory
  {

  protected:

    //! Trajectory starting date.
    Date Date_min;
    //! Trajectory time step in seconds.
    T Delta_t;
    //! Number of time steps in the trajectory file.
    int Nt;

    //! File name storing the trajectory.
    string file_name;

    //! Is backward loading? Initially true.
    bool backward;
    //! Data stored on trajectory file at the step before \a load_date.
    Data<T, N> FileData_i;
    //! Data stored on trajectory file at the step after \a load_date.
    Data<T, N> FileData_f;
    /*! \brief Step value for FileData_i in trajectory file. Its valid range
      is from zero to \a Nt - 2. */
    int step_i;

  public:

    /*** Constructor and destructor ***/

    Trajectory();
    Trajectory(const TinyVector<int, N>& shape, Date Date_min_, T Delta_t_,
               string file_name_, bool backward_ = true);
    ~Trajectory();

    /*** Access ***/

    bool IsBackward();
    void SetBackward(bool flag);

    /*** Methods ***/

    void Init(const TinyVector<int, N>& shape, Date Date_min_, T Delta_t_,
              string file_name_, bool backward_ = true);
    void InitLoad();

    template <class T_in>
    void Append(Data<T_in, N>& data, Date date);

    template <class T_out>
    void Load(Date date, Data<T_out, N>& data, bool interpolated = false);

  protected:

    template <class T_out>
    void Interpolate(Date date_min_file, T delta_t_file, int step_i,
                     Data<T, N>& file_data_i, Data<T, N>& file_data_f,
                     Date date, Data<T_out, N>& interpolated_data);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_TRAJECTORY_HXX
#endif
