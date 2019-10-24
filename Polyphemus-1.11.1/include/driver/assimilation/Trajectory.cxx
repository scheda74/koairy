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


#ifndef POLYPHEMUS_FILE_DRIVER_TRAJECTORY_CXX


#include "Trajectory.hxx"


namespace Polyphemus
{


  //////////////////////////////
  // CONSTRUCTOR & DESTRUCTOR //
  //////////////////////////////


  //! Constructor.
  template<class T, int N>
  Trajectory<T, N>::Trajectory()
    : Delta_t(0.), Nt(0), backward(true)
  {
    file_name = "";
  }

  //! Constructor.
  /*!
    \param shape vector of lengths of the data arrays for trajectory manager.
    \param Date_min_ trajectory starting date.
    \param Delta_t_ trajectory time step in seconds.
    \param file_name_ file name storing the trajectory.
    \param backward_ flag for backward loading.
  */
  template<class T, int N>
  Trajectory<T, N>::Trajectory(const TinyVector<int, N>& shape,
                               Date Date_min_, T Delta_t_,
                               string file_name_, bool backward_):
    Date_min(Date_min_), Delta_t(Delta_t_), Nt(0)
  {
    // Set variables.
    file_name = file_name_;
    backward = backward_;

    // Empties trajectory file.
    ofstream(file_name.c_str()).close();

    // Allocates data arrays.
    FileData_i.Resize(shape);
    FileData_f.Resize(shape);
  }


  //! Destructor.
  template<class T, int N>
  Trajectory<T, N>::~Trajectory()
  {
    // Empties trajectory file.
    ofstream(file_name.c_str()).close();
  }


  ////////////
  // ACCESS //
  ////////////


  //! Is backward loading?
  /*!
    \return True for backward loading; false for forward loading.
  */
  template<class T, int N>
  bool Trajectory<T, N>::IsBackward()
  {
    return backward;
  }


  //! Sets flag for loading direction.
  /*!
    \param flag true for backward loading; false for forward loading.
  */
  template<class T, int N>
  void Trajectory<T, N>::SetBackward(bool flag)
  {
    backward = flag;
  }


  /////////////
  // METHODS //
  /////////////


  //! Initializations.
  /*!
    \param shape vector of lengths of the data arrays for trajectory manager.
    \param Date_min_ trajectory starting date.
    \param Delta_t_ trajectory time step in seconds.
    \param file_name_ file name storing the trajectory.
    \param backward_ flag for backward loading.
  */
  template<class T, int N>
  void Trajectory<T, N>::Init(const TinyVector<int, N>& shape,
                              Date Date_min_, T Delta_t_,
                              string file_name_, bool backward_)
  {
    // Set variables.
    Date_min = Date_min_;
    Delta_t = Delta_t_;
    file_name = file_name_;
    backward = backward_;

    Nt = 0;

    // Empties trajectory file.
    ofstream(file_name.c_str()).close();

    // Allocates data arrays.
    FileData_i.Resize(shape);
    FileData_f.Resize(shape);
  }


  //! It appends data at a given date to the trajectory file.
  /*!
    \param data Data instance to be appended to the trajectory file.
    \param date the date of the Data instance \a data.
    \warning the caller is responsible for the correctness of date.
  */
  template<class T, int N>
  template<class T_in>
  void Trajectory<T, N>::Append(Data<T_in, N>& data, Date date)
  {
    int step = int(T(date.GetSecondsFrom(Date_min)) / Delta_t);
    if (step != Nt)
      {
        Date close_date = Date_min;
        close_date.AddSeconds(double(Delta_t * T(Nt)));
        throw string("ERROR! Appending date \"") +
          date.GetDate("%y-%m-%d_%hh%i") + "\" too far, the closest is \""
          + close_date.GetDate("%y-%m-%d_%hh%i") + "\"";
      }
    Nt++;
    FormatBinary<T>().Append(data, file_name);
  }


  //! Initialization for loading.
  /*! It initializes data arrays surrounding the loading date. */
  template<class T, int N>
  void Trajectory<T, N>::InitLoad()
  {
    if (Nt < 2)
      throw "Not enought data in trajectory file for loading.";

    if (backward)
      step_i = Nt - 2;
    else
      step_i = 0;
    FormatBinary<T>().ReadRecord(file_name, step_i, FileData_i);
    FormatBinary<T>().ReadRecord(file_name, step_i + 1, FileData_f);
  }


  /*! \brief It loads data for a given date to the Data instance from the
    trajectory file. */
  /*!
    \param date the loading date for the Data instance \a data.
    \param data Data instance to which the data from the trajectory is loaded.
    \param interpolated flag indicating whether interpolation is performed.
  */
  template<class T, int N>
  template<class T_out>
  void Trajectory<T, N>::Load(Date date, Data<T_out, N>& data,
                              bool interpolated)
  {
    T diff = date.GetSecondsFrom(Date_min);
    int record = int(diff / Delta_t);

    // Reads boundary values.
    if (record == 0)
      {
        FormatBinary<T>().ReadRecord(file_name, 0, data);
        return;
      }
    else if (record == Nt - 1)
      {
        FormatBinary<T>().ReadRecord(file_name, Nt - 1, data);
        return;
      }
    // Out of range.
    else if (record < 0 || record > Nt - 1)
      throw string("Record index \"") + to_str(record) +
        string("\" out of range; valid range in trajectory file ") +
        file_name + string(" is from 0 to ") + to_str(Nt - 1);

    // Updates FileData_i and FileData_f when necessary.
    if (record != step_i)
      if (backward)
        {
          if (record == step_i - 1)
            {
              FileData_f.GetArray() = FileData_i.GetArray();
              FormatBinary<T>().ReadRecord(file_name, record, FileData_i);
            }
          else
            {
              FormatBinary<T>().ReadRecord(file_name, record, FileData_i);
              FormatBinary<T>().ReadRecord(file_name, record + 1, FileData_f);
            }
        }
      else // case of forward loading.
        {
          if (record == step_i + 1)
            {
              FileData_i.GetArray() = FileData_f.GetArray();
              FormatBinary<T>().ReadRecord(file_name, record + 1, FileData_f);
            }
          else
            {
              FormatBinary<T>().ReadRecord(file_name, record, FileData_i);
              FormatBinary<T>().ReadRecord(file_name, record + 1, FileData_f);
            }
        }

    // Updates step index for FileData_i.
    step_i = record;

    if (interpolated)
      Interpolate(Date_min, Delta_t, step_i,
                  FileData_i, FileData_f, date, data);
    else
      {
        for (int i = 0; i < (int) data.GetArray().size(); i++)
          {
            TinyVector<int, N> index = 0;
            T* d_i = &FileData_i.GetArray()(index);
            T_out* d = &data.GetArray()(index);
            *(d + i) = T_out(*(d_i + i));
          }
      }
  }


  //! Linear interpolation in time.
  /*! It interpolates the data from trajectory file, and copies the
    interpolation result to Data instance.
    \param date_min_file starting date of the trajectory file.
    \param Delta_t_file trajectory file time-step.
    \param step_i step index of \a file_data_i in the trajectory file.
    \param file_data_i data read from the trajectory file at record \a step_i.
    \param file_data_f data read from the trajectory file at record \a step_i
    + 1.
    \param date date at which data is interpolated.
    \param InterpolatedData (output) data linearly interpolated between
    \a FileData_i and \a FileData_f.
  */
  template<class T, int N>
  template<class T_out>
  void Trajectory<T, N>::Interpolate(Date date_min_file, T delta_t_file,
                                     int step_i, Data<T, N>& file_data_i,
                                     Data<T, N>& file_data_f, Date date,
                                     Data<T_out, N>& interpolated_data)
  {
    // Pointers to the first elements of data arrays.
    TinyVector<int, N> index = 0;
    T* d_i = &file_data_i.GetArray()(index);
    T* d_f = &file_data_f.GetArray()(index);
    T_out* d = &interpolated_data.GetArray()(index);

    T time_distance = date.GetSecondsFrom(date_min_file);
    T weight_f = (time_distance - T(step_i) * delta_t_file) / delta_t_file;
    T weight_i = 1. - weight_f;

    // Interpolations.
    for (int i = 0; i < (int) file_data_i.GetArray().size(); i++)
      *(d + i) = T_out(*(d_i + i) * weight_i + * (d_f + i) * weight_f);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_TRAJECTORY_CXX
#endif
