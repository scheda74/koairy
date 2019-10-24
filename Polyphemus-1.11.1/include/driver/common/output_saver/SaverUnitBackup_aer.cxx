// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Edouard Debry
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITBACKUP_AER_CXX


#include "SaverUnitBackup_aer.hxx"

#include "BaseSaverUnit.cxx"


namespace Polyphemus
{


  /////////////////////////
  // SAVERUNITBACKUP_AER //
  /////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitBackup_aer<T, ClassModel>
  ::SaverUnitBackup_aer(): BaseSaverUnit<T, ClassModel>()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitBackup_aer<T, ClassModel>::~SaverUnitBackup_aer()
  {
  }


  //! Type of saver.
  /*!
    \return The string "backup_aer".
  */
  template<class T, class ClassModel>
  string SaverUnitBackup_aer<T, ClassModel>::GetType()  const
  {
    return "backup_aer";
  }


  //! First initialization.
  /*!
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetNx()
    <li> GetNy()
    <li> GetNz()
    <li> GetNs_aer()
    <li> GetNbin_aer()
    <li> GetSpeciesList_aer()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitBackup_aer<T, ClassModel>::Init(ConfigStream& config_stream,
                                                ClassModel& Model)
  {
    config_stream.PeekValue("Type", type);

    string element;
    // Interval length.
    config_stream.PeekValue("Interval_length", element);
    if (is_num(element))
      this->interval_length = to_num<int>(element);
    else
      throw string("Value of \"Interval_length\" must be a")
        + string(" number of iterations.");

    // Output filename.
    string filename = config_stream.GetValue("Output_file");

    vector<string> list_aer = Model.GetSpeciesList_aer();
    Ns_aer = Model.GetNs_aer();
    Nbin_aer = Model.GetNbin_aer();

    string file;
    for (int s = 0; s < Ns_aer; s++)
      {
        output_file.push_back(vector<string>());
        for (int b = 0; b < Nbin_aer; b++)
          {
            file = find_replace(filename, "&f", list_aer[s]);
            file = find_replace(file, "&n", to_str(b));
            output_file[s].push_back(file);
          }
      }

    // Date file.
    date_file = config_stream.GetValue("Date_file");

    // Dimensions of the underlying model.
    this->base_Nx = Model.GetNx();
    this->base_Ny = Model.GetNy();
    this->base_Nz = Model.GetNz();

    this->counter = 0;
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitBackup_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration_aer()
    <li> GetCurrentDate()
    <li> GetCurrentTime()
    <li> GetDelta_t()
    <li> GetNt()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitBackup_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int this_iteration = (int)((Model.GetCurrentTime() * 1.0001)
                               / Model.GetDelta_t());
    int last_iteration = Model.GetNt();
    Date this_date = Model.GetCurrentDate();
    int s, b;
    string file, file_buf, date_file_buf;
    ifstream date_stream_in;
    ofstream date_stream_out, date_stream_buf;


    // Update buffer BEFORE.
    if (this->interval_length == this->counter + 2 &&
        this_iteration > this->interval_length)
      {
        date_file_buf = date_file + ".buf";

        // Remove previous date_file.
        date_stream_buf.open(date_file_buf.c_str(), ios::trunc);
        date_stream_buf << "!! BUFFER SAVING NOT FINISHED !!" << endl;
        date_stream_buf.close();

        for (s = 0; s < Ns_aer; s++)
          for (b = 0; b < Nbin_aer; b++)
            {
              file = output_file[s][b];
              file_buf = output_file[s][b] + ".buf";
              Data<T, 3> Concentration_tmp;
              Concentration_tmp.Resize(this->base_Nz, this->base_Ny,
                                       this->base_Nx);
              FormatBinary<float>().Read(file, Concentration_tmp);
              FormatBinary<float>().Write(Concentration_tmp, file_buf);
            }

        // Update the date_file.
        date_stream_in.open(date_file.c_str(), ios::in);
        date_stream_buf.open(date_file_buf.c_str(), ios::trunc);
        string tmp;
        while (getline(date_stream_in, tmp))
          date_stream_buf << tmp << endl;
        date_stream_in.close();
        date_stream_buf.close();
      }


    // Backup saving.
    if (this->counter % this->interval_length == 0)
      {
        // Remove previous date_file.
        date_stream_out.open(date_file.c_str(), ios::trunc);
        date_stream_out << "!! BACKUP SAVING NOT FINISHED !!" << endl;
        date_stream_out.close();

        // Save all species at given time.
        for (s = 0; s < Ns_aer; s++)
          for (b = 0; b < Nbin_aer; b++)
            {
              file = output_file[s][b];
              Data<T, 3>
                Concentration_tmp(&Model.GetConcentration_aer()
                                  (s, b, 0, 0, 0),
                                  shape(this->base_Nz, this->base_Ny,
                                        this->base_Nx));
              FormatBinary<float>().Write(Concentration_tmp, file);
            }

        // Update the date_file.
        date_stream_out.open(date_file.c_str(), ios::trunc);
        date_stream_out << "Iteration saved: #" << this_iteration << endl;
        date_stream_out << "Date of saved iteration: "
                        << this_date.GetDate("%y-%m-%d+%h-%i") << endl;
        date_stream_out.close();

        this->counter = 0;
      }


    // When last iteration remove buffers.
    if (this_iteration == last_iteration)
      {
        for (s = 0; s < Ns_aer; s++)
          for (b = 0; b < Nbin_aer; b++)
            {
              file = output_file[s][b] + ".buf";
              remove(file.c_str());
            }
        file = date_file + ".buf";
        remove(file.c_str());
      }
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITBACKUP_AER_CXX
#endif
