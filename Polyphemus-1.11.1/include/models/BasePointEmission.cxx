// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Yelva Roustan, Régis Briant.
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


#ifndef POLYPHEMUS_FILE_MODELS_BASEPOINTEMISSION_CXX


#include "BasePointEmission.hxx"

#include "ContinuousEmission.cxx"
#include "ContinuousEmission_aer.cxx"
#include "ContinuousLineEmission.cxx"
#include "ContinuousLineEmission_aer.cxx"
#include "PuffEmission.cxx"
#include "PuffEmission_aer.cxx"
#include "TemporalEmission.cxx"
#include "TemporalEmission_aer.cxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  BasePointEmission<T>::BasePointEmission()
  {
  }


  //! Destructor.
  template<class T>
  BasePointEmission<T>::~BasePointEmission()
  {
    for (unsigned int i = 0; i < Emission.size(); i++)
      delete Emission[i];
  }


  //! Model initialization.
  /*! It reads the configuration and allocates memory.
   */
  template<class T>
  void BasePointEmission<T>::Init(string config_file,
                                  vector<string> species_list)
  {
    ConfigStream config_stream(config_file);

    // Browses all sections "[source]" that define emissions.
    string line, type;
    int i = 0;
    while (!config_stream.IsEmpty())
      {
        line = config_stream.GetLine();
        if (split(line)[0] == "[source]")
          {
            config_stream.PeekValue("Type", "continuous | puff "
                                    "| temporal | continuous_line", type);
            if (type == "continuous")
              Emission.push_back(new ContinuousEmission<T>());
            else if (type == "puff")
              Emission.push_back(new PuffEmission<T>());
            else if (type == "temporal")
              Emission.push_back(new TemporalEmission<T>());
            else if (type == "continuous_line")
              Emission.push_back(new ContinuousLineEmission<T>);
            Emission[i]->Init(config_stream, species_list);
            i++;
          }
      }
  }

  //! Model initialization for aerosol emissions.
  /*! It reads the configuration and allocates memory.
   */
  template<class T>
  void BasePointEmission<T>::Init(string config_file,
                                  vector<string> species_list,
                                  vector<string> species_list_aer,
                                  int Nbin_aer)
  {
    ConfigStream config_stream(config_file);

    // Browses all sections "[source]" that define emissions.
    string line, type;
    int i = 0;
    while (!config_stream.IsEmpty())
      {
        line = config_stream.GetLine();
        if (split(line)[0] == "[source]")
          {
            config_stream.PeekValue("Type", "continuous | puff "
                                    "| temporal | continuous_line", type);
            if (type == "continuous")
              Emission.push_back(new ContinuousEmission_aer<T>());
            else if (type == "puff")
              Emission.push_back(new PuffEmission_aer<T>());
            else if (type == "temporal")
              Emission.push_back(new TemporalEmission_aer<T>());
            else if (type == "continuous_line")
              Emission.push_back(new ContinuousLineEmission_aer<T>);
            Emission[i]->Init(config_stream, species_list, species_list_aer, Nbin_aer);
            i++;
          }
      }
  }

  //! Gets the number of emissions.
  /*!
    \return The number of emissions.
  */
  template<class T>
  int BasePointEmission<T>::GetNumberEmission()
  {
    return Emission.size();
  }


  //! Gets the emission coordinates.
  /*!
    Gets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
    \param emission the emission index.
  */
  template<class T>
  void BasePointEmission<T>::GetEmissionCoordinates(T& abscissa, T& ordinate,
                                                    T& height,
                                                    int emission) const
  {
    Emission[emission]->GetCoordinates(abscissa, ordinate, height);
  }

  //Gets the source ID
  template<class T>
  string BasePointEmission<T>::GetEmissionSourceId(int emission)
  {
    return Emission[emission]->GetSourceId();
  }

  //Gets the source water
  template<class T>
  T BasePointEmission<T>::GetEmissionSourceWater(int emission) const
  {
    return Emission[emission]->GetSourceWater();
  }


  //! Gets the list of emission coordinates.
  /*!
    Gets the list of emission coordinates.
    \param coordinate_list the list of coordinates.
  */
  template<class T>
  void BasePointEmission<T>::GetEmissionCoordinates(list<Array<T, 1> >&
                                                    coordinate_list,
                                                    int emission)
  {
    Emission[emission]->GetCoordinates(coordinate_list);
  }


  //! Sets the emission coordinates.
  /*!
    Sets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
    \param emission the emission index.
  */
  template<class T>
  void BasePointEmission<T>::SetEmissionCoordinates(T abscissa, T ordinate,
                                                    T height, int emission)
  {
    Emission[emission]->SetCoordinates(abscissa, ordinate, height);
  }


  //! Sets the list of emission coordinates.
  /*!
    \param[in] coordinate_list the list of coordinates.
  */
  template<class T>
  void BasePointEmission<T>::SetEmissionCoordinates(const list<Array<T, 1> >&
                                                    coordinate_list,
                                                    int emission)
  {
    Emission[emission]->SetCoordinates(coordinate_list);
  }


  //! Gets the type of an emission.
  /*!
    \return The type of an emission.
  */
  template<class T>
  string BasePointEmission<T>::GetEmissionType(int emission)
  {
    return Emission[emission]->GetType();
  }


  //! Gets the beginning date of an emission.
  /*!
    \return The beginning date of an emission.
  */
  template<class T>
  Date BasePointEmission<T>::GetDateBeg(int emission)
  {
    return Emission[emission]->GetDateBeg();
  }


  //! Returns the ending date of the given emission.
  /*!
    \param[in] emission the emission index.
    \return The ending date of the emission.
  */
  template<class T>
  Date BasePointEmission<T>::GetDateEnd(int emission)
  {
    return Emission[emission]->GetDateEnd();
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool BasePointEmission<T>::HasBegun(Date date, int emission)
  {
    return Emission[emission]->HasBegun(date);
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool BasePointEmission<T>::HasEnded(Date date, int emission)
  {
    return Emission[emission]->HasEnded(date);
  }


  //! Checks if the source is emitting during a time interval.
  /*!
    \return True if source is emitting.
  */
  template<class T>
  bool BasePointEmission<T>::IsEmitting(Date current_date,
                                        Date next_date,
                                        int emission)
  {
    return Emission[emission]->HasBegun(next_date)
      && !Emission[emission]->HasEnded(current_date);
  }


  //! Gets the list of emitted species.
  /*!
    \return The list of emitted species.
  */
  template<class T>
  vector<int> BasePointEmission<T>::GetEmittedSpeciesIndex(int emission) const
  {
    return Emission[emission]->GetEmittedSpeciesIndex();
  }

  //! Gets the list of emitted aerosol species.
  /*!
    \return The list of emitted aerosol species.
  */
  template<class T>
  vector<int> BasePointEmission<T>::GetEmittedSpeciesIndex_aer(int emission) const
  {
    return Emission[emission]->GetEmittedSpeciesIndex_aer();
  }

  //! Gets the map-type list of emitted aerosol species.
  /*!
    \return The map-type list of emitted aerosol species .
  */
  template<class T>
  vector<map <string, string> > BasePointEmission<T>::GetEmittedAerosolSpecies(int emission) const
  {
    return Emission[emission]->GetEmittedAerosolSpecies();
  }


  //! Gets the species rate.
  /*!
    \return The species rate.
  */
  template<class T>
  T BasePointEmission<T>::GetRate(int emission, int species) const
  {
    return Emission[emission]->GetRate(species);
  }


  //! Gets the species rate for a given source section.
  /*!
    \return The species rate.
  */
  template<class T>
  T BasePointEmission<T>::GetRate(int emission, int species,
                                  int id_section) const
  {
    return Emission[emission]->GetRate(species, id_section);
  }


  //! Multiplies all emission rates by some value.
  /*!
    \param[in] factor the factor by which to multiply the rates.
  */
  template<class T>
  void BasePointEmission<T>::MultiplyRate(T factor)
  {
    for (int emission = 0; emission < int(Emission.size()); emission++)
      Emission[emission]->MultiplyRate(factor);
  }


  //! Multiplies all emission rates of a given species by some value.
  /*!
    \param[in] species index of the species.
    \param[in] factor the factor by which to multiply the rates.
  */
  template<class T>
  void BasePointEmission<T>::MultiplyRate(int species, T factor)
  {
    for (int emission = 0; emission < int(Emission.size()); emission++)
      Emission[emission]->MultiplyRate(species, factor);
  }


  //! Shifts the emission times by some value.
  /*! The initial and final release times are both shifted.
    \param[in] shift the time shift in seconds (a positive value will postpone
    the release).
  */
  template<class T>
  void BasePointEmission<T>::ApplyTimeShift(T shift)
  {
    for (int emission = 0; emission < int(Emission.size()); emission++)
      Emission[emission]->ApplyTimeShift(shift);
  }


  //! Shifts the emission altitudes by some value.
  /*!
    \param[in] shift the altitude shift in meters.
  */
  template<class T>
  void BasePointEmission<T>::ApplyAltitudeShift(T shift)
  {
    for (int emission = 0; emission < int(Emission.size()); emission++)
      Emission[emission]->ApplyAltitudeShift(shift);
  }


  //! Gets emitted quantities for one species.
  /*!
    \return The emitted quantities for one species.
  */
  template<class T>
  void BasePointEmission<T>::GetEmission(Date date_beg, Date date_end,
                                         int species, int emission,
                                         Array<T, 2>& point_emission)
  {
    Emission[emission]->GetEmission(date_beg, date_end, species,
                                    point_emission);
  }

  //! Gets emitted quantities for one aerosol species.
  /*!
    \return The emitted quantities for one aerosol species.
  */
  template<class T>
  void BasePointEmission<T>::GetEmission_aer(Date date_beg, Date date_end,
                                             int species, int emission,
                                             T& quantity_aer_bin)
  {
    Emission[emission]->GetEmission_aer(date_beg, date_end, species,
                                        quantity_aer_bin);
  }


  //! Checks if plume rise parameters are given for this source.
  /*!
    \return True if there is plume rise.
  */
  template<class T>
  bool BasePointEmission<T>::HasPlumeRise(int emission)
  {
    return Emission[emission]->HasPlumeRise();
  }


  //! Gets the source parameters needed for plume rise computation.
  /*!
    \return The source parameters needed for plume rise computation.
  */
  template<class T>
  void BasePointEmission<T>::GetPlumeRiseParam(T& velocity, T& temperature,
                                               T& diameter, int emission)
  {
    Emission[emission]->GetPlumeRiseParam(velocity, temperature, diameter);
  }


  //! Gets the emission timespan in a time interval.
  /*!
    \return The emission timespan.
  */
  template<class T>
  T BasePointEmission<T>::GetEmissionTimespan(const Date& current_date,
                                              const Date& next_date,
                                              int emission) const
  {
    return Emission[emission]->GetEmissionTimespan(current_date, next_date);
  }


  //! Gets the width of the line source (works for line sources only).
  /*!
    \param[in] emission the emission index.
    \param[in] id_section the line source section index.
    \return The line source width.
    \warning Only works for line sources.
  */
  template<class T>
  T BasePointEmission<T>::GetWidth(int emission, int id_section)
  {
    return Emission[emission]->GetWidth(id_section);
  }


  //! Gets the Vehicle Velocity in case of a line source (m/s).
  /*!
    \param emission point emission index.
    \param id_section identifier of the section to be considered.
    \returns The line source Vehicle Velocity.
    \warning Only valid for line sources.
  */
  template<class T>
  T BasePointEmission<T>::GetVehicleVelocity(int emission, int id_section)
  {
    return Emission[emission]->GetVehicleVelocity(id_section);
  }


  //! Gets the area in case of a line source (m²).
  /*!
    \param emission point emission index.
    \param id_section identifier of the section to be considered.
    \returns The line source area.
    \warning Only valid for line sources.
  */
  template<class T>
  T BasePointEmission<T>::GetArea(int emission, int id_section)
  {
    return Emission[emission]->GetArea(id_section);
  }


  //! Gets the Density in case of a line source (veh/m).
  /*!
    \param emission point emission index.
    \param id_section identifier of the section to be considered.
    \returns The line source density.
    \warning Only valid for line sources.
  */
  template<class T>
  T BasePointEmission<T>::GetDensity(int emission, int id_section)
  {
    return Emission[emission]->GetDensity(id_section);
  }


  //! Adds a new continuous emission.
  /*!
    \param[in] abscissa: abscissa (m).
    \param[in] ordinate: ordinate (m).
    \param[in] height: height (m).
    \param[in] rate: the list of release rates (mass unit per seconds).
    \param[in] date beg: the release date.
    \param[in] date end: the release ending date.
    \param[in] species: list of species names.
    \param[in] diameter: the source diameter (m).
    \param[in] velocity: the efflux speed (m/s).
    \param[in] temperature: the temperature of emissions (Celsius degrees).
    \return The line source width.
  */
  template<class T>
  void BasePointEmission<T>::AddContinuousEmission(T abscissa, T ordinate,
                                                   T height,
                                                   const vector<T>& rate,
                                                   Date date_beg,
                                                   Date date_end,
                                                   const vector<int>& species,
                                                   T diameter, T velocity,
                                                   T temperature)
  {
    Emission.push_back(new ContinuousEmission<T>());
    Emission.back()->Init(abscissa, ordinate, height, rate,
                          date_beg, date_end, species, diameter,
                          velocity, temperature);
  }


  //! Checks whether the source is a volume source or a point source.
  /*!
    \return True if the source is a volume source.
  */
  template<class T>
  bool BasePointEmission<T>::IsVolumeSource(int emission)
  {
    return Emission[emission]->IsVolumeSource();
  }


  //! Gets the source parameters needed for volume source.
  /*!
    \return The source parameters needed for volume source.
  */
  template<class T>
  void BasePointEmission<T>::GetVolumeSource(T& width, T& length, int emission)
  {
    Emission[emission]->GetVolumeSource(width, length);
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_BASEPOINTEMISSION_CXX
#endif
