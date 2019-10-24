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


#ifndef POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_CXX


#include "PointEmissionUnit.hxx"
#include "AtmoData.hxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  PointEmissionUnit<T>::PointEmissionUnit()
  {
    velocity_ = 0.;
    temperature_ = 0.;
    diameter_ = 0.;
    width_ = 0.;
    length_ = 0.;
    is_volume_source_ = false;
  }


  //! Destructor.
  template<class T>
  PointEmissionUnit<T>::~PointEmissionUnit()
  {
  }


  //! Model initialization.
  /*! It reads the configuration and allocates memory.
    \param[in] config the ConfigStream instance containing the sources parameters.
    \param[in] species_list the list of species names.
  */
  template<class T>
  void PointEmissionUnit<T>::Init(ConfigStream& config,
                                  vector<string> species_list)
  {
    config.PeekValue("ID", source_id);
    config.PeekValue("Source_water", source_water);
    config.PeekValue("Type", type);
    date_beg = config.PeekValue("Date_beg");
    config.Find("Species");
    vector<string> emitted_species_list =  split(config.GetLine());
    Ns_emis = emitted_species_list.size();
    for (int s = 0; s < Ns_emis; s++)
      emitted_species_index.push_back(GetSpeciesIndex(emitted_species_list[s],
                                                      species_list));
  }

  //! Model initialization with aerosol emissions.
  /*! It reads the configuration and allocates memory.
   */
  template<class T>
  void PointEmissionUnit<T>::Init(ConfigStream& config,
                                  vector<string> species_list,
                                  vector<string> species_list_aer,
                                  int Nbin_aer)
  {
    config.Find("Species_aer");
    vector<string> emitted_species_list_aer =  split(config.GetLine());
    vector<string> str_split;

    Ns_emis_aer = emitted_species_list_aer.size();
    for (int s = 0; s < Ns_emis_aer; s++)
      {
        emitted_species_list_aer_bin.push_back(map<string, string>());
        str_split = split(emitted_species_list_aer[s], "_");
        emitted_species_list_aer_bin[s]["species"] = str_split[0];
        emitted_species_list_aer_bin[s]["bin"] = str_split[1];

        emitted_species_index_aer.push_back(GetSpeciesIndex_aer(str_split[0], species_list_aer));
        if (convert<int>(str_split[1]) >= Nbin_aer)
          throw string("Given aerosol bin no \"") + str_split[1] + "\" is wrong. It should be less than \"" + to_str(Nbin_aer) + "\" in \"PointEmissionUnit<T>::Init(ConfigStream&, vector<string>, vector<string>, int)\"";
      }
  }


  //! Model initialization.
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
  */
  template<class T>
  void PointEmissionUnit<T>::Init(T abscissa, T ordinate, T height,
                                  const vector<T>& rate, Date date_beg,
                                  Date date_end, const vector<int>& species,
                                  T diameter, T velocity, T temperature)
  {
    throw "'PointEmissionUnit<T>::Init(...)' is undefined.";
  }


  //! Sets the type of emission.
  /*!
    \param type_name The emission type.
  */
  template<class T>
  void PointEmissionUnit<T>::SetType(string type_name)
  {
    type = type_name;
  }


  //! Sets the emission coordinates.
  /*!
    Sets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
  */
  template<class T>
  void PointEmissionUnit<T>::SetCoordinates(T abscissa, T ordinate, T height)
  {
    throw string("\"PointEmissionUnit<T>::SetCoordinates(T&, T&, T&)\"")
      + " is not defined.";
  }


  //! Sets the emission coordinates.
  /*!
    \param coordinate_list the list of coordinates.
  */
  template<class T>
  void PointEmissionUnit<T>::SetCoordinates(const list<Array<T, 1> >&
                                            coordinate_list)
  {
    throw "'PointEmissionUnit<T>::SetCoordinates(list<Array<T, 1> >&"
      " coordinate_list)' is undefined.";
  }


  //! Gets the emission type.
  /*!
    \return The emission type.
  */
  template<class T>
  string PointEmissionUnit<T>::GetType() const
  {
    return type;
  }


  // Get source Id
  template<class T>
  string PointEmissionUnit<T>::GetSourceId() const
  {
    return source_id;
  }

  // Get source Water
  template<class T>
  T PointEmissionUnit<T>::GetSourceWater() const
  {
    return source_water;
  }

  //! Gets the emission coordinates.
  /*!
    Gets the emission coordinates.
    \param abscissa the source abscissa.
    \param ordinate the source ordinate.
    \param height the source height.
  */
  template<class T>
  void PointEmissionUnit<T>::GetCoordinates(T& abscissa, T& ordinate,
                                            T& height) const
  {
    throw string("\"PointEmissionUnit<T>::GetCoordinates(T&, T&, T&)\"")
      + " is not defined.";
  }


  //! Gets the list of emission coordinates.
  /*!
    Gets the list of emission coordinates.
    \param coordinate_list the list of coordinates.
  */
  template<class T>
  void PointEmissionUnit<T>::GetCoordinates(list<Array<T, 1> >&
                                            coordinate_list) const
  {
    throw string("\"PointEmissionUnit<T>::GetCoordinates")
      + "(list<Array<T ,6> >&)\" is not defined.";
  }


  //! Gets the emission beginning date.
  /*!
    \return The emission beginning date.
  */
  template<class T>
  Date PointEmissionUnit<T>::GetDateBeg() const
  {
    return date_beg;
  }

  //! Gets the emission ending date.
  /*!
    \return The emission ending date.
  */
  template<class T>
  Date PointEmissionUnit<T>::GetDateEnd() const
  {
    return date_end;
  }


  //! Returns the index in a given vector of a species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \param ref_species_list a vector of species names.
    \return The index in 'ref_species_list' of the species named \a species.
  */
  template<class T>
  int PointEmissionUnit<T>::GetSpeciesIndex(string species,
                                            vector<string>& ref_species_list)
    const
  {
    int index = 0;
    int length = ref_species_list.size();
    while (index < length && ref_species_list[index] != species)
      index++;
    if (index == length)
      throw string("Species \"") + species + "\" unknown.";
    return index;
  }

  //! Returns the index of an aerosol species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \return The index (in the model) of the species named \a species.
  */
  template<class T>
  int PointEmissionUnit<T>::GetSpeciesIndex_aer(string species, vector<string> species_list_aer) const
  {
    return this->GetSpeciesIndex(species, species_list_aer);
  }

  //! Gets the list of emitted species.
  /*!
    \return The list of emitted species.
  */
  template<class T>
  vector<int> PointEmissionUnit<T>::GetEmittedSpeciesIndex() const
  {
    return emitted_species_index;
  }

  //! Gets the list of emitted species.
  /*!
    \return The list of emitted species.
  */
  template<class T>
  vector<int> PointEmissionUnit<T>::GetEmittedSpeciesIndex_aer() const
  {
    return emitted_species_index_aer;
  }

  //! Gets the map-type list of emitted aerosol species.
  /*!
    \return The map-type list of emitted aerosol species .
  */
  template<class T>
  vector<map <string, string> > PointEmissionUnit<T>::GetEmittedAerosolSpecies() const
  {
    return emitted_species_list_aer_bin;
  }



  //! Gets the timespan of the emission.
  /*!
    \return The emission timespan.
  */
  template<class T>
  T PointEmissionUnit<T>::GetEmissionTimespan(const Date& current_date,
                                              const Date& next_date) const
  {
    throw string("\"PointEmissionUnit<T>::GetEmissionTimespan")
      + "(const Date&, const Date&)\" is not defined.";
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool PointEmissionUnit<T>::HasBegun(Date date)
  {
    return date > date_beg;
  }


  //! Checks if the emission has begun at the given date.
  /*!
    \return True if the emission has begun at date date.
  */
  template<class T>
  bool PointEmissionUnit<T>::HasEnded(Date date)
  {
    throw string("\"PointEmissionUnit<T>::HasEnded(Date)\"")
      + " is not defined.";
  }


  //! Gets the species rate.
  /*!
    \return The species rate (mass unit per seconds).
  */
  template<class T>
  T PointEmissionUnit<T>::GetRate(int species) const
  {
    throw string("\"PointEmissionUnit<T>::GetRate(int)\"")
      + " is not defined.";
  }


  //! Gets the species rate for a given source section.
  /*!
    \return The species rate (mass unit per seconds).
  */
  template<class T>
  T PointEmissionUnit<T>::GetRate(int species, int id_section) const
  {
    throw "'PointEmissionUnit<T>::GetRate(int, int)'"
      " is not defined.";
  }


  //! Multiplies all rates by some factor.
  /*!
    \param[in] factor the factor by which to multiply all emission rates.
  */
  template<class T>
  void PointEmissionUnit<T>::MultiplyRate(T factor)
  {
    throw string("\"PointEmissionUnit<T>::MultiplyRate(T)\"")
      + " is not defined.";
  }


  //! Multiplies all rates of a given species by some factor.
  /*!
    \param[in] species index of the species.
    \param[in] factor the factor by which to multiply all emission rates.
  */
  template<class T>
  void PointEmissionUnit<T>::MultiplyRate(int species, T factor)
  {
    throw string("\"PointEmissionUnit<T>::MultiplyRate(int, T)\"")
      + " is not defined.";
  }


  //! Shifts the emission times by some value.
  /*! The initial and final release times are both shifted.
    \param[in] shift the time shift in seconds (a positive value will postpone
    the release).
  */
  template<class T>
  void PointEmissionUnit<T>::ApplyTimeShift(T shift)
  {
    throw string("\"PointEmissionUnit<T>::ApplyTimeShift(T)\"")
      + " is not defined.";
  }


  //! Shifts the emission altitudes by some value.
  /*!
    \param[in] shift the altitude shift in meters.
  */
  template<class T>
  void PointEmissionUnit<T>::ApplyAltitudeShift(T shift)
  {
    throw string("\"PointEmissionUnit<T>::ApplyAltitudeShift(T)\"")
      + " is not defined.";
  }


  //! Gets emitted quantity for one species during a given time interval.
  /*!
    \param date_beg beginning date of the current time interval
    \param date_end ending date of the current time interval
    \param s species index
    \param point_emission a 2D array containing the source coordinates and the
    species quantity.
  */
  template<class T>
  void PointEmissionUnit<T>::GetEmission(Date date_beg, Date date_end, int s,
                                         Array<T, 2>& point_emission)
  {
    throw "'PointEmissionUnit<T>::GetEmission(Date, Date, int, Array<T, 2>)'"
      " is undefined.";
  }

  //! Gets emitted quantity for one species during a given time interval.
  /*!
    \param date_beg beginning date of the current time interval
    \param date_end ending date of the current time interval
    \param s species index
    \param quantity species quantity.
  */
  template<class T>
  void PointEmissionUnit<T>::GetEmission_aer(Date date_beg, Date date_end, int s, T& quantity)
  {
    throw string("\"PointEmissionUnit<T>::GetEmission_aer(Date, Date, int)\"")
      + " is not defined.";
  }


  //! Gets the source parameters needed for plume rise computation.
  /*!
    \param velocity the source ejection rate (m/s)
    \param temperature the source temperature (K)
    \param diameter the source diameter (m)
  */
  template<class T>
  void PointEmissionUnit<T>::GetPlumeRiseParam(T& velocity, T& temperature,
                                               T& diameter)
  {
    velocity = velocity_;
    temperature = temperature_;
    diameter = diameter_;
  }


  //! Checks if plume rise parameters are given for this source.
  /*!
    \return True if there is plume rise.
  */
  template<class T>
  bool PointEmissionUnit<T>::HasPlumeRise()
  {
    return (velocity_ != 0. && temperature_ != 0. && diameter_ != 0.);
  }


  //! Gets the width of a line source at the given section (m).
  /*!
    \param The identifier of the line source section to consider.
    \return The source with (m).
    \warning Is only defined for lines sources.
  */
  template<class T>
  T PointEmissionUnit<T>::GetWidth(int id_section) const
  {
    throw "'PointEmissionUnit<T>::GetWidth(int)' is undefined.";
  }


  //! Returns the source identifier.
  /*!
    \return  the source identifier.
  */
  template<class T>
  int PointEmissionUnit<T>::GetIdSource() const
  {
    throw "'PointEmissionUnit<T>::GetIdSource()' is undefined.";
  }


  //! Gets the id of the source section.
  /*!
    \param Id_section.
  */
  template<class T>
  int PointEmissionUnit<T>::GetIdSection() const
  {
    throw "'PointEmissionUnit<T>::GetIdSection()' is undefined.";
  }


  //! Gets the VehicleVelocity of a line source at the given section (m/s).
  /*!
    \param The identifier of the source section to consider.
    \returns VehicleVelocity of a line source.
    \warning Is only defined for lines sources.
  */
  template<class T>
  T PointEmissionUnit<T>::GetVehicleVelocity(int id_section) const
  {
    throw "'PointEmissionUnit<T>::GetVehicleVelocity(int)' is undefined.";
  }


  //! Gets the Area of a line source at the given section (m²).
  /*!
    \param The identifier of the source section to consider.
    \returns area of a line source.
    \warning Is only defined for lines sources.
  */
  template<class T>
  T PointEmissionUnit<T>::GetArea(int id_section) const
  {
    throw "'PointEmissionUnit<T>::GetArea(int)' is undefined.";
  }


  //! Gets the Density of a line source at the given section (veh/m).
  /*!
    \param The identifier of the source section to consider.
    \returns Density of a line source.
    \warning Is only defined for lines sources.
  */
  template<class T>
  T PointEmissionUnit<T>::GetDensity(int id_section) const
  {
    throw "'PointEmissionUnit<T>::GetDensity(int)' is undefined.";
  }


  //! Checks whether the source is a volume source or a point source.
  /*!
    \return True if the source is a volume source.
  */
  template<class T>
  bool PointEmissionUnit<T>::IsVolumeSource()
  {
    return (is_volume_source_);
  }


  //! Gets the source parameters needed for volume source.
  /*!
    \param width the source width (m)
    \param length the source lenth (m)
  */
  template<class T>
  void PointEmissionUnit<T>::GetVolumeSource(T& width, T& length)
  {
    width = width_;
    length = length_;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_CXX
#endif
