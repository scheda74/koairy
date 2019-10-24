// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Ir�ne Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_AER_CXX


#include "PointEmissionUnit_aer.hxx"


namespace Polyphemus
{


  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! Main constructor.
  template<class T>
  PointEmissionUnit_aer<T>::PointEmissionUnit_aer():
    PointEmissionUnit<T>()
  {
    // velocity_ = 0.;
    // temperature_ = 0.;
    // diameter_ = 0.;
  }


  //! Destructor.
  template<class T>
  PointEmissionUnit_aer<T>::~PointEmissionUnit_aer()
  {
  }


  //! Model initialization.
  /*! It reads the configuration and allocates memory.
   */
  template<class T>
  void PointEmissionUnit_aer<T>::Init(ConfigStream& config,
                                      vector<string> species_list,
                                      vector<string> species_list_aer)
  {
    //    PointEmissionUnit<T>::Init(config, species_list);
    // config.PeekValue("Type", type);
    // date_beg = config.PeekValue("Date_beg");
    //    cout << "Enter PointEmissionUnit_aer::Init !! " << endl;
    config.Find("Species_aer");
    vector<string> emitted_species_list_aer =  split(config.GetLine());
    vector<string> str_split;

    Ns_emis_aer = emitted_species_list_aer.size();
    for (int s = 0; s < Ns_emis_aer; s++)
      {
        emitted_species_list_aer_bin.push_back(map<string, string>());
        str_split = split(emitted_species_list_aer[s], "_");
        emitted_species_list_aer_bin[s]["species"] = str_split[0];
        emitted_species_list_aer_bin[s]["bin"] = str_split[0];

        emitted_species_index_aer.push_back(GetSpeciesIndex_aer(str_split[0], species_list_aer));
        cout << str_split[0] << ", " << emitted_species_index_aer[s] << endl;
        //        throw string("\"PointEmissionUnit_aer<T>::Init()\"")
        //       + " is not defined.";
      }
  }


  // //! Sets the type of emission.
  // /*!
  //   \param type_name The emission type.
  // */
  // template<class T>
  // void PointEmissionUnit<T>::SetType(string type_name)
  // {
  //   type = type_name;
  // }


  // //! Sets the emission coordinates.
  // /*!
  //   Sets the emission coordinates.
  //   \param abscissa the source abscissa.
  //   \param ordinate the source ordinate.
  //   \param height the source height.
  // */
  // template<class T>
  // void PointEmissionUnit<T>::SetCoordinates(T abscissa, T ordinate, T height)
  // {
  //   throw string("\"PointEmissionUnit<T>::SetCoordinates(T&, T&, T&)\"")
  //     + " is not defined.";
  // }


  // //! Gets the emission type.
  // /*!
  //   \return The emission type.
  // */
  // template<class T>
  // string PointEmissionUnit<T>::GetType() const
  // {
  //   return type;
  // }


  // //! Gets the emission coordinates.
  // /*!
  //   Gets the emission coordinates.
  //   \param abscissa the source abscissa.
  //   \param ordinate the source ordinate.
  //   \param height the source height.
  // */
  // template<class T>
  // void PointEmissionUnit<T>::GetCoordinates(T& abscissa, T& ordinate,
  //                  T& height) const
  // {
  //   throw string("\"PointEmissionUnit<T>::GetCoordinates(T&, T&, T&)\"")
  //     + " is not defined.";
  // }


  // //! Gets the list of emission coordinates.
  // /*!
  //   Gets the list of emission coordinates.
  //   \param coordinate_list the list of coordinates.
  // */
  // template<class T>
  // void PointEmissionUnit<T>::GetCoordinates(list<Array<T, 1> >&
  //                                           coordinate_list) const
  // {
  //   throw string("\"PointEmissionUnit<T>::GetCoordinates")
  //     + "(list<Array<T ,6> >&)\" is not defined.";
  // }


  // //! Gets the emission beginning date.
  // /*!
  //   \return The emission beginning date.
  // */
  // template<class T>
  // Date PointEmissionUnit<T>::GetDateBeg() const
  // {
  //   return date_beg;
  // }


  // //! Returns the index in a given vector of a species.
  // /*! If the species is not found, an exception is thrown.
  //   \param species the species name.
  //   \param ref_species_list a vector of species names.
  //   \return The index in 'ref_species_list' of the species named \a species.
  // */
  // template<class T>
  // int PointEmissionUnit<T>::GetSpeciesIndex(string species,
  //                  vector<string>& ref_species_list)
  //   const
  // {
  //   int index = 0;
  //   int length = ref_species_list.size();
  //   while (index < length && ref_species_list[index] != species)
  //     index++;
  //   //    cout << "PointEmissionUnit.cxx: " << species << " " << index << endl;
  //   if (index == length)
  //     throw string("Species \"") + species + "\" unknown.";

  //   return index;
  // }

  //! Returns the index of an aerosol species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \return The index (in the model) of the species named \a species.
  */
  template<class T>
  int PointEmissionUnit_aer<T>::GetSpeciesIndex_aer(string species, vector<string> species_list_aer) const
  {
    return this->GetSpeciesIndex(species, species_list_aer);
  }


  // //! Gets the list of emitted species.
  // /*!
  //   \return The list of emitted species.
  // */
  // template<class T>
  // vector<int> PointEmissionUnit<T>::GetEmittedSpeciesIndex() const
  // {
  //   return emitted_species_index;
  // }


  // //! Gets the timespan of the emission.
  // /*!
  //   \return The emission timespan.
  // */
  // template<class T>
  // T PointEmissionUnit<T>::GetEmissionTimespan(const Date& current_date,
  //                    const Date& next_date) const
  // {
  //   throw string("\"PointEmissionUnit<T>::GetEmissionTimespan")
  //     + "(const Date&, const Date&)\" is not defined.";
  // }


  // //! Checks if the emission has begun at the given date.
  // /*!
  //   \return True if the emission has begun at date date.
  // */
  // template<class T>
  // bool PointEmissionUnit<T>::HasBegun(Date date)
  // {
  //   return date > date_beg;
  // }


  // //! Checks if the emission has begun at the given date.
  // /*!
  //   \return True if the emission has begun at date date.
  // */
  // template<class T>
  // bool PointEmissionUnit<T>::HasEnded(Date date)
  // {
  //   throw string("\"PointEmissionUnit<T>::HasEnded(Date)\"")
  //     + " is not defined.";
  // }


  // //! Gets the species rate.
  // /*!
  //   \return The species rate (mass unit per seconds).
  // */
  // template<class T>
  // T PointEmissionUnit<T>::GetRate(int species) const
  // {
  //   throw string("\"PointEmissionUnit<T>::GetRate(int)\"")
  //     + " is not defined.";
  // }


  // //! Gets emitted quantity for one species during a given time interval.
  // /*!
  //   \param date_beg beginning date of the current time interval
  //   \param date_end ending date of the current time interval
  //   \param s species index
  //   \param point_emission a 2D array containing the source coordinates and the
  //   species quantity.
  // */
  // template<class T>
  // void PointEmissionUnit<T>::GetEmission(Date date_beg, Date date_end, int s,
  //               Array<T, 2>& point_emission)
  // {
  //   throw string("\"PointEmissionUnit<T>::GetEmission(Date, Date, int)\"")
  //     + " is not defined.";
  // }


  // //! Gets the source parameters needed for plume rise computation.
  // /*!
  //   \param velocity the source ejection rate (m/s)
  //   \param temperature the source temperature (K)
  //   \param temperature the source diameter (m)
  // */
  // template<class T>
  // void PointEmissionUnit<T>::GetPlumeRiseParam(T& velocity, T& temperature,
  //                     T& diameter)
  // {
  //   velocity = velocity_;
  //   temperature = temperature_;
  //   diameter = diameter_;
  // }


  // //! Checks if plume rise parameters are given for this source.
  // /*!
  //   \return True if there is plume rise.
  // */
  // template<class T>
  // bool PointEmissionUnit<T>::HasPlumeRise()
  // {
  //   return (velocity_ != 0. && temperature_ != 0. && diameter_ != 0.);
  // }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_AER_CXX
#endif
